#include "babelsim/parallel_writer.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace babelsim {
namespace {

std::filesystem::path rankDirectory(
    const std::filesystem::path& time_directory,
    int rank)
{
    std::ostringstream name;
    name << "rank-" << std::setfill('0') << std::setw(4) << rank;
    return time_directory / name.str();
}

void ensureDirectory(const std::filesystem::path& time_directory, const ParallelContext& parallel) {
    int error = 0;
    if (parallel.rank == 0) {
        std::error_code code;
        std::filesystem::create_directories(time_directory, code);
        error = code ? 1 : 0;
    }
    if (parallel.maximum(error) != 0) {
        throw std::runtime_error("cannot create result directory: " + time_directory.string());
    }
    parallel.barrier();
    std::error_code code;
    std::filesystem::create_directories(rankDirectory(time_directory, parallel.rank), code);
    if (parallel.maximum(code ? 1 : 0) != 0) {
        throw std::runtime_error("cannot create rank result directory: " + time_directory.string());
    }
}

template <typename T, typename WriteValue>
void writeField(
    const std::filesystem::path& time_directory,
    const Field<T>& field,
    const ParallelContext& parallel,
    int components,
    WriteValue write_value)
{
    const Mesh& mesh = field.mesh();
    if (field.location() != FieldLocation::Cell || field.name().empty()) {
        throw std::invalid_argument("parallel result output requires a named cell field");
    }
    ensureDirectory(time_directory, parallel);
    const auto path = rankDirectory(time_directory, parallel.rank) / (field.name() + ".csv");
    std::ofstream output(path);
    if (parallel.maximum(output ? 0 : 1) != 0) {
        throw std::runtime_error("cannot create field output: " + path.string());
    }
    output << "global_id,x,y,z";
    for (int component = 0; component < components; ++component) {
        output << ",value" << component;
    }
    output << '\n' << std::setprecision(17);
    for (Index cell : mesh.owned_cells) {
        const Vec3& centre = mesh.cellCentre(cell);
        output << mesh.globalCellId(cell) << ',' << centre.x << ',' << centre.y << ',' << centre.z;
        write_value(output, field[cell]);
        output << '\n';
    }
    output.close();
    if (parallel.maximum(output ? 0 : 1) != 0) {
        throw std::runtime_error("cannot flush field output: " + path.string());
    }
    parallel.barrier();
}

const char* locationName(FieldLocation location) {
    switch (location) {
        case FieldLocation::Cell: return "cell";
        case FieldLocation::Face: return "face";
        case FieldLocation::Vertex: return "vertex";
    }
    throw std::invalid_argument("unknown field location");
}

}  // 匿名命名空间

void writeOwnedFieldCsv(
    const std::filesystem::path& time_directory,
    const ScalarField& field,
    const ParallelContext& parallel)
{
    writeField(time_directory, field, parallel, 1, [](std::ostream& output, double value) {
        output << ',' << value;
    });
}

void writeOwnedFieldCsv(
    const std::filesystem::path& time_directory,
    const VectorField& field,
    const ParallelContext& parallel)
{
    writeField(time_directory, field, parallel, 3, [](std::ostream& output, const Vec3& value) {
        output << ',' << value.x << ',' << value.y << ',' << value.z;
    });
}

void writeOwnedFieldCsv(
    const std::filesystem::path& time_directory,
    const TensorField& field,
    const ParallelContext& parallel)
{
    writeField(time_directory, field, parallel, 9, [](std::ostream& output, const Tensor3& value) {
        for (const Vec3& row : value.rows) output << ',' << row.x << ',' << row.y << ',' << row.z;
    });
}

void writeOwnedResultMetadata(
    const std::filesystem::path& time_directory,
    const Mesh& mesh,
    const ParallelContext& parallel,
    const std::string& time_name,
    const std::vector<FieldOutputInfo>& fields)
{
    if (time_name.empty() || fields.empty()) {
        throw std::invalid_argument("result metadata needs a time name and fields");
    }
    ensureDirectory(time_directory, parallel);
    const auto path = rankDirectory(time_directory, parallel.rank) / "metadata.bs";
    std::ofstream output(path);
    if (parallel.maximum(output ? 0 : 1) != 0) {
        throw std::runtime_error("cannot create result metadata: " + path.string());
    }
    output << "format babelsim_result 1\n"
           << "time " << time_name << '\n'
           << "rank " << parallel.rank << '\n'
           << "ranks " << parallel.size << '\n'
           << "global_dimensions " << mesh.global_dimensions[0] << ' '
           << mesh.global_dimensions[1] << ' ' << mesh.global_dimensions[2] << '\n'
           << "owned_cells " << mesh.ownedCellCount() << '\n';
    for (const FieldOutputInfo& field : fields) {
        if (field.name.empty() || field.type.empty()) {
            throw std::invalid_argument("result metadata field is incomplete");
        }
        output << "field " << field.name << ' ' << field.type << ' '
               << locationName(field.location) << '\n';
    }
    output.close();
    if (parallel.maximum(output ? 0 : 1) != 0) {
        throw std::runtime_error("cannot flush result metadata: " + path.string());
    }
    parallel.barrier();
}

}  // babelsim 命名空间
