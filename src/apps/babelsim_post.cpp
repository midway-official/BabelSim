#include "babelsim/case.h"
#include "babelsim/mesh_io.h"
#include "babelsim/result_reader.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace babelsim {
namespace {

struct Arguments {
    std::filesystem::path case_directory;
    std::string time_name;
    std::vector<std::string> formats;
};

Arguments parseArguments(int argc, char* argv[]) {
    Arguments result;
    for (int index = 1; index < argc;) {
        const std::string option = argv[index++];
        if (option == "-case" || option == "-time") {
            if (index == argc) throw std::invalid_argument("missing value for " + option);
            if (option == "-case") result.case_directory = argv[index];
            else result.time_name = argv[index];
            ++index;
        } else if (option == "-format") {
            while (index < argc && argv[index][0] != '-') result.formats.emplace_back(argv[index++]);
        } else {
            throw std::invalid_argument("unknown option " + option);
        }
    }
    if (result.case_directory.empty() || result.formats.empty()) {
        throw std::invalid_argument(
            "usage: babelsim-post -case <case-directory> [-time <name|latest>] -format <vtk|tecplot> [...]");
    }
    return result;
}

std::filesystem::path resolveTimeDirectory(
    const CaseDefinition& definition,
    const OutputControl& output,
    const std::string& requested)
{
    const auto result_root = definition.root / output.directory;
    if (requested.empty()) return result_root / output.time_name;
    if (requested != "latest") return result_root / requested;
    std::filesystem::path latest;
    for (const auto& entry : std::filesystem::directory_iterator(result_root)) {
        if (entry.is_directory() && (latest.empty() || entry.path().filename() > latest.filename())) {
            latest = entry.path();
        }
    }
    if (latest.empty()) throw std::runtime_error("case has no saved result times");
    return latest;
}

std::string componentName(const ResultField& field, int component) {
    if (field.components == 1) return field.info.name;
    static constexpr const char* suffix[] = {"x", "y", "z", "xx", "xy", "xz", "yx", "yy", "yz"};
    return field.info.name + "_" + suffix[component];
}

void writeVtk(const std::filesystem::path& path, const Mesh& mesh, const ResultData& data) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create VTK output: " + path.string());
    output << "# vtk DataFile Version 3.0\nBabelSim " << data.time_name
           << "\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << mesh.vertexCount()
           << " double\n" << std::setprecision(17);
    for (const Vec3& point : mesh.vertices) output << point.x << ' ' << point.y << ' ' << point.z << '\n';
    output << "CELLS " << mesh.cellCount() << ' ' << mesh.cellCount() * 9 << '\n';
    for (Index k = 0; k < mesh.dimensions[2]; ++k) {
        for (Index j = 0; j < mesh.dimensions[1]; ++j) {
            for (Index i = 0; i < mesh.dimensions[0]; ++i) {
                output << "8 " << mesh.vertexId(i, j, k) << ' ' << mesh.vertexId(i + 1, j, k) << ' '
                       << mesh.vertexId(i + 1, j + 1, k) << ' ' << mesh.vertexId(i, j + 1, k) << ' '
                       << mesh.vertexId(i, j, k + 1) << ' ' << mesh.vertexId(i + 1, j, k + 1) << ' '
                       << mesh.vertexId(i + 1, j + 1, k + 1) << ' ' << mesh.vertexId(i, j + 1, k + 1) << '\n';
            }
        }
    }
    output << "CELL_TYPES " << mesh.cellCount() << '\n';
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) output << "12\n";
    output << "CELL_DATA " << mesh.cellCount() << '\n';
    for (const ResultField& field : data.fields) {
        if (field.components == 1) {
            output << "SCALARS " << field.info.name << " double 1\nLOOKUP_TABLE default\n";
            for (Index cell = 0; cell < mesh.cellCount(); ++cell) output << field.values[cell] << '\n';
        } else if (field.components == 3) {
            output << "VECTORS " << field.info.name << " double\n";
            for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
                const std::size_t begin = static_cast<std::size_t>(cell) * 3;
                output << field.values[begin] << ' ' << field.values[begin + 1] << ' '
                       << field.values[begin + 2] << '\n';
            }
        } else {
            output << "TENSORS " << field.info.name << " double\n";
            for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
                const std::size_t begin = static_cast<std::size_t>(cell) * 9;
                for (int row = 0; row < 3; ++row) {
                    output << field.values[begin + 3 * row] << ' '
                           << field.values[begin + 3 * row + 1] << ' '
                           << field.values[begin + 3 * row + 2] << '\n';
                }
            }
        }
    }
}

void writeTecplot(const std::filesystem::path& path, const Mesh& mesh, const ResultData& data) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create Tecplot output: " + path.string());
    output << "TITLE = \"BabelSim " << data.time_name << "\"\nVARIABLES = \"X\" \"Y\" \"Z\"";
    int variables = 3;
    for (const ResultField& field : data.fields) {
        for (int component = 0; component < field.components; ++component) {
            output << " \"" << componentName(field, component) << "\"";
            ++variables;
        }
    }
    output << "\nZONE N=" << mesh.vertexCount() << ", E=" << mesh.cellCount()
           << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK, VARLOCATION=([4-" << variables
           << "]=CELLCENTERED)\n" << std::setprecision(17);
    for (int coordinate = 0; coordinate < 3; ++coordinate) {
        for (const Vec3& point : mesh.vertices) output << point[coordinate] << '\n';
    }
    for (const ResultField& field : data.fields) {
        for (int component = 0; component < field.components; ++component) {
            for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
                output << field.values[static_cast<std::size_t>(cell) * field.components + component] << '\n';
            }
        }
    }
    for (Index k = 0; k < mesh.dimensions[2]; ++k) {
        for (Index j = 0; j < mesh.dimensions[1]; ++j) {
            for (Index i = 0; i < mesh.dimensions[0]; ++i) {
                output << mesh.vertexId(i, j, k) + 1 << ' ' << mesh.vertexId(i + 1, j, k) + 1 << ' '
                       << mesh.vertexId(i + 1, j + 1, k) + 1 << ' ' << mesh.vertexId(i, j + 1, k) + 1 << ' '
                       << mesh.vertexId(i, j, k + 1) + 1 << ' ' << mesh.vertexId(i + 1, j, k + 1) + 1 << ' '
                       << mesh.vertexId(i + 1, j + 1, k + 1) + 1 << ' ' << mesh.vertexId(i, j + 1, k + 1) + 1 << '\n';
            }
        }
    }
}

int run(const Arguments& arguments) {
    const CaseDefinition definition = readCase(arguments.case_directory);
    const OutputControl output = readOutputControl(definition);
    const auto time_directory = resolveTimeDirectory(definition, output, arguments.time_name);
    const Mesh mesh = readMeshFile(definition.mesh_file);
    const ResultData results = readParallelResults(time_directory, mesh.cellCount());
    if (results.global_dimensions != mesh.dimensions) {
        throw std::runtime_error("result global dimensions do not match the case mesh");
    }
    const auto post_directory = definition.root / "post";
    std::filesystem::create_directories(post_directory);
    for (const std::string& format : arguments.formats) {
        if (format == "vtk") writeVtk(post_directory / (results.time_name + ".vtk"), mesh, results);
        else if (format == "tecplot") writeTecplot(post_directory / (results.time_name + ".dat"), mesh, results);
        else throw std::invalid_argument("unknown output format " + format);
    }
    std::cout << "BabelSim postprocessed " << time_directory << " into " << post_directory << '\n';
    return 0;
}

}  // 匿名命名空间
}  // babelsim 命名空间

int main(int argc, char* argv[]) {
    try {
        return babelsim::run(babelsim::parseArguments(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "babelsim-post: " << error.what() << '\n';
        return 1;
    }
}
