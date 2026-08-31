#include "internal/mesh_access.h"
#include "babelsim/postprocess.h"
#include "babelsim/case.h"
#include "babelsim/mesh_io.h"
#include "babelsim/result_reader.h"

#include <algorithm>
#include <cmath>
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
            "usage: babelsim-post -case <case-directory> [-time <name|latest|all>] -format <vtk|tecplot> [...]");
    }
    return result;
}

bool numericalTime(const std::string& name, double& time) {
    try {
        std::size_t consumed = 0;
        time = std::stod(name, &consumed);
        return consumed == name.size() && std::isfinite(time);
    } catch (const std::exception&) {
        return false;
    }
}

std::vector<std::filesystem::path> resolveTimeDirectories(
    const CaseDefinition& definition, const OutputControl& output, const std::string& requested)
{
    const std::filesystem::path selection = requested.empty() ? output.time_name : requested;
    if (selection.is_absolute()) throw std::invalid_argument("result selection must be relative");
    for (const auto& part : selection) {
        if (part == ".." || part == ".") throw std::invalid_argument("invalid result selection");
    }
    const auto root = definition.root / output.directory;
    const std::string mode = selection.filename().string();
    if (mode != "latest" && mode != "all") return {root / selection};
    std::vector<std::pair<double, std::filesystem::path>> sorted;
    for (const auto& entry : std::filesystem::directory_iterator(root / selection.parent_path())) {
        double time;
        if (entry.is_directory() && numericalTime(entry.path().filename().string(), time))
            sorted.emplace_back(time, entry.path());
    }
    std::sort(sorted.begin(), sorted.end());
    if (sorted.empty()) throw std::runtime_error("case has no numerical result times");
    std::vector<std::filesystem::path> times;
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        if (i > 0 && sorted[i].first == sorted[i - 1].first)
            throw std::runtime_error("duplicate physical time in results");
        times.push_back(sorted[i].second);
    }
    if (mode == "latest") return {times.back()};
    return times;
}

std::string xml(const std::string& text) {
    std::string result;
    for (char c : text) {
        if (c == '&') result += "&amp;";
        else if (c == '<') result += "&lt;";
        else if (c == '>') result += "&gt;";
        else if (c == '"') result += "&quot;";
        else result += c;
    }
    return result;
}

void writePvd(const std::filesystem::path& path, const std::vector<std::string>& times) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create PVD index: " + path.string());
    output << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" "
              "byte_order=\"LittleEndian\">\n<Collection>\n";
    for (const std::string& name : times) {
        double time;
        if (!numericalTime(name, time)) throw std::runtime_error("PVD needs physical times");
        output << "  <DataSet timestep=\"" << std::setprecision(17) << time
               << "\" group=\"\" part=\"0\" file=\"" << xml(name) << ".vtu\"/>\n";
    }
    output << "</Collection>\n</VTKFile>\n";
    output.close();
    if (!output) throw std::runtime_error("cannot flush PVD index");
}

std::string componentName(const ResultField& field, int component) {
    if (field.components == 1) return field.info.name;
    static constexpr const char* vector_suffix[] = {"x", "y", "z"};
    static constexpr const char* tensor_suffix[] = {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"};
    return field.info.name + "_" +
        (field.components == 3 ? vector_suffix[component] : tensor_suffix[component]);
}

void writeVtk(const std::filesystem::path& path, const Mesh& mesh, const ResultData& data) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create VTU output: " + path.string());
    output << "<?xml version=\"1.0\"?>\n"
              "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
              "<UnstructuredGrid><Piece NumberOfPoints=\"" << mesh.vertexCount()
           << "\" NumberOfCells=\"" << mesh.cellCount() << "\">\n"
              "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n"
           << std::setprecision(17);
    for (const Vec3& point : detail::meshData(mesh).vertices) output << point.x << ' ' << point.y << ' ' << point.z << '\n';
    output << "</DataArray></Points>\n<Cells>\n"
              "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (Index k = 0; k < detail::meshData(mesh).dimensions[2]; ++k) {
        for (Index j = 0; j < detail::meshData(mesh).dimensions[1]; ++j) {
            for (Index i = 0; i < detail::meshData(mesh).dimensions[0]; ++i) {
                output << mesh.vertexId(i, j, k) << ' ' << mesh.vertexId(i + 1, j, k) << ' '
                       << mesh.vertexId(i + 1, j + 1, k) << ' ' << mesh.vertexId(i, j + 1, k) << ' '
                       << mesh.vertexId(i, j, k + 1) << ' ' << mesh.vertexId(i + 1, j, k + 1) << ' '
                       << mesh.vertexId(i + 1, j + 1, k + 1) << ' ' << mesh.vertexId(i, j + 1, k + 1) << '\n';
            }
        }
    }
    output << "</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) output << 8LL * (cell + 1) << '\n';
    output << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) output << "12\n";
    output << "</DataArray></Cells>\n<CellData>\n";
    for (const ResultField& field : data.fields) {
        output << "<DataArray type=\"Float64\" Name=\"" << xml(field.info.name)
               << "\" NumberOfComponents=\"" << field.components << "\" format=\"ascii\">\n";
        for (double value : field.values) output << value << '\n';
        output << "</DataArray>\n";
    }
    output << "</CellData>\n</Piece></UnstructuredGrid>\n</VTKFile>\n";
    output.close();
    if (!output) throw std::runtime_error("cannot flush VTU output");
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
        for (const Vec3& point : detail::meshData(mesh).vertices) output << point[coordinate] << '\n';
    }
    for (const ResultField& field : data.fields) {
        for (int component = 0; component < field.components; ++component) {
            for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
                output << field.values[static_cast<std::size_t>(cell) * field.components + component] << '\n';
            }
        }
    }
    for (Index k = 0; k < detail::meshData(mesh).dimensions[2]; ++k) {
        for (Index j = 0; j < detail::meshData(mesh).dimensions[1]; ++j) {
            for (Index i = 0; i < detail::meshData(mesh).dimensions[0]; ++i) {
                output << mesh.vertexId(i, j, k) + 1 << ' ' << mesh.vertexId(i + 1, j, k) + 1 << ' '
                       << mesh.vertexId(i + 1, j + 1, k) + 1 << ' ' << mesh.vertexId(i, j + 1, k) + 1 << ' '
                       << mesh.vertexId(i, j, k + 1) + 1 << ' ' << mesh.vertexId(i + 1, j, k + 1) + 1 << ' '
                       << mesh.vertexId(i + 1, j + 1, k + 1) + 1 << ' ' << mesh.vertexId(i, j + 1, k + 1) + 1 << '\n';
            }
        }
    }
    output.close();
    if (!output) throw std::runtime_error("cannot flush Tecplot output");
}

int run(const Arguments& arguments) {
    for (const std::string& format : arguments.formats) {
        if (format != "vtk" && format != "tecplot")
            throw std::invalid_argument("unknown output format " + format);
    }
    const CaseDefinition definition = readCase(arguments.case_directory);
    const OutputControl output = readOutputControl(definition);
    const Mesh mesh = readMeshFile(definition.mesh_file);
    const auto time_directories = resolveTimeDirectories(definition, output, arguments.time_name);
    const auto selection = std::filesystem::path(arguments.time_name);
    const auto post_directory = definition.root / "post" / selection.parent_path();
    std::filesystem::create_directories(post_directory);
    const bool all = selection.filename() == "all";
    std::vector<std::string> times;
    bool writes_vtk = false;
    for (const auto& directory : time_directories) {
        const ResultData results = readParallelResults(directory, mesh.cellCount());
        if (results.global_dimensions != detail::meshData(mesh).dimensions)
            throw std::runtime_error("result global dimensions do not match the case mesh");
        const std::string name = directory.filename().string();
        double actual, expected;
        if (all && (!numericalTime(results.time_name, actual) || !numericalTime(name, expected) ||
                    actual != expected))
            throw std::runtime_error("result metadata time does not match its directory");
        for (const std::string& format : arguments.formats) {
            if (format == "vtk") {
                writeVtk(post_directory / (name + ".vtu"), mesh, results);
                writes_vtk = true;
            } else {
                writeTecplot(post_directory / (name + ".dat"), mesh, results);
            }
        }
        times.push_back(name);
    }
    if (all && writes_vtk) writePvd(post_directory / "series.pvd", times);
    std::cout << "BabelSim postprocessed " << times.size()
              << " time set(s) into " << post_directory << '\n';
    return 0;
}

}  // 匿名命名空间
}  // babelsim 命名空间

int babelsim::runPostprocess(int argc, char* argv[]) {
    try {
        return babelsim::run(babelsim::parseArguments(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "babelsim-post: " << error.what() << '\n';
        return 1;
    }
}
