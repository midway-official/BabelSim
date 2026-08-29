#include "babelsim/parallel_writer.h"
#include "babelsim/result_reader.h"

#include "test_util.h"

#include <filesystem>
#include <iostream>

using namespace babelsim;

namespace {

const ResultField& field(const ResultData& data, const std::string& name) {
    for (const ResultField& result : data.fields) {
        if (result.info.name == name) return result;
    }
    throw std::runtime_error("missing output field " + name);
}

}  // 匿名命名空间

int main() {
    const Mesh mesh = Mesh::cartesian({2, 2, 1}, {0, 0, 0}, {2, 2, 1});
    ScalarField pressure(mesh, FieldLocation::Cell, "p");
    VectorField velocity(mesh, FieldLocation::Cell, "U");
    TensorField tensor(mesh, FieldLocation::Cell, "gradU");
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        pressure[cell] = 10.0 + cell;
        velocity[cell] = {double(cell), double(2 * cell), -double(cell)};
        tensor[cell].rows[0] = {double(cell), 1.0, 2.0};
        tensor[cell].rows[1] = {3.0, 4.0, 5.0};
        tensor[cell].rows[2] = {6.0, 7.0, 8.0};
    }
    const std::filesystem::path output = "build/field-writer-test/final";
    std::filesystem::remove_all(output.parent_path());
    const ParallelContext serial;
    writeOwnedFieldCsv(output, pressure, serial);
    writeOwnedFieldCsv(output, velocity, serial);
    writeOwnedFieldCsv(output, tensor, serial);
    writeOwnedResultMetadata(output, mesh, serial, "final", {
        {"p", "scalar", FieldLocation::Cell},
        {"U", "vector", FieldLocation::Cell},
        {"gradU", "tensor", FieldLocation::Cell},
    });
    const ResultData result = readParallelResults(output, mesh.cellCount());
    require(result.time_name == "final" && result.fields.size() == 3,
            "parallel output metadata was not reconstructed");
    const ResultField& p = field(result, "p");
    const ResultField& u = field(result, "U");
    const ResultField& grad = field(result, "gradU");
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        require(near(p.values[static_cast<std::size_t>(cell)], pressure[cell]),
                "scalar output value changed");
        const std::size_t vector = static_cast<std::size_t>(cell) * 3;
        require(near(u.values[vector], velocity[cell].x) &&
                near(u.values[vector + 1], velocity[cell].y) &&
                near(u.values[vector + 2], velocity[cell].z),
                "vector output value changed");
        const std::size_t tensor_index = static_cast<std::size_t>(cell) * 9;
        require(near(grad.values[tensor_index], tensor[cell].rows[0].x) &&
                near(grad.values[tensor_index + 8], tensor[cell].rows[2].z),
                "tensor output value changed");
    }
    std::cout << "field_writer_test: scalar/vector/tensor owned output passed\n";
}
