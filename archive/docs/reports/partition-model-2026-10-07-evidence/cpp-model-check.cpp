#include "babelsim/mesh_io.h"
#include "internal/mesh_partition.h"
#include <iostream>
#include <string>

// Run outside Physics. Independently read native face rings, then expose the
// framework ownership for comparison with the selected Python graph model.
int main(int argc, char** argv) {
    if (argc < 3) return 2;
    const int parts = std::stoi(argv[1]);
    for (int file = 2; file < argc; ++file) {
        const auto mesh = babelsim::readMeshFile(argv[file]);
        mesh.validate();
        const auto owners = babelsim::detail::partitionOwners(mesh, parts);
        std::cout << "{\"mesh\":\"" << argv[file] << "\",\"parts\":" << parts << ",\"owners\":[";
        for (std::size_t cell = 0; cell < owners.size(); ++cell) {
            if (cell) std::cout << ',';
            std::cout << owners[cell];
        }
        std::cout << "]}\n";
    }
}
