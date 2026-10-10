#include "babelsim/mesh_io.h"
#include "internal/mesh_partition.h"
#include <metis.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

// A standalone audit, linked with --wrap=METIS_PartGraphKway. This captures
// the production partitioner's actual CSR and floating-point parameters.
// No solver or framework implementation is changed by this observation.
namespace {
std::string prefix;
template<class T> void dump(const char* suffix, const T* values, std::size_t size) {
    std::ofstream file(prefix + suffix, std::ios::binary);
    file.write(reinterpret_cast<const char*>(values), size * sizeof(T));
    if (!file) throw std::runtime_error("cannot write partition audit data");
}
}

extern "C" int __real_METIS_PartGraphKway(idx_t*, idx_t*, idx_t*, idx_t*, idx_t*,
    idx_t*, idx_t*, idx_t*, real_t*, real_t*, idx_t*, idx_t*, idx_t*);
extern "C" int __wrap_METIS_PartGraphKway(idx_t* cells, idx_t* constraints,
    idx_t* offsets, idx_t* neighbours, idx_t* weights, idx_t* sizes,
    idx_t* edge_weights, idx_t* parts, real_t* targets, real_t* imbalance,
    idx_t* options, idx_t* objective, idx_t* owners) {
    dump("-offsets.bin", offsets, *cells + 1);
    dump("-neighbours.bin", neighbours, offsets[*cells]);
    dump("-weights.bin", weights, *cells * *constraints);
    dump("-targets.bin", targets, *parts * *constraints);
    dump("-imbalance.bin", imbalance, *constraints);
    dump("-options.bin", options, METIS_NOPTIONS);
    const int status = __real_METIS_PartGraphKway(cells, constraints, offsets,
        neighbours, weights, sizes, edge_weights, parts, targets, imbalance,
        options, objective, owners);
    dump("-raw-owners.bin", owners, *cells);
    return status;
}

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    prefix = argv[3];
    const auto mesh = babelsim::readMeshFile(argv[1]);
    mesh.validate();
    const int parts = std::stoi(argv[2]);
    const auto owners = babelsim::detail::partitionOwners(mesh, parts);
    dump("-owners.bin", owners.data(), owners.size());
    std::cout << "{\"cells\":" << mesh.cellCount() << ",\"parts\":" << parts
        << ",\"indexBytes\":" << sizeof(idx_t) << ",\"realBytes\":" << sizeof(real_t)
        << ",\"prefix\":\"" << prefix << "\"}\n";
}
