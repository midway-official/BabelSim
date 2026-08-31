#include "babelsim/runtime.h"
#include "babelsim/mpi_support.h"
#include "babelsim/parallel.h"
#include "babelsim/operators.h"
#include "internal/field_access.h"
#include "internal/mesh_access.h"
#include "test_util.h"

#include <iostream>

using namespace babelsim;

namespace {
double magnitude(double value) { return std::abs(value); }
double magnitude(Vec3 value) { return norm(value); }
double difference(double a, double b) { return magnitude(a-b); }
double difference(Vec3 a, Vec3 b) { return magnitude(a-b); }
double difference(const Tensor3& a, const Tensor3& b) {
    return difference(a[0], b[0]) + difference(a[1], b[1]) + difference(a[2], b[2]);
}

// 测试维护接口故意破坏 ghost/重复面值，确保每个公开 fvc 操作独立履行同步契约。
template <typename T>
void poison(Field<T>& field) {
    const Mesh& mesh = field.mesh();
    if (field.location() == FieldLocation::Cell) {
        for (Index cell = 0; cell < mesh.cellCount(); ++cell)
            if (!detail::isOwned(mesh, cell)) detail::fieldData(field)[cell] = T{};
    } else if (detail::meshData(mesh).global_i_offset > 0) {
        for (Index cell : detail::meshData(mesh).owned_cells) {
            if (cell % detail::meshData(mesh).dimensions[0] == detail::meshData(mesh).owned_i_begin)
                detail::fieldData(field)[detail::meshData(mesh).cell_faces[cell][0]] = T{};
        }
    }
}

template <typename T>
double compare(const Field<T>& local, const Field<T>& serial) {
    const Mesh& mesh = local.mesh();
    double error = 0.0;
    if (local.location() == FieldLocation::Cell) {
        for (Index cell = 0; cell < mesh.cellCount(); ++cell)
            error = std::max(error, difference(detail::fieldData(local)[cell],
                detail::fieldData(serial)[detail::globalCellId(mesh, cell)]));
    } else {
        for (Index face : detail::meshData(mesh).owned_faces) {
            bool found = false;
            for (Index other = 0; other < serial.mesh().faceCount(); ++other) {
                if (norm(mesh.faceCentre(face) - serial.mesh().faceCentre(other)) > 1e-12) continue;
                error = std::max(error, difference(detail::fieldData(local)[face], detail::fieldData(serial)[other]));
                found = true;
                break;
            }
            require(found, "serial face mapping missing");
        }
    }
    return error;
}

struct Fields {
    explicit Fields(const Mesh& mesh)
        : p(mesh, FieldLocation::Cell, "p"), k(mesh, FieldLocation::Cell, "k"),
          U(mesh, FieldLocation::Cell, "U"), gradP(mesh, FieldLocation::Cell),
          gradU(mesh, FieldLocation::Cell), phi(mesh, FieldLocation::Face),
          scalar(mesh, FieldLocation::Cell), vector(mesh, FieldLocation::Cell),
          faceScalar(mesh, FieldLocation::Face), faceVector(mesh, FieldLocation::Face) {
        p.evaluate([](Vec3 x) { return 1 + 2*x.x - 0.7*x.y + 0.3*x.z; });
        k.evaluate([](Vec3 x) { return 2 + 0.2*x.x; });
        U.evaluate([](Vec3 x) { return Vec3{x.x + 0.3*x.y, -x.y + 0.1*x.z, 0.5*x.z}; });
        for (Index patch = 0; patch < mesh.patchCount(); ++patch) {
            if (mesh.patchKind(patch) == PatchKind::Processor) continue;
            if (detail::meshData(mesh).patches[patch].faces.empty()) continue;
            const Vec3 n = mesh.faceNormal(detail::meshData(mesh).patches[patch].faces.front());
            p.boundary(patch) = fixedGradient(dot(Vec3{2, -0.7, 0.3}, n));
            k.boundary(patch) = fixedGradient(0.2*n.x);
            U.boundary(patch) = fixedGradient(Vec3{n.x+0.3*n.y, -n.y+0.1*n.z, 0.5*n.z});
        }
    }
    ScalarField p, k;
    VectorField U, gradP;
    TensorField gradU;
    ScalarField phi, scalar;
    VectorField vector;
    ScalarField faceScalar;
    VectorField faceVector;
};

// 两个运行阶段按同一操作序列记录/比较，串行参考在 MPI_Init 之前计算。
struct Answers {
    std::vector<ScalarField> scalars;
    std::vector<VectorField> vectors;
    std::vector<TensorField> tensors;
    std::size_t si = 0, vi = 0, ti = 0;
    double error = 0;
    void check(const ScalarField& value, bool record) {
        if (record) scalars.push_back(value); else error = std::max(error, compare(value, scalars.at(si++)));
    }
    void check(const VectorField& value, bool record) {
        if (record) vectors.push_back(value); else error = std::max(error, compare(value, vectors.at(vi++)));
    }
    void check(const TensorField& value, bool record) {
        if (record) tensors.push_back(value); else error = std::max(error, compare(value, tensors.at(ti++)));
    }
};

void exercise(Fields& f, Answers& answers, bool record) {
    const auto check = [&](const auto& field) { answers.check(field, record); };
    poison(f.p); fvc::evaluate(fvc::grad(f.p), f.gradP); check(f.gradP);
    poison(f.U); fvc::evaluate(fvc::grad(f.U), f.gradU); check(f.gradU);
    poison(f.U); fvc::evaluate(fvc::flux(f.U), f.phi); check(f.phi);
    poison(f.phi); fvc::evaluate(fvc::div(f.phi), f.scalar); check(f.scalar);
    poison(f.U); fvc::evaluate(fvc::div(f.U), f.scalar); check(f.scalar);
    poison(f.p); fvc::evaluate(fvc::interpolate(f.p), f.faceScalar); check(f.faceScalar);
    poison(f.U); fvc::evaluate(fvc::interpolate(f.U), f.faceVector); check(f.faceVector);
    poison(f.p); poison(f.gradP);
    fvc::evaluate(fvc::reconstruct(f.p, f.gradP), f.faceScalar); check(f.faceScalar);
    poison(f.U); poison(f.gradU);
    fvc::evaluate(fvc::reconstruct(f.U, f.gradU), f.faceVector); check(f.faceVector);
    poison(f.phi); poison(f.p); fvc::evaluate(fvc::div(f.phi, f.p), f.scalar); check(f.scalar);
    poison(f.phi); poison(f.U); fvc::evaluate(fvc::div(f.phi, f.U), f.vector); check(f.vector);
    poison(f.p); fvc::evaluate(fvc::laplacian(f.p), f.scalar); check(f.scalar);
    poison(f.p); poison(f.k); fvc::evaluate(fvc::laplacian(f.k, f.p), f.scalar); check(f.scalar);
    poison(f.p); poison(f.k); fvc::evaluate(fvc::flux(f.k, f.p), f.faceScalar); check(f.faceScalar);
    poison(f.p); fvc::evaluate(fvc::normalGradient(f.p), f.faceScalar); check(f.faceScalar);
    f.vector.fill({1, 2, 3}); poison(f.k); poison(f.p);
    fvc::subtract(f.k, fvc::grad(f.p), f.vector); check(f.vector);
    f.faceScalar.fill(3); poison(f.k); poison(f.p);
    fvc::subtract(fvc::flux(f.k, f.p), f.faceScalar); check(f.faceScalar);

    poison(f.faceVector);
    fvc::evaluate(fvc::flux(f.faceVector), f.phi); check(f.phi);
    f.phi.fill(3); poison(f.phi); poison(f.faceVector);
    fvc::add(fvc::flux(f.faceVector), f.phi); check(f.phi);
    f.phi.fill(3); poison(f.phi); poison(f.faceVector);
    fvc::add(fvc::flux(f.faceVector), f.phi, fvc::FaceRegion::Interior); check(f.phi);
    // 给定梯度有意不等于 grad(p)，确保执行层确实使用传入的重构量。
    f.gradP.fill({0.7, -0.2, 0.4});
    f.faceScalar.evaluate([](Vec3 x) { return 1.7 + 0.1*x.x; });
    poison(f.faceScalar); poison(f.p); poison(f.gradP); poison(f.phi);
    fvc::subtract(fvc::flux(f.faceScalar, fvc::reconstruct(f.p, f.gradP)),
                  f.phi, fvc::FaceRegion::Interior); check(f.phi);
    // 与原 Rhie--Chow 逐面数学式直接对照，尤其检查分区界面与物理边界的区别。
    const Mesh& mesh = f.p.mesh();
    for (Index face : detail::meshData(mesh).owned_faces) {
        double expected = 3.0;
        if (detail::meshData(mesh).face_neighbour[face] != invalid_index)
            expected += dot(detail::fieldData(f.faceVector)[face], mesh.faceAreaVector(face)) -
                detail::fieldData(f.faceScalar)[face] * integratedNormalGradient(
                    f.p, f.gradP, face, numericalMethods().diffusion);
        require(near(detail::fieldData(f.phi)[face], expected, 1e-12), "composed flux differs from the original face formula");
    }
    poison(f.faceScalar); poison(f.p); poison(f.gradP);
    fvc::evaluate(fvc::flux(f.faceScalar, fvc::reconstruct(f.p, f.gradP)), f.phi); check(f.phi);
}
}  // 匿名命名空间

int main(int argc, char* argv[]) {
    // 仿射倾斜三维网格同时覆盖非正交面、第二层 halo 和分区交界。
    std::vector<Vec3> points;
    for (int k = 0; k <= 3; ++k)
        for (int j = 0; j <= 4; ++j)
            for (int i = 0; i <= 16; ++i) {
                const double x = i/16.0, y = j/4.0, z = k/3.0;
                points.push_back({x+0.2*y, y+0.15*z, z+0.1*x});
            }
    const Mesh global = Mesh::structured({16, 4, 3}, std::move(points));
    Answers answers;
    for (DiffusionMethod method : {DiffusionMethod::Orthogonal, DiffusionMethod::Corrected,
                                   DiffusionMethod::LimitedCorrected}) {
        RuntimeControl control;
        control.methods.diffusion = method;
        RunTime time = RunTime::forMesh(global, control);
        Fields serial(global);
        exercise(serial, answers, true);
    }
    detail::checkMpi(MPI_Init(&argc, &argv), "MPI_Init");
    try {
        const ParallelContext parallel = ParallelContext::world();
        const Mesh local = decompose(global, parallel);
        for (DiffusionMethod method : {DiffusionMethod::Orthogonal, DiffusionMethod::Corrected,
                                       DiffusionMethod::LimitedCorrected}) {
            RuntimeControl control;
            control.methods.diffusion = method;
            RunTime time = RunTime::forMesh(local, control);
            Fields fields(local);
            exercise(fields, answers, false);
        }
        double maximum = 0;
        parallel.maximum(&answers.error, &maximum, 1);
        require(maximum < 1e-9, "public fvc synchronization differs from serial oracle");
        if (parallel.rank == 0) std::cout << "parallel_fvc_test: 22 operations x 3 diffusion methods, poisoned halos, maxError=" << maximum << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        detail::checkMpi(MPI_Abort(MPI_COMM_WORLD, 1), "MPI_Abort");
    }
    detail::checkMpi(MPI_Finalize(), "MPI_Finalize");
}
