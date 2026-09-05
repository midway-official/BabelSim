# 默认只构建框架库和正式程序；测试仅由显式 test*/validate* 目标触发。
.DEFAULT_GOAL := all

CXX := mpic++
AR := gcc-ar

# 面向本机计算速度：跨文件优化、矢量化、浮点重结合/倒数优化与融合乘加。
# 保留 NaN/Inf 检查，不能让非法输入或发散结果被视为正常收敛。
# fat LTO 同时保存机器码，允许外部 Solver 不启用 LTO 时链接静态库。
OPTFLAGS ?= -O3 -march=native -mtune=native -flto=auto -ffat-lto-objects \
            -ffast-math -fno-finite-math-only -ffp-contract=fast -DNDEBUG
CXXFLAGS ?= -std=c++17 $(OPTFLAGS) -Wall -Wextra -Wpedantic -Wshadow \
            -DOMPI_SKIP_MPICXX=1 -DMPICH_SKIP_MPICXX=1
CPPFLAGS ?= -Iinclude -Isrc -I/usr/include/eigen3

BUILD := build
LIB := $(BUILD)/libbabelsim.a
# 计算后端只需实现 internal/compute_backend.h 的 makeComputeBackend()。
# 整组替换可同时移除默认 Eigen 装配/求解实现，数值前端和 Physics 无需修改。
COMPUTE_BACKEND_SOURCES ?= src/backend/eigen_mpi.cpp \
                          src/backend/eigen_assembly.cpp \
                          src/backend/algebraic_multigrid.cpp \
                          src/algebra/linear_solver.cpp \
                          src/algebra/distributed_solver.cpp
SOURCES := src/core/mesh.cpp \
           src/io/config.cpp \
           src/io/case_reader.cpp \
           src/io/case.cpp \
           src/io/field_reader.cpp \
           src/io/numerics_reader.cpp \
           src/io/mesh_reader.cpp \
           src/io/result_reader.cpp \
           src/io/postprocess.cpp \
           src/discretization/operators.cpp \
           src/discretization/equation_expression.cpp \
           src/discretization/fvm_execution.cpp \
           $(COMPUTE_BACKEND_SOURCES) \
           src/parallel/parallel_context.cpp \
           src/parallel/parallel_writer.cpp \
           src/runtime/runtime.cpp \
           src/runtime/solver_api.cpp \
           src/runtime/application.cpp \
           $(filter-out $(wildcard src/physics/*/main.cpp),$(wildcard src/physics/*/*.cpp))
OBJECTS := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SOURCES))
SOLVER_SOURCES := $(wildcard src/physics/*/main.cpp)
SOLVER_OBJECTS := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SOLVER_SOURCES))
HEADERS := $(wildcard include/babelsim/*.h)

TEST_SOURCES := tests/mesh_geometry_test.cpp \
                tests/field_boundary_test.cpp \
                tests/operators_test.cpp \
                tests/math_runtime_test.cpp \
                tests/public_equation_test.cpp \
                tests/backend_interface_test.cpp \
                tests/assembly_solver_test.cpp \
                tests/simple_solver_test.cpp \
                tests/mesh_file_test.cpp \
                tests/case_lifecycle_test.cpp \
                tests/case_io_test.cpp \
                tests/field_writer_test.cpp \
                tests/heat_solver_test.cpp \
                tests/time_history_test.cpp \
                tests/transport_solver_test.cpp \
                tests/cavity_regression_test.cpp \
                tests/cavity_3d_test.cpp \
                tests/nonorthogonal_cavity_test.cpp \
                tests/nonorthogonal_cavity_3d_test.cpp
TESTS := $(patsubst tests/%.cpp,$(BUILD)/%,$(TEST_SOURCES))
MPI_TESTS := $(BUILD)/parallel_domain_test $(BUILD)/parallel_simple_test \
             $(BUILD)/parallel_math_test \
             $(BUILD)/parallel_channel_test $(BUILD)/parallel_cavity_3d_test \
             $(BUILD)/parallel_transport_test
APPS := $(BUILD)/babelsim-solve $(BUILD)/babelsim-post

all: $(LIB) $(APPS)

# 目录依赖使删除/新增 Solver 文件后也会重新归档，避免残留旧模块。
PHYSICS_DIRECTORIES := src/physics $(wildcard src/physics/*/)
$(LIB): $(OBJECTS) $(PHYSICS_DIRECTORIES) Makefile
	@mkdir -p $(dir $@)
	$(RM) $@
	$(AR) rcs $@ $(OBJECTS)

$(BUILD)/%.o: src/%.cpp Makefile
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD)/%: tests/%.cpp tests/test_util.h $(HEADERS) $(LIB)
	@mkdir -p $(dir $@)
	+$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $< $(LIB) $(LDLIBS) -o $@

# 注册源文件作为目标文件直接链接，不能仅藏在静态库中等待按需抽取。
# 链接命令前的 + 让 LTO 使用 Make 的并行作业配额，避免另起不限额编译进程。
$(BUILD)/babelsim-solve: src/apps/babelsim_solve.cpp $(HEADERS) $(SOLVER_OBJECTS) $(LIB) $(PHYSICS_DIRECTORIES)
	@mkdir -p $(dir $@)
	+$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $< $(SOLVER_OBJECTS) $(LIB) $(LDLIBS) -o $@

$(BUILD)/babelsim-post: src/apps/babelsim_post.cpp $(HEADERS) $(LIB)
	@mkdir -p $(dir $@)
	+$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/parallel_channel_test: tests/parallel_channel_test.cpp tests/test_util.h $(HEADERS) $(LIB)
	@mkdir -p $(dir $@)
	+$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/parallel_cavity_3d_test: tests/parallel_cavity_3d_test.cpp tests/test_util.h $(HEADERS) $(LIB)
	@mkdir -p $(dir $@)
	+$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $< $(LIB) $(LDLIBS) -o $@

test-architecture:
	python3 tests/architecture_test.py

test-external: $(LIB)
	python3 tests/external_solver_test.py

test: test-architecture $(TESTS)
	@set -e; for test in $(TESTS); do $$test; done

$(BUILD)/case_programming_test: tests/case_programming_test.cpp tests/examples/coupled_scalar.cpp $(HEADERS) $(LIB)
	+$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) tests/case_programming_test.cpp tests/examples/coupled_scalar.cpp $(LIB) $(LDLIBS) -o $@

test-workflow: test-architecture $(APPS) $(BUILD)/case_programming_test $(BUILD)/time_history_test
	$(BUILD)/time_history_test
	python3 tests/solver_workflow_test.py

test-mpi: $(MPI_TESTS)
	TMPDIR=/tmp mpirun -np 1 $(BUILD)/parallel_math_test
	TMPDIR=/tmp mpirun -np 2 $(BUILD)/parallel_math_test
	TMPDIR=/tmp mpirun -np 4 $(BUILD)/parallel_math_test
	TMPDIR=/tmp mpirun -np 2 $(BUILD)/parallel_domain_test $(BUILD)/mpi-output
	TMPDIR=/tmp mpirun -np 2 $(BUILD)/parallel_channel_test \
		cases/poiseuille/mesh/poiseuille.mesh $(BUILD)/mpi-output
	TMPDIR=/tmp mpirun -np 1 $(BUILD)/parallel_simple_test $(BUILD)/mpi-output
	TMPDIR=/tmp mpirun -np 2 $(BUILD)/parallel_simple_test $(BUILD)/mpi-output
	TMPDIR=/tmp mpirun -np 4 $(BUILD)/parallel_simple_test $(BUILD)/mpi-output
	TMPDIR=/tmp mpirun -np 2 $(BUILD)/parallel_cavity_3d_test $(BUILD)/mpi-output
	TMPDIR=/tmp mpirun -np 2 $(BUILD)/parallel_transport_test
	$(MAKE) test-mpi-heat

test-mpi-heat: $(BUILD)/babelsim-solve
	TMPDIR=/tmp mpirun -np 1 $(BUILD)/babelsim-solve \
		-case cases/heat -time mpi-np1
	TMPDIR=/tmp mpirun -np 2 $(BUILD)/babelsim-solve \
		-case cases/heat -time mpi-np2
	python3 tools/compare_parallel_results.py \
		cases/heat/results/mpi-np1 cases/heat/results/mpi-np2 \
		--atol 1e-9 --rtol 1e-9

test-mpi-poiseuille: $(BUILD)/babelsim-solve $(BUILD)/babelsim-post
	TMPDIR=/tmp mpirun -np 1 $(BUILD)/babelsim-solve \
		-case cases/poiseuille -time mpi-np1
	TMPDIR=/tmp mpirun -np 2 $(BUILD)/babelsim-solve \
		-case cases/poiseuille -time mpi-np2
	TMPDIR=/tmp mpirun -np 4 $(BUILD)/babelsim-solve \
		-case cases/poiseuille -time mpi-np4
	python3 tools/compare_parallel_results.py \
		cases/poiseuille/results/mpi-np1 cases/poiseuille/results/mpi-np2 \
		--atol 5e-6 --rtol 5e-6
	python3 tools/compare_parallel_results.py \
		cases/poiseuille/results/mpi-np1 cases/poiseuille/results/mpi-np4 \
		--atol 5e-6 --rtol 5e-6
	$(BUILD)/babelsim-post -case cases/poiseuille -time mpi-np4 -format vtk tecplot

validate-cavity: $(BUILD)/babelsim-solve
	TMPDIR=/tmp mpirun -np 4 $(BUILD)/babelsim-solve \
		-case cases/cavity -time validation-mpi4
	python3 cases/cavity/validation/validate_cavity.py \
		cases/cavity/results/validation-mpi4

postprocess-mpi-poiseuille:
	@test -d cases/poiseuille/results/mpi-np4 || \
		(echo "run make test-mpi-poiseuille first" && exit 1)
	$(BUILD)/babelsim-post -case cases/poiseuille -time mpi-np4 -format vtk tecplot

validate-poiseuille: $(BUILD)/babelsim-solve
	$(BUILD)/babelsim-solve -case cases/poiseuille -time validation
	python3 tools/validate_poiseuille.py cases/poiseuille/results/validation/rank-0000/U.csv \
		--y-min 0.0 --y-max 1.0

validate: test validate-cavity validate-poiseuille

clean:
	$(RM) -r $(BUILD)

.PHONY: all test test-architecture test-external test-workflow test-mpi test-mpi-heat test-mpi-poiseuille postprocess-mpi-poiseuille \
	validate validate-cavity validate-poiseuille clean

-include $(OBJECTS:.o=.d) $(SOLVER_OBJECTS:.o=.d)
