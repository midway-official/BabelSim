CXX := mpic++
AR ?= ar

CXXFLAGS ?= -std=c++17 -O3 -march=native -Wall -Wextra -Wpedantic -Wshadow \
            -DOMPI_SKIP_MPICXX=1 -DMPICH_SKIP_MPICXX=1
CPPFLAGS ?= -Iinclude -Isrc -I/usr/include/eigen3

BUILD := build
LIB := $(BUILD)/libbabelsim.a
SOURCES := src/core/mesh.cpp \
           src/apps/solver_selection.cpp \
           src/io/config.cpp \
           src/io/case_reader.cpp \
           src/io/case.cpp \
           src/io/field_reader.cpp \
           src/io/numerics_reader.cpp \
           src/io/mesh_reader.cpp \
           src/io/result_reader.cpp \
           src/discretization/operators.cpp \
           src/discretization/simple_discretization.cpp \
           src/discretization/fvm_expression.cpp \
           src/discretization/assembly.cpp \
           src/algebra/linear_solver.cpp \
           src/algebra/distributed_solver.cpp \
           src/parallel/parallel_context.cpp \
           src/parallel/parallel_writer.cpp \
           src/runtime/runtime.cpp \
           src/physics/heat/transient_heat_solver.cpp \
           src/physics/transport/transient_scalar_transport_solver.cpp \
           src/physics/simple/simple_solver.cpp \
           src/physics/simple/create_fields.cpp \
           src/physics/simple/convergence.cpp \
           src/physics/simple/momentum.cpp \
           src/physics/simple/pressure.cpp \
           $(wildcard src/physics/*/main.cpp)
OBJECTS := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SOURCES))
HEADERS := $(wildcard include/babelsim/*.h)

TEST_SOURCES := tests/mesh_geometry_test.cpp \
                tests/field_boundary_test.cpp \
                tests/operators_test.cpp \
                tests/fvc_runtime_test.cpp \
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
             $(BUILD)/parallel_channel_test $(BUILD)/parallel_cavity_3d_test \
             $(BUILD)/parallel_transport_test
APPS := $(BUILD)/babelsim-solve $(BUILD)/babelsim-post

.DEFAULT_GOAL := all

all: $(LIB) $(APPS)

$(LIB): $(OBJECTS)
	@mkdir -p $(dir $@)
	$(RM) $@
	$(AR) rcs $@ $^

$(BUILD)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD)/%: tests/%.cpp tests/test_util.h $(HEADERS) $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(BUILD)/babelsim-solve: src/apps/babelsim_solve.cpp $(HEADERS) $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(BUILD)/babelsim-post: src/apps/babelsim_post.cpp $(HEADERS) $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(BUILD)/parallel_channel_test: tests/parallel_channel_test.cpp tests/test_util.h $(HEADERS) $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(BUILD)/parallel_cavity_3d_test: tests/parallel_cavity_3d_test.cpp tests/test_util.h $(HEADERS) $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

test: $(TESTS)
	@set -e; for test in $(TESTS); do $$test; done

$(BUILD)/case_programming_test: tests/case_programming_test.cpp tests/examples/coupled_scalar.cpp $(HEADERS) $(LIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) tests/case_programming_test.cpp tests/examples/coupled_scalar.cpp $(LIB) -o $@

test-workflow: $(APPS) $(BUILD)/case_programming_test $(BUILD)/time_history_test
	$(BUILD)/time_history_test
	python3 tests/solver_workflow_test.py

test-mpi: $(MPI_TESTS)
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
	$(BUILD)/babelsim-solve -case cases/cavity -time validation
	python3 tools/validate_cavity.py cases/cavity/results/validation/rank-0000/U.csv

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

.PHONY: all test test-workflow test-mpi test-mpi-heat test-mpi-poiseuille postprocess-mpi-poiseuille \
	validate validate-cavity validate-poiseuille clean

-include $(OBJECTS:.o=.d)
