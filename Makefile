CXX := mpic++
AR ?= ar

CXXFLAGS ?= -std=c++17 -O3 -march=native -Wall -Wextra -Wpedantic -Wshadow \
            -DOMPI_SKIP_MPICXX=1 -DMPICH_SKIP_MPICXX=1
CPPFLAGS ?= -Iinclude -I/usr/include/eigen3

BUILD := build
LIB := $(BUILD)/libbabelsim.a
TAIHO_POISEUILLE ?= /home/midway/TaihoCFD/examples/meshes/poiseuille

SOURCES := src/mesh.cpp \
           src/legacy_taiho.cpp \
           src/operators.cpp \
           src/incompressible.cpp \
           src/assembly.cpp \
           src/linear_solver.cpp \
           src/distributed_solver.cpp \
           src/parallel.cpp
OBJECTS := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SOURCES))
HEADERS := $(wildcard include/babelsim/*.h)

TEST_SOURCES := tests/mesh_geometry_test.cpp \
                tests/field_boundary_test.cpp \
                tests/legacy_import_test.cpp \
                tests/operators_test.cpp \
                tests/assembly_solver_test.cpp \
                tests/simple_solver_test.cpp \
                tests/cavity_regression_test.cpp \
                tests/cavity_3d_test.cpp \
                tests/nonorthogonal_cavity_test.cpp
TESTS := $(patsubst tests/%.cpp,$(BUILD)/%,$(TEST_SOURCES))
MPI_TESTS := $(BUILD)/parallel_domain_test $(BUILD)/parallel_simple_test \
             $(BUILD)/parallel_imported_test $(BUILD)/parallel_cavity_3d_test
EXAMPLES := $(BUILD)/cavity $(BUILD)/imported_simple

.DEFAULT_GOAL := all

all: $(LIB) $(EXAMPLES)

$(LIB): $(OBJECTS)
	@mkdir -p $(dir $@)
	$(AR) rcs $@ $^

$(BUILD)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD)/%: tests/%.cpp tests/test_util.h $(HEADERS) $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(BUILD)/cavity: examples/cavity.cpp $(HEADERS) $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(BUILD)/imported_simple: examples/imported_simple.cpp $(HEADERS) $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(BUILD)/parallel_cavity_3d_test: tests/parallel_cavity_3d_test.cpp tests/test_util.h $(HEADERS) $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

test: $(TESTS)
	@set -e; for test in $(TESTS); do $$test; done

test-mpi: $(MPI_TESTS)
	TMPDIR=/tmp mpirun -np 2 $(BUILD)/parallel_domain_test $(BUILD)/mpi-output
	TMPDIR=/tmp mpirun -np 1 $(BUILD)/parallel_simple_test $(BUILD)/mpi-output
	TMPDIR=/tmp mpirun -np 2 $(BUILD)/parallel_simple_test $(BUILD)/mpi-output
	TMPDIR=/tmp mpirun -np 4 $(BUILD)/parallel_simple_test $(BUILD)/mpi-output
	TMPDIR=/tmp mpirun -np 2 $(BUILD)/parallel_cavity_3d_test $(BUILD)/mpi-output

test-mpi-poiseuille: $(BUILD)/parallel_imported_test
	@rm -rf $(BUILD)/mpi-poiseuille-np1 $(BUILD)/mpi-poiseuille-np2
	TMPDIR=/tmp mpirun -np 1 $(BUILD)/parallel_imported_test \
		$(TAIHO_POISEUILLE) $(BUILD)/mpi-poiseuille-np1
	TMPDIR=/tmp mpirun -np 2 $(BUILD)/parallel_imported_test \
		$(TAIHO_POISEUILLE) $(BUILD)/mpi-poiseuille-np2
	python3 tools/compare_parallel_csv.py \
		$(BUILD)/mpi-poiseuille-np1 $(BUILD)/mpi-poiseuille-np2 \
		--reference-ranks 1 --candidate-ranks 2 --prefix poiseuille \
		--atol 5e-6 --rtol 5e-6

validate-cavity: $(BUILD)/cavity
	@mkdir -p $(BUILD)/validation
	$(BUILD)/cavity 64 5000 0.01 $(BUILD)/validation/cavity64.csv
	python3 tools/validate_cavity.py $(BUILD)/validation/cavity64.csv

validate-poiseuille: $(BUILD)/imported_simple
	@mkdir -p $(BUILD)/validation
	$(BUILD)/imported_simple $(TAIHO_POISEUILLE) 1200 1.0 0.01 \
		$(BUILD)/validation/poiseuille.csv
	python3 tools/validate_poiseuille.py $(BUILD)/validation/poiseuille.csv \
		--y-min 0.0032154057 --y-max 0.9967845943

validate: test validate-cavity validate-poiseuille

clean:
	$(RM) -r $(BUILD)

.PHONY: all test test-mpi test-mpi-poiseuille validate validate-cavity validate-poiseuille clean

-include $(OBJECTS:.o=.d)
