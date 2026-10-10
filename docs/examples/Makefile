# Run from the BabelSim root: make -f docs/examples/Makefile
include Makefile
.DEFAULT_GOAL := examples

EXAMPLE_DIR := $(BUILD)/examples
REGISTERED_SOLVERS := $(filter-out $(BUILD)/physics/custom_heat/main.o,$(SOLVER_OBJECTS))

.PHONY: examples sdk-example
examples: $(EXAMPLE_DIR)/custom-heat $(EXAMPLE_DIR)/coupled-scalar $(EXAMPLE_DIR)/registered-heat

# Compile example sources against public headers only, using standard g++.
$(EXAMPLE_DIR)/%.o: docs/examples/%.cpp $(HEADERS)
	@mkdir -p $(dir $@)
	g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -Iinclude -MMD -MP -c $< -o $@

$(EXAMPLE_DIR)/custom-heat: $(EXAMPLE_DIR)/external_heat.o $(EXAMPLE_DIR)/external_main.o $(LIB) docs/examples/Makefile
	+$(CXX) $(LDFLAGS) $(filter %.o,$^) $(LIB) $(LDLIBS) -o $@

$(EXAMPLE_DIR)/coupled-scalar: $(EXAMPLE_DIR)/coupled_scalar.o $(EXAMPLE_DIR)/external_main.o $(LIB) docs/examples/Makefile
	+$(CXX) $(LDFLAGS) $(filter %.o,$^) $(LIB) $(LDLIBS) -o $@

# Demonstrate registration beside the shipped solvers without editing src/physics.
$(EXAMPLE_DIR)/registered-heat: src/apps/babelsim_solve.cpp $(EXAMPLE_DIR)/external_heat.o $(REGISTERED_SOLVERS) $(LIB) docs/examples/Makefile
	+$(CXX) -Iinclude $(CXXFLAGS) $(LDFLAGS) $< $(filter %.o,$^) $(LIB) $(LDLIBS) -o $@

sdk-example: $(LIB) $(HEADERS)
	@mkdir -p $(EXAMPLE_DIR)/sdk-project/sdk/lib
	cp -R include $(EXAMPLE_DIR)/sdk-project/sdk/
	cp $(LIB) $(EXAMPLE_DIR)/sdk-project/sdk/lib/
	cp docs/examples/external_heat.cpp docs/examples/external_main.cpp $(EXAMPLE_DIR)/sdk-project/
	cp docs/examples/external_project.mk $(EXAMPLE_DIR)/sdk-project/Makefile
	+$(MAKE) -C $(EXAMPLE_DIR)/sdk-project PETSC_DIR='$(abspath $(PETSC_DIR))' PARTITION_LIBS='$(PARTITION_LIBS)'

-include $(wildcard $(EXAMPLE_DIR)/*.d)
