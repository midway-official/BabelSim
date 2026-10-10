# Standalone project: needs only public headers, libbabelsim.a and runtime libraries.
CXX := g++
MPI_CXX ?= mpic++
BABELSIM_SDK ?= sdk
PETSC_DIR ?= ../../petsc
PARTITION_LIBS ?=
PETSC_LIBS = $(shell PKG_CONFIG_PATH="$(abspath $(PETSC_DIR))/lib/pkgconfig" pkg-config --libs PETSc)
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror
OBJECTS := external_heat.o external_main.o

.DEFAULT_GOAL := custom-heat
.PHONY: clean
custom-heat: $(OBJECTS) $(BABELSIM_SDK)/lib/libbabelsim.a
	$(MPI_CXX) $(OBJECTS) $(BABELSIM_SDK)/lib/libbabelsim.a $(PETSC_LIBS) \
		-Wl,-rpath,$(abspath $(PETSC_DIR))/lib $(PARTITION_LIBS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(BABELSIM_SDK)/include -MMD -MP -c $< -o $@

clean:
	$(RM) $(OBJECTS) $(OBJECTS:.o=.d) custom-heat

-include $(OBJECTS:.o=.d)
