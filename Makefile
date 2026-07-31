# ==== Toolchain detection ====
# The platform selects the toolchain automatically:
#   Linux  -> intel       (icpx, cluster modules; the default on main)
#   Darwin -> appleclang  (Apple clang + Homebrew netcdf/boost/libomp)
#
# Overrides, all settable on the command line:
#   make TOOLCHAIN=intel|appleclang   force a toolchain
#   make CXX=g++-15                   force a compiler (env CXX works too)
#   make NETCDF_PATH=/path/to/netcdf  skip nc-config discovery
#   make BOOST_PATH=... LIBOMP_PATH=... BREW_PREFIX=...   (appleclang only)
#   make ARCH_FLAGS=-xHost            build for a non-Sapphire-Rapids Intel node
# Run `make config` to print what was detected.

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  TOOLCHAIN ?= appleclang
else ifeq ($(UNAME_S),Linux)
  TOOLCHAIN ?= intel
else
  $(error Unsupported platform '$(UNAME_S)'. Run: make TOOLCHAIN=appleclang  or  make TOOLCHAIN=intel)
endif

# ==== Compiler ====
# `CXX ?=` cannot be used here: make predefines CXX, so its origin is "default"
# and ?= would never assign. Testing the origin also lets an environment or
# command-line CXX win.
ifeq ($(origin CXX),default)
  ifeq ($(TOOLCHAIN),appleclang)
    CXX := clang++
  else
    CXX := icpx
  endif
endif

# ==== Per-toolchain flags ====
ifeq ($(TOOLCHAIN),intel)

  ARCH_FLAGS ?= -xSapphireRapids
  DIAG_FLAGS ?= -qopt-report-phase=vec -qopt-prefetch -Rpass=loop-vectorize -Rpass=inline
  CXXFLAGS := -O3 -ipo -fp-model fast=2 -qopenmp -fma $(ARCH_FLAGS) $(DIAG_FLAGS) -DNDEBUG -std=c++17

  # -qopenmp above covers both compiling and linking (the link rule uses CXXFLAGS).
  OMP_CPPFLAGS :=
  OMP_LDFLAGS  :=

  # Boost headers come from the module environment (CPATH), so no -I is needed.
  # Boost.odeint is header-only; -lboost_system is kept to match the cluster
  # build that this branch inherited.
  BOOST_CPPFLAGS :=
  BOOST_LDFLAGS  := -lboost_system

else ifeq ($(TOOLCHAIN),appleclang)

  # Homebrew prefix (arm64: /opt/homebrew, Intel Mac: /usr/local).
  BREW_PREFIX ?= $(shell brew --prefix 2>/dev/null)
  BOOST_PATH  ?= $(BREW_PREFIX)/opt/boost
  LIBOMP_PATH ?= $(BREW_PREFIX)/opt/libomp

  CXXFLAGS := -O3 -ffast-math -std=c++17 -DNDEBUG

  # Apple clang needs libomp from Homebrew.
  OMP_CPPFLAGS := -Xpreprocessor -fopenmp -I$(LIBOMP_PATH)/include
  OMP_LDFLAGS  := -L$(LIBOMP_PATH)/lib -lomp

  BOOST_CPPFLAGS := -I$(BOOST_PATH)/include
  BOOST_LDFLAGS  :=

else
  $(error Unknown TOOLCHAIN '$(TOOLCHAIN)'. Valid values: intel, appleclang)
endif

# ==== NetCDF ====
# An explicit NETCDF_PATH wins (the cluster netcdf module sets it); otherwise
# ask nc-config, which tracks Homebrew version bumps automatically.
ifneq ($(NETCDF_PATH),)
  NETCDF_LIBDIR   := $(if $(wildcard $(NETCDF_PATH)/lib64),$(NETCDF_PATH)/lib64,$(NETCDF_PATH)/lib)
  NETCDF_CPPFLAGS := -I$(NETCDF_PATH)/include
  NETCDF_LDFLAGS  := -L$(NETCDF_LIBDIR) -lnetcdf
else
  NC_CONFIG := $(shell command -v nc-config 2>/dev/null)
  ifeq ($(NC_CONFIG),)
    $(error Cannot locate NetCDF. Load the netcdf module, put nc-config on PATH, or run: make NETCDF_PATH=/path/to/netcdf)
  endif
  NETCDF_CPPFLAGS := $(shell $(NC_CONFIG) --cflags)
  NETCDF_LDFLAGS  := $(shell $(NC_CONFIG) --libs)
endif

# ==== Composed flags ====
CPPFLAGS := $(NETCDF_CPPFLAGS) $(BOOST_CPPFLAGS) $(OMP_CPPFLAGS)
LDFLAGS  := $(BOOST_LDFLAGS) $(NETCDF_LDFLAGS) $(OMP_LDFLAGS)

# ==== Source Files ====
SRC := src/main.cpp \
       src/build_info.cpp \
       src/omp_info.cpp \
       src/model_setup.cpp \
       src/routing.cpp \
       src/end_info.cpp \
       src/I_O/node_info.cpp \
       src/I_O/output_series.cpp \
       src/I_O/inputs.cpp \
       src/I_O/config_loader.cpp \
       src/utils/time.cpp

# ==== Build and Binary Directories ====
BUILD_DIR := build
BIN_DIR := bin

# ==== Object files in build dir ====
OBJ := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SRC))

# ==== Executable ====
BIN := $(BIN_DIR)/routing

# ==== Default Target ====
all: $(BIN)
	@echo "Build successful: $(BIN)"

# ==== Link executable ====
$(BIN): $(OBJ)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# ==== Compile source files ====
$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# ==== Clean ====
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# ==== Show detected build configuration ====
config:
	@echo "platform  : $(UNAME_S)"
	@echo "toolchain : $(TOOLCHAIN)"
	@echo "CXX       : $(CXX)  (origin: $(origin CXX))"
	@echo "CXXFLAGS  : $(CXXFLAGS)"
	@echo "CPPFLAGS  : $(CPPFLAGS)"
	@echo "LDFLAGS   : $(LDFLAGS)"

.PHONY: all clean config
