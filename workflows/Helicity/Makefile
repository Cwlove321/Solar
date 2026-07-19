FC = ifort
FFLAGS = -O3
BUILD_DIR ?= build
TARGET ?= helicity

MKL_INCLUDE := $(MKLROOT)/include
MKL_LIBDIR := $(MKLROOT)/lib/intel64
MODULE_FLAGS := -module $(BUILD_DIR) -I$(BUILD_DIR) -I$(MKL_INCLUDE)
MKL_LIBS := -L$(MKL_LIBDIR) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

OBJECTS := \
	$(BUILD_DIR)/mkl_dfti.o \
	$(BUILD_DIR)/mkl_poisson.o \
	$(BUILD_DIR)/helicity_numerics.o \
	$(BUILD_DIR)/helicity_diagnostics.o \
	$(BUILD_DIR)/helicity.o

.PHONY: all check-env clean

all: check-env $(TARGET)

check-env:
	@test -n "$(MKLROOT)" || (echo "MKLROOT 未设置；请先加载 Intel oneAPI/MKL 环境。" && exit 1)
	@test -f "$(MKL_INCLUDE)/mkl_dfti.f90" || (echo "找不到 $(MKL_INCLUDE)/mkl_dfti.f90" && exit 1)
	@test -f "$(MKL_INCLUDE)/mkl_poisson.f90" || (echo "找不到 $(MKL_INCLUDE)/mkl_poisson.f90" && exit 1)
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/mkl_dfti.o: | check-env
	$(FC) $(FFLAGS) -c $(MKL_INCLUDE)/mkl_dfti.f90 $(MODULE_FLAGS) -o $@

$(BUILD_DIR)/mkl_poisson.o: $(BUILD_DIR)/mkl_dfti.o
	$(FC) $(FFLAGS) -c $(MKL_INCLUDE)/mkl_poisson.f90 $(MODULE_FLAGS) -o $@

$(BUILD_DIR)/helicity_numerics.o: helicity_numerics.f90 $(BUILD_DIR)/mkl_poisson.o
	$(FC) $(FFLAGS) -c $< $(MODULE_FLAGS) -o $@

$(BUILD_DIR)/helicity_diagnostics.o: helicity_diagnostics.f90 $(BUILD_DIR)/helicity_numerics.o
	$(FC) $(FFLAGS) -c $< $(MODULE_FLAGS) -o $@

$(BUILD_DIR)/helicity.o: helicity.f90 $(BUILD_DIR)/helicity_diagnostics.o
	$(FC) $(FFLAGS) -c $< $(MODULE_FLAGS) -o $@

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) $(MKL_LIBS) -o $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET)
