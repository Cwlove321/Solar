GLOBAL_FMHD_LIBDIR ?= /share/home/guoyang1/apps/global_fmhd_deps/chenpf_libs
INTEL_COMPILER_LIBDIR ?= /share/apps/intel/oneapi_2023.1.0/compiler/2023.1.0/linux/compiler/lib/intel64_lin

LIBS += pdfi_ss
LIBS += fishpack
LIBS += cfitsio
LIB_DIRS += $(GLOBAL_FMHD_LIBDIR)

LIBS += ifcore
LIBS += ifport
LIBS += imf
LIBS += svml
LIBS += intlc
LIBS += irc
LIB_DIRS += $(INTEL_COMPILER_LIBDIR)
