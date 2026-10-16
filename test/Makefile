# PhyML test suite.
#
#   make -C test            build the binaries used by the tests
#   make -C test check      build, then run the unit and functional tests
#   make -C test check-unit / check-functional / check-asan
#
# The suite builds its own copies of PhyML from ../src (it does not need
# autotools), so it always tests the sources in the working tree:
#
#   build/phyml          default build (SIMD kernels when the target has them)
#   build/phyml-generic  -DDISABLE_NATIVE: scalar likelihood kernels only
#   build/phyml-asan     AddressSanitizer + UBSan, used by memory-safety tests
#   build/phyml-ndebug   -DNDEBUG, as a distribution build might be compiled
#   build/phyml-omp      OpenMP (--threads), when the compiler supports it
#
# Variables: CC, OPTFLAGS (default -O2), ARCHFLAGS (default -march=native on
# x86-64, empty elsewhere),
# PYTHON (default python3), TESTS (unittest pattern, e.g. TESTS='test_cli*'),
# SANITIZE=0 to skip the sanitizer build when the toolchain lacks ASan,
# OPENMP=0 to skip the OpenMP build (default: built when OPENMP_CFLAGS and
# OPENMP_LIBS compile a test program; on macOS they point to Homebrew libomp).

ROOT      := $(abspath ..)
SRC       := $(ROOT)/src
BUILD     := build
CC        ?= cc
PYTHON    ?= python3
OPTFLAGS  ?= -O2
# configure defaults to -march=native. On x86-64 the SSE kernels need it (or
# at least -msse3): they are enabled by __SSE2__ but call SSE3 intrinsics.
ifeq ($(shell uname -m),x86_64)
ARCHFLAGS ?= -march=native
else
ARCHFLAGS ?=
endif
TESTS     ?= test_*.py
SANITIZE  ?= 1

ifeq ($(shell uname -s),Darwin)
LIBOMP        ?= $(shell brew --prefix libomp 2>/dev/null)
OPENMP_CFLAGS ?= -Xpreprocessor -fopenmp -I$(LIBOMP)/include
OPENMP_LIBS   ?= -L$(LIBOMP)/lib -lomp
else
OPENMP_CFLAGS ?= -fopenmp
OPENMP_LIBS   ?= -fopenmp
endif
ifndef OPENMP
OPENMP := $(shell printf '\043include <omp.h>\nint main(void){return omp_get_max_threads() < 1;}\n' | \
            $(CC) $(OPENMP_CFLAGS) -x c - $(OPENMP_LIBS) -o /dev/null >/dev/null 2>&1 && echo 1 || echo 0)
endif

# Same sources as phyml_SOURCES in src/Makefile.am, minus main.c.
LIB_SRCS := utilities optimiz lk bionj models free help simu eigen pars alrt \
            interface cl spr draw stats io make nexus init xml sse avx tbe   \
            ancestral mixt evolve rates date times m4 cv sms structure lg4

DEFS   := -DUNIX -DPHYML -DHAVE_CONFIG_H -I$(ROOT) -I$(SRC)
BASE   := -std=gnu99 -w $(DEFS)
LIBS   := -lm

VARIANTS := default generic asan ndebug omp
CFLAGS_default := $(OPTFLAGS) $(ARCHFLAGS) -g
CFLAGS_generic := $(OPTFLAGS) -g -DDISABLE_NATIVE
CFLAGS_asan    := -O1 -g $(ARCHFLAGS) -fno-omit-frame-pointer -fsanitize=address,undefined \
                  -fsanitize-recover=address,undefined
CFLAGS_ndebug  := $(OPTFLAGS) $(ARCHFLAGS) -DNDEBUG
CFLAGS_omp     := $(OPTFLAGS) $(ARCHFLAGS) -g $(OPENMP_CFLAGS)
LDFLAGS_asan   := -fsanitize=address,undefined -fsanitize-recover=address,undefined

UNIT_SRCS := $(wildcard unit/test_*.c)
UNIT_BIN  := $(BUILD)/unit_tests

.PHONY: all check check-unit check-functional check-asan check-mpi clean data
.SECONDARY:

BINS := $(BUILD)/phyml $(BUILD)/phyml-generic $(BUILD)/phyml-ndebug
ifeq ($(SANITIZE),1)
BINS += $(BUILD)/phyml-asan
endif
ifeq ($(OPENMP),1)
BINS += $(BUILD)/phyml-omp
endif

all: $(BINS) $(UNIT_BIN)

# ---- object files, one directory per variant ---------------------------
define VARIANT_RULES
$(BUILD)/obj-$(1)/%.o: $(SRC)/%.c $(wildcard $(SRC)/*.h) | $(BUILD)/obj-$(1)
	$$(CC) $$(BASE) $$(CFLAGS_$(1)) -c $$< -o $$@

$(BUILD)/obj-$(1):
	mkdir -p $$@

$(BUILD)/libphyml-$(1).a: $(addprefix $(BUILD)/obj-$(1)/,$(addsuffix .o,$(LIB_SRCS)))
	rm -f $$@
	ar rcs $$@ $$^
endef
$(foreach v,$(VARIANTS),$(eval $(call VARIANT_RULES,$(v))))

$(BUILD)/phyml: $(BUILD)/obj-default/main.o $(BUILD)/libphyml-default.a
	$(CC) $(CFLAGS_default) $^ $(LIBS) -o $@

$(BUILD)/phyml-generic: $(BUILD)/obj-generic/main.o $(BUILD)/libphyml-generic.a
	$(CC) $(CFLAGS_generic) $^ $(LIBS) -o $@

$(BUILD)/phyml-asan: $(BUILD)/obj-asan/main.o $(BUILD)/libphyml-asan.a
	$(CC) $(CFLAGS_asan) $(LDFLAGS_asan) $^ $(LIBS) -o $@

$(BUILD)/phyml-ndebug: $(BUILD)/obj-ndebug/main.o $(BUILD)/libphyml-ndebug.a
	$(CC) $(CFLAGS_ndebug) $^ $(LIBS) -o $@

$(BUILD)/phyml-omp: $(BUILD)/obj-omp/main.o $(BUILD)/libphyml-omp.a
	$(CC) $(CFLAGS_omp) $^ $(OPENMP_LIBS) $(LIBS) -o $@

# ---- unit tests -----------------------------------------------------------
$(UNIT_BIN): $(UNIT_SRCS) unit/tinytest.h unit/tinytest.c $(BUILD)/libphyml-default.a
	$(CC) $(BASE) $(CFLAGS_default) -Iunit unit/tinytest.c $(UNIT_SRCS) \
	  $(BUILD)/libphyml-default.a $(LIBS) -o $@

# ---- running --------------------------------------------------------------
TEST_ENV = PHYML=$(abspath $(BUILD)/phyml) \
           PHYML_GENERIC=$(abspath $(BUILD)/phyml-generic) \
           PHYML_NDEBUG=$(abspath $(BUILD)/phyml-ndebug) \
           $(if $(filter 1,$(OPENMP)),PHYML_OMP=$(abspath $(BUILD)/phyml-omp)) \
           $(if $(filter 1,$(SANITIZE)),PHYML_ASAN=$(abspath $(BUILD)/phyml-asan) \
           ASAN_OPTIONS=detect_leaks=0:detect_stack_use_after_return=0)

check: all
	@status=0; \
	$(MAKE) --no-print-directory check-unit || status=1; \
	$(MAKE) --no-print-directory check-functional || status=1; \
	exit $$status

check-unit: $(UNIT_BIN)
	./$(UNIT_BIN)

check-functional: $(BINS)
	cd functional && $(TEST_ENV) $(PYTHON) -m unittest discover -p '$(TESTS)' -v

# Runs the whole functional suite (or TESTS=...) with the sanitizer build as
# the main binary. Slow; not part of `check`, which already runs the
# memory-safety tests against build/phyml-asan.
check-asan: $(BUILD)/phyml-asan
	cd functional && PHYML=$(abspath $(BUILD)/phyml-asan) \
	  PHYML_ASAN=$(abspath $(BUILD)/phyml-asan) \
	  ASAN_OPTIONS=detect_leaks=0:abort_on_error=0:detect_stack_use_after_return=0 \
	  $(PYTHON) -m unittest discover -p '$(TESTS)' -v

# Optional MPI build and tests (needs an MPI compiler and mpirun):
#   make -C test check-mpi [MPICC=mpicc] [MPIRUN="mpirun --oversubscribe"]
MPICC  ?= mpicc
MPIRUN ?= mpirun
MPI_SRCS := main $(LIB_SRCS) mpi_boot

$(BUILD)/obj-mpi/%.o: $(SRC)/%.c $(wildcard $(SRC)/*.h) | $(BUILD)/obj-mpi
	$(MPICC) $(BASE) -DMPI $(OPTFLAGS) $(ARCHFLAGS) -g -c $< -o $@

$(BUILD)/obj-mpi:
	mkdir -p $@

$(BUILD)/phyml-mpi: $(addprefix $(BUILD)/obj-mpi/,$(addsuffix .o,$(MPI_SRCS)))
	$(MPICC) $^ $(LIBS) -o $@

check-mpi: $(BUILD)/phyml $(BUILD)/phyml-mpi
	cd functional && PHYML=$(abspath $(BUILD)/phyml) PHYML_MPI=$(abspath $(BUILD)/phyml-mpi) \
	  MPIRUN="$(MPIRUN)" $(PYTHON) -m unittest test_mpi -v

data:
	$(PYTHON) tools/make_data.py

clean:
	rm -rf $(BUILD)
