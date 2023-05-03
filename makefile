CXX = g++
C99 = gcc -std=c99
LINK = g++
AR = ar
#DEBUG_FLAG=-g
CXXFLAGS = -O1 -Wall -fPIC $(DEBUG_FLAG)
CFLAGS = $(CXXFLAGS)
ARFLAGS = -rv
OUT_DIR = ./build
OBJS = $(OUT_DIR)/objs/world/cheaptrick.o $(OUT_DIR)/objs/world/common.o $(OUT_DIR)/objs/world/d4c.o $(OUT_DIR)/objs/world/dio.o $(OUT_DIR)/objs/world/fft.o $(OUT_DIR)/objs/world/harvest.o $(OUT_DIR)/objs/world/matlabfunctions.o $(OUT_DIR)/objs/world/stonemask.o $(OUT_DIR)/objs/world/synthesis.o $(OUT_DIR)/objs/world/synthesisrealtime.o
LIBS =
MKDIR = mkdir -p $(1)
ifeq ($(shell echo "check_quotes"),"check_quotes")
	# Windows
	MKDIR = mkdir $(subst /,\,$(1)) > nul 2>&1 || (exit 0)
endif

all: default test

###############################################################################################################
### Tests
###############################################################################################################
test: $(OUT_DIR)/test $(OUT_DIR)/ctest

test_OBJS=$(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/test/test.o
$(OUT_DIR)/test: $(OUT_DIR)/libworld.a $(test_OBJS)
	$(LINK) $(CXXFLAGS) -o $(OUT_DIR)/test $(test_OBJS) $(OUT_DIR)/libworld.a -lm

ctest_OBJS=$(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/test/ctest.o
$(OUT_DIR)/ctest: $(OUT_DIR)/libworld.a $(ctest_OBJS)
	$(LINK) $(CXXFLAGS) -o $(OUT_DIR)/ctest $(ctest_OBJS) $(OUT_DIR)/libworld.a -lm

$(OUT_DIR)/objs/test/test.o : tools/audioio.h src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld
$(OUT_DIR)/objs/test/ctest.o : tools/audioio.h src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld

###############################################################################################################
### Library
###############################################################################################################
default: $(OUT_DIR)/libworld.a

$(OUT_DIR)/libworld.a: $(OBJS)
	$(AR) $(ARFLAGS) $(OUT_DIR)/libworld.a $(OBJS) $(LIBS)
	@echo Done.

$(OUT_DIR)/objs/world/cheaptrick.o : src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld
$(OUT_DIR)/objs/world/common.o : src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld
$(OUT_DIR)/objs/world/d4c.o : src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld
$(OUT_DIR)/objs/world/dio.o : src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld
$(OUT_DIR)/objs/world/fft.o : src/TorchWorld src/TorchWorld
$(OUT_DIR)/objs/world/harvest.o : src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld
$(OUT_DIR)/objs/world/matlabfunctions.o : src/TorchWorld src/TorchWorld src/TorchWorld
$(OUT_DIR)/objs/world/stonemask.o : src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld
$(OUT_DIR)/objs/world/synthesis.o : src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld
$(OUT_DIR)/objs/world/synthesisrealtime.o : src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld src/TorchWorld


###############################################################################################################
### Global rules
###############################################################################################################
$(OUT_DIR)/objs/test/%.o : test/%.c
	$(call MKDIR,$(OUT_DIR)/objs/test)
	$(C99) $(CFLAGS) -Isrc -Itools -o "$@" -c "$<"

$(OUT_DIR)/objs/test/%.o : test/%.cpp
	$(call MKDIR,$(OUT_DIR)/objs/test)
	$(CXX) $(CXXFLAGS) -Isrc -Itools -o "$@" -c "$<"

$(OUT_DIR)/objs/tools/%.o : tools/%.cpp
	$(call MKDIR,$(OUT_DIR)/objs/tools)
	$(CXX) $(CXXFLAGS) -Isrc -o "$@" -c "$<"

$(OUT_DIR)/objs/world/%.o : src/%.cpp
	$(call MKDIR,$(OUT_DIR)/objs/world)
	$(CXX) $(CXXFLAGS) -Isrc -o "$@" -c "$<"

clean:
	@echo 'Removing all temporary binaries... '
	@$(RM) $(OUT_DIR)/libworld.a $(OBJS)
	@$(RM) $(test_OBJS) $(ctest_OBJS) $(OUT_DIR)/test $(OUT_DIR)/ctest
	@echo Done.

clear: clean

.PHONY: clean clear test default
.DELETE_ON_ERRORS:
