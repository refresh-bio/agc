# all: agc libagc py_agc_api
all: agc libagc py_agc_api

# *** REFRESH makefile utils
include refresh.mk

$(call INIT_SUBMODULES)
$(call INIT_GLOBALS)
$(call CHECK_OS_ARCH, $(PLATFORM))

# *** Project directories
$(call SET_SRC_OBJ_BIN,src,obj,bin)
3RD_PARTY_DIR := ./3rd_party

# *** Project configuration
$(call CHECK_NASM)
$(call ADD_MIMALLOC, $(3RD_PARTY_DIR)/mimalloc)
$(call PROPOSE_ISAL, $(3RD_PARTY_DIR)/isa-l)
$(call PROPOSE_ZLIB_NG, $(3RD_PARTY_DIR)/zlib-ng)
$(call CHOOSE_GZIP_DECOMPRESSION)
$(call ADD_LIBDEFLATE, $(3RD_PARTY_DIR)/libdeflate)
$(call ADD_LIBZSTD, $(3RD_PARTY_DIR)/zstd)
$(call ADD_RADULS_INPLACE,$(3RD_PARTY_DIR)/raduls-inplace)
$(call ADD_PYBIND11,$(3RD_PARTY_DIR)/pybind11/include)
$(call SET_STATIC, $(STATIC_LINK))

$(call SET_C_CPP_STANDARDS, c11, c++20)
$(call ADD_REFRESH_LIB, $(3RD_PARTY_DIR))

$(call SET_GIT_COMMIT)

$(call SET_FLAGS, $(TYPE))

$(call SET_COMPILER_VERSION_ALLOWED, GCC, Linux_x86_64, 10, 20)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Linux_aarch64, 11, 20)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Darwin_x86_64, 11, 13)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Darwin_arm64, 11, 13)

ifneq ($(MAKECMDGOALS),clean)
$(call CHECK_COMPILER_VERSION)
endif

# *** Source files and rules
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,APP,app))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,CORE,core))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,COMMON,common))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,EXAMPLES,examples))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,LIB_CXX,lib-cxx))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,PY_AGC_API,py_agc_api,$(PY_FLAGS)))


# *** Targets
agc: $(OUT_BIN_DIR)/agc
$(OUT_BIN_DIR)/agc: \
	$(OBJ_APP) $(OBJ_CORE) $(OBJ_COMMON)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) -o $@  \
	$(MIMALLOC_OBJ) \
	$(OBJ_APP) $(OBJ_CORE) $(OBJ_COMMON) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

libagc: $(OUT_BIN_DIR)/libagc
$(OUT_BIN_DIR)/libagc: \
	$(OBJ_LIB_CXX) $(OBJ_COMMON)
	-mkdir -p $(OUT_BIN_DIR)
	$(AR) $(AR_OPT) $@.a  \
	$(OBJ_LIB_CXX) $(OBJ_COMMON)


#.PHONY:py_agc_api
py_agc_api: $(OUT_BIN_DIR)/py_agc_api
$(OUT_BIN_DIR)/py_agc_api: \
	$(OBJ_PY_AGC_API) $(OBJ_LIB_CXX) $(OBJ_COMMON)
	-mkdir -p $(OUT_BIN_DIR)
	$(CXX) $(PY_FLAGS) $(INCLUDE_DIRS) \
	$(OBJ_PY_AGC_API) $(OBJ_LIB_CXX) $(OBJ_COMMON) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS) \
	-o $@$(PY_EXTENSION_SUFFIX)


# *** Cleaning
.PHONY: clean init
clean: clean-libzstd clean-zlib-ng clean-isa-l clean-libdeflate clean-mimalloc_obj
	-rm -r $(OBJ_DIR)
	-rm -r $(OUT_BIN_DIR)

init:
	$(call INIT_SUBMODULES)
