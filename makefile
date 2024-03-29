AGC_ROOT_DIR = .
AGC_MAIN_DIR = src
AGC_EXAMPLES_DIR = src/examples
AGC_CORE_DIR = src/core
AGC_APP_DIR = src/app
AGC_CXX_DIR = src/lib-cxx
AGC_LIBS_DIR = libs
LIBS_DIR = . #/usr/local/lib
INC_DIRS =. /usr/local/include 3rd_party/mimalloc/include 3rd_party/zstd/lib 3rd_party/zlib-ng/ 3rd_party/raduls-inplace/Raduls 3rd_party/isa-l/include 3rd_party/libdeflate
INCLUDE_DIR=$(foreach d, $(INC_DIRS), -I$d)
PY_AGC_API_DIR = py_agc_api
PYBIND11_LIB = $(PY_AGC_API_DIR)/pybind11-2.11.1/pybind11/include

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
	uname_M := "x86_64"
else                          # If uname not available => 'not'
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
	uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
endif

NASM_V := $(shell nasm --version 2>/dev/null)

RADULS_DIR = 3rd_party/raduls-inplace/Raduls
ZSTD_DIR = 3rd_party/zstd
ZLIB_DIR = 3rd_party/zlib-ng
ISAL_DIR = 3rd_party/isa-l
LIBDEFLATE_DIR = 3rd_party/libdeflate

MIMALLOC_INLUCDE_DIR = 3rd_party/mimalloc/include
MIMALLOC_OBJ=libs/mimalloc.o


ifeq ($(PLATFORM), arm8)
$(info *** ARMv8 with NEON extensions ***)
	ARCH_FLAGS := -march=armv8-a  -DARCH_ARM
else ifeq ($(PLATFORM), m1)
$(info *** Apple M1(or never) with NEON extensions ***)
	ARCH_FLAGS := -march=armv8.4-a  -DARCH_ARM
else ifeq ($(PLATFORM), sse2)
$(info *** x86-64 with SSE2 extensions ***)
	ARCH_FLAGS := -msse2 -m64 -DARCH_X64 
else ifeq ($(PLATFORM), avx)
$(info *** x86-64 with AVX extensions ***)
	ARCH_FLAGS := -mavx -m64  -DARCH_X64
else ifeq ($(PLATFORM), avx2)
$(info *** x86-64 with AVX2 extensions ***)
	ARCH_FLAGS := -mavx2 -m64  -DARCH_X64
else
$(info *** Unspecified platform - use native compilation)
	ifeq ($(uname_M),x86_64)
		ARCH_FLAGS := -march=native -DARCH_X64
	else
		ARCH_FLAGS := -march=native -DARCH_ARM
	endif	
endif


#CXX = g++ #(by default)

AR 	= ar
CFLAGS	= -fPIC -Wall -g -O3 $(ARCH_FLAGS) -std=c++17 -pthread -I $(INCLUDE_DIR) -fpermissive
#CLINK	= -lm -lz -lpthread -std=c++17
#CLINK	= -lm -lpthread -std=c++17 -lc
PY_CFLAGS = -Wl,-undefined,dynamic_lookup -fPIC -Wall -shared -std=c++17  -O3 -I $(INCLUDE_DIR)

ifeq ($(uname_S),Linux)
	CLINK	= -lm -static -O -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++17 -lc -fabi-version=6
	AR_OPT=rcs -o
	PY_AGC_API_CFLAGS = -fPIC -Wall -shared -std=c++14 -O3
endif

ifeq ($(uname_S),Darwin)
	AR_OPT=-rcs
	PY_AGC_API_CFLAGS = -Wl,-undefined,dynamic_lookup -fPIC -Wall -shared -std=c++14 -O3
	CLINK	= -lm -lpthread -std=c++17 -lc -static-libgcc
endif

LIB_ZSTD=libzstd.a
LIB_RADULS=libraduls.a
LIB_DEFLATE=libdeflate.a

ifeq ($(uname_M),x86_64)
	ifdef NASM_V
		GZ_LIB = isa-l.a
		gz_target = isa-l
		CFLAGS+=-DREFRESH_USE_IGZIP
	else
		GZ_LIB = libz.a
		gz_target = ng_zlib
		CFLAGS+=-DREFRESH_USE_ZLIB
	endif
else
	GZ_LIB = libz.a
	gz_target = ng_zlib
	CFLAGS+=-DREFRESH_USE_ZLIB
endif


# default install location (binary placed in the /bin folder)
prefix      = /usr/local

# optional install location
exec_prefix = $(prefix)


all: agc libagc py_agc_api raduls zstd $(gz_target) libdeflate

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

agc: raduls zstd $(gz_target) libdeflate $(MIMALLOC_OBJ) \
	$(AGC_APP_DIR)/main.o \
	$(AGC_CORE_DIR)/agc_basic.o \
	$(AGC_CORE_DIR)/agc_compressor.o \
	$(AGC_CORE_DIR)/agc_decompressor.o \
	$(AGC_CORE_DIR)/agc_decompressor_lib.o \
	$(AGC_APP_DIR)/application.o \
	$(AGC_CORE_DIR)/archive.o \
	$(AGC_CORE_DIR)/collection.o \
	$(AGC_CORE_DIR)/collection_v1.o \
	$(AGC_CORE_DIR)/collection_v2.o \
	$(AGC_CORE_DIR)/collection_v3.o \
	$(AGC_CORE_DIR)/genome_io.o \
	$(AGC_CORE_DIR)/lz_diff.o \
	$(AGC_CORE_DIR)/segment.o \
	$(AGC_CORE_DIR)/utils.o
	$(CXX) -o $(AGC_ROOT_DIR)/$@  \
	$(MIMALLOC_OBJ) \
	$(AGC_APP_DIR)/main.o \
	$(AGC_CORE_DIR)/agc_basic.o \
	$(AGC_CORE_DIR)/agc_compressor.o \
	$(AGC_CORE_DIR)/agc_decompressor.o \
	$(AGC_CORE_DIR)/agc_decompressor_lib.o \
	$(AGC_APP_DIR)/application.o \
	$(AGC_CORE_DIR)/archive.o \
	$(AGC_CORE_DIR)/collection.o \
	$(AGC_CORE_DIR)/collection_v1.o \
	$(AGC_CORE_DIR)/collection_v2.o \
	$(AGC_CORE_DIR)/collection_v3.o \
	$(AGC_CORE_DIR)/genome_io.o \
	$(AGC_CORE_DIR)/lz_diff.o \
	$(AGC_CORE_DIR)/segment.o \
	$(AGC_CORE_DIR)/utils.o \
	$(AGC_LIBS_DIR)/$(LIB_ZSTD) \
	$(AGC_LIBS_DIR)/$(GZ_LIB) \
	$(AGC_LIBS_DIR)/$(LIB_RADULS) \
	$(AGC_LIBS_DIR)/$(LIB_DEFLATE) \
	$(CLINK)

libagc: zstd \
	$(AGC_CXX_DIR)/lib-cxx.o \
	$(AGC_CORE_DIR)/agc_basic.o \
	$(AGC_CORE_DIR)/agc_decompressor_lib.o \
	$(AGC_CORE_DIR)/archive.o \
	$(AGC_CORE_DIR)/collection.o \
	$(AGC_CORE_DIR)/collection_v1.o \
	$(AGC_CORE_DIR)/collection_v2.o \
	$(AGC_CORE_DIR)/collection_v3.o \
	$(AGC_CORE_DIR)/genome_io.o \
	$(AGC_CORE_DIR)/lz_diff.o \
	$(AGC_CORE_DIR)/segment.o \
	$(AGC_CORE_DIR)/utils.o
	$(AR) $(AR_OPT)  $(AGC_ROOT_DIR)/$@.a  \
	$(AGC_CXX_DIR)/lib-cxx.o \
	$(AGC_CORE_DIR)/agc_basic.o \
	$(AGC_CORE_DIR)/agc_decompressor_lib.o \
	$(AGC_CORE_DIR)/archive.o \
	$(AGC_CORE_DIR)/collection.o \
	$(AGC_CORE_DIR)/collection_v1.o \
	$(AGC_CORE_DIR)/collection_v2.o \
	$(AGC_CORE_DIR)/collection_v3.o \
	$(AGC_CORE_DIR)/genome_io.o \
	$(AGC_CORE_DIR)/lz_diff.o \
	$(AGC_CORE_DIR)/segment.o \
	$(AGC_CORE_DIR)/utils.o 

.PHONY:py_agc_api
py_agc_api: zstd \
    $(PY_AGC_API_DIR)/py_agc_api.cpp $(AGC_CXX_DIR)/lib-cxx.o \
	$(AGC_CORE_DIR)/agc_basic.o \
	$(AGC_CORE_DIR)/agc_decompressor_lib.o \
	$(AGC_CORE_DIR)/archive.o \
	$(AGC_CORE_DIR)/collection.o \
	$(AGC_CORE_DIR)/collection_v1.o \
	$(AGC_CORE_DIR)/collection_v2.o \
	$(AGC_CORE_DIR)/collection_v3.o \
	$(AGC_CORE_DIR)/lz_diff.o \
	$(AGC_CORE_DIR)/segment.o \
	$(AGC_CORE_DIR)/utils.o
	$(CXX) $(PY_CFLAGS)  \
	$(PY_AGC_API_DIR)/py_agc_api.cpp \
	$(AGC_CXX_DIR)/lib-cxx.o \
	$(AGC_CORE_DIR)/agc_basic.o \
	$(AGC_CORE_DIR)/agc_decompressor_lib.o \
	$(AGC_CORE_DIR)/archive.o \
	$(AGC_CORE_DIR)/collection.o \
	$(AGC_CORE_DIR)/collection_v1.o \
	$(AGC_CORE_DIR)/collection_v2.o \
	$(AGC_CORE_DIR)/collection_v3.o \
	$(AGC_CORE_DIR)/lz_diff.o \
	$(AGC_CORE_DIR)/segment.o \
	$(AGC_CORE_DIR)/utils.o \
	$(AGC_LIBS_DIR)/$(LIB_ZSTD) \
	-I $(AGC_MAIN_DIR) \
	-I $(PYBIND11_LIB) \
	-I `python3 -c "import sysconfig;print(sysconfig.get_paths()['include'])"` \
	-o $@`python3-config --extension-suffix`

raduls:
	cd $(RADULS_DIR) && $(MAKE)
	cp $(RADULS_DIR)/libraduls.a $(AGC_LIBS_DIR)

zstd:
	cd $(ZSTD_DIR) && $(MAKE) lib
	cp $(ZSTD_DIR)/lib/libzstd.* $(AGC_LIBS_DIR)
	
ng_zlib:
	cd $(ZLIB_DIR) && ./configure --zlib-compat && $(MAKE)
	cp $(ZLIB_DIR)/libz.* $(AGC_LIBS_DIR)

libdeflate:
	cd $(LIBDEFLATE_DIR) && cmake -B build && cmake --build build
	cp $(LIBDEFLATE_DIR)/build/libdeflate.* $(AGC_LIBS_DIR)

isa-l:
	cd $(ISAL_DIR) && $(MAKE) -f Makefile.unx
	cp $(ISAL_DIR)/bin/isa-l.a $(AGC_LIBS_DIR)
	cp $(ISAL_DIR)/bin/libisal.* $(AGC_LIBS_DIR)

$(MIMALLOC_OBJ):
	$(CC) -DMI_MALLOC_OVERRIDE -O3 -DNDEBUG -fPIC -Wall -Wextra -Wno-unknown-pragmas -fvisibility=hidden -Wstrict-prototypes -ftls-model=initial-exec -fno-builtin-malloc -std=gnu11 -c -I 3rd_party/mimalloc/include 3rd_party/mimalloc/src/static.c -o $(MIMALLOC_OBJ)

clean:
	-rm $(AGC_EXAMPLES_DIR)/*.o
	-rm $(AGC_APP_DIR)/*.o
	-rm $(AGC_CORE_DIR)/*.o
	-rm $(AGC_CXX_DIR)/*.o
	-rm agc
	-rm libagc.a
	-rm -f $(PY_AGC_API_DIR)/*.o
	-rm *.so
	-rm $(AGC_LIBS_DIR)/libraduls.*
	-rm $(AGC_LIBS_DIR)/libzstd.*
	-rm $(AGC_LIBS_DIR)/libz.*
	-rm $(AGC_LIBS_DIR)/isa-l.*
	-rm $(AGC_LIBS_DIR)/libisal.*
	-rm $(AGC_LIBS_DIR)/mimalloc.*
	-rm $(AGC_LIBS_DIR)/libdeflate.*
	cd $(RADULS_DIR) && $(MAKE) clean
	cd $(ZSTD_DIR) && $(MAKE) clean
	cd $(ZLIB_DIR) && $(MAKE) -f Makefile.in clean
	cd $(ISAL_DIR) && $(MAKE) -f Makefile.unx clean
	-cd $(LIBDEFLATE_DIR) && rm -r build
