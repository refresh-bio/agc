all: agc libagc py_agc_api

AGC_ROOT_DIR = .
AGC_MAIN_DIR = src
AGC_EXAMPLES_DIR = src/examples
AGC_CORE_DIR = src/core
AGC_APP_DIR = src/app
AGC_CXX_DIR = src/lib-cxx
AGC_LIBS_DIR = libs
LIBS_DIR = . #/usr/local/lib
INCLUDE_DIR= . #/usr/local/include
PY_AGC_API_DIR = py_agc_api
PYBIND11_LIB = $(PY_AGC_API_DIR)/pybind11-2.8.1

CC 	= g++
AR 	= ar
CFLAGS	= -fPIC -Wall -g -O3 -m64 -std=c++17 -pthread -mavx -I $(INCLUDE_DIR) -fpermissive
#CLINK	= -lm -static -O -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++17 
#CLINK	= -lm -lz -lpthread -std=c++17
CLINK	= -lm -lpthread -std=c++17
PY_CFLAGS = -Wl,-undefined,dynamic_lookup -fPIC -Wall -shared -std=c++17  -O3 

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
else                          # If uname not available => 'not' 
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
endif
ifeq ($(uname_S),Linux)
	CLINK+=-fabi-version=6
	LIB_ALLOC=libmimalloc.a
	LIB_ZSTD=libzstd.a
	LIB_ZLIB=cloudflare-zlib/libz.a
	LIB_RADULS=libraduls.a
	AR_OPT=rcs -o
	PY_AGC_API_CFLAGS = -fPIC -Wall -shared -std=c++14 -O3
endif

ifeq ($(uname_S),Darwin)
	LIB_ALLOC=libmimalloc.mac.a
	LIB_ZSTD=libzstd.mac.a
	LIB_ZLIB=cloudflare-zlib/libz.mac.a
	LIB_RADULS=libraduls.mac.a
	AR_OPT=-rcs
	CC = g++-11
	PY_AGC_API_CFLAGS = -Wl,-undefined,dynamic_lookup -fPIC -Wall -shared -std=c++14 -O3
endif

# default install location (binary placed in the /bin folder)
prefix      = /usr/local

# optional install location
exec_prefix = $(prefix)


%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

agc: $(AGC_APP_DIR)/main.o \
	$(AGC_CORE_DIR)/agc_basic.o \
	$(AGC_CORE_DIR)/agc_compressor.o \
	$(AGC_CORE_DIR)/agc_decompressor.o \
	$(AGC_CORE_DIR)/agc_decompressor_lib.o \
	$(AGC_APP_DIR)/application.o \
	$(AGC_CORE_DIR)/archive.o \
	$(AGC_CORE_DIR)/collection.o \
	$(AGC_CORE_DIR)/genome_io.o \
	$(AGC_CORE_DIR)/lz_diff.o \
	$(AGC_CORE_DIR)/segment.o
	$(CC) -o $(AGC_ROOT_DIR)/$@  \
	$(AGC_APP_DIR)/main.o \
	$(AGC_CORE_DIR)/agc_basic.o \
	$(AGC_CORE_DIR)/agc_compressor.o \
	$(AGC_CORE_DIR)/agc_decompressor.o \
	$(AGC_CORE_DIR)/agc_decompressor_lib.o \
	$(AGC_APP_DIR)/application.o \
	$(AGC_CORE_DIR)/archive.o \
	$(AGC_CORE_DIR)/collection.o \
	$(AGC_CORE_DIR)/genome_io.o \
	$(AGC_CORE_DIR)/lz_diff.o \
	$(AGC_CORE_DIR)/segment.o \
	$(AGC_LIBS_DIR)/$(LIB_ZSTD) \
	$(AGC_LIBS_DIR)/$(LIB_ZLIB) \
	$(AGC_LIBS_DIR)/$(LIB_RADULS) \
	$(AGC_LIBS_DIR)/mimalloc/$(LIB_ALLOC) \
	$(CLINK)

libagc: $(AGC_CXX_DIR)/lib-cxx.o \
	$(AGC_CORE_DIR)/agc_basic.o \
	$(AGC_CORE_DIR)/agc_decompressor_lib.o \
	$(AGC_CORE_DIR)/archive.o \
	$(AGC_CORE_DIR)/collection.o \
	$(AGC_CORE_DIR)/genome_io.o \
	$(AGC_CORE_DIR)/lz_diff.o \
	$(AGC_CORE_DIR)/segment.o
	$(AR) $(AR_OPT)  $(AGC_ROOT_DIR)/$@.a  \
	$(AGC_CXX_DIR)/lib-cxx.o \
	$(AGC_CORE_DIR)/agc_basic.o \
	$(AGC_CORE_DIR)/agc_decompressor_lib.o \
	$(AGC_CORE_DIR)/archive.o \
	$(AGC_CORE_DIR)/collection.o \
	$(AGC_CORE_DIR)/genome_io.o \
	$(AGC_CORE_DIR)/lz_diff.o \
	$(AGC_CORE_DIR)/segment.o \
	$(AGC_LIBS_DIR)/$(LIB_ZSTD)

.PHONY:py_agc_api
py_agc_api: $(PY_AGC_API_DIR)/py_agc_api.cpp $(AGC_CXX_DIR)/lib-cxx.o \
	$(AGC_CORE_DIR)/agc_basic.o \
	$(AGC_CORE_DIR)/agc_decompressor_lib.o \
	$(AGC_CORE_DIR)/archive.o \
	$(AGC_CORE_DIR)/collection.o \
	$(AGC_CORE_DIR)/lz_diff.o \
	$(AGC_CORE_DIR)/segment.o
	$(CC) $(PY_CFLAGS)  \
	$(PY_AGC_API_DIR)/py_agc_api.cpp \
	$(AGC_CXX_DIR)/lib-cxx.o \
	$(AGC_CORE_DIR)/agc_basic.o \
	$(AGC_CORE_DIR)/agc_decompressor_lib.o \
	$(AGC_CORE_DIR)/archive.o \
	$(AGC_CORE_DIR)/collection.o \
	$(AGC_CORE_DIR)/lz_diff.o \
	$(AGC_CORE_DIR)/segment.o \
	$(AGC_LIBS_DIR)/$(LIB_ZSTD) \
	-I $(AGC_MAIN_DIR) \
	-I $(PYBIND11_LIB)/include \
	-I `python3 -c "import sysconfig;print(sysconfig.get_paths()['include'])"` \
	-o $@`python3-config --extension-suffix`

clean:
	-rm $(AGC_EXAMPLES_DIR)/*.o
	-rm $(AGC_APP_DIR)/*.o
	-rm $(AGC_CORE_DIR)/*.o
	-rm $(AGC_CXX_DIR)/*.o
	-rm agc
	-rm libagc.a
	-rm -f $(PY_AGC_API_DIR)/*.o
	-rm *.so

