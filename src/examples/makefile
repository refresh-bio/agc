all: example-agc-lib-c example-agc-lib-cpp

AGC_ROOT_DIR = ../..
EXAMPLE_APP_DIR = .
AGC_CXX_DIR = ../lib-cxx
AGC_LIBS_DIR = ../../libs
LIBS_DIR = . #/usr/local/lib
INCLUDE_DIR= . #/usr/local/include

CXX = g++
CC 	= gcc
CFLAGS	= -Wall -O3 -m64 -std=c11 -pthread -mavx -I $(INCLUDE_DIR) -I $(AGC_CXX_DIR)
CXXFLAGS	= -Wall -O3 -m64 -std=c++17 -pthread -mavx -I $(INCLUDE_DIR) -I $(AGC_CXX_DIR) -fpermissive
CLINK	= -lm -pthread -std=c++17 

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
else                          # If uname not available => 'not' 
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
endif
ifeq ($(uname_S),Linux)
	CLINK+=-fabi-version=6
	LIB_ZSTD=libzstd.a
	LIB_AGC=libagc.a
endif

ifeq ($(uname_S),Darwin)
	LIB_ZSTD=libzstd.mac.a
	LIB_AGC=libagc.a
	CXX=g++-11	
endif
# default install location (binary placed in the /bin folder)
prefix      = /usr/local

# optional install location
exec_prefix = $(prefix)


%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

example-agc-lib-c: $(EXAMPLE_APP_DIR)/example-agc-lib-c.o
	$(CXX) -o $(EXAMPLE_APP_DIR)/$@  \
	$(EXAMPLE_APP_DIR)/example-agc-lib-c.o \
	$(AGC_ROOT_DIR)/$(LIB_AGC) \
	$(AGC_LIBS_DIR)/$(LIB_ZSTD) \
	$(CLINK)

example-agc-lib-cpp: $(EXAMPLE_APP_DIR)/example-agc-lib-cpp.o
	$(CXX) -o $(EXAMPLE_APP_DIR)/$@  \
	$(EXAMPLE_APP_DIR)/example-agc-lib-cpp.o \
	$(AGC_ROOT_DIR)/$(LIB_AGC) \
	$(AGC_LIBS_DIR)/$(LIB_ZSTD) \
	$(CLINK)

clean:
	-rm $(EXAMPLE_APP_DIR)/*.o
	-rm example-agc-lib-c
	-rm example-agc-lib-cpp
