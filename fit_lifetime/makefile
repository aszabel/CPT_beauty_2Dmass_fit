#compiler definition
#CC=clang++
CC=g++

#flags declaration
CFLAGS= -std=c++17 `root-config --cflags --libs` -fopenmp -Wall -O3 -march=haswell -fPIE -lMathMore -I include -I ../common/include
CFLAGS_COMP= -std=c++17 `root-config --cflags` -fopenmp -Wall -O3 -march=haswell -I include -I ../common/include
#MPICFLAGS= -fopenmp  -std=c++20 --include-directory /mnt/opt/spack-0.17/opt/spack/linux-centos7-ivybridge/gcc-8.3.0/root-6.24.06-f2vmdt3rp2x3mqemsppuq6tvxnztd3sf/include

SRCS = src/fit_lifetime.cpp src/pdfs_cpt_model.cpp src/config.cpp
COMMON = ../common/src/D_M_fit_shape.cpp ../common/src/M_B_2missPT_fit.cpp ../common/src/ChebyshevPDF.cpp ../common/src/sweights.cpp
COMMON_NAMES = $(notdir ${COMMON})
OBJS = $(SRCS:.cpp=.o)
OBJS_COMMON = $(COMMON_NAMES:%.cpp=common/%.o)

TARGET = fit2time

all: $(TARGET)

# Linking
$(TARGET): $(OBJS) $(OBJS_COMMON)
	$(CC) $(CFLAGS)  -o $@ $^

# Compilation
src/%.o: src/%.cpp
	$(CC) $(CFLAGS_COMP) -c $< -o $@

common/%.o: ../common/src/%.cpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS_COMP) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) $(OBJS_COMMON) $(TARGET)

# Phony targets
.PHONY: all clean


