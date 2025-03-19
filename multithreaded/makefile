#compiler definition
#CC=clang++
CC=g++

#flags declaration
CFLAGS= -std=c++17 `root-config --cflags --libs` -fopenmp -Wall -O3 -march=haswell -I ../common/include -I include 
CFLAGS_COMP= -std=c++17 `root-config --cflags` -fopenmp -Wall -O3 -march=haswell -I ../common/include -I include
#MPICFLAGS= -fopenmp  -std=c++20 --include-directory /mnt/opt/spack-0.17/opt/spack/linux-centos7-ivybridge/gcc-8.3.0/root-6.24.06-f2vmdt3rp2x3mqemsppuq6tvxnztd3sf/include

SRCS1 = src/Dfit_2D_main.cpp 
SRCS2 = src/Dfit_1D_sidebands.cpp 
COMMON = ../common/src/D_M_fit_shape.cpp ../common/src/M_B_2missPT_fit.cpp ../common/src/ChebyshevPDF.cpp ../common/src/sweights.cpp
CONFIG = src/config.cpp
OBJS1 = $(SRCS1:.cpp=.o)
OBJS2 = $(SRCS2:.cpp=.o)
OBJS_CONFIG = $(CONFIG:.cpp=.o)
OBJS_COMMON = $(COMMON:.cpp=.o)

TARGET1 = fit2D_mass 
TARGET2 = fit1D_mass

all: $(TARGET1) $(TARGET2)

# Linking
$(TARGET1): $(OBJS1) $(OBJS_COMMON) $(CONFIG)
	$(CC) $(CFLAGS)  -o $@ $^

$(TARGET2): $(OBJS2) $(OBJS_COMMON) $(CONFIG)
	$(CC) $(CFLAGS)  -o $@ $^

# Compilation
%.o: %.cpp
	$(CC) $(CFLAGS_COMP) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS1) $(OBJS2) $(OBJS_COMMON) $(OBJS_CONFIG) $(TARGET1) $(TARGET2)

# Phony targets
.PHONY: all clean


