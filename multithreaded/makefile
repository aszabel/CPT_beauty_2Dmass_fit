#compiler definition
#CC=clang++
CC=g++

#flags declaration
CFLAGS= `root-config --cflags --libs` -fopenmp -Wall -O3 -march=haswell
CFLAGS_COMP= `root-config --cflags` -fopenmp -Wall -O3 -march=haswell
#MPICFLAGS= -fopenmp  -std=c++20 --include-directory /mnt/opt/spack-0.17/opt/spack/linux-centos7-ivybridge/gcc-8.3.0/root-6.24.06-f2vmdt3rp2x3mqemsppuq6tvxnztd3sf/include

SRCS = Dfit_2D_main.cpp D_M_fit_shape.cpp M_B_2missPT_fit.cpp ChebyshevPDF.cpp
OBJS = $(SRCS:.cpp=.o)

TARGET = fit2D_mass

all: $(TARGET)


# Linking
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS)  -o $@ $^

# Compilation
%.o: %.cpp
	$(CC) $(CFLAGS_COMP) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean


