#compiler definition
#CC=clang++
CC=g++

#flags declaration
CFLAGS= -std=c++14 `root-config --cflags --libs` -fopenmp -Wall -O3 -march=haswell -fPIE -lMathMore
CFLAGS_COMP= `root-config --cflags` -fopenmp -Wall -O3 -march=haswell 
#MPICFLAGS= -fopenmp  -std=c++20 --include-directory /mnt/opt/spack-0.17/opt/spack/linux-centos7-ivybridge/gcc-8.3.0/root-6.24.06-f2vmdt3rp2x3mqemsppuq6tvxnztd3sf/include

SRCS = fit_lifetime.cpp pdfs_cpt_model.cpp  ../multithreaded/D_M_fit_shape.cpp ../multithreaded/M_B_2missPT_fit.cpp ../multithreaded/ChebyshevPDF.cpp config_time.cpp sweights.cpp
OBJS = $(SRCS:.cpp=.o)

TARGET = fit2time

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


