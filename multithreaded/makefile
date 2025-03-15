#compiler definition
#CC=clang++
CC=g++

#flags declaration
CFLAGS= `root-config --cflags --libs` -fopenmp -Wall -O3 -march=haswell
CFLAGS_COMP= `root-config --cflags` -fopenmp -Wall -O3 -march=haswell
#MPICFLAGS= -fopenmp  -std=c++20 --include-directory /mnt/opt/spack-0.17/opt/spack/linux-centos7-ivybridge/gcc-8.3.0/root-6.24.06-f2vmdt3rp2x3mqemsppuq6tvxnztd3sf/include

SRCS1 = Dfit_2D_main.cpp 
SRCS2 = Dfit_1D_sidebands.cpp 
COMMON =  D_M_fit_shape.cpp M_B_2missPT_fit.cpp ChebyshevPDF.cpp config.cpp sweights.cpp
OBJS1 = $(SRCS1:.cpp=.o) $(COMMON)
OBJS2 = $(SRCS2:.cpp=.o) $(COMMON)

TARGET1 = fit2D_mass 
TARGET2 = fit1D_mass

all: $(TARGET1) $(TARGET2)


# Linking
$(TARGET1): $(OBJS1)
	$(CC) $(CFLAGS)  -o $@ $^

$(TARGET2): $(OBJS2)
	$(CC) $(CFLAGS)  -o $@ $^

# Compilation
%.o: %.cpp
	$(CC) $(CFLAGS_COMP) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) $(TARGET1) $(TARGET2)

# Phony targets
.PHONY: all clean


