INC=-I/usr/local/cuda/include
NVCC=/usr/local/cuda/bin/nvcc
NVCC_OPT=-std=c++11

all:
    $(NVCC) $(NVCC_OPT) gpu-example.cpp -o gpu-example

clean:
    -rm -f gpu-example