#
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# Compiler definition
CC = g++

# Compiler flags
CFLAGS = -Wall -O3 -fopenmp -std=c++11 -pthread

# External libraries
LIBS = -lm -lfftw3 -lfftw3_omp

# Source files
SRCS = main.cpp utils.cpp evolution.cpp potential.cpp sim.cpp
OBJS = $(SRCS:.c=.o)

# Executable file
MAIN = tfm

all: $(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(MAIN) $(OBJS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN) Data/PhiG/*.dat Data/PhiE/*.dat Data/PhiT/*.dat \
		 Images/PhiG/*.png Images/PhiE/*.png Images/PhiT/*.png
