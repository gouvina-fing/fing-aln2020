CC = gcc

# Compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS = -g -O3

BIN_DIR = bin

all: iteration multiplication

iteration: benchmark_matrix_iteration.c
	$(CC) $(CFLAGS) -o $(BIN_DIR)/benchmark_matrix_iteration.x benchmark_matrix_iteration.c

multiplication:
	$(CC) $(CFLAGS) -o $(BIN_DIR)/benchmark_matrix_multiplication.x benchmark_matrix_multiplication.c

clean: 
	$(RM) *.o *.x