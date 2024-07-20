# Compiler and flags
CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LDFLAGS = -lm

# Target executable
TARGET = symnmf

# Source files
SRCS = symnmf.c

# Header files
HDRS = symnmf.h

# Object files
OBJS = $(SRCS:.c=.o)

# Default target to build the executable
all: $(TARGET)

# Build target executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compile source files to object files
%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up build artifacts
clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: all clean
