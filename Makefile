# Compiler configuration
CC = mpicc
CFLAGS = -Wall -O3

# Target executable name
TARGET = sandpile

OBJS = sandpile.o sandpile_init.o sandpile_simulation.o sandpile_io.o arraymalloc.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET)

# --- Compilation Rules (Dependencies) ---

# 1. Main (sandpile.c)
# Depends on its source and the shared context header.
sandpile.o: sandpile.c sandpile_context.h
	$(CC) $(CFLAGS) -c sandpile.c

# 2. Initialization Module
# Depends on context, shared sandpile definitions, and its source.
sandpile_init.o: sandpile_init.c sandpile_context.h 
	$(CC) $(CFLAGS) -c sandpile_init.c

# 3. Simulation Module
# Depends on context, shared sandpile definitions, and its source.
sandpile_simulation.o: sandpile_simulation.c sandpile_context.h 
	$(CC) $(CFLAGS) -c sandpile_simulation.c

# 4. I/O and Cleanup Module (Replaces sandio.o)
# Depends on context, shared sandpile definitions, and its source.
sandpile_io.o: sandpile_io.c sandpile_context.h 
	$(CC) $(CFLAGS) -c sandpile_io.c

# 5. Utility Module
arraymalloc.o: arraymalloc.c arraymalloc.h
	$(CC) $(CFLAGS) -c arraymalloc.c 

# Cleanup rule
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
