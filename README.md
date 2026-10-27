# 2D Sandpile MPI Implementation

## Compilation

```bash
make -f Makefile
```

## Usage

```bash
# Default size (L=512), centered offset
mpirun -np <P> ./sandpile <offset>

# Custom size and offset
mpirun -np <P> ./sandpile <size> <offset_x> <offset_y>
```

A batch job script `sandpile.job`.
For running multiple jobs `sandpile3.job` uses 3 nodes.

To increase the number of iteration change the `MAX_LOOPS` macro in sandpile_init.c

## Program Structure

| File | Description |
|------|-------------|
| `sandpile.c` | Main driver: initialisation, simulation loop, output |
| `sandpile_init.c` | Input parsing, Cartesian topology setup, domain decomposition |
| `sandpile_simulation.c` | Single simulation step with halo exchange |
| `sandpile_io.c` | Output writing (.ppm visualisation, .mat final state) |
| `sandpile_context.h` | Data structures and function declarations |
| `arraymalloc.h` | Dynamic array allocation utilities |
| `arraymalloc.c` | Dynamic array allocation implementation |

## Output Files

| File | Description |
|------|-------------|
| `.mat` | Final state of sandpile after simulation |
| `.ppm` | Visualisation of the final sandpile |

## Testing

```bash
# Compile test utility
make -f Makefile_test

# Compare two output matrices
./sandpile_test <file_1> <file_2>
```
