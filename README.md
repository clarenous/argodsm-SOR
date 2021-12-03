# Comparison of the acceleration of SOR algorithm

Comparison of the acceleration of Successive Over-Relaxation algorithm using `C`, `MPI` and `ArgoDSM`.

The implementations are located in the following directories:

| Source Code Path | Binary Name |
|---|---|
| SOR-C-reference | `sor_c` |
| SOR-MPI | `sor_mpi` |
| SOR-DSM-1 | `sor_argo_1` |
| SOR-DSM-2 | `sor_argo_2` |
| SOR-DSM-3 | `sor_argo_3` |

## Parameters

For each implementation, modify the corresponding `sor_params.h` file to build & test for different parameters.

## Build

To build all implementations, run `./build.sh`.

The binaries (`sor_c`, `sor_mpi`, `sor_argo_1`, `sor_argo_2`, `sor_argo_3`) will be generated in `./bin/`.

**Or**, to build one certain implementation, follow the next steps:

```
# IMPL_PATH indicates the implementation source code path stated above.

mkdir -p build/$IMPL_PATH && cd build/$IMPL_PATH

cmake ../../$IMPL_PATH

# The corresponding binary will be generated.
make
```

# Run

Run each implementation with `mpirun -n <node-count>`.

# Test

Run 900 set of tests according to the following parameters:

| Parameter | Values |
|---|---|
| im | `16`, `32`, `64`, `128`, `256`, `512` |
| jm | `16`, `32`, `64`, `128`, `256`, `512` |
| km | `16`, `32`, `64`, `128`, `256` |
| node_count | `1`, `2`, `4`, `8`, `16` |

For each `node_count`, there are more than one pair of `dsmNX` and `dsmNY`:

| Node Count | (dsmNX, dsmNY) |
|---|---|
| 1 | `(1, 1)` |
| 2 | `(1, 2)`, `(2, 1)` |
| 4 | `(1, 4)`, `(2, 2)`, `(4, 1)` |
| 8 | `(1, 8)`, `(2, 4)`, `(4, 2)`, `(8, 1)` |
| 16 | `(1, 16)`, `(2, 8)`, `(4, 4)`, `(8, 2)`, `(16, 1)` |

Considering all possible parameters, we will launch
9180 (`180 * 5 * 1 + 180 * 1 * 1 + 180 * (1 + 2 + 3 + 4 + 5) * 3 = 9180`) single tests.

```bash
nohup python testing.py >> nohup.out 2>&1 &
tail -f nohup.out
```
