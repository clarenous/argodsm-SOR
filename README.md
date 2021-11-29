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
