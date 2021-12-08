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

Run 81 set of tests according to the following parameters:

| Parameter | Values |
|---|---|
| im | `128`, `256`, `512` |
| jm | `128`, `256`, `512` |
| km | `64`, `128`, `256` |
| node_count | `2`, `8`, `16` |

For each `node_count`, there are more than one pair of `dsmNX` and `dsmNY`:

| Node Count | (dsmNX, dsmNY) |
|---|---|
| 2 | `(1, 2)`, `(2, 1)` |
| 8 | `(1, 8)`, `(2, 4)`, `(4, 2)`, `(8, 1)` |
| 16 | `(1, 16)`, `(2, 8)`, `(4, 4)`, `(8, 2)`, `(16, 1)` |

Considering all possible parameters, we will launch
999 (`27 * 3 * 1 + 27 * 1 * 1 + 27 * (2 + 4 + 5) * 3 = 999`) single tests.

```bash
nohup python testing.py >> nohup.out 2>&1 &
tail -f nohup.out
```

Additionally, whether using the compiler `-O3` optimization flag may also affect the performance. So there would be actually `1998 = 999 * 2` single tests to run.

# Results

## Manifest

The test result can be obtained from the `misc` directory.

There are 5 files within `misc` directory:

- `cvm-info.txt` is the test environment on Google Cloud Platform.
- `results.tar.gz` is the original test logs collected from `main` branch.
- `results-flags.tar.gz` is the original test logs collected from `flags` branch.
- `test-results.csv` is the test results in CSV format, containing 1998 records.
- `test-results.json` is the test results in JSON format, containing 1998 records.

## How to read the results

For each test record, `test-results.json` stores the SOR times from all computing nodes, while `test-results.csv` only stores the worst SOR time among all computing nodes.

The explanation for `test-results.csv` headers is as follows:

- `im`: im from array size.
- `jm`: jm from array size.
- `km`: km from array size.
- `method`: numbers from `0` to `4` represents `Reference C Code`, `ArgoDSM V1 Code`, `ArgoDSM V2 Code`, `ArgoDSM V3 Code`, `MPI Code`, respectively.
- `optimized`: whether the `-O3` optimization flag has been set within `CMAKE_CXX_FLAGS`. `0` means `Not Set`, and `1` means `Set`.
- `node_count`: number of computing nodes.
- `dsm_nx`: `dsmNX` parameter for ArgoDSM.
- `dsm_ny`: `dsmNY` parameter for ArgoDSM.
- `time`: SOR time in seconds.

It is highly recommended reading the results via `Python`:

```python
>>> import pandas as pd
>>> df = pd.read_csv('test-results.csv')
>>> df
       im   jm   km  method  optimized  node_count  dsm_nx  dsm_ny      time
0     128  128   64       0          0           1       0       0  0.412469
1     128  128   64       0          1           1       0       0  0.079119
2     128  128   64       1          0           2       1       2  7.100646
3     128  128   64       1          0           2       2       1  6.763574
4     128  128   64       1          0           8       1       8  0.426822
...   ...  ...  ...     ...        ...         ...     ...     ...       ...
1993  512  512  256       4          0           8       0       0  6.401710
1994  512  512  256       4          0          16       0       0  4.882370
1995  512  512  256       4          1           2       0       0  8.194650
1996  512  512  256       4          1           8       0       0  1.832680
1997  512  512  256       4          1          16       0       0  1.667990

[1998 rows x 9 columns]
```

## How to regenerate CSV and JSON result

To regenerate `test-results.csv` and `test-results.json` from original tset logs, follow the next steps:

```bash
cd misc

tar -xzf results.tar.gz
tar -xzf results-flags.tar.gz

# Golang is required.
go run ../collector/main.go
```
