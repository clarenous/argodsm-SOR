#include <iostream>
#include <mpi.h>
#include "array_index_f2c1d.h"
#include "sor_params.h"

void divide_params(int world_rank, int world_size, int *base_k_size, int *actual_k_size) {
    *base_k_size = (km + 2) / world_size;
    int k_rem = (km + 2) % world_size;
    if (world_rank == 0) {
        *base_k_size += k_rem;
    }
    *actual_k_size = *base_k_size + 2;
}

void divide_location(int world_rank, int world_size, int k, int *target_world_rank, int *target_k) {
    int base_k_size = (km + 2) / world_size;
    *target_world_rank = k / base_k_size;
    if (*target_world_rank >= world_size) {
        *target_world_rank = world_size - 1;
    }
    *target_k = k - (*target_world_rank) * base_k_size;
}

void size_params(int base_k_size, int actual_k_size,
                 int64_t *base_2d_size, int64_t *base_3d_size, int64_t *actual_3d_size) {
    *base_2d_size = (int64_t) (im + 2) * (int64_t) (jm + 2);
    *base_3d_size = *base_2d_size * (int64_t) base_k_size;
    *actual_3d_size = *base_2d_size * (int64_t) actual_k_size;
}

void init_array(float **p0_ptr, float **p1_ptr, float **rhs_ptr, int actual_k_size, int64_t actual_3d_size) {
    *p0_ptr = new float[actual_3d_size];
    *p1_ptr = new float[actual_3d_size];
    *rhs_ptr = new float[actual_3d_size];
    // initialize
    for (int k = 0; k < actual_k_size; k += 1) {
        for (int j = 0; j < jm + 2; j += 1) {
            for (int i = 0; i < im + 2; i += 1) {
                (*rhs_ptr)[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i)] = 1.0;
                (*p0_ptr)[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i)] = 1.0;
            }
        }
    }
    printf("Allocated float slice, actual_k_size = %d, actual_3d_size = %lld\n", actual_k_size, actual_3d_size);
}

void find_neighbor_proc(int world_rank, int world_size, int *left, int *right) {
    if (world_rank > 0) {
        *left = world_rank - 1;
    } else {
        *left = MPI_PROC_NULL;
    }
    if (world_rank < world_size - 1) {
        *right = world_rank + 1;
    } else {
        *right = MPI_PROC_NULL;
    }
}

void sor(float *p0, float *p1, float *rhs, int world_rank, int world_size, int actual_k_size, int64_t base_2d_size) {
    // define constants
    const float cn1 = 1.0 / 3.0;
    const float cn2l = 0.5;
    const float cn2s = 0.5;
    const float cn3l = 0.5;
    const float cn3s = 0.5;
    const float cn4l = 0.5;
    const float cn4s = 0.5;
    const float omega = 1.0;

    // find neighbor process
    int down_neighbor, up_neighbor;
    find_neighbor_proc(world_rank, world_size, &down_neighbor, &up_neighbor);
//    printf("find neighbor process: world_rank = %d, world_size = %d, down = %d, up = %d\n", world_rank, world_size, down_neighbor, up_neighbor);

    // transfer p0 boundary data
//    printf("transfer p0 boundary data: world_rank = %d, world_size = %d\n", world_rank, world_size);
    int tag_up = 1023, tag_down = 1024;
    MPI_Status status;
    MPI_Sendrecv(p0 + base_2d_size * (actual_k_size - 2), int(base_2d_size), MPI_FLOAT, up_neighbor, tag_up,
                 p0, int(base_2d_size), MPI_FLOAT, down_neighbor, tag_up,
                 MPI_COMM_WORLD, &status);
    MPI_Sendrecv(p0 + base_2d_size, int(base_2d_size), MPI_FLOAT, down_neighbor, tag_down,
                 p0 + base_2d_size * (actual_k_size - 1), int(base_2d_size), MPI_FLOAT, up_neighbor, tag_down,
                 MPI_COMM_WORLD, &status);

    // iterations
//    printf("computation: world_rank = %d, world_size = %d\n", world_rank, world_size);
    float rel_tmp;
    for (int k = 1; k < actual_k_size - 1; k += 1) {
        for (int j = 0; j < jm + 2; j += 1) {
            for (int i = 0; i < im + 2; i += 1) {
                if (i == im + 1) {
                    p1[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i)] =
                            p0[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i - im)];
                } else if (i == 0) {
                    p1[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i)] =
                            p0[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i + im)];
                } else if (j == jm + 1) {
                    p1[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i)] =
                            p0[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i - 1)];
                } else if (j == 0) {
                    //  We keep the original values
                    p1[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i)] =
                            p0[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i)];
                } else if ((world_rank > 0 && world_rank < world_size - 1) || (world_rank == 0 && k > 1) ||
                           (world_rank == world_size - 1 && k < actual_k_size - 2)) {
                    //  The actual SOR expression
                    rel_tmp = omega * (cn1 * (
                            cn2l * p0[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i + 1)] +
                            cn2s * p0[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i - 1)] +
                            cn3l * p0[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j + 1, i)] +
                            cn3s * p0[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j - 1, i)] +
                            cn4l * p0[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k + 1, j, i)] +
                            cn4s * p0[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k - 1, j, i)] -
                            rhs[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i)]) -
                                       p0[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i)]);
                    p1[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i)] =
                            p0[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, k, j, i)] + rel_tmp;
                }
            }
        }
    }
}

void print_number_from_array(float *arr, int world_rank, int world_size, int actual_k_size, int k, int j, int i) {
    int target_world_rank, target_k;
    divide_location(world_rank, world_size, k, &target_world_rank, &target_k);
    if (world_rank == target_world_rank) {
        std::cout << "" << arr[F3D2C_kji(actual_k_size, jm + 2, 0, 0, 0, target_k + 1, j, i)] << "\n";
    }
}

int main(int argc, char **argv) {
    MPI_Init(NULL, NULL);

    clock_t total_start = clock();

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int base_k_size, actual_k_size;
    divide_params(world_rank, world_size, &base_k_size, &actual_k_size);

    int64_t base_2d_size, base_3d_size, actual_3d_size;
    size_params(base_k_size, actual_k_size, &base_2d_size, &base_3d_size, &actual_3d_size);

    float *p0, *p1, *rhs = nullptr;
    init_array(&p0, &p1, &rhs, actual_k_size, actual_3d_size);

    const int num_iterations = 5;
    for (int iter = 1; iter <= num_iterations; iter += 1) {
        printf("num_iteration: %d / %d (world_rank = %d, world_size = %d)\n", iter, num_iterations,
               world_rank, world_size);
        if (iter % 2 == 0) {
            sor(p1, p0, rhs, world_rank, world_size, actual_k_size, base_2d_size);
        } else {
            sor(p0, p1, rhs, world_rank, world_size, actual_k_size, base_2d_size);
        }
    }

    clock_t total_end = clock();
    double total_time = (double) (total_end - total_start) / CLOCKS_PER_SEC;
    print_number_from_array(p0, world_rank, world_size, actual_k_size, km / 2, jm / 2, im / 2);
    std::cout << "Total: " << total_time << std::endl;
    delete[] p0;
    delete[] p1;
    delete[] rhs;

    MPI_Finalize();
}
