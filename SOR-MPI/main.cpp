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
    for (int k = 1; k < actual_k_size - 1; k += 1) {
        for (int i = 0; i < im + 2; i += 1) {
            for (int j = 0; j < jm + 2; j += 1) {
                (*rhs_ptr)[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)] = 1.0;
                (*p0_ptr)[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)] = 1.0;
            }
        }
    }
    printf("Allocated slice, actual_k_size = %d, actual_3d_size = %lld\n", actual_k_size, actual_3d_size);
}

void sor(float *p0, float *p1, float *rhs, int world_rand, int world_size) {

}

int main(int argc, char **argv) {
    MPI_Init(NULL, NULL);

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

    MPI_Finalize();
}
