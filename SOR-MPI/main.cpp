#include <iostream>
#include <mpi.h>
#include "array_index_f2c1d.h"
#include "sor_params.h"
#include <unistd.h>

void divide_params(int world_rank, int world_size, int *base_k_size, int *base_k_rem, int *actual_k_size) {
    *base_k_size = (km + 2) / world_size;
    *base_k_rem = (km + 2) % world_size;
    if (world_rank == 0) {
        *base_k_size += *base_k_rem;
    }
    *actual_k_size = *base_k_size + 2;
}

void size_params(int base_k_size, int actual_k_size,
                 int64_t *base_2d_size, int64_t *base_3d_size, int64_t *actual_3d_size) {
    *base_2d_size = (int64_t) (im + 2) * (int64_t) (jm + 2);
    *base_3d_size = *base_2d_size * (int64_t) base_k_size;
    *actual_3d_size = *base_2d_size * (int64_t) actual_k_size;
}

struct Context {
    int world_size;
    int world_rank;
    int base_k_size;
    int base_k_rem;
    int actual_k_size;
    int64_t base_2d_size;
    int64_t base_3d_size;
    int64_t actual_3d_size;
};

void init_context(Context *ctx) {
    MPI_Comm_size(MPI_COMM_WORLD, &ctx->world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &ctx->world_rank);
    divide_params(ctx->world_rank, ctx->world_size, &ctx->base_k_size, &ctx->base_k_rem, &ctx->actual_k_size);
    size_params(ctx->base_k_size, ctx->actual_k_size, &ctx->base_2d_size, &ctx->base_3d_size, &ctx->actual_3d_size);
}

void init_array(Context *ctx, float **p0_ptr, float **p1_ptr, float **rhs_ptr) {
    int actual_k_size = ctx->actual_k_size;
    int64_t actual_3d_size = ctx->actual_3d_size;
    *p0_ptr = new float[actual_3d_size];
    *p1_ptr = new float[actual_3d_size];
    *rhs_ptr = new float[actual_3d_size];
    // initialize
    for (int i = 0; i < im + 2; i += 1) {
        for (int j = 0; j < jm + 2; j += 1) {
            for (int k = 0; k < actual_k_size; k += 1) {
                (*rhs_ptr)[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)] = 1.0;
                (*p0_ptr)[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)] = 1.0;
                (*p1_ptr)[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)] = 0.0;
            }
        }
    }
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

void sor(Context *ctx, float *p0, float *p1, const float *rhs) {
    // define variables
    int world_rank = ctx->world_rank;
    int world_size = ctx->world_size;
    int actual_k_size = ctx->actual_k_size;
    int64_t base_2d_size = ctx->base_2d_size;

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

    // transfer p0 boundary data
    int tag_up = 1023, tag_down = 1024;
    MPI_Status status;
    MPI_Sendrecv(p0 + base_2d_size * (actual_k_size - 2), int(base_2d_size), MPI_FLOAT, up_neighbor, tag_up,
                 p0, int(base_2d_size), MPI_FLOAT, down_neighbor, tag_up,
                 MPI_COMM_WORLD, &status);
    MPI_Sendrecv(p0 + base_2d_size, int(base_2d_size), MPI_FLOAT, down_neighbor, tag_down,
                 p0 + base_2d_size * (actual_k_size - 1), int(base_2d_size), MPI_FLOAT, up_neighbor, tag_down,
                 MPI_COMM_WORLD, &status);

    // computation
    float rel_tmp;
    for (int i = 0; i < im + 2; i += 1) {
        for (int j = 0; j < jm + 2; j += 1) {
            for (int k = 1; k < actual_k_size - 1; k += 1) {
                if (i == im + 1) {
                    p1[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)] =
                            p0[F3D2C(im + 2, jm + 2, 0, 0, 0, i - im, j, k)];
                } else if (i == 0) {
                    p1[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)] =
                            p0[F3D2C(im + 2, jm + 2, 0, 0, 0, i + im, j, k)];
                } else if (j == jm + 1) {
                    p1[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)] =
                            p0[F3D2C(im + 2, jm + 2, 0, 0, 0, i - 1, j, k)];
                } else if (j == 0) {
                    //  We keep the original values
                    p1[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)] =
                            p0[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)];
                } else if ((world_size == 1 && k > 1 && k < actual_k_size - 2) ||
                           (world_size > 1 &&
                            ((world_rank > 0 && world_rank < world_size - 1) || (world_rank == 0 && k > 1) ||
                             (world_rank == world_size - 1 && k < actual_k_size - 2)))) {
                    //  The actual SOR expression
                    rel_tmp = omega * (cn1 * (
                            cn2l * p0[F3D2C(im + 2, jm + 2, 0, 0, 0, i + 1, j, k)] +
                            cn2s * p0[F3D2C(im + 2, jm + 2, 0, 0, 0, i - 1, j, k)] +
                            cn3l * p0[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j + 1, k)] +
                            cn3s * p0[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j - 1, k)] +
                            cn4l * p0[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k + 1)] +
                            cn4s * p0[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k - 1)] -
                            rhs[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)]) -
                                       p0[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)]);
                    p1[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)] = p0[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)] + rel_tmp;
                }
            }
        }
    }
}

void divide_location(Context *ctx, int k, int *target_world_rank, int *target_k) {
    k = k - ctx->base_k_rem;
    *target_world_rank = k / ctx->base_k_size;
    *target_k = k - (*target_world_rank) * (ctx->base_k_size);
    if (*target_world_rank == 0) {
        *target_k += ctx->base_k_rem;
    }
}

void print_number_from_array(Context *ctx, float *arr, int i, int j, int k) {
    int target_world_rank, target_k;
    divide_location(ctx, k, &target_world_rank, &target_k);
    if (ctx->world_rank == target_world_rank) {
        std::cout << "" << arr[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, target_k + 1)] << "\n";
    }
}

__attribute__((unused)) void print_array(float *arr, int actual_k_size) {
    std::cout << "print array start" << std::endl;
    for (int k = 1; k < actual_k_size - 1; k++) {
        for (int j = 0; j < jm + 2; j++) {
            for (int i = 0; i < im + 2; i++) {
                std::cout << arr[F3D2C(im + 2, jm + 2, 0, 0, 0, i, j, k)] << std::endl;
            }
        }
    }
    std::cout << "print array end" << std::endl;
}

int main(__attribute__((unused)) int argc, __attribute__((unused)) char **argv) {
    MPI_Init(nullptr, nullptr);

    clock_t total_start = clock();

    Context ctx{};
    init_context(&ctx);

    float *p0, *p1, *rhs = nullptr;
    init_array(&ctx, &p0, &p1, &rhs);

    const int num_iterations = 5;
    for (int iter = 1; iter <= num_iterations; iter += 1) {
        printf("num_iteration: %d / %d (world_rank = %d, world_size = %d)\n", iter, num_iterations,
               ctx.world_rank, ctx.world_size);
        if (iter % 2 == 0) {
            sor(&ctx, p1, p0, rhs);
        } else {
            sor(&ctx, p0, p1, rhs);
        }
    }

    clock_t total_end = clock();
    double total_time = (double) (total_end - total_start) / CLOCKS_PER_SEC;
    print_number_from_array(&ctx, p0, im / 2, jm / 2, km / 2);
    std::cout << "Total: " << total_time << std::endl;
//    sleep(world_rank+2);
//    print_array(p0, actual_k_size);
    delete[] p0;
    delete[] p1;
    delete[] rhs;

    MPI_Finalize();
}
