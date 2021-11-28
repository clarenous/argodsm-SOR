#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include "array_index_f2c1d.h"
#include "../src/dsm-api.h"
#include "sor_params.h"

#ifndef _SOR_DSM_1
#define _SOR_DSM_1
#endif

typedef struct par_params {
    int i_start;
    int i_stop;
    int j_start;
    int j_stop;
} par_params;

void sor(float *p0,float *p1,float *rhs);

void sor(DSMArrayPtr p0, DSMArrayPtr p1, DSMArrayPtr rhs) {

    int i;
    int j;
    int k;
    const float cn1 = 1.0 / 3.0;
    const float cn2l = 0.5;
    const float cn2s = 0.5;
    const float cn3l = 0.5;
    const float cn3s = 0.5;
    const float cn4l = 0.5;
    const float cn4s = 0.5;
    const float omega = 1.0;
    float reltmp;
    for (i = 0;i <= im+1;i += 1) {
        for (j = 0;j <= jm+1;j += 1) {
            for (k = 0;k <= km+1;k += 1) {
                //std::cout << "("<<i <<","<< j <<","<< k << ");"<<F3D2C(im+2,jm+2,0,0,0,i,j,k)<<"\n";
                //  assume i=x =  west to east, y=j = south to north, k=z = vertical
                if (i==im+1) {
                    //std::cout <<"i==im+1\n";
                //  circular
                //  i=im+1
                    int64_t idx_p1 = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    int64_t idx_p0 = F3D2C(im+2,jm+2,0,0,0,i-im,j,k);
                    float val;
                    dsmreadfloatarrayf_(p0, &idx_p0, &val);
                    dsmwritefloatarrayf_(p1, &idx_p1, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i-im,j,k)];
                } else if (i==0) {
                    //std::cout <<"i==0\n";
                //  i=0
                //  circular
                    int64_t idx_p1 = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    int64_t idx_p0 = F3D2C(im+2,jm+2,0,0,0,i+im,j,k);
                    float val;
                    dsmreadfloatarrayf_(p0, &idx_p0, &val);
                    dsmwritefloatarrayf_(p1, &idx_p1, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i+im,j,k)];
                } else if (j==jm+1) {
                    //std::cout <<"j==jm+1\n";
                //  open
                //  j = jm+1
                    int64_t idx_p1 = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    int64_t idx_p0 = F3D2C(im+2,jm+2,0,0,0,i-1,j,k);
                    float val;
                    dsmreadfloatarrayf_(p0, &idx_p0, &val);
                    dsmwritefloatarrayf_(p1, &idx_p1, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i-1,j,k)];
                } else if (j==0) {
                    //std::cout <<"j==0\n";
                //  fixed
                //  j = 0
                //  We keep the original values
                    int64_t idx_p1 = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    int64_t idx_p0 = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    float val;
                    dsmreadfloatarrayf_(p0, &idx_p0, &val);
                    dsmwritefloatarrayf_(p1, &idx_p1, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i,j,k)];
                 } else if (i>0 && i<im+1 && j>0 && j<jm+1 && k>0 && k<km+1) {
                    //std::cout <<"core\n";
                    //  the core
                    //  The actual SOR expression
                    int64_t idx_cn2l = F3D2C(im+2,jm+2,0,0,0,i+1,j,k);
                    float val_cn2l = 0;
                    dsmreadfloatarrayf_(p0, &idx_cn2l, &val_cn2l);
                    int64_t idx_cn2s = F3D2C(im+2,jm+2,0,0,0,i-1,j,k);
                    float val_cn2s = 0;
                    dsmreadfloatarrayf_(p0, &idx_cn2s, &val_cn2s);
                    int64_t idx_cn3l = F3D2C(im+2,jm+2,0,0,0,i,j+1,k);
                    float val_cn3l = 0;
                    dsmreadfloatarrayf_(p0, &idx_cn3l, &val_cn3l);
                    int64_t idx_cn3s = F3D2C(im+2,jm+2,0,0,0,i,j-1,k);
                    float val_cn3s = 0;
                    dsmreadfloatarrayf_(p0, &idx_cn3s, &val_cn3s);
                    int64_t idx_cn4l = F3D2C(im+2,jm+2,0,0,0,i,j,k+1);
                    float val_cn4l = 0;
                    dsmreadfloatarrayf_(p0, &idx_cn4l, &val_cn4l);
                    int64_t idx_cn4s = F3D2C(im+2,jm+2,0,0,0,i,j,k-1);
                    float val_cn4s = 0;
                    dsmreadfloatarrayf_(p0, &idx_cn4s, &val_cn4s);
                    int64_t idx_rhs = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    float val_rhs = 0;
                    dsmreadfloatarrayf_(rhs, &idx_rhs, &val_rhs);
                    int64_t idx_core = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    float val_p0_core = 0;
                    dsmreadfloatarrayf_(p0, &idx_core, &val_p0_core);

                    reltmp = omega*(cn1*(
                                cn2l*val_cn2l+
                                cn2s*val_cn2s+
                                cn3l*val_cn3l+
                                cn3s*val_cn3s+
                                cn4l*val_cn4l+
                                cn4s*val_cn4s-
                                val_rhs)-
                            val_p0_core);
                    float val_new = val_p0_core+reltmp;
                    dsmwritefloatarrayf_(p1, &idx_core, &val_new);
                }
             }
         }
    }
}
    
void sor_par(DSMArrayPtr p0, DSMArrayPtr p1, DSMArrayPtr rhs, par_params params) {
    int i;
    int j;
    int k;

    // printf("h: %d\n", h);
    // printf("i_span: %d, i_start: %d, i_stop: %d\n", i_span, i_start, i_stop);
    // printf("j_span: %d, j_start: %d, j_stop: %d\n", j_span, j_start, j_stop);

    const float cn1 = 1.0 / 3.0;
    const float cn2l = 0.5;
    const float cn2s = 0.5;
    const float cn3l = 0.5;
    const float cn3s = 0.5;
    const float cn4l = 0.5;
    const float cn4s = 0.5;
    const float omega = 1.0;
    float reltmp;
    for (i = params.i_start; i <= params.i_stop; i += 1) {
        for (j = params.j_start;j <= params.j_stop;j += 1) {
            for (k = 0;k <= km+1;k += 1) {
                //std::cout << "("<<i <<","<< j <<","<< k << ");"<<F3D2C(im+2,jm+2,0,0,0,i,j,k)<<"\n";
                //  assume i=x =  west to east, y=j = south to north, k=z = vertical
                if (i==im+1) {
                    //std::cout <<"i==im+1\n";
                //  circular
                //  i=im+1
                    int64_t idx_p1 = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    int64_t idx_p0 = F3D2C(im+2,jm+2,0,0,0,i-im,j,k);
                    float val;
                    dsmreadfloatarrayf_(p0, &idx_p0, &val);
                    dsmwritefloatarrayf_(p1, &idx_p1, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i-im,j,k)];
                } else if (i==0) {
                    //std::cout <<"i==0\n";
                //  i=0
                //  circular
                    int64_t idx_p1 = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    int64_t idx_p0 = F3D2C(im+2,jm+2,0,0,0,i+im,j,k);
                    float val;
                    dsmreadfloatarrayf_(p0, &idx_p0, &val);
                    dsmwritefloatarrayf_(p1, &idx_p1, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i+im,j,k)];
                } else if (j==jm+1) {
                    //std::cout <<"j==jm+1\n";
                //  open
                //  j = jm+1
                    int64_t idx_p1 = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    int64_t idx_p0 = F3D2C(im+2,jm+2,0,0,0,i-1,j,k);
                    float val;
                    dsmreadfloatarrayf_(p0, &idx_p0, &val);
                    dsmwritefloatarrayf_(p1, &idx_p1, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i-1,j,k)];
                } else if (j==0) {
                    //std::cout <<"j==0\n";
                //  fixed
                //  j = 0
                //  We keep the original values
                    int64_t idx_p1 = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    int64_t idx_p0 = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    float val;
                    dsmreadfloatarrayf_(p0, &idx_p0, &val);
                    dsmwritefloatarrayf_(p1, &idx_p1, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i,j,k)];
                 } else if (i>0 && i<im+1 && j>0 && j<jm+1 && k>0 && k<km+1) {
                    //std::cout <<"core\n";
                    //  the core
                    //  The actual SOR expression
                    int64_t idx_cn2l = F3D2C(im+2,jm+2,0,0,0,i+1,j,k);
                    float val_cn2l = 0;
                    dsmreadfloatarrayf_(p0, &idx_cn2l, &val_cn2l);
                    int64_t idx_cn2s = F3D2C(im+2,jm+2,0,0,0,i-1,j,k);
                    float val_cn2s = 0;
                    dsmreadfloatarrayf_(p0, &idx_cn2s, &val_cn2s);
                    int64_t idx_cn3l = F3D2C(im+2,jm+2,0,0,0,i,j+1,k);
                    float val_cn3l = 0;
                    dsmreadfloatarrayf_(p0, &idx_cn3l, &val_cn3l);
                    int64_t idx_cn3s = F3D2C(im+2,jm+2,0,0,0,i,j-1,k);
                    float val_cn3s = 0;
                    dsmreadfloatarrayf_(p0, &idx_cn3s, &val_cn3s);
                    int64_t idx_cn4l = F3D2C(im+2,jm+2,0,0,0,i,j,k+1);
                    float val_cn4l = 0;
                    dsmreadfloatarrayf_(p0, &idx_cn4l, &val_cn4l);
                    int64_t idx_cn4s = F3D2C(im+2,jm+2,0,0,0,i,j,k-1);
                    float val_cn4s = 0;
                    dsmreadfloatarrayf_(p0, &idx_cn4s, &val_cn4s);
                    int64_t idx_rhs = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    float val_rhs = 0;
                    dsmreadfloatarrayf_(rhs, &idx_rhs, &val_rhs);
                    int64_t idx_core = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    float val_p0_core = 0;
                    dsmreadfloatarrayf_(p0, &idx_core, &val_p0_core);

                    reltmp = omega*(cn1*(
                                cn2l*val_cn2l+
                                cn2s*val_cn2s+
                                cn3l*val_cn3l+
                                cn3s*val_cn3s+
                                cn4l*val_cn4l+
                                cn4s*val_cn4s-
                                val_rhs)-
                            val_p0_core);
                    float val_new = val_p0_core+reltmp;
                    dsmwritefloatarrayf_(p1, &idx_core, &val_new);
                }
             }
         }
    }
}

int main() {
    // Starting barrier test!
    clock_t total_start = clock();
    clock_t init_start = clock();

    assert(im % dsmNX == 0); // im must be divisible by dsmNX.
    assert(jm % dsmNY == 0); // jm must be divisible by dsmNY.
    
    // Size of a single matrix.
    int64_t mat_size = (int64_t)(im+2)*(int64_t)(jm+2)*(int64_t)(km+2);
    // Mumber of matrices.
    int64_t mat_num = 3;
    // Allocation in bytes, but we're allocating for floats & doubles.
    int64_t bytes_per_float = 4;
    int64_t bytes_per_double = 8;

    // Scale by number of nodes.
    int64_t node_mult = dsmNX*dsmNY;

    // Account for the array of pointers.
    int64_t pa_size = mat_num*bytes_per_double*node_mult;

    int64_t mem_size;
    int64_t cache_size;
    int64_t gb = 1024UL*1024UL*1024UL;
    if (mode == 1) {
        if (node_mult == 1) {
            mat_num++; // Avoid bad alloc by allocating extra memory.
        }
        mem_size = node_mult*bytes_per_float*mat_num*mat_size;
        mem_size += pa_size;
        cache_size = (mat_num*bytes_per_float*mat_size) + pa_size;
    } else if (mode == 2) {
        mat_num++; // Avoid bad alloc by allocating extra memory.
        node_mult = 1; // Static memory size.
        mem_size = node_mult*bytes_per_float*mat_num*mat_size;
        mem_size += pa_size;
        cache_size = (mat_num*bytes_per_float*mat_size) + pa_size;
    } else if (mode == 3) {
        mem_size = 20UL*gb;
        cache_size = mem_size / 2UL;
        mem_size = mem_size*node_mult;
    } else if (mode == 4) {
        mem_size = 27UL*gb;
        cache_size = mem_size / 8UL;
        mem_size = mem_size*node_mult;
    } else if (mode == 5) {
        mem_size = 30UL*gb;
        cache_size = mem_size / 32UL;
        mem_size = mem_size*node_mult;
    }
    // Fit the equivalent of 3 LOCAL matrices in the cache.
    // int64_t cache_size = (mat_num*bytes_per_float*mat_size) / (dsmNX * dsmNY);
    // Testing vals for memory:
    
    #ifdef VERBOSE
        printf("mem_size: %ld\n", mem_size);
        printf("cache_size: %ld\n", cache_size);
    #endif // VERBOSE

    dsminitf_(&mem_size, &cache_size);

    #ifdef VERBOSE
        printf("Memory Initialised\n");
    #endif // VERBOSE

    int num_nodes;
    dsmgetnnodesf_(&num_nodes);
    assert(num_nodes == dsmNX*dsmNY); // There must be an MPI node for each sub-array.

    int64_t node_id;
    dsmgetnidf_(&node_id);
    // printf("Number of nodes: %d. This node: %d\n", num_nodes, node_id);

    clock_t init_end = clock();
    double init_time = (double)(init_end - init_start) / CLOCKS_PER_SEC;

    clock_t alloc_start = clock();

    int64_t p0 = 0;
    int64_t p0_size = mat_size;
    dsmallocfloatarrayf_(&p0, &p0_size);
    // argo::barrier();

    int64_t p1 = 0;
    int64_t p1_size = mat_size;
    dsmallocfloatarrayf_(&p1, &p1_size);
    // argo::barrier();

    int64_t rhs = 0;
    int64_t rhs_size = mat_size;
    dsmallocfloatarrayf_(&rhs, &rhs_size);
    // argo::barrier();

    int ilh = 1;
    int jlh = 1;
    
    // int h = (int)sqrt(num_nodes);
    int i_span = im / dsmNX;
    int i_start;
    if (node_id % dsmNX == 0) {
        i_span += ilh;
        i_start = 0;
    } else if (node_id % dsmNX == dsmNX-1) {
        i_start = i_span*(dsmNX-1) + ilh;
        i_span += ilh;
    } else {
        i_start = i_span*(node_id % dsmNX) + ilh;
    }
    if (num_nodes == 1) {
        i_span += ilh;
    }
    int i_stop = i_start + i_span - 1;

    int j_span = jm / dsmNY;
    int j_start;
    if (0 <= node_id && node_id < dsmNX) {
        j_span += jlh;
        j_start = 0;
    } else if (dsmNX*(dsmNY-1) <= node_id && node_id < dsmNX*dsmNY) {
        j_start = j_span*(dsmNY-1) + jlh;
        j_span += jlh;
    } else {
        j_start = j_span*(node_id / dsmNX) + jlh;
    }
    if (num_nodes < 3) {
        j_span += jlh;
    }
    int j_stop = j_start + j_span - 1;
    
    par_params params = {
        i_start,
        i_stop,
        j_start,
        j_stop
    };
    argo::barrier();
    clock_t alloc_end = clock();
    double alloc_time = (double)(alloc_end - alloc_start) / CLOCKS_PER_SEC;

    clock_t fill_start = clock();
    
    #ifdef VERBOSE
        if (node_id == 0) {
            printf("N[%ld] i_start: %d, i_stop: %d, j_start: %d, j_stop: %d\n",
                node_id, i_start, i_stop, j_start, j_stop
            );
        }
    #endif // VERBOSE

    if (node_id == 0) {
        int64_t i;
        float val1 = 1.0;
        for (i = 0; i < mat_size; i++) {
            dsmwritefloatarrayf_(&rhs, &i, &val1);
            dsmwritefloatarrayf_(&p0, &i, &val1);
        }
    }

    /*
    int i;
    int j;
    int k;
    for (i = params.i_start; i <= params.i_stop; i += 1) {
        for (j = params.j_start; j <= params.j_stop; j += 1) {
            for (k = 0;k <= km+1;k += 1) {
                int64_t idx = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                float val = 1.0;
                dsmwritefloatarrayf_(&rhs, &idx, &val);
                dsmwritefloatarrayf_(&p0, &idx, &val);
                // Avoid possibility of undefined values in p1
                // dsmwritefloatarrayf_(&p1, &idx, &val);
            }
        }
    }
    */

    /*
    if (node_id == 0) {
        for (i = 0; i <= im+1; i++) {
            for (j = 0; j <= jm+1; j++) {
                for (k = 0; k <= km+1; k++) {
                    int64_t idx = F3D2C(im+2,jm+2,0,0,0,i,j,k);
                    float val1 = 0;
                    float val2 = 0;
                    dsmreadfloatarrayf_(&p0, &idx, &val1);
                    dsmreadfloatarrayf_(&rhs, &idx, &val2);
                    if (val1 != 1.0F) {
                        printf("Failed at p0, val: %f, i: %d, j: %d, k: %d\n\n", val1, i,j,k);
                        exit(1);
                    }
                    // assert(val1 == 1.0F);
                    if (val2 != 1.0F) {
                        printf("Failed at rhs, val: %f, i: %d, j: %d, k: %d\n\n", val2, i,j,k);
                        exit(1);
                    }
                    // assert(val2 == 1.0F);
                }
            }
        }
    }
    */
    argo::barrier();

    clock_t fill_end = clock();
    double fill_time = (double)(fill_end - fill_start) / CLOCKS_PER_SEC;
    
    #ifdef VERBOSE
        if (node_id == 0) {
            printf("Array(s) filled\n");
        }
    #endif // VERBOSE

    clock_t sor_start = clock();
    
    int iter;
    for (iter = 1;iter <= niters;iter += 1) {
        if (iter % 10 == 0 && node_id == 0) {
            printf("Iter %d\n", iter);
        }
        if (iter % 2 == 0) {
            sor_par(&p1, &p0, &rhs, params);
        } else {
            sor_par(&p0, &p1, &rhs, params);
        }
        argo::barrier();
    }

    clock_t sor_end = clock();
    double sor_time = (double)(sor_end - sor_start) / CLOCKS_PER_SEC;
    clock_t total_end = clock();
    double total_time = (double)(total_end - total_start) / CLOCKS_PER_SEC;

    if (node_id == 0) {
        int64_t idx_final = F3D2C(im+2,jm+2, 0,0,0,im/2,jm/2,km/2);
        float val_final = 0;
        dsmreadfloatarrayf_(&p0, &idx_final, &val_final);
        printf("\nDimensions: (%d, %d, %d)\nMemory: %ld\nCache: %ld\nMode: %d\nNumber of nodes: %d\n",
            im, jm, km, mem_size, cache_size, mode, num_nodes);
        printf("DSM version: 1\n");
        printf("Result: %f\nInit: %lf\nAlloc: %lf\nFill: %lf\nSOR: %lf\nTotal: %lf\n\n", 
            val_final, init_time, alloc_time, fill_time, sor_time, total_time);
        // std::cout << p0[F3D2C(im+2,jm+2, 0,0,0,im/2,jm/2,km/2)] << "\n";
        const int str_len = 50;

        char dimensions_s[str_len],
            mem_size_s[str_len],
            cache_size_s[str_len],
            dsm_version_s[str_len],
            mode_s[str_len],
            num_nodes_s[str_len],
            val_final_s[str_len], 
            init_time_s[str_len], 
            alloc_time_s[str_len], 
            fill_time_s[str_len], 
            sor_time_s[str_len], 
            total_time_s[str_len];

        snprintf(dimensions_s, str_len, "Dimensions: (%d, %d, %d)\n", im, jm, km);
        snprintf(mem_size_s, str_len, "Memory: %ld\n", mem_size);
        snprintf(cache_size_s, str_len, "Cache: %ld\n", cache_size);
        snprintf(dsm_version_s, str_len, "DSM version: 1\n");
        snprintf(mode_s, str_len, "Mode: %f\n", val_final);
        snprintf(num_nodes_s, str_len, "Number of nodes: %d\n", num_nodes);
        snprintf(val_final_s, str_len, "Result: %f\n", val_final);
        snprintf(init_time_s, str_len, "Init: %f\n", init_time);
        snprintf(alloc_time_s, str_len, "Alloc: %f\n", alloc_time);
        snprintf(fill_time_s, str_len, "Fill: %f\n", fill_time);
        snprintf(sor_time_s, str_len, "SOR: %f\n", sor_time);
        snprintf(total_time_s, str_len, "Total: %f\n", total_time);
        
        FILE* fPtr;
        std::string filename = "./results/result";
        filename += std::to_string(num_nodes);
        filename += ".txt";
        fPtr = fopen(filename.c_str(), "w");
        if (fPtr == NULL) {
            printf("Unable to create file.\n");
            exit(EXIT_FAILURE);
        }

        fputs(dimensions_s, fPtr);
        fputs(mem_size_s, fPtr);
        fputs(cache_size_s, fPtr);
        fputs(mode_s, fPtr);
        fputs(num_nodes_s, fPtr);
        fputs(val_final_s, fPtr);
        fputs(init_time_s, fPtr);
        fputs(alloc_time_s, fPtr);
        fputs(fill_time_s, fPtr);
        fputs(sor_time_s, fPtr);
        fputs(total_time_s, fPtr);

        fclose(fPtr);
    }

    argo::barrier();
    dsmdeletefloatarrayf_(&p0);
    dsmdeletefloatarrayf_(&p1);
    dsmdeletefloatarrayf_(&rhs);
    dsmfinalisef_();
}
    
