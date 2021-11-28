#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "array_index_f2c1d.h"
#include "sor_params.h"
#include "../src/dsm-api.h"
#include "../src/dsm-C-2D-grid-halos.cc"

#ifndef _SOR_DSM_2
#define _SOR_DSM_2
#endif

typedef struct par_params {
    int64_t i_start;
    int64_t i_stop;
    int64_t j_start;
    int64_t j_stop;
} par_params;

// void sor(float *p0,float *p1,float *rhs);
    
void sor_par(DSM3DArray& p0, DSM3DArray& p1, DSM3DArray& rhs, par_params params) {

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
        for (j = params.j_start; j <= params.j_stop; j += 1) {
            for (k = 0; k <= km+1; k += 1) {
                //std::cout << "("<<i <<","<< j <<","<< k << ");"<<F3D2C(im+2,jm+2,0,0,0,i,j,k)<<"\n";
                //  assume i=x =  west to east, y=j = south to north, k=z = vertical
                if (i==im+1) {
                    //std::cout <<"i==im+1\n";
                //  circular
                //  i=im+1
                    float val = dsmRead3DReal4Array(p0, i-im, j, k);
                    dsmWrite3DReal4Array(p1, i,j,k, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i-im,j,k)];
                } else if (i==0) {
                    //std::cout <<"i==0\n";
                //  i=0
                //  circular
                    // printf("i==0\n");
                    float val = dsmRead3DReal4Array(p0, i+im, j, k);
                    // printf("Read success! Val: %f", val);
                    dsmWrite3DReal4Array(p1, i,j,k, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i+im,j,k)];
                } else if (j==jm+1) {
                    //std::cout <<"j==jm+1\n";
                    // printf("j==jm+1\n");
                //  open
                //  j = jm+1
                    float val = dsmRead3DReal4Array(p0, i-1, j, k);
                    dsmWrite3DReal4Array(p1, i,j,k, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i-1,j,k)];
                } else if (j==0) {
                    //std::cout <<"j==0\n";
                //  fixed
                //  j = 0
                    // printf("j==0\n");
                //  We keep the original values
                    float val = dsmRead3DReal4Array(p0, i, j, k);
                    dsmWrite3DReal4Array(p1, i,j,k, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i,j,k)];
                 } else if (i>0 && i<im+1 && j>0 && j<jm+1 && k>0 && k<km+1) {
                    //std::cout <<"core\n";
                    //  the core
                    //  The actual SOR expression
                    // printf("Core start\n");
                    // Right
                    float val_cn2l = dsmRead3DReal4Array(p0, i+1, j, k);

                    // Left
                    float val_cn2s = dsmRead3DReal4Array(p0, i-1, j, k);

                    // Up
                    float val_cn3l = dsmRead3DReal4Array(p0, i, j+1, k);

                    // Down
                    float val_cn3s = dsmRead3DReal4Array(p0, i, j-1, k);

                    // top
                    float val_cn4l = dsmRead3DReal4Array(p0, i, j, k+1);

                    // Bottom
                    float val_cn4s = dsmRead3DReal4Array(p0, i, j, k-1);

                    float val_rhs = dsmRead3DReal4Array(rhs, i, j, k);

                    float val_p0_core = dsmRead3DReal4Array(p0, i, j, k);

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
                    dsmWrite3DReal4Array(p1, i, j, k, &val_new);
                    // printf("Core end\n");
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
    printf("Memory Initialised\n");

    int num_nodes;
    dsmgetnnodesf_(&num_nodes);
    assert(num_nodes == dsmNX*dsmNY); // There must be an MPI node for each sub-array.

    int64_t node_id;
    dsmgetnidf_(&node_id);
    // printf("Number of nodes: %d. This node: %d\n", num_nodes, node_id);
    // sleep(1);

    clock_t init_end = clock();
    double init_time = (double)(init_end - init_start) / CLOCKS_PER_SEC;

    clock_t alloc_start = clock();
    // reversed p1 and p0 to see if it still segfaults. 
    // It doesn't. Does the other case now segfault?
    // It's a Heisenbug. Either segfaults or gives the "backing node" error.
    // I wonder what would happen if I created a single array for p0 and p1?
    int64_t p0 = 0;
    int64_t p0_size = dsmNX*dsmNY;
    dsmallocpointerarrayf_(&p0, &p0_size);
    DSM3DArray a0;
    a0.isz=im+2; a0.jsz=jm+2; a0.ksz=km+2;
    a0.ioff=0; a0.joff=0; a0.koff=0;
    a0.ilh=1; a0.jlh=1; a0.klh=1;
    a0.ihh=1; a0.jhh=1; a0.khh=1;
    a0.ptr=0; a0.shared=true;
    a0.coll_ptr = &p0;

    int64_t p1 = 0; 
    int64_t p1_size = dsmNX*dsmNY;
    dsmallocpointerarrayf_(&p1, &p1_size);
    DSM3DArray a1;
    a1.isz=im+2; a1.jsz=jm+2; a1.ksz=km+2;
    a1.ioff=0; a1.joff=0; a1.koff=0;
    a1.ilh=1; a1.jlh=1; a1.klh=1;
    a1.ihh=1; a1.jhh=1; a1.khh=1;
    a1.ptr=0; a1.shared=true;
    a1.coll_ptr = &p1;

    int64_t rhs = 0;
    int64_t rhs_size = dsmNX*dsmNY;
    dsmallocpointerarrayf_(&rhs, &rhs_size);
    DSM3DArray arhs;
        arhs.isz=im+2; arhs.jsz=jm+2; arhs.ksz=km+2;
        arhs.ioff=0; arhs.joff=0; arhs.koff=0;
        arhs.ilh=1; arhs.jlh=1; arhs.klh=1;
        arhs.ihh=1; arhs.jhh=1; arhs.khh=1;
        arhs.ptr=0; arhs.shared=true;
    arhs.coll_ptr = &rhs;
    argo::barrier();

    #ifdef VERBOSE
    if (node_id == 0) {
        printf("Global arrays allocated\n");
    }
    #endif // VERBOSE
    
    int64_t i_span = im / dsmNX;
    int64_t i_start;
    if (node_id % dsmNX == 0) {
        i_span += a0.ilh;
        i_start = 0;
    } else if (node_id % dsmNX == dsmNX-1) {
        i_start = i_span*(dsmNX-1) + a0.ilh;
        i_span += a0.ilh;
    } else {
        i_start = i_span*(node_id % dsmNX) + a0.ilh;
    }
    if (num_nodes == 1) {
        i_span += a0.ilh;
    }
    int64_t i_stop = i_start + i_span - 1;

    int64_t j_span = jm / dsmNY;
    int64_t j_start;
    if (0 <= node_id && node_id < dsmNX) {
        j_span += a0.jlh;
        j_start = 0;
    } else if (dsmNX*(dsmNY-1) <= node_id && node_id < dsmNX*dsmNY) {
        j_start = j_span*(dsmNY-1) + a0.jlh;
        j_span += a0.jlh;
    } else {
        j_start = j_span*(node_id / dsmNX) + a0.jlh;
    }
    if (num_nodes < 3) {
        j_span += a0.jlh;
    }
    int64_t j_stop = j_start + j_span - 1;
    
    par_params params = {
        i_start,
        i_stop,
        j_start,
        j_stop
    };
    
    #ifdef VERBOSE
    if (node_id == 0) {
        printf("node: %ld, i_span: %ld, i_start: %ld, i_stop: %ld, j_span: %ld, j_start: %ld, j_stop: %ld\n",
            node_id, i_span, i_start, i_stop, j_span, j_start, j_stop
        );
    }
    #endif // VERBOSE

    // printf("Allocating first local array\n");

    
    int64_t a0_local = 0;
    int64_t a0_local_size = (int64_t)i_span*(int64_t)j_span*(int64_t)(km+2);
    dsmallocfloatlocalarrayf_(&a0_local, &a0_local_size);
    argo::barrier();
    // printf("Allocating second local array\n");
    int64_t a1_local = 0;
    int64_t a1_local_size = (int64_t)i_span*(int64_t)j_span*(int64_t)(km+2);
    dsmallocfloatlocalarrayf_(&a1_local, &a1_local_size);
    argo::barrier();
    // printf("Allocating third local array\n");
    int64_t arhs_local = 0;
    int64_t arhs_local_size = (int64_t)i_span*(int64_t)j_span*(int64_t)(km+2);
    dsmallocfloatlocalarrayf_(&arhs_local, &arhs_local_size);
    argo::barrier();

    #ifdef VERBOSE
    if (node_id == 0) {
        printf("Local arrays allocated\n");
    }
    #endif // VERBOSE

    // dsmwritepointerarrayf_(a0.ptr, &node_id, &a0_local);
    // dsmwritepointerarrayf_(a1.ptr, &node_id, &a1_local);
    // dsmwritepointerarrayf_(arhs.ptr, &node_id, &arhs_local);
    // argo::barrier();

    // #ifdef VERBOSE
    //     printf("Global arrays written to\n");
    // #endif // VERBOSE


    clock_t alloc_end = clock();
    double alloc_time = (double)(alloc_end - alloc_start) / CLOCKS_PER_SEC;

    clock_t fill_start = clock();

    int64_t idx;
    float val1 = 1.0;
    for (idx = 0; idx < a0_local_size; idx++) {
        dsmwritefloatlocalarrayf_(&arhs_local, &idx, &val1);
        dsmwritefloatlocalarrayf_(&a0_local, &idx, &val1);
    }
    argo::barrier();

    /*
    int i;
    int j;
    int k;
    for (i = 0; i < i_span; i += 1) {
        for (j = 0; j < j_span; j += 1) {
            for (k = 0; k <= km+1; k += 1) {
                    int64_t idx = F3D2C(i_span,j_span,0,0,0,i,j,k);
                    // dsmwritefloatarrayf_(&arhs_local, &idx, &val1);
                    // dsmwritefloatarrayf_(&a0_local, &idx, &val1);
                    dsmwritefloatarrayf_(&a1_local, &idx, &val1);
            }
        }
    }
    argo::barrier();       
    */ 

    #ifdef VERBOSE
    if (node_id == 0) {
        printf("Local arrays filled\n");
    }
    #endif // VERBOSE

    // WV: This was a0 first, then a1; but it still segfaults with this order.
    dsmwritepointerarrayf_(a0.coll_ptr, &node_id, &a0_local);
    argo::barrier();
    dsmwritepointerarrayf_(a1.coll_ptr, &node_id, &a1_local);
    argo::barrier();
    dsmwritepointerarrayf_(arhs.coll_ptr, &node_id, &arhs_local);
    argo::barrier();

    #ifdef VERBOSE
    if (node_id == 0) {
        printf("Global arrays written to\n");
    }
    #endif // VERBOSE

    // Copy collective array of arrays to local array of arrays.
    a0.ptr = new int64_t*[num_nodes];
    a1.ptr = new int64_t*[num_nodes];
    arhs.ptr = new int64_t*[num_nodes];
    for (int64_t i = 0; i < num_nodes; i++) {
        int64_t val1 = 0;
        int64_t val2 = 0;
        int64_t val3 = 0;
        dsmreadpointerarrayf_(a0.coll_ptr, &i, &val1);
        dsmreadpointerarrayf_(a1.coll_ptr, &i, &val2);
        dsmreadpointerarrayf_(arhs.coll_ptr, &i, &val3);
        a0.ptr[i] = (int64_t*)val1;
        a1.ptr[i] = (int64_t*)val2;
        arhs.ptr[i] = (int64_t*)val3;
    }

    clock_t fill_end = clock();
    double fill_time = (double)(fill_end - fill_start) / CLOCKS_PER_SEC;
    /*
    argo::barrier();

    if (node_id == 0) {
        for (i = 0; i <= im+1; i++) {
            for (j = 0; j <= jm+1; j++) {
                for (k = 0; k <= km+1; k++) {
                    float val1 = dsmRead3DReal4Array(a0, i, j, k);
                    float val2 = dsmRead3DReal4Array(arhs, i, j, k);
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
    /*
    for (int64_t i = 0; i < num_nodes; i++) {
        int64_t tt1;
        dsmreadpointerarrayf_(a0.ptr, &i, &tt1);
        assert(tt1 != 0);
        // printf("N[%ld] a0[%ld] : %ld\n", node_id, i, tt);
        int64_t tt2;
        dsmreadpointerarrayf_(a1.ptr, &i, &tt2);
        assert(tt2 != 0);
        assert(tt2 != tt1);
        // printf("N[%ld] a1[%ld] : %ld\n", node_id, i, tt);
        int64_t tt3;
        dsmreadpointerarrayf_(arhs.ptr, &i, &tt3);
        assert(tt3 != 0);
        assert(tt3 != tt1);
        assert(tt3 != tt2);
        // printf("N[%ld] rhs[%ld] : %ld\n", node_id, i, tt);
    }
    argo::barrier();
    argo::barrier();
    */

   #ifdef VERBOSE
    if (node_id == 0) {
            printf("Entering SOR loop\n");
    }
    #endif // VERBOSE

    clock_t sor_start = clock();
    
    int iter;
    for (iter = 1;iter <= niters;iter += 1) {
        if (
            iter % 10 == 1 && node_id == 0) {
            printf("Iter %d\n", iter);
        }
        if (iter % 2 == 0) {
            sor_par(a1, a0, arhs, params); // This segfaults. We read from a1 and write to a0. And if it does not segfault, if gives
//             ArgoDSM failed to fetch a valid backing node. Please report a bug.
// terminate called after throwing an instance of 'std::system_error'
//   what():  ArgoDSM failed to fetch a valid backing node. Please report a bug.: No such file or directory

        }
        else {
            // sor_par(a0, a1, arhs, params); // This works fine. We read from a0 and write to a1
            // If I change this to sor_par(a0, a0,...) it gives the "backing node" error
            // 
            sor_par(a0, a1, arhs, params);
        }
        argo::barrier();
    }

    clock_t sor_end = clock();
    double sor_time = (double)(sor_end - sor_start) / CLOCKS_PER_SEC;
    clock_t total_end = clock();
    double total_time = (double)(total_end - total_start) / CLOCKS_PER_SEC;

    if (node_id == 0) {
        float val_final = dsmRead3DReal4Array(a0, im/2, jm/2, km/2);
        printf("\nDimensions: (%d, %d, %d)\nMemory: %ld\nCache: %ld\nMode: %d\nNumber of nodes: %d\n",
            im, jm, km, mem_size, cache_size, mode, num_nodes);
        printf("DSM version: 2\n");
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
        snprintf(dsm_version_s, str_len, "DSM version: 2\n");
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
    dsmdeletepointerarrayf_(a0.coll_ptr);
    dsmdeletepointerarrayf_(a1.coll_ptr);
    dsmdeletepointerarrayf_(arhs.coll_ptr);
    delete(a0.ptr);
    delete(a1.ptr);
    delete(arhs.ptr);
    dsmdeletefloatlocalarrayf_(&p0);
    dsmdeletefloatlocalarrayf_(&p1);
    dsmdeletefloatlocalarrayf_(&rhs);
    dsmfinalisef_();
}
    
