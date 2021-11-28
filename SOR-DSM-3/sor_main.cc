#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <thread>
#include <vector>
#include "array_index_f2c1d.h"
#include "sor_params.h"
#include "../src/dsm-api.h"
#include "../src/dsm-C-2D-grid-halos.cc"

#ifndef _SOR_DSM_3
#define _SOR_DSM_3
#endif

typedef struct par_params {
    int i_start;
    int i_stop;
    int j_start;
    int j_stop;
} par_params;

DSM3DArrayH a0;
DSM3DArrayH a1;
DSM3DArrayH arhs;

// void sor(float *p0,float *p1,float *rhs);
    
void sor_par(DSM3DArrayH& p0, DSM3DArrayH& p1, DSM3DArrayH& rhs, par_params params) {

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
                    float val = dsmReadArrayH(p0, i-im, j, k);
                    dsmWriteArrayH(p1, i,j,k, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i-im,j,k)];
                } else if (i==0) {
                    //std::cout <<"i==0\n";
                //  i=0
                //  circular
                    // printf("i==0\n");
                    float val = dsmReadArrayH(p0, i+im, j, k);
                    // printf("Read success! Val: %f", val);
                    dsmWriteArrayH(p1, i,j,k, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i+im,j,k)];
                } else if (j==jm+1) {
                    //std::cout <<"j==jm+1\n";
                    // printf("j==jm+1\n");
                //  open
                //  j = jm+1
                    float val = dsmReadArrayH(p0, i-1, j, k);
                    dsmWriteArrayH(p1, i,j,k, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i-1,j,k)];
                } else if (j==0) {
                    //std::cout <<"j==0\n";
                //  fixed
                //  j = 0
                    // printf("j==0\n");
                //  We keep the original values
                    float val = dsmReadArrayH(p0, i, j, k);
                    dsmWriteArrayH(p1, i,j,k, &val);
                    // p1[F3D2C(im+2,jm+2,0,0,0,i,j,k)] = p0[F3D2C(im+2,jm+2,0,0,0,i,j,k)];
                 } else if (i>0 && i<im+1 && j>0 && j<jm+1 && k>0 && k<km+1) {
                    //std::cout <<"core\n";
                    //  the core
                    //  The actual SOR expression
                    // printf("Core start\n");
                    // Right
                    float val_cn2l = dsmReadArrayH(p0, i+1, j, k);

                    // Left
                    float val_cn2s = dsmReadArrayH(p0, i-1, j, k);

                    // Up
                    float val_cn3l = dsmReadArrayH(p0, i, j+1, k);

                    // Down
                    float val_cn3s = dsmReadArrayH(p0, i, j-1, k);

                    // top
                    float val_cn4l = dsmReadArrayH(p0, i, j, k+1);

                    // Bottom
                    float val_cn4s = dsmReadArrayH(p0, i, j, k-1);

                    float val_rhs = dsmReadArrayH(rhs, i, j, k);

                    float val_p0_core = dsmReadArrayH(p0, i, j, k);

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
                    dsmWriteArrayH(p1, i, j, k, &val_new);
                    // printf("Core end\n");
                }
             }
         }
    }
}

void threaded_sor_loop(par_params& sub_params) {
    int node_id = argo::node_id();
    for (int iter = 1; iter <= niters; iter += 1) {
        if (iter % 2 == 0) {
            sor_par(a1, a0, arhs, sub_params);
        }
        else {
            sor_par(a0, a1, arhs, sub_params); 
        }
        // printf("N[%d] on iter %d\n", node_id, iter);
        argo::barrier(16);
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
        mem_size = 2*mat_num*bytes_per_float*2*(dsmNX+dsmNY)*(im+2)*(km+2) + pa_size;
        cache_size = (mem_size / node_mult) / 2UL;
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
    // DSM3DArrayH a0;
        a0.isz=im+2; a0.jsz=jm+2; a0.ksz=km+2;
        a0.ioff=0; a0.joff=0; a0.koff=0;
        a0.ilh=1; a0.jlh=1; a0.klh=1;
        a0.ihh=1; a0.jhh=1; a0.khh=1;
        a0.shared=true; a0.ptr = 0;
    a0.coll_ptr = &p0;
#ifdef NEW_DSM3        
        a0.core_ptr=0;
#endif        

    int64_t p1 = 0; 
    int64_t p1_size = dsmNX*dsmNY;
    dsmallocpointerarrayf_(&p1, &p1_size);
    // DSM3DArrayH a1;
        a1.isz=im+2; a1.jsz=jm+2; a1.ksz=km+2;
        a1.ioff=0; a1.joff=0; a1.koff=0;
        a1.ilh=1; a1.jlh=1; a1.klh=1;
        a1.ihh=1; a1.jhh=1; a1.khh=1;
        a1.shared=true; a1.ptr = 0;
    a1.coll_ptr = &p1;
#ifdef NEW_DSM3        
        a1.core_ptr=0;
#endif        

    int64_t rhs = 0;
    int64_t rhs_size = dsmNX*dsmNY;
    dsmallocpointerarrayf_(&rhs, &rhs_size);
    // DSM3DArrayH arhs;
        arhs.isz=im+2; arhs.jsz=jm+2; arhs.ksz=km+2;
        arhs.ioff=0; arhs.joff=0; arhs.koff=0;
        arhs.ilh=1; arhs.jlh=1; arhs.klh=1;
        arhs.ihh=1; arhs.jhh=1; arhs.khh=1;
        arhs.shared=true; arhs.ptr = 0;
    arhs.coll_ptr = &rhs;
#ifdef NEW_DSM3        
        arhs.core_ptr=0;
#endif        

    argo::barrier();

    #ifdef VERBOSE
    if (node_id == 0) {
        printf("Global arrays allocated\n");
    }
    #endif // VERBOSE
    
    // Work out start index, stop index and span for i & j
    int i_span = im / dsmNX;
    int i_start;
    int ilhe = a0.ilh;
    int ihhe = a0.ihh;
    if (node_id % dsmNX == 0) {
        ilhe = 2*a0.ilh;
        i_span += a0.ilh;
        i_start = 0;
    } else if (node_id % dsmNX == dsmNX-1) {
        i_start = i_span*(dsmNX-1) + a0.ilh;
        ihhe = 2*a0.ihh;
        i_span += a0.ihh;
    } else {
        i_start = i_span*(node_id % dsmNX) + a0.ilh;
    }
    // If num nodes is 1, i-dimension will run edge to edge.
    if (num_nodes == 1) {
        ilhe = 2*a0.ilh;
        ihhe = 2*a0.ihh;
        i_span += a0.ihh;
    }
    int i_stop = i_start + i_span - 1;
    

    int j_span = jm / dsmNY;
    int j_start;
    int jlhe = a0.jlh;
    int jhhe = a0.jhh;
    if (0 <= node_id && node_id < dsmNX) {
        jlhe = 2*a0.jlh;
        j_span += a0.jlh;
        j_start = 0;
    } else if (dsmNX*(dsmNY-1) <= node_id && node_id < dsmNX*dsmNY) {
        j_start = j_span*(dsmNY-1) + a0.jlh;
        jhhe = 2*a0.jhh;
        j_span += a0.jhh;
    } else {
        j_start = j_span*(node_id / dsmNX) + a0.jlh;
    }
    // If num nodes is 1 or 2, j-dimension will run edge to edge.
    if (num_nodes < 3) {
        jlhe = 2*a0.jlh;
        jhhe = 2*a0.jhh;
        j_span += a0.jhh;
    }
    int j_stop = j_start + j_span - 1;
    
    #ifndef THREADED
    par_params params = {
        i_start,
        i_stop,
        j_start,
        j_stop
    };
    #endif // THREADED
    
    #ifdef VERBOSE
        printf("N[%ld], i_span: %d, i_start: %d, i_stop: %d\nj_span: %d, j_start: %d, j_stop: %d\n",
            node_id, i_span, i_start, i_stop, j_span, j_start, j_stop
        );
        printf("node: %ld, ilhe: %d, ihhe: %d, jlhe: %d, jhhe: %d\n",
            node_id, ilhe, ihhe, jlhe, jhhe
        );
    #endif // VERBOSE

#ifndef NEW_DSM3
    // Matrix block consists of 5 segments: core, top, bottom, left, right.
    int64_t local_size = 5;

    // These variables will be filled with pointers to arrays of arrays and segments.
    int64_t a0_local = 0;
    int64_t a0_local_core = 0;
    int64_t a0_local_top = 0;
    int64_t a0_local_bottom = 0;
    int64_t a0_local_left = 0;
    int64_t a0_local_right = 0;

    int64_t a1_local = 0;
    int64_t a1_local_core = 0;
    int64_t a1_local_top = 0;
    int64_t a1_local_bottom = 0;
    int64_t a1_local_left = 0;
    int64_t a1_local_right = 0;

    int64_t arhs_local = 0;
    int64_t arhs_local_core = 0;
    int64_t arhs_local_top = 0;
    int64_t arhs_local_bottom = 0;
    int64_t arhs_local_left = 0;
    int64_t arhs_local_right = 0;

    // Allocate array of arrays.
    dsmallocpointerlocalarrayf_(&a0_local, &local_size);
    dsmallocpointerlocalarrayf_(&a1_local, &local_size);
    dsmallocpointerlocalarrayf_(&arhs_local, &local_size);
#endif
    // Work out segment sizes.
    int64_t local_core_size = (int64_t)(i_span-ilhe-ihhe)*(int64_t)(j_span-jlhe-jhhe)*(int64_t)(km+2);
    int64_t local_top_size = jlhe*(i_span)*(km+2);
    int64_t local_bottom_size = jhhe*(i_span)*(km+2);
    int64_t local_left_size = ilhe*(j_span-jlhe-jhhe)*(km+2);
    int64_t local_right_size = ihhe*(j_span-jlhe-jhhe)*(km+2);
#ifdef NEW_DSM3
    // core is a local new() array, to be stored in .local_ptr
    float* a0_local_core_f = new float[local_core_size];
    a0.core_ptr = (int64_t*)toWord<float*>(a0_local_core_f);
    // There is a single argo::new_array combining all edges
    int64_t local_edges_size = local_top_size + local_bottom_size + local_left_size + local_right_size;
    int64_t a0_local_edges = 0;    
    dsmallocfloatlocalarrayf_(&a0_local_edges, &local_edges_size);
    argo::barrier();
#else
    // Allocate segment arrays for a0, a1 & arhs.
    dsmallocfloatlocalarrayf_(&a0_local_core, &local_core_size);
    dsmallocfloatlocalarrayf_(&a0_local_top, &local_top_size);
    dsmallocfloatlocalarrayf_(&a0_local_bottom, &local_bottom_size);
    dsmallocfloatlocalarrayf_(&a0_local_left, &local_left_size);
    dsmallocfloatlocalarrayf_(&a0_local_right, &local_right_size);
    argo::barrier();
#endif
#ifdef NEW_DSM3
    // core is a local new() array, to be stored in .local_ptr
    float* a1_local_core_f = new float[local_core_size];
    a1.core_ptr = (int64_t*)toWord<float*>(a1_local_core_f);
    int64_t a1_local_edges = 0;    
    dsmallocfloatlocalarrayf_(&a1_local_edges, &local_edges_size);
    argo::barrier();
#else
    dsmallocfloatlocalarrayf_(&a1_local_core, &local_core_size);
    dsmallocfloatlocalarrayf_(&a1_local_top, &local_top_size);
    dsmallocfloatlocalarrayf_(&a1_local_bottom, &local_bottom_size);
    dsmallocfloatlocalarrayf_(&a1_local_left, &local_left_size);
    dsmallocfloatlocalarrayf_(&a1_local_right, &local_right_size);
    argo::barrier();
#endif
#ifdef NEW_DSM3
    // core is a local new() array, to be stored in .local_ptr
    float* arhs_local_core_f = new float[local_core_size];
    arhs.core_ptr = (int64_t*)toWord<float*>(arhs_local_core_f);
    int64_t arhs_local_edges = 0;    
    dsmallocfloatlocalarrayf_(&arhs_local_edges, &local_edges_size);
    argo::barrier();
#else
    dsmallocfloatlocalarrayf_(&arhs_local_core, &local_core_size);
    dsmallocfloatlocalarrayf_(&arhs_local_top, &local_top_size);
    dsmallocfloatlocalarrayf_(&arhs_local_bottom, &local_bottom_size);
    dsmallocfloatlocalarrayf_(&arhs_local_left, &local_left_size);
    dsmallocfloatlocalarrayf_(&arhs_local_right, &local_right_size);
    argo::barrier();
#endif
#ifndef NEW_DSM3
    // Fill collective array of arrays.
    int64_t segment = 0;
    dsmwritepointerlocalarrayf_(&a0_local, &segment, &a0_local_core);
    dsmwritepointerlocalarrayf_(&a1_local, &segment, &a1_local_core);
    dsmwritepointerlocalarrayf_(&arhs_local, &segment, &arhs_local_core);
    segment = 1;
    dsmwritepointerlocalarrayf_(&a0_local, &segment, &a0_local_top);
    dsmwritepointerlocalarrayf_(&a1_local, &segment, &a1_local_top);
    dsmwritepointerlocalarrayf_(&arhs_local, &segment, &arhs_local_top);
    segment = 2;
    dsmwritepointerlocalarrayf_(&a0_local, &segment, &a0_local_bottom);
    dsmwritepointerlocalarrayf_(&a1_local, &segment, &a1_local_bottom);
    dsmwritepointerlocalarrayf_(&arhs_local, &segment, &arhs_local_bottom);
    segment = 3;
    dsmwritepointerlocalarrayf_(&a0_local, &segment, &a0_local_left);
    dsmwritepointerlocalarrayf_(&a1_local, &segment, &a1_local_left);
    dsmwritepointerlocalarrayf_(&arhs_local, &segment, &arhs_local_left);
    segment = 4;
    dsmwritepointerlocalarrayf_(&a0_local, &segment, &a0_local_right);
    dsmwritepointerlocalarrayf_(&a1_local, &segment, &a1_local_right);
    dsmwritepointerlocalarrayf_(&arhs_local, &segment, &arhs_local_right);

    argo::barrier();
#endif
    #ifdef VERBOSE
    if (node_id == 0) {
        printf("Local arrays allocated\n");
    }
    #endif // VERBOSE
#ifdef NEW_DSM3
    // Set collective array of edge arrays pointer.
    dsmwritepointerarrayf_(a0.coll_ptr, &node_id, &a0_local_edges);
    dsmwritepointerarrayf_(a1.coll_ptr, &node_id, &a1_local_edges);
    dsmwritepointerarrayf_(arhs.coll_ptr, &node_id, &arhs_local_edges);
#else
    // Set collective array of arrays pointer.
    dsmwritepointerarrayf_(a0.coll_ptr, &node_id, &a0_local);
    dsmwritepointerarrayf_(a1.coll_ptr, &node_id, &a1_local);
    dsmwritepointerarrayf_(arhs.coll_ptr, &node_id, &arhs_local);
#endif
    argo::barrier();

    clock_t alloc_end = clock();
    double alloc_time = (double)(alloc_end - alloc_start) / CLOCKS_PER_SEC;

    #ifdef VERBOSE
        printf("Global arrays written to\n");
    #endif // VERBOSE

    // Copy collective array of arrays to local array of arrays.
    #ifdef NEW_DSM3
    a0.ptr = new int64_t[num_nodes];
    a1.ptr = new int64_t[num_nodes];
    arhs.ptr = new int64_t[num_nodes];
    #else
    a0.ptr = new int64_t*[num_nodes];
    a1.ptr = new int64_t*[num_nodes];
    arhs.ptr = new int64_t*[num_nodes];
    #endif
    for (int64_t i = 0; i < num_nodes; i++) {
        int64_t val1 = 0;
        dsmreadpointerarrayf_(a0.coll_ptr, &i, &val1);
        int64_t val2 = 0;
        dsmreadpointerarrayf_(a1.coll_ptr, &i, &val2);
        int64_t val3 = 0;
        dsmreadpointerarrayf_(arhs.coll_ptr, &i, &val3);
#ifdef NEW_DSM3
        a0.ptr[i] = val1;
        a1.ptr[i] = val2;
        arhs.ptr[i] = val3;
#else
        a0.ptr[i] = new int64_t[5];
        a1.ptr[i] = new int64_t[5];
        arhs.ptr[i] = new int64_t[5];
        for (int64_t j = 0; j < 5; j++) {
            int64_t val = 0;
            dsmreadpointerarrayf_(&val1, &j, &val);
            a0.ptr[i][j] = val;
            dsmreadpointerarrayf_(&val2, &j, &val);
            a1.ptr[i][j] = val;
            dsmreadpointerarrayf_(&val3, &j, &val);
            arhs.ptr[i][j] = val;
        }
#endif        
    }

    #ifdef VERBOSE
        printf("Collective array fetched\n");
    #endif // VERBOSE

    argo::barrier();

    clock_t fill_start = clock();
    
    int ten_percent_point = i_span / 10;
    // Prevent division by 0 when checking modulo.
    if (ten_percent_point == 0) {
        ten_percent_point = 1;
    }
    // Fill 3D arrays a0 & arhs with values (1.0).
    for (int i = i_start; i <= i_stop; i++) {
        #ifdef VERBOSE
            if (i % ten_percent_point == 0) {
                printf("Node[%ld] filled to %d\n", node_id, i);
            }
        #endif // VERBOSE
        float val = 1.0;
        for (int j = j_start; j <= j_stop; j++) {
            for (int k = 0; k <= km+1; k++) {
                dsmWriteArrayH(a0, i, j, k, &val);
                dsmWriteArrayH(arhs, i, j, k, &val);
            }
        }
    }
    /*
    if (node_id == 0) {
        // Fill 3D arrays a0 & arhs with values (1.0).
        for (int i = 0; i <= im+1; i++) {
            float val = 1.0;
            for (int j = 0; j <= jm+2; j++) {
                for (int k = 0; k <= km+1; k++) {
                    dsmWriteArrayH(a0, i, j, k, &val);
                    dsmWriteArrayH(arhs, i, j, k, &val);
                }
            }
        }
    }
    */
    argo::barrier();

    #ifdef VERBOSE
    if (node_id == 0) {
        printf("Local arrays filled\n");
    }
    #endif // VERBOSE

    clock_t fill_end = clock();
    double fill_time = (double)(fill_end - fill_start) / CLOCKS_PER_SEC;
#ifndef NEW_DSM3
    #ifdef DEBUG
        // Verify local array of arrays correctness.
        for (int64_t i = 0; i < num_nodes; i++) {
            int64_t val1 = 0;
            dsmreadpointerarrayf_(a0.coll_ptr, &i, &val1);
            int64_t val2 = 0;
            dsmreadpointerarrayf_(a1.coll_ptr, &i, &val2);
            int64_t val3 = 0;
            dsmreadpointerarrayf_(arhs.coll_ptr, &i, &val3);

            for (int64_t j = 0; j < 5; j++) {
                int64_t val = 0;
                dsmreadpointerarrayf_(&val1, &j, &val);
                if (a0.ptr[i][j] != val) {
                    printf("Pointer failure at node[%ld], array[%ld]\n", i, j);
                }
                // assert(a0.ptr[i][j] == val);
                dsmreadpointerarrayf_(&val2, &j, &val);
                if (a1.ptr[i][j] != val) {
                    printf("Pointer failure at node[%ld], array[%ld]\n", i, j);
                }
                // assert(a1.ptr[i][j] == val);
                dsmreadpointerarrayf_(&val3, &j, &val);
                if (arhs.ptr[i][j] != val) {
                    printf("Pointer failure at node[%ld], array[%ld]\n", i, j);
                }
                // assert(arhs.ptr[i][j] == val);
            }
        }
    #endif // DEBUG
#endif    
    /*
    // Verify all points are correctly written to.
    if (node_id == 0) {
        for (int i = 0; i < im+2; i++) {
            float tval = 1.0;
            for (int j = 0; j < jm+2; j++) {
                for (int k = 0; k < km+2; k++) {
                    float rval = 0.0;
                    rval = dsmReadArrayH(a0, i, j, k);
                    assert(rval == tval);
                    rval = dsmReadArrayH(arhs, i, j, k);
                    assert(rval == tval);
                }
            }
        }
    }
    */

    argo::barrier();
#ifndef NEW_DSM3    
    // Check local segment arrays contain expected values.
    #ifdef DEBUG
        // To keep all these temporary values off the stack frame.
        if (true) {
            float val1 = 0.0;
            float val2 = 0.0;
            int core_fails = 0;
            int top_fails = 0;
            int bottom_fails = 0;
            int left_fails = 0;
            int right_fails = 0;
            int core_fails2 = 0;
            int top_fails2 = 0;
            int bottom_fails2 = 0;
            int left_fails2 = 0;
            int right_fails2 = 0;
            for (int64_t i = 0; i < local_core_size; i++) {
                dsmreadfloatlocalarrayf_(&a0_local_core, &i, &val1);
                dsmreadfloatlocalarrayf_(&arhs_local_core, &i, &val2);
                if (val1 != 1.0) {
                    // printf("Incorrect a0 value at core[%ld]\n", i);
                    core_fails++;
                }
                if (val2 != 1.0) {
                    // printf("Incorrect arhs value at core[%ld]\n", i);
                    core_fails2++;
                }
            }
            for (int64_t i = 0; i < local_top_size; i++) {
                dsmreadfloatlocalarrayf_(&a0_local_top, &i, &val1);
                dsmreadfloatlocalarrayf_(&arhs_local_top, &i, &val2);
                if (val1 != 1.0) {
                    // printf("Incorrect a0 value at top[%ld]\n", i);
                    top_fails++;
                }
                if (val2 != 1.0) {
                    // printf("Incorrect arhs value at top[%ld]\n", i);
                    top_fails2++;
                }
            }
            for (int64_t i = 0; i < local_bottom_size; i++) {
                dsmreadfloatlocalarrayf_(&a0_local_bottom, &i, &val1);
                dsmreadfloatlocalarrayf_(&arhs_local_bottom, &i, &val2);
                if (val1 != 1.0) {
                    // printf("Incorrect a0 value at bottom[%ld]\n", i);
                    bottom_fails++;
                }
                if (val2 != 1.0) {
                    // printf("Incorrect arhs value at bottom[%ld]\n", i);
                    bottom_fails2++;
                }
            }
            for (int64_t i = 0; i < local_left_size; i++) {
                dsmreadfloatlocalarrayf_(&a0_local_left, &i, &val1);
                dsmreadfloatlocalarrayf_(&arhs_local_left, &i, &val2);
                if (val1 != 1.0) {
                    // printf("Incorrect a0 value at left[%ld]\n", i);
                    left_fails++;
                }
                if (val2 != 1.0) {
                    // printf("Incorrect arhs value at left[%ld]\n", i);
                    left_fails2++;
                }
            }
            for (int64_t i = 0; i < local_right_size; i++) {
                dsmreadfloatlocalarrayf_(&a0_local_right, &i, &val1);
                dsmreadfloatlocalarrayf_(&arhs_local_right, &i, &val2);
                if (val1 != 1.0) {
                    // printf("Incorrect a0 value at right[%ld]\n", i);
                    right_fails++;
                }
                if (val2 != 1.0) {
                    // printf("Incorrect arhs value at right[%ld]\n", i);
                    right_fails2++;
                }
            }
            printf("N[%ld] a0: Core fails: %d/%ld\nTop fails: %d/%ld, Bottom fails: %d/%ld\nLeft fails: %d/%ld, Right fails: %d/%ld\n", 
                node_id,
                core_fails, local_core_size, 
                top_fails, local_top_size, 
                bottom_fails, local_bottom_size, 
                left_fails, local_left_size, 
                right_fails, local_right_size
            );
            printf("N[%ld] arhs: Core fails: %d/%ld\nTop fails: %d/%ld, Bottom fails: %d/%ld\nLeft fails: %d/%ld, Right fails: %d/%ld\n", 
                node_id,
                core_fails2, local_core_size, 
                top_fails2, local_top_size, 
                bottom_fails2, local_bottom_size, 
                left_fails2, local_left_size, 
                right_fails2, local_right_size
            );
        }
    #endif //DEBUG
#endif
   #ifdef VERBOSE
    if (node_id == 0) {
            printf("Entering SOR loop\n");
    }
    #endif // VERBOSE

    clock_t sor_start = clock();

    #ifdef THREADED
        int sub_i_span = i_span / 4;
        int sub_j_span = j_span / 4;
        std::vector<std::thread> threads;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int sub_i_start = i_start + i*sub_i_span;
                int sub_i_stop = i == 3 ? i_stop : i_start + (i+1)*sub_i_span - 1;
                int sub_j_start = j_start + j*sub_j_span;
                int sub_j_stop = j == 3 ? j_stop : j_start + (j+1)*sub_j_span - 1;
                par_params sub_params = {sub_i_start, sub_i_stop, sub_j_start, sub_j_stop};
                printf("N[%ld] T[%d] - i_start: %d, i_stop: %d, j_start: %d, j_stop: %d\n",
                    node_id, 4*i + j, sub_i_start, sub_i_stop, sub_j_start, sub_j_stop);
                threads.push_back(std::thread(threaded_sor_loop, std::ref(sub_params)));
            }
        }
        for (int i = 0; i < 16; i++) {
            threads[i].join();
        }
    #else
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
    #endif // THREADED

    clock_t sor_end = clock();
    double sor_time = (double)(sor_end - sor_start) / CLOCKS_PER_SEC;
    clock_t total_end = clock();
    double total_time = (double)(total_end - total_start) / CLOCKS_PER_SEC;

    if (node_id == 0) {
        float val_final = dsmReadArrayH(a0, im/2, jm/2, km/2);
        printf("\nDimensions: (%d, %d, %d)\nMemory: %ld\nCache: %ld\nMode: %d\nNumber of nodes: %d\n",
            im, jm, km, mem_size, cache_size, mode, num_nodes);
        printf("DSM version: 3\n");
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
        snprintf(dsm_version_s, str_len, "DSM version: 3\n");
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
#ifdef NEW_DSM3
    delete[] a0.core_ptr;
    delete[] a1.core_ptr;
    delete[] arhs.core_ptr;
#endif    
    // dsmdeletepointerarrayf_(a0.ptr);
    // dsmdeletepointerarrayf_(a1.ptr);
    // dsmdeletepointerarrayf_(arhs.ptr);
    dsmdeletefloatlocalarrayf_(&p0);
    dsmdeletefloatlocalarrayf_(&p1);
    dsmdeletefloatlocalarrayf_(&rhs);
    dsmfinalisef_();
}
    
