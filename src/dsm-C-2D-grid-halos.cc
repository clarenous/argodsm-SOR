#include "dsm-api.h"

// For devel
// #ifdef _SOR_PARAMS_H_
//     #define dsmNX 4
//     #define dsmNY 2
// #endif // _SOR_PARAMS_H_
#include <iostream>

#ifdef _SOR_DSM_1
#include "../SOR-DSM-1/sor_params.h"
#endif // _SOR_DSM_1

#ifdef _SOR_DSM_2
#include "../SOR-DSM-2/sor_params.h"
#endif // _SOR_DSM_2

#ifdef _SOR_DSM_3
#include "../SOR-DSM-3/sor_params.h"
#endif // _SOR_DSM_3


void indexing_case3(DSM3DArrayH& a, int i, int j, int k, int64_t* part_idx, int64_t* sel, int64_t* comb_lin_idx) {
    // icsz = a%isz - a%ilh - a%ihh
    int icsz = a.isz - a.ilh - a.ihh;
    // jcsz = a%jsz - a%jlh - a%jhh
    int jcsz = a.jsz - a.jlh - a.jhh;

    // nicsz = icsz / dsmNX
    int nicsz = icsz / dsmNX;
    // njcsz = jcsz / dsmNY
    int njcsz = jcsz / dsmNY;

    // nx = (i-1-ioff) / nicsz

    int nx = i-1-a.ioff < dsmNX*nicsz ? (i-1-a.ioff) / nicsz : dsmNX-1;
    // printf("(i:%d-1-a.ioff:%ld) / nicsz:%d = nx:%d\n",i,a.ioff,nicsz,nx);
    // ny = (j-1-joff) / njcsz
    int ny = (j-1-a.joff) < dsmNY*njcsz ? (j-1-a.joff) / njcsz : dsmNY-1;
    // printf("(j:%d-1-a.joff:%ld) / njcsz:%d = ny:%d\n",j,a.joff,njcsz,ny);
    int in=-1; // for error checking at it should never be <0
    int jn=-1;

    // if (nx>0 .and. (nx < dsmNX-1 .or. (nx == dsmNX-1 .and. i<=icsz)) )
    if (nx > 0 && (nx < dsmNX-1 || (nx == dsmNX-1 && i <= icsz))) {
        // in = mod(i-1-ioff, nicsz)
        in = (i-1-a.ioff) % nicsz;
    // else if (nx == 0)
    } else if (nx == 0) {
        // in = i-1-ioff+ilh
        in = i-1-a.ioff+a.ilh;
    // else if (nx == dsmNX-1 .and. i>icsz)
    } else if (nx == dsmNX-1 && i > icsz) {
        // in = i-1-ioff-icsz+nicsz
        in = i-1-a.ioff-icsz+nicsz;
    }

    //if(ny> 0.and.(ny< dsmNY-1 .or. (ny == dsmNY-1 .and. j<=jcsz)) )
    if (ny > 0 && (ny < dsmNY-1 || (ny == dsmNY-1 && j <= jcsz))) {
        // jn = mod(j-1-joff, njcsz)
        jn = (j-1-a.joff) % njcsz;
    // else if (ny == 0)
    } else if (ny == 0) {
        //jn=j-1-  joff+  jlh
        jn = j-1-a.joff+a.jlh;
    // else if(ny == dsmNY-1.and.j> jcsz)
    } else if (ny == dsmNY-1 && j > jcsz) {
        //jn=j-1-  joff-jcsz+njcsz
        jn = j-1-a.joff-jcsz+njcsz;
    }

    // part_id =nx+ny*dsmNX
    *part_idx = nx+ny*dsmNX;
    int64_t lin_idx = 0;
    // if(nx>0 .and. nx<dsmNX-1 )
    if (nx > 0 && nx < dsmNX-1) {
        // lin_idx = k-1-koff+klh+ksz*(in+nicsz*jn) ! for the core
        lin_idx = k-1-a.koff+a.klh+a.ksz*(in+nicsz*jn);
    // else if (nx==0)
    } else if (nx == 0) {
        // lin_idx = k-1-koff+klh+ksz*(in+(nicsz+ilh)*jn) ! left
        lin_idx = k-1-a.koff+a.klh+a.ksz*(in+(nicsz+a.ilh)*jn);
    } else {
        // lin_idx = k-1-koff+klh+ksz*(in+(nicsz+ihh)*jn) ! right
        lin_idx = k-1-a.koff+a.klh+a.ksz*(in+(nicsz+a.ihh)*jn);
    }
    /*
    if (*lin_idx < 0 || *lin_idx >= (icsz+2)*(jcsz+2)*(km+2)) {
        printf("i: %d, j: %d, k: %d, dsmNX: %d, dsmNY: %d, \n", 
            i, j, k, dsmNX, dsmNY
            );
        printf("nicsz: %d, njcsz: %d, nx: %d, ny: %d, in: %d, jn: %d, lin_idx: %ld, part_idx: %ld\n",
            nicsz, njcsz, nx, ny, in, jn, *lin_idx, *part_idx
        );
    }
    */

    // ! * Case 3
    // ! ========

    // ! The calculation of nx,ny and in, jn is the same as for Case 2

    // ! Then we need to work out in which of the 5 arrays ni,nj is located
    // ! It means that we pass 0,1,2,3,4 on to the access function and select the pointer based on that.

    // ! We need to redefine the halos for the edges, because they twice as large

     // ilhe =   ilh
    int ilhe = a.ilh;
    //  nicsze = nicsz
    int nicsze = nicsz;
 // if (nx == 0) then
    if (nx == 0) {
     // ilhe = 2*  ilh
        ilhe = 2*a.ilh;
    //  nicsze = nicsz-  ilh+ilhe
        nicsze = nicsz-a.ilh+ilhe;
    }
    // end if

     // ihhe =   ihh
    int ihhe = a.ihh;
 // if (nx == dsmNX-1) then
    if (nx == dsmNX-1) {
     // ihhe = 2*  ihh
        ihhe = 2*a.ihh;
    //  nicsze = nicsz-  ihh+ihhe
        nicsze = nicsz-a.ihh+ihhe;
    }
    // end if

     // jlhe =   jlh
    int jlhe = a.jlh;
    //  njcsze = njcsz
    int njcsze = njcsz;
 // if (ny == 0) then
    if (ny == 0) {
     // jlhe = 2*  jlh
        jlhe = 2*a.jlh;
    //  njcsze = njcsz-  jlh+jlhe
        njcsze = njcsz-a.jlh+jlhe;
    }
    // end if

     // jhhe =   jhh
    int jhhe = a.jhh;
 // if (ny == dsmNY-1) then
    if (ny == dsmNY-1) {
     // jhhe = 2*  jhh
        jhhe = 2*a.jhh;
    //  njcsze = njcsz-  jhh+jhhe
        njcsze = njcsz-a.jhh+jhhe; 
    }
    // end if

    int in3, jn3;
    // ! Because all indices start at 0, this should be correct:
 // if (jn < jlhe) then

 //(0,1,0) 0;1;90 val:11
    if (jn < jlhe) {
    //   sel = 1
        *sel = 1;
    //     !top
    //  in3 = in
        in3 = in;
    //  jn3 = jn
        jn3 = jn;
    //   lin_idx = k-1-  koff+  klh+  ksz*(in3+nicsze*jn3)
        // *lin_idx = k-1-a.koff+a.klh+a.ksz*(in3+nicsze*jn3);
        lin_idx = nicsze*jlhe*(k-1-a.koff+a.klh)+nicsze*jn3+in3;
        // std::cout <<"SEL=1 ("<< i<<","<<j<<","<<k<<") "<<"("<< in<<","<<jn<<") "<<"("<< in3<<","<<jn3<<") "<<a.ksz<<";"<<nicsze << ";"<<*lin_idx  << "\n";        
    //else if (jn > njcsze-jhhe-1) then
    } else if (jn > njcsze-jhhe-1) {
    //   sel = 2
        *sel = 2;
    //  in3 = in
        in3 = in;
    //  jn3 = jn - njcsz + jhhe
        jn3 = jn - njcsze + jhhe;
    //   lin_idx = k-1-  koff+  klh+  ksz*(in3+nicsze*jn3)
    
        // *lin_idx = k-1-a.koff+a.klh+a.ksz*(in3+nicsze*jn3);
        lin_idx = nicsze*jhhe*(k-1-a.koff+a.klh)+nicsze*jn3+in3;
        // std::cout <<"SEL=2 ("<< i<<","<<j<<","<<k<<") "<<"("<< in<<","<<jn<<") "<<"("<< in3<<","<<jn3<<") "<<a.ksz<<";"<<nicsze << ";"<<njcsz<<";" <<*lin_idx  << "\n";
    //     !bottom
   // else
    } else {
    //  if (in < ilhe) then
        if (in < ilhe) {
    //       sel = 3
            *sel = 3;
    //         !left
    //      in3 = in
            in3 = in;
    //      jn3 = jn - jlhe
            jn3 = jn - jlhe;
        //   lin_idx = k-1-  koff+  klh+  ksz*(in3+ilhe*jn3)
            // *lin_idx = k-1-a.koff+a.klh+a.ksz*(in3+ilhe*jn3);
            lin_idx = (k-1-a.koff+a.klh)*ilhe*(njcsze-jlhe-jhhe)+(in3+ilhe*jn3);
    //    else if (in > nicsze-ihhe-1) then
        } else if (in > nicsze-ihhe-1) {
    //       sel = 4
            *sel = 4;
    //         !right
    //      in3 = in - nicsz + ihhe
            in3 = in - nicsze + ihhe;
    //      jn3 = jn - jlhe
            jn3 = jn - jlhe;
        //   lin_idx = k-1-  koff+  klh+  ksz*(in3+ihhe*jn3)
            // *lin_idx = k-1-a.koff+a.klh+a.ksz*(in3+ihhe*jn3);
            lin_idx = (k-1-a.koff+a.klh)*ihhe*(njcsze-jlhe-jhhe)+(in3+ihhe*jn3);
        // std::cout <<"SEL=4 ("<< i<<","<<j<<","<<k<<") "<<"("<< in<<","<<jn<<") "<<"("<< in3<<","<<jn3<<") "<<ihhe<<";"<<njcsze << ";"<<njcsz<<";" <<*lin_idx  << "\n";

    //    else
        } else {
    //       sel = 0
            *sel = 0;
    //         !core
    //      in3 = in - ilhe
            in3 = in - ilhe;
    //      jn3 = jn - jlhe
            jn3 = jn - jlhe;
        //   lin_idx = k-1-koff+klh+ksz*(in3+(nicsze-ilhe-ihhe)*jn3)
            // *lin_idx = k-1-a.koff+a.klh+a.ksz*(in3+(nicsze-ilhe-ihhe)*jn3);
            lin_idx = (k-1-a.koff+a.klh)*(nicsze-ilhe-ihhe)*(njcsze-jlhe-jhhe)+
            (in3+(nicsze-ilhe-ihhe)*jn3);
        }
    // end if 
    }
    // end if

    // The core and the top start at 0
    if (*sel==0 || *sel==1) {
        *comb_lin_idx=lin_idx;
    } else if (*sel==2) { // offset is the size of the top
        int64_t local_top_sz = nicsze*jlhe*a.ksz;
        *comb_lin_idx=local_top_sz+lin_idx;
    } else if (*sel==3) { // plus the size of the bottom
        int64_t local_top_sz = nicsze*jlhe*a.ksz;
        int64_t local_bottom_sz = local_top_sz + nicsze*jhhe*a.ksz;
        *comb_lin_idx=local_bottom_sz+lin_idx;
    } else { // plus the size of the left
        int64_t local_top_sz = nicsze*jlhe*a.ksz;
        int64_t local_bottom_sz = local_top_sz + nicsze*jhhe*a.ksz;
        int64_t local_left_sz = local_bottom_sz +a.ksz*ilhe*(njcsze-jlhe-jhhe);
        *comb_lin_idx=local_left_sz+lin_idx;
    }
    
}
#ifdef NEW_DSM3
void dsmWriteArrayH(DSM3DArrayH& a, int64_t i, int64_t j, int64_t k, float* val) {
    int64_t lin_idx = 0;
    int64_t part_idx = 0;
    int64_t sel = 0;
    indexing_case3(a,i,j,k,&part_idx, &sel, &lin_idx);
    int64_t dsm_local_ptr = sel==0 ? (int64_t)a.core_ptr : a.ptr[part_idx];
    dsmwritefloatarrayf_(&dsm_local_ptr, &lin_idx, val);
}

float dsmReadArrayH(DSM3DArrayH& a, int64_t i, int64_t j, int64_t k) {
    int64_t lin_idx = 0;
    int64_t part_idx = 0;
    int64_t sel = 0;
    indexing_case3(a,i,j,k,&part_idx, &sel, &lin_idx);
    int64_t dsm_local_ptr = sel==0 ? (int64_t)a.core_ptr : a.ptr[part_idx];
    float val = 0.0;
    dsmreadfloatarrayf_(&dsm_local_ptr, &lin_idx, &val);
    return val;
}
#else
void dsmWriteArrayH(DSM3DArrayH& a, int64_t i, int64_t j, int64_t k, float* val) {
    int64_t lin_idx = 0;
    int64_t part_idx = 0;
    int64_t sel = 0;
    indexing_case3(a,i,j,k,&part_idx, &sel, &lin_idx);
    int64_t dsm_local_ptr = a.ptr[part_idx][sel];
    // printf("Pidx: %ld, Lidx: %ld, Sel: %d, Val: %f\n", part_idx, lin_idx, sel, *val);
    dsmwritefloatarrayf_(&dsm_local_ptr, &lin_idx, val);
}

float dsmReadArrayH(DSM3DArrayH& a, int64_t i, int64_t j, int64_t k) {
    int64_t lin_idx = 0;
    int64_t part_idx = 0;
    int64_t sel = 0;
    indexing_case3(a,i,j,k,&part_idx, &sel, &lin_idx);
    // Using pre-fetched pointer:
    int64_t dsm_local_ptr = a.ptr[part_idx][sel];
    // Using collective pointer:
    /*
    int64_t dsm_local_ptr = 0;
    int64_t part_array_ptr = 0;
    dsmreadpointerarrayf_(a.coll_ptr, &part_idx, &part_array_ptr);
    int64_t part_ptr = 0;
    dsmreadpointerarrayf_(&part_array_ptr, &sel, &part_ptr);
    */
    float val = 0.0;
    dsmreadfloatarrayf_(&dsm_local_ptr, &lin_idx, &val);
    return val;
}
#endif

void indexing_case2(DSM3DArray& a, int i, int j, int k, int64_t* part_idx, int64_t* lin_idx) {
    // icsz = a%isz - a%ilh - a%ihh
    int icsz = a.isz - a.ilh - a.ihh;
    // jcsz = a%jsz - a%jlh - a%jhh
    int jcsz = a.jsz - a.jlh - a.jhh;

    // nicsz = icsz / dsmNX
    int nicsz = icsz / dsmNX;
    // njcsz = jcsz / dsmNY
    int njcsz = jcsz / dsmNY;

    // nx = (i-1-ioff) / nicsz

    int nx = i-1-a.ioff < dsmNX*nicsz ? (i-1-a.ioff) / nicsz : dsmNX-1;
    // printf("(i:%d-1-a.ioff:%ld) / nicsz:%d = nx:%d\n",i,a.ioff,nicsz,nx);
    // ny = (j-1-joff) / njcsz
    int ny = (j-1-a.joff) < dsmNY*njcsz ? (j-1-a.joff) / njcsz : dsmNY-1;
    // printf("(j:%d-1-a.joff:%ld) / njcsz:%d = ny:%d\n",j,a.joff,njcsz,ny);
    int in=-1; // for error checking at it should never be <0
    int jn=-1;

    // if (nx>0 .and. (nx < dsmNX-1 .or. (nx == dsmNX-1 .and. i<=icsz)) )
    if (nx > 0 && (nx < dsmNX-1 || (nx == dsmNX-1 && i <= icsz))) {
        // in = mod(i-1-ioff, nicsz)
        in = (i-1-a.ioff) % nicsz;
    // else if (nx == 0)
    } else if (nx == 0) {
        // in = i-1-ioff+ilh
        in = i-1-a.ioff+a.ilh;
    // else if (nx == dsmNX-1 .and. i>icsz)
    } else if (nx == dsmNX-1 && i > icsz) {
        // in = i-1-ioff-icsz+nicsz
        in = i-1-a.ioff-icsz+nicsz;
    }

    // if (ny>0 .and. (ny < dsmNY-1 .or. (ny == dsmNY-1 .and. j<=jcsz)) )
    if (ny > 0 && (ny < dsmNY-1 || (ny == dsmNY-1 && j <= jcsz))) {
        // jn = mod(j-1-joff, njcsz)
        jn = (j-1-a.joff) % njcsz;
    // else if (ny == 0)
    } else if (ny == 0) {
        // jn = j-1-joff+jlh
        jn = j-1-a.joff+a.jlh;
    // else if (ny == dsmNY-1 .and. j>jcsz)
    } else if (ny == dsmNY-1 && j > jcsz) {
        // jn = j-1-joff-jcsz+njcsz
        jn = j-1-a.joff-jcsz+njcsz;
    }

    // part_id = nx+ny*dsmNX
    *part_idx = nx+ny*dsmNX;

    // if(nx>0 .and. nx<dsmNX-1 )
    if (nx > 0 && nx < dsmNX-1) {
        // lin_idx = k-1-koff+klh+ksz*(in+nicsz*jn) ! for the core
        *lin_idx = (int64_t)(k-1-a.koff+a.klh+a.ksz)*(int64_t)(in+nicsz*jn);
    // else if (nx==0)
    } else if (nx == 0) {
        // lin_idx = k-1-koff+klh+ksz*(in+(nicsz+ilh)*jn) ! left
        *lin_idx = (int64_t)(k-1-a.koff+a.klh+a.ksz)*(int64_t)(in+(nicsz+a.ilh)*jn);
    } else {
        // lin_idx = k-1-koff+klh+ksz*(in+(nicsz+ihh)*jn) ! right
        *lin_idx = (int64_t)(k-1-a.koff+a.klh+a.ksz)*(int64_t)(in+(nicsz+a.ihh)*jn);
    }
    /*
    if (*lin_idx < 0 || *lin_idx >= (icsz+2)*(jcsz+2)*(km+2)) {
        printf("i: %d, j: %d, k: %d, dsmNX: %d, dsmNY: %d, \n", 
            i, j, k, dsmNX, dsmNY
            );
        printf("nicsz: %d, njcsz: %d, nx: %d, ny: %d, in: %d, jn: %d, lin_idx: %ld, part_idx: %ld\n",
            nicsz, njcsz, nx, ny, in, jn, *lin_idx, *part_idx
        );
    }
    */
}

// These are the routines to write to and read from a Case-2 3-D DSM array
void dsmWrite3DReal4Array(DSM3DArray& a, int64_t i, int64_t j, int64_t k, float* val) {
    int64_t lin_idx;
    int64_t part_idx;
    // Calculate lin_idx and part_idx
    indexing_case2(a,i,j,k,&part_idx,&lin_idx);
    // assert(part_idx >= 0);
    // assert(part_idx < dsmNX*dsmNY);
    // assert(lin_idx >= 0);
    // Get the pointer from the DSM3DArray
    /*
    DSMArrayPtr dsm_a_ptr = a.ptr;
	float** dsm_a_ptr_f = fromWord<float**>(*dsm_a_ptr);
	float* dsm_part_ptr_f = dsm_a_ptr_f[part_idx];
    dsm_part_ptr_f[lin_idx] = *val;
    */
    int64_t dsm_local_ptr = (int64_t)a.ptr[part_idx];
    // dsmreadpointerarrayf_(a.ptr, &part_idx, &dsm_local_ptr);
    // assert(dsm_local_ptr != 0);
    /*
    int64_t node_id;
    dsmgetnidf_(&node_id);
    if (node_id != part_idx) {
        printf("N[%ld] WTA, i: %ld, j: %ld, k: %ld, part_idx: %ld, lin_idx: %ld, ptr: %ld\n",
            node_id, i, j, k, part_idx, lin_idx, dsm_local_ptr
        );
    }
    */
    dsmwritefloatarrayf_(&dsm_local_ptr, &lin_idx, val);
    // argo::barrier();
    // dsmwritepointerarrayf_(a.ptr, &part_idx, &dsm_local_ptr);
    // argo::barrier();
    // printf("dsmWrite done\n");
}  

float dsmRead3DReal4Array(DSM3DArray& a, int64_t i, int64_t j, int64_t k) {
    // int64_t node_id;
    // dsmgetnidf_(&node_id);
    // int num_nodes;
    // dsmgetnnodesf_(&num_nodes);
    /*
    int64_t lin_idx;
    int64_t part_idx;
    // Calculate lin_idx and part_idx
    indexing_case2(a,i,j,k,&part_idx,&lin_idx);
    // Get the pointer from the DSM3DArray
    DSMArrayPtr dsm_a_ptr = a.ptr;
	float** dsm_a_ptr_f = fromWord<float**>(*dsm_a_ptr);
	float* dsm_part_ptr_f = dsm_a_ptr_f[part_idx];
    float val = dsm_part_ptr_f[lin_idx];
    printf("After array index retrieval2\n");
    return val;
    */
    int64_t lin_idx = 0;
    int64_t part_idx = 0;//num_nodes - 1 - node_id ;
    // Calculate lin_idx and part_idx
    indexing_case2(a,i,j,k,&part_idx,&lin_idx);
    // assert(part_idx >= 0);
    // assert(part_idx < dsmNX*dsmNY);
    // assert(lin_idx >= 0);
    // if (i==4 && j == 4 && k == 4) {
        // std::cerr << "NId "<<node_id<<": part_idx: " << part_idx << "\n";
    // printf("Indexing_case2: part idx: %ld, lin_idx: %ld\n", (long)part_idx, (long)lin_idx);
    // }
    // Get the pointer from the DSM3DArray
    int64_t dsm_local_ptr = (int64_t)a.ptr[part_idx];
    // assert(dsm_local_ptr != 0);
    // if (i==4 && j == 4 && k == 4) {
        // std::cerr<< "NId "<<node_id<< " dsm_local_ptr: "<< std::hex <<dsm_local_ptr << std::dec << "\n";
    // printf("Retrieved pointer: %ld\n", (long)dsm_local_ptr);
    // }
    /*
    int64_t node_id;
    dsmgetnidf_(&node_id);
    if (node_id != part_idx) {
        printf("N[%ld] RFA, i: %ld, j: %ld, k: %ld, part_idx: %ld, lin_idx: %ld, ptr: %ld\n",
            node_id, i, j, k, part_idx, lin_idx, dsm_local_ptr
        );
    }
    */
//    std::cerr << "NId "<<node_id<<": lin_idx: "<<lin_idx<<"\n";
    float val = 0.0;
    dsmreadfloatarrayf_(&dsm_local_ptr, &lin_idx, &val);
    // std::cerr<< "NId "<< node_id <<  "dsmRead done for "<< std::hex << dsm_local_ptr << std::dec << "lin_idx "<< lin_idx <<"\n";
    // printf("dsmRead done\n");
    return val;
}  




void indexing_case1(DSM3DArray& a, int* i,int* j,int* k,int* lin_idx) {

    // Properly handle halos and offsets so that i1,j1,k1 start at 0
    int k1 = *k-1-a.koff+a.klh;
    int i1 = *i-1-a.ioff+a.ilh;
    int j1 = *j-1-a.joff+a.jlh;

    *lin_idx = k1+a.ksz*(i1+a.isz*j1);
    
} // end subroutine indexing_case1
