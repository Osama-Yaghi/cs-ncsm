#include "cgeev_wrapper.h"
#include "mplapack__Float128.h"
#include "mpblas__Float128.h"
#include <stdio.h>
#include <cstring>
#include <algorithm>
#include <complex>
using namespace std;
    // Call the Cgeev function from mplapack_qd
//    Cgeev(jobvl, jobvr, *n, a, *lda, w, vl, *ldvl, vr, *ldvr, work, *lwork, rwork, info);
extern "C" void cgeev_wr_(const int* l,_Float128* a_r,_Float128* a_i,_Float128* w_r,_Float128*w_i,_Float128* vr_r,_Float128* vr_i, int* ierr) {
   cgeev_wrapper_c(l,a_r,a_i,w_r,w_i,vr_r,vr_i,ierr);
  return;
}
