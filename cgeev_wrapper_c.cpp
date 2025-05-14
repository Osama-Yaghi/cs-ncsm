//#include"iostream"
#include "mplapack__Float128.h"
#include "mpblas__Float128.h"
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#define BUFLEN 1024
void printnum(_Float128 rtmp)
{
    int width = 42;
    char buf[BUFLEN];
    strfromf128(buf, sizeof(buf), "%.35e", rtmp);
    if (rtmp >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
    return;
}
void printnum(std::complex<_Float128> rtmp)
{
    int width = 42;
    char buf[BUFLEN];
    strfromf128(buf, sizeof(buf), "%.35e", rtmp.real());
    if (rtmp.real() >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
    strfromf128(buf, sizeof(buf), "%.35e", rtmp.imag());
    if (rtmp.imag() >= 0.0)
        printf ("+%si", buf);
    else
        printf ("%si", buf);
    return;
}
template <class X> void printvec(X *a, int len) {
    X tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}
void cgeev_wrapper_c(const int* l,_Float128* a_r,_Float128* a_i,_Float128* w_r,_Float128*w_i,_Float128* vr_r,_Float128* vr_i,int* ierr) {
  int i,j,n,ind;
  n=*l;
    std::complex<_Float128> *vl = new std::complex<_Float128>[n * n];
    std::complex<_Float128> *w = new std::complex<_Float128>[n];
    std::complex<_Float128> *vr = new std::complex<_Float128>[n * n];
    mplapackint lwork = 4 * n;
    std::complex<_Float128> *work = new std::complex<_Float128>[lwork];    
    _Float128 *rwork = new _Float128[lwork];
    mplapackint info;
/*  for(i=0;i<4;i++){
     for (j=0;j<4;j++)
     {
        printnum(a_r[i*n+j]);
              }
              printf("\n");
  }
  for(i=0;i<4;i++){
     for (j=0;j<4;j++)
     {
        printnum(a_i[i*n+j]);
              }
              printf("\n");
  }
*/
  std::complex<_Float128>* a=new std::complex<_Float128>[n * n];
  for(i=0;i<n;i++)
     for(j=0;j<n;j++){
        ind=i*n+j;
        a[ind]=std::complex<_Float128>(a_r[ind],a_i[ind]);
              }
    Cgeev("N", "V", n, a, n, w, vl, n, vr, n, work, lwork, rwork, info);
  for(i=0;i<n;i++) {
     w_r[i]=w[i].real();
     w_i[i]=w[i].imag();
  }
  for(i=0;i<n;i++)
     for(j=0;j<n;j++){
        ind=i*n+j;
        vr_r[ind]=vr[ind].real();
        vr_i[ind]=vr[ind].imag();
     }
    *ierr=info;
//    printf("lambda ="); printvec(w,n); printf("\n");    
    delete[] rwork;
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] w;
    delete[] a;
 return;
}
