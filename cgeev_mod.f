      module cgeev_mod
      use iso_c_binding
!      interface 
!         subroutine ftoctest_wr(n,a_r,a_i) bind(c)
!      use iso_fortran_env
!         integer:: n
!         real(real128):: a_r(n,n),a_i(n,n)
!         end subroutine
!      end interface 
      contains
      subroutine cgeev_wrapper(jobvl,jobvr,n,a,w,vl,ldvl,vr,ldvr
     &,info) !bind(c)
      use iso_fortran_env
         implicit none
         !character(kind=c_char,len=1):: jobvl,jobvr
         character(len=1):: jobvl,jobvr
         integer:: n,ldvl,ldvr,info,i,j
         complex(real128)::a(n,n),w(n),vl(n,n),vr(n,n)
         real(real128)::a_r(n,n),a_i(n,n),vr_r(n*n),vr_i(n*n),
     &  w_r(n),w_i(n)
         a_r(:,:)=real(a(:,:))
         a_i(:,:)=aimag(a(:,:))
         !convert from fortran complex to c qd_complex

         !call cgeev_wrapper_c()
         call cgeev_wr(n,a_r,a_i,w_r,w_i,vr_r,vr_i,info)
         do i=1,n
            w(i)=complex(w_r(i),w_i(i))
         enddo
         do i=1,n
            do j=1,n
               vr(i,j)= complex(vr_r(i+n*(j-1)),vr_i(i+(j-1)*n))
            enddo
         enddo
         endsubroutine 
       end module
