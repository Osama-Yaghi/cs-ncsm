      module HO_wf
       contains

      subroutine laguerre_rootfinder_kd(ngr,l,roots)
      use iso_fortran_env
         implicit none
         integer,intent(in) :: ngr
         integer,parameter :: kd=real128
         real(kd),intent(in) :: l
         real(kd),intent(inout) :: roots(ngr,ngr)
         integer :: n,i,j,m,np,nm,np_nm,n1
         real(kd) :: x1,x2,y1,t1,y2,t2,x3,y3,
     &                                      t3
         real(kd),allocatable ::zero(:)
         character(len=4) :: string1,string2
         character(len=80) :: file_name='LaguerrePoly_'
         logical :: filetest=.false.

         if(ngr==0) return
         write(string2,'(i4)') int(l)
         string2=adjustl(string2)
         do i=1500,ngr+1,-1 
            write(string1,'(i4)') i
            string1=adjustl(string1)
            inquire(file=trim(file_name)//trim(string1)//'_leq'//
     &             trim(string2)//'_Roots_kd.bin',exist=filetest)
            if(filetest) exit
         enddo
         if(filetest) then
            open(3,file=trim(file_name)//trim(string1)//'_leq'//
     &          trim(string2)//'_Roots_kd.bin',status='old'
     &,action='read')
            roots=0._kd
            do n=1,ngr
               read(3,*) roots(:n,n)
            enddo
            close(3)
            print*,"Laguerre roots read from files"
            return
         else
            print*,"calculating Laguerre roots"
            open(4,file=trim(file_name)//trim(string1)//'_leq'//
     &          trim(string2)//'_Roots_kd.bin',status='new'
     &,action='write')
         endif

         roots=0._kd
         if(ngr==1) then
            roots(1,1)=1._kd+l
            return
         endif
         if(allocated(zero)) deallocate(zero)
         allocate(zero(ngr+1))
         zero(1)=2._kd+l
! we iteratively look for the roots of the n-order laguerre polynomial
         do n=2,ngr
            zero(n+1)=zero(n-1)*2._kd
!            zero(:n )=eoshift(zero(:n),-1,boundary=1.d0/real(n,kd))
            do m=n,1,-1
               if(m>1) then
                  zero(m)=zero(m-1)
               else
                  zero(m)=1._kd/real(n,kd)
               endif
            end do
! roots finding loop
            do m=1,n
               np=0
               nm=0
               n1=2
               x1=zero(m)
               x2=zero(m+1)
               call glaguer_poly_kd(n,l,x1,y1,t1)
               call glaguer_poly_kd(n,l,x2,y2,t2)
               if (y1*y2.gt.0._kd) then
                  print*,'error in laguerre algo'
                  stop
               endif
               do while(abs(x2-x1)>10._kd**-30._kd)
                  select case (np+nm)
                  case(:4)
                     x3=x1+(x2-x1)*y1/(y1-y2)
                  case default
                     x3=0.5_kd*(x1+x2)
                  end select
                  if ((n1.lt.9).and.(x3.eq.x1)) x3=0.2_kd*(4._kd*x1+x2)
                  call glaguer_poly_kd(n,l,x3,y3,t3)
                  n1=n1+1
                  if (y3.eq.0._kd) exit
                  select case(y1*y3.lt.0._kd)
                  case (.false.)
                     if (x1.eq.x3) exit
                     nm=nm+1
                     np=0
                     x1=x3
                     y1=y3
                  case (.true.)
                     if (x2.eq.x3) exit
                     np=np+1
                     nm=0
                     x2=x3
                     y2=y3
                  end select
               end do
               zero(m)=x3
            end do
            roots(:n,n)=zero(:n)
         end do
         if(ngr>1) then
            roots(:ngr,ngr)=zero(:ngr)
         endif
         do n=1,ngr
            write(4,*) roots(:n,n)
         end do
         close(4)
         return
      end subroutine laguerre_rootfinder_kd

      subroutine glaguer_poly_kd(n,l,x,y,dy)
! y is the n-th order (l) laguerre polynomial at x and dy its derivative
! w is integration weight
      use iso_fortran_env
         implicit none
         integer,intent(in) :: n
         integer,parameter :: kd=real128
         real(kd),intent(in) :: l
         real(kd),intent(in) :: x
         real(kd),intent(out) :: y,dy
         integer :: k
         real(kd) :: ln2,ln1

         ln2 =1._kd
         ln1 =1._kd+l-x
         do k=2,n
            y =((real(2*k-1,kd)+l-x)*ln1-(real(k-1,kd)+l)*ln2)
     &         /real(k,kd)
            ln2 =ln1
            ln1 =y
         end do
         dy=(real(n,kd)*ln1-(real(n,kd)+l)*ln2)/x
         return
      end subroutine glaguer_poly_kd




!      subroutine wf_coulomb_old(n,l,anu,r,wave,roots)
!         use iso_fortran_env
!         use paramdef
!         use gamalog
!         implicit real(real128) (a-h,o-z)
!!         implicit double precision (a-h,o-z)
!         integer,intent(in) :: n,l
!         integer,parameter :: kd=real128
!         real(kd),intent(in) :: r,anu
!!        the radial part wave function (wave=rr(r)) int wave(r)**2 dr =1)
!!        anu is the size parameter of nuclear well
!!        anu=mass times omega over h bar
!!         integer,parameter :: kd=kind(0.d0)
!         real(   kd),intent(in) :: roots(n)
!         real(   kd) :: logdfacq
!         real(kd) :: aq,bq,guerpq
!         real(kd) :: zz,z,wavel,wave,a,b,guerp
!
!         dlanu=log(anu)
!         wave =0._kd
!         zz=r
!         z =sqrt(r/anu)
!         wavel=0.25_kd*dlanu
!     &        +(real(l,kd)+1._kd)*(0.5d0*dlanu+log(z))
!         if (n==0) then
!            guerp=exp(0.5d0*(log(2.d0)-gamal(2*l+3)))
!            wave =exp(wavel)*guerp
!         elseif (n==1) then
!            guerp=exp(0.5d0*(log(2.d0)-gamal(2*l+5)))
!     &                     *(real(l,kd)+1.5d0-zz)
!            wave =exp(wavel)*guerp
!         else
!            guerp=1._kd
!            do nnn=1,n
!               guerp=guerp*(zz-roots(nnn))
!            end do
!            if(mod(n,2)==0) then
!               wave= exp(wavel+0.5_kd*(log(2._kd)-gamal(2*l+2*n+3)
!     &                  -logdfacq(n)))*guerp
!            else
!               wave=-exp(wavel+0.5_kd*(log(2._kd)-gamal(2*l+2*n+3)
!     &                  -logdfacq(n)))*guerp
!            endif
!         endif
!      end subroutine wf_coulomb_old

      subroutine wf_coulomb(n,l,anu,r,wave,roots)
         use iso_fortran_env
         use paramdef
         use gamalog
         implicit real(real128) (a-h,o-z)
!         implicit double precision (a-h,o-z)
         integer,intent(in) :: n,l
         integer,parameter :: kd=real128
         real(kd),intent(in) :: r,anu
!        the radial part wave function (wave=rr(r)) int wave(r)**2 dr =1)
!        anu is the size parameter of nuclear well
!        anu=mass times omega over h bar
!         integer,parameter :: kd=kind(0.d0)
         real(   kd),intent(in) :: roots(max(n,1))
         real(kd) :: aq,bq,guerpq
         real(kd) :: zz,z,wavel,wave,a,b,guerp,log_wave,sgn

         dlanu=log(anu)
         wave =0._kd

         zz=r*2
         z =sqrt(2*r/anu)
         wavel=0.25_kd*dlanu
     &        +(real(l,kd)+1._kd)*(0.5d0*dlanu+log(z))
     & +(r-zz)/2._kd
         if (n==0) then
            guerp=exp(0.5d0*(log(2.d0)-gamal(2*l+3)))
            wave =exp(wavel)*guerp
         elseif (n==1) then
            guerp=exp(0.5d0*(log(2.d0)-gamal(2*l+5)))
     &                     *(real(l,kd)+1.5d0-zz)
            wave =exp(wavel)*guerp
         else
            guerp=0._kd
            sgn=1._kd
            do nnn=1,n
               guerp=guerp+log(abs(zz-roots(nnn)))
               if(zz<roots(nnn)) sgn=sgn*-1._kd
            end do
            log_wave=wavel+0.5_kd*(log(2._kd)-gamal(2*l+2*n+3)
     &                  -logdfacq(n))+guerp
            if(mod(n,2)==0) then
               wave= sgn*exp(log_wave)
            else
               wave=-sgn*exp(log_wave)
            endif
         endif
         contains
!**************
      real(real128) function logdfacq(n)
         use iso_fortran_env
         implicit none
         integer,intent(in) :: n
         integer :: i
         integer,save :: i_save=1
         integer,parameter :: kd=real128
         real(kd),save :: log_fact(0:1000)
         log_fact(0)=0._kd
         log_fact(1)=0._kd
         if(n<=i_save) then
            logdfacq=log_fact(n)
         elseif(n<=1000) then
            do i=i_save+1,n
               log_fact(i)=log_fact(i-1)+log(real(i,kd))
            end do
            i_save =n
            logdfacq=log_fact(n)
         else
            print*,'***error: log factorial function out-of-range'
            stop
         end if
      end function logdfacq
! *************
      end subroutine wf_coulomb

      subroutine wave_func_c(n,l,anu,r,wave,roots,laguerr,theta_in)
      use paramdef
      use gamalog
      use iso_fortran_env
      implicit real(real128) (a-h,o-z)
!      implicit double precision (a-h,o-z)
      integer,intent(in) :: n,l
      double precision,intent(in),optional:: theta_in
      integer,parameter :: kd=real128
      real(kd),intent(in) :: r,anu
      complex(kd),intent(inout) :: wave
      logical,optional:: laguerr
!     the radial part wave function (wave=rr(r)) int wave(r)**2 dr =1)
!     anu is the size parameter of nuclear well
!     1/anu=mass times omega over h bar
!      integer,parameter :: kd=kind(0.d0)
      real(   kd),intent(in) :: roots(max(n,1))
      real(   kd) :: rr,theta_t
      complex(kd) :: aq,bq,guerpq,lag_shift
      complex(kd) :: zz,z,wavel,a,b,guerp,waved
      complex(kd):: a1,a2,a3,a4,a5,a6,a7
      complex(kd):: expo
      integer:: nnn
      if(present(theta_in)) then
         theta_t=theta_in
      else
         theta_t=theta
      endif
      wave=complex(0._kd,0._kd)
      z=r*complex(cos(theta_t),sin(theta_t))
      zz=z*z*anu ! zz= hbar/(omega*m)*z^2
!      zz=r/cos(2._kd*theta)/2.d0
!     &   *cmplx(cos(2._kd*theta),sin(2._kd*theta))
!     ! rr=r
!      z =sqrt(r/2.d0/anu/cos(2._kd*theta))
!     &   *cmplx(cos(      theta),sin(      theta))
! from the addition of the exp(+ r)exp(-r) to the integrand to fit laguerre quadrat
! -ure:
!        integral ( f(x) )= integral( exp(-x)exp(+x)f(x))=sum[ w_i f(x_i)exp(+x_i)]
      wavel=(anu**(0.25_kd))*sqrt(2._kd)
      a5=acos(-1._kd)**(0.25_kd)
      expo=-zz/2._kd
      lag_shift=real(r*r*(anu*cos(2*theta_t)/4._kd),kd)
      if(laguerr) expo=expo+lag_shift

      do nnn=1,2*l+2*n+1,2
         a5=a5*sqrt(real(nnn,kd)/2)
      enddo
!      print*,"wavel=",wavel
!      print*,"gamal",a5
      if (n==0) then
         wave =wavel*(sqrt(anu)*z)**(l+1._kd)*exp(expo)/a5
      elseif (n==1) then
         guerp=(real(l,kd)+1.5_kd-zz)
         wave =wavel*(sqrt(anu)*z)**(l+1._kd)*exp(expo)/a5*guerp
      else
         guerp=(sqrt(anu)*z)**(l+1._kd)
         a2=exp(expo/n)
         do nnn=1,n
            a1=zz-roots(nnn)
            a3=a1*a2/sqrt(real(nnn,kd))
            guerp=guerp*a3
 !           print*,"step ",nnn,a1,a2,a3
         end do
         wave=wavel*guerp/a5*(-1)**n
      endif

      end subroutine wave_func_c
      end module HO_wf
!***************************
