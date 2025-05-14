      module v2b_quadrature
      use iso_fortran_env
      contains

!**********************
      subroutine inter_2b_quadrature(ipotin)
      use constants
      !$ use OMP_LIB
      use paramdef
      use harmosc
      use SRG_module
      use iso_fortran_env
      use nninteraction
      use pless_eft_param, only: coulomb_switch
      !use trep_integration
      use omp_lib
      use H_matrices
      use quadratures
      use HO_wf
      implicit none
      integer,parameter :: kd=real128
      integer,parameter :: dp=kind(0.d0)
!*      parameter (nstep=500,pc=30.d0,pinf=0.000000001)
      complex(kd) :: vnp(6),vpp(6),vnn(6),sumv(6)  
      !file variables
      character(len=5):: string1,string2,string3,string4
     & ,string5,string6
      character(len=100):: v_file_name
      logical:: filetest=.false.
      !timer variables
      real:: time1,time2
      integer:: ir,time3,time4,time5,time6
      !parameters 
      integer:: Nf,Nmax,Nd,v_n,n_threads
      real(kd):: fmbyh2,hbc3,pfmev,rdmavg,bosc,th_kd,prange
      !temporary variables 
      integer:: i,k,nq,ia,ib,ind,nra,nrb,lra,lrb,jrel
      real(kd) ::pi,pk,rrel
      real(dp)::pii,pkk
      real(kd)::e_kd
      complex(kd),allocatable,dimension(:,:,:,:,:):: H0
      complex(kd),allocatable,dimension(:,:,:):: vint

      integer:: ipotin,icceff,inn,inp,ipp,ispin,iuc,iucmax,lbmax
      integer::lramin,mtot2,nrbmax,nna,nreffmax,sra
      double precision::cnp,cpp,cnn
      
      !quadrature 
      real(kd),allocatable::mesh_range(:),m_temp(:,:,:),mesh_we(:)
     & ,v_x_mesh(:)
      real(kd)::theta2,interr,regcut,regpow,temp_kd
      complex(kd),allocatable::u4(:,:,:)
      complex(kd):: wave,integ,reg,vrel
      integer,allocatable:: N_arr(:,:)
      integer::mesh_s,NN,totn,qn
      real(kd),allocatable:: gau_leg_poi(:)
     & ,gau_leg_we(:)
      real(kd),allocatable:: roots(:,:,:)
      real(real128),allocatable:: mesh(:,:,:),step_arr(:,:,:)
      real(kd):: v_scale,v_theta

      ! eigen solver variables
      complex(kd),allocatable:: Heff_kd(:,:),eigr_kd(:),eigv_kd(:,:)
     &, Heff_kd2(:,:)
      complex(kind(0.d0)),allocatable::Heff(:,:),eigr(:),eigv(:,:)

      
      complex(kd),allocatable,dimension(:,:,:,:):: sumvi
!----------------------------------------------------------------------------
      !common block for n3lo subroutine
      real(dp):: KREAD,KWRITE,KPUNCH,XMEV,YMEV
     &,ipot,V(6),KDA(9)
      LOGICAL HEFORM,SING,TRIP,COUP,ENDEP
      character*4 label
      integer::J
!----------------------------------------------------------------------------
!        THE USER NEEDS TO BE FAMILIAR WITH THE FOLLOWING
!        THREE COMMON BLOCKS:
!
      COMMON /CRDWRT/ KREAD,KWRITE,KPUNCH,KDA
!
!        ARGUMENTS AND VALUES OF CDBONN:
      COMMON /CPOT/   V,XMEV,YMEV
      COMMON /CSTATE/ J,HEFORM,SING,TRIP,COUP,ENDEP,LABEL
      common /cnn/ ipot
      interface
         subroutine mkl_zgeev(A,n,m,W,Z)
            integer, intent(in) :: n,m
            complex(kind(0.d0)),intent(in) :: A(n,n)
            complex(kind(0.d0)), intent(out) :: W(n),Z(n,n)
         end subroutine mkl_zgeev
      end interface
!----------------------------------------------------------------------------
!
!        ONCE FOR EVER, SET THE FOLLOWING 6 PARAMETERS AND FORGET THEM:
!        --------------------------------------------------------------
      KREAD=5
      KWRITE=6
      HEFORM=.FALSE.
      SING=.TRUE.
      TRIP=.TRUE.
      COUP=.TRUE.
!        IGNORE KPUNCH, KDA(9), ENDEP, and LABEL.
!        ----------------------------------------
      open(223,file="inter_2b_fit.out")
      call system_clock(count=time5, count_rate=ir)
!***********************
      !initialize variables
      ipot=ipotin
      prange=6.00
      Nf=9
      v_n=100
      v_scale=1._kd
      v_theta=0._kd
      e_kd=exp(1._kd)
      ipot=ipotin
      if (ipot==1) then
         ipp=1
         cpp=1.d0
         print *,' ipp=',ipp
         print *,' cpp=',cpp
      else
         ipp=0
      endif
      rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
      fmbyh2=rdmavg/hbc**2
      anu_kd=fmbyh2*hbo ! anu=(m*h_bar*omega)/ (h_bar*c)^2
      pfmev=hbc/dsqrt(2.d0)
      bsquare_kd=1.d0/anu_kd
      bosc=sqrt(real(bsquare_kd))
      hbc3=hbc*hbc*hbc/dsqrt(8.d0)
      th_kd=theta
!***********************
      !printing a description
      print*,"called inter_2b_fit"
      print*,"$$$$$ theta=",th_kd
!$OMP parallel
!$       n_threads=OMP_GET_MAX_THREADS()
!$OMP end parallel
      print*,"The number of threads being used with Open MP:"
     &                                            ,n_threads
      print *,'  iucdim=',iucdim

      if (icutnum==100) then
         write(2,*)
     &        'N3LO NN potential in momentum space'
      elseif (icutnum==102) then
         write(2,*)
     &        'NNLO(POUNDerS) NN potential in momentum space'
      elseif (icutnum==105) then
         write(2,*)
     &        'NNLOsat NN potential in momentum space'
      elseif (icutnum==106) then
         write(2,*)
     &        'N4LO500 NN potential in momentum space'
      elseif (icutnum==107) then
         write(2,*)
     &        'N3LO500new NN potential in momentum space'
      elseif (icutnum==108) then
         write(2,*)
     &      'N2LO500 NN potential in momentum space'
      elseif (icutnum==109) then
         write(2,*)
     &        'NLO500 NN potential in momentum space'
      elseif (icutnum==110) then
         write(2,*)
     &        'LO500 NN potential in momentum space'
      elseif (icutnum==115) then
         write(2,*)
     &        'chiral NN potential in momentum space from inputs'
         if(.not.coulomb_switch.and.ipp/=0) ipp=0
      endif

      if (ipp/=0) then
      endif
!*****************************
      theta2=theta
      !theta=0
      allocate(gau_leg_poi(100))
      allocate(gau_leg_we(100))
      if (allocated(roots)) deallocate(roots)
      allocate(roots(nrelm,nrelm,0:lrelm))
      do lra=0,lrelm
         call laguerre_rootfinder_kd(nrelm,lra+0.5_kd
     &                                 ,roots(1:nrelm,1:nrelm,lra))
      end do
      print*,"laguerre roots allocated"
      allocate(m_temp(0:nrelm,0:lrelm,50))
      m_temp=0
      NN=0
      allocate(mesh(0:5000,0:nrelm,0:lrelm),N_arr(0:nrelm,0:lrelm))
      allocate(step_arr(0:5000,0:nrelm,0:lrelm))
      call gauleg2(50,0._kd,1._kd,gau_leg_poi,gau_leg_we)
         call cpu_time(time1)
         call system_clock(count=time3, count_rate=ir)
!$OMP parallel do
!$OMP& private(nra,lra,totn,interr,sumv)
!$OMP& schedule(dynamic)
      do nra=0,nrelm
         do lra=0,lrelm
            totn=0
            interr=0
            wave=adapt_integ_leg(nra,nra,lra,0._kd,
     &   15._kd,-100._kd,50,totn,interr)
            N_arr(nra,lra)=totn
            NN=max(NN,totn)
          if((mod(nra,20)==0 .or. nra==10) .and. lra<3) then
            print*,"nra,:",nra," integral: ", real(wave,8),
     &          " points ",totn,"error",real(interr,8)
          endif
         enddo
      enddo
!$OMP end parallel do
         call cpu_time(time2)
         call system_clock(count=time4, count_rate=ir)
      print*,'time', 
     &      real(time4 - time3,kind=8)/real(ir,kind=8),time2-time1
      print*,"NN= ",NN
      allocate(mesh_range(120))
      mesh_range=0
      mesh_s=0
      do nra=0,nrelm
         do lra=0,2
            do k=1,N_arr(nra,lra)/50+1
         pi=m_temp(nra,lra,k)
         i=1
         do while( mesh_range(i)<pi .and. i<=mesh_s)
            i=i+1
         enddo
         if(i<=mesh_s .and. pi==mesh_range(i)) cycle
         do j=mesh_s+1,i+1,-1
            mesh_range(j)=mesh_range(j-1)
         enddo
         mesh_range(i)=pi
         mesh_s=mesh_s+1
            enddo
         enddo
      enddo
      print*,mesh_s
      do i=1,mesh_s
         print*,mesh_range(i)
      enddo
      allocate(v_x_mesh(50*(mesh_s-1)),mesh_we(50*(mesh_s-1)))
      do i=1,mesh_s-1
         temp_kd=mesh_range(i+1)-mesh_range(i)
         do k=1,50
            v_x_mesh((i-1)*50+k)=mesh_range(i)+gau_leg_poi(k)
     &             *temp_kd
            mesh_we((i-1)*50+k)=temp_kd*gau_leg_we(k)
         enddo
      enddo
         qn=(mesh_s-1)*50
         print*,"quadrature points:",qn

!
!      qn=2000
!      deallocate(gau_leg_poi,gau_leg_we,v_x_mesh,mesh_we)
!      allocate(gau_leg_poi(qn),v_x_mesh(qn),mesh_we(qn))
!      allocate(gau_leg_we(qn))
!      call gauleg2(qn,0._kd,15._kd,gau_leg_poi,gau_leg_we)
!      mesh_s=2
!      v_x_mesh(:)=gau_leg_poi
!      mesh_we=gau_leg_we
!
         theta=theta2
         allocate(u4(0:qn,0:nrelm,0:lrelm))
         print*,"bsquare", bsquare_kd
         do lra=0,lrelm
            do nra=0,nrelm
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,pi,wave)
!$OMP& SCHEDULE(DYNAMIC)
               do i=1,qn
!                  pp=gau_leg_poi(i)
                  pi=v_x_mesh(i)
                  !call waver(nr,lr,bsquare,pp,wave)
                  call wave_func_c(nra,lra,bsquare_kd,pi,wave
     &,roots(:nra,nra,lra),.false.)
                  u4(i,nra,lra)=wave*dble((-1)**nra)
     & *complex(cos(theta/2),sin(theta/2))
               end do
!$OMP end parallel do
            end do
         end do
!**************************
!**************************
         theta=theta2
         print*,"wave function created"
!*****************************
      print*,"hbc3",hbc3
      print*,"bsquare", bsquare_kd
!***********************
!calculating v_n3lo
      allocate(vint(0:qn,0:qn,6))
      allocate(sumvi(0:nrelm,0:2,0:qn,6))
      allocate(H0(0:nrelm,0:lrelm,0:nrelm,0:lrelm,6))
      do jrel=0,min(jrelm,mxeng+1)
         J=jrel
         print *,' doing J=',J

         do i=1,qn
            pi=v_x_mesh(i)
            xmev=pi*pfmev
!**            do k=0,nstep
            do k=1,qn
!**               pk=pinf+k*pstep
               pk=v_x_mesh(k)
               ymev=pk*pfmev
               if (pk<=25.d0.and.pk>1.d-9.and.
     &             pi<=25.d0.and.pi>1.d-9     ) then
               if (icutnum==100) then
                  call n3lo
               elseif (icutnum==102) then
                  call nnlo
               elseif (icutnum==105) then
                  call nnlo_sat
               elseif (icutnum==106) then
                  call n4lo500
               elseif (icutnum==107) then
                  call n3lo500new
               elseif (icutnum==108) then
                  call n2lo500
               elseif (icutnum==109) then
                  call nlo500
               elseif (icutnum==110) then
                  call lo500
               elseif (icutnum==115) then
                  call nnlo_input
               endif
               vint(i,k,:)=v(:)
!               vint(i,k,:)=-20._kd*exp(
!     &         -(pi**2._kd+pk**2._kd)/2._kd)/hbc3
!               pii=pi
!               pkk=pk
!               vint(i,k,:)=-20._rp*exp(
!     &         -(pii**2._rp+pkk**2._rp)/2._rp)/hbc3
!     &    *exp(-(pi/6._kd)**7._kd -(pk/6._kd)**7._kd) ! regulator
               else
               vint(i,k,:)=0.d0
               endif
               vint(i,k,:)=vint(i,k,:)
            end do
         end do

         print *,' potential calculated'
!*****************************************
       
       sumvi=complex(0._kd,0._kd)
       H0=complex(0._kd,0._kd)
       lramin=max(0,jrel-1)
!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(nra,lra,i,k,pk)
!$OMP& SCHEDULE(DYNAMIC)
         do nra=0,nrelm
            do lra=lramin,jrel+1
               if (2*nra+lra>2*nrelm) cycle
               do i=1,qn
                  do k=1,qn
                     pk=v_x_mesh(k)
                     sumvi(nra,lra-lramin,i,:)=
     &               sumvi(nra,lra-lramin,i,:)+mesh_we(k)*vint(i,k,:)
     &                 *u4(k,nra,lra)*pk
                  end do
               end do
            end do
         end do
!$OMP END PARALLEL DO
         print *,' initial integral calculated'

!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(nra,lra,nrb,lrb,sumv,i,pi)
!$OMP& SCHEDULE(DYNAMIC)
         do nra=0,nrelm
            do lra=lramin,jrel+1
               if (2*nra+lra>2*nrelm) cycle
               do nrb=nra,nrelm
                  do lrb=lramin,jrel+1
                     if (2*nrb+lrb>2*nrelm) cycle
                     if ((-1)**(lra+lrb)/=1) cycle

                     sumv=complex(0._kd,0._kd)
                     do i=1,qn
                        pi=v_x_mesh(i)
                        sumv(:)=sumv(:)
     &                    +sumvi(nrb,lrb-lramin,i,:)
     &                    *mesh_we(i)*u4(i,nra,lra)*pi
                     end do
                     sumv(:)=sumv(:)*real(hbc3,kd)
                     H0(nra,lra,nrb,lrb,:)=sumv(:)
                     ! fill the second half of H0
                     H0(nrb,lrb,nra,lra,:)=sumv(:)
                     H0(nrb,lrb,nra,lra,6)=sumv(5)
                     H0(nrb,lrb,nra,lra,5)=sumv(6)
                  end do
               end do
            end do
         end do
!$OMP END PARALLEL DO
         print *,' matrix elements calculated'
         do nra=0,nrelm
            do nrb=0,nrelm
               do lra=max(0,jrel-1),jrel+1
                  do lrb=max(0,lra-2),min(lra+2,jrel+1)
                     if((-1)**(lra+lrb)/=1) cycle
                     if(2*nra+lra>2*nrelm) cycle
                     if(2*nrb+lrb>2*nrelm) cycle

                    sumv(:)=H0(nra,lra,nrb,lrb,:)
                    ia=nra*2+(lra-jrel+1)/2
                    ib=nrb*2+(lrb-jrel+1)/2
                    !****************
                     if (jrel==0) then
                        if (lra==0) then
                           v_uc(nra,nrb,1)=v_uc(nra,nrb,1)
     &                                     +sumv(1)
                        elseif(lra==1) then
                           v_uc(nra,nrb,2)=v_uc(nra,nrb,2)
     &                                     +sumv(3)
                        endif
                     else
                        if (lra==jrel.and.lrb==jrel) then
                           v_uc(nra,nrb,2*jrel+1)=v_uc(nra,nrb,2*jrel+1
     &                                            )+sumv(1)
                           v_uc(nra,nrb,2*jrel+2)=v_uc(nra,nrb,2*jrel+2
     &                                           )+sumv(2)
                        elseif (lra==jrel+1.and.lrb==jrel+1) then
                           v_cc(ia,ib,jrel)=v_cc(ia,ib,jrel) 
     &                                      +sumv(3)
                        elseif (lra==jrel-1.and.lrb==jrel-1) then
                           v_cc(ia,ib,jrel)=v_cc(ia,ib,jrel)
     &                                      +sumv(4)
                        elseif (lra==jrel+1.and.lrb==jrel-1) then
                           v_cc(ia,ib,jrel)=v_cc(ia,ib,jrel)
     &                                      +sumv(5)
                        elseif (lra==jrel-1.and.lrb==jrel+1) then
                           v_cc(ia,ib,jrel)=v_cc(ia,ib,jrel)
     &                                      +sumv(6)
                        endif
                     endif
                    !****************

                  enddo
               enddo
            enddo
         enddo
      

      end do



!**********************************************
      close(223)
      return


!      stop
! end write
      close(223)

!************
      contains

      recursive function adapt_integ_leg(nra,nrb,lr,
     &    p_b,p_e,area,n,totn,interr) result(res)
         implicit none
         !integer,parameter:: kd=real128
         real(kd)::p_b,p_e,area,interr
         integer::n,nra,nrb,lr,totn

         real(kd):: res,tol,u1,u2
         real(kd)::p_m,area1,area2,step,pi,cf,tmp2,length,vrel
         complex(kd)::wave
         integer:: i,reste,tmp1
         area1=0._kd
         area2=0._kd
         tol=10._kd**(-25._kd)
         if(nra<70) tol=10._kd**(-26._kd)
         if(nra<40) tol=10._kd**(-28._kd)
         if(nra<10) tol=10._kd**(-30._kd)
         tol=10._kd**(-15._kd)
         p_m=(p_b+p_e)/2._kd
!******************************
         length=p_m-p_b
         do i=1,n
            pi=gau_leg_poi(i)*(length)+p_b
            call wave_func_c(nra,lr,bsquare_kd,pi,wave,roots(:nra,
     &               nra,lr),.false.)
            u1=wave*(-1)**nra
            area1=area1+u1*pi*gau_leg_we(i)*length!*vrel
!     & *exp(-(pi/10)**6)
         enddo
!******************************
         length=p_e-p_m
         do i=1,n
            pi=gau_leg_poi(i)*(length)+p_m
            call wave_func_c(nra,lr,bsquare_kd,pi,wave,roots(:nra,
     &               nra,lr),.false.)
            u1=wave*(-1)**nra
            area2=area2+u1*pi*gau_leg_we(i)*length!*vrel
!     & *exp(-(pi/10)**6)
         enddo
!******************************
!         if(abs(area1+area2-area)>tol*(p_e-p_b)/20._kd) then
         if(abs(area1+area2-area)>tol) then
            tmp2=0
            res=adapt_integ_leg(nra,nrb,lr,p_b,p_m,area1,n,totn,
     &             interr)
            res=res +adapt_integ_leg(nra,nrb,lr,p_m,p_e,area2,n,totn
     &             ,interr)
            interr=interr+tmp2
         else
            length=p_e-p_b
            if(totn==0) then
               m_temp(nra,lr,1)=p_b
            endif
            m_temp(nra,lr,int(totn/n)+2)=p_e
            do i=1,n
               if(p_b/=0 .and. i==0) then
                  step_arr(totn+i,nra,lr)=(step_arr(totn+i,nra,lr)
     &                   +step)/2._kd
               else
                  mesh(totn+i,nra,lr)=p_b+gau_leg_poi(i)*length
                  step_arr(totn+i,nra,lr)=gau_leg_we(i)*length
               endif
            enddo
            totn=totn+n
            interr=interr+abs(area1+area2-area)
            res=area
         endif
      end function adapt_integ_leg
!**********************
      end subroutine inter_2b_quadrature

!**********************
!**********************
      subroutine inter_2b_quadrature_qp(ipotin)
      use constants
      !$ use OMP_LIB
      use paramdef
      use harmosc
      use SRG_module
      use iso_fortran_env
      use nninteraction
      use pless_eft_param, only: coulomb_switch
      !use trep_integration
      use omp_lib
      use H_matrices
      use quadratures
      use HO_wf
      use n3lo_qp_module
      implicit none
      integer,parameter :: kd=real128
      integer,parameter :: dp=kind(0.d0)
!*      parameter (nstep=500,pc=30.d0,pinf=0.000000001)
      complex(kd) :: vnp(6),vpp(6),vnn(6),sumv(6)  
      !file variables
      character(len=5):: string1,string2,string3,string4
     & ,string5,string6
      character(len=100):: v_file_name
      logical:: filetest=.false.
      !timer variables
      real:: time1,time2
      integer:: ir,time3,time4,time5,time6
      !parameters 
      integer:: Nf,Nmax,Nd,v_n,n_threads
      real(kd):: fmbyh2,hbc3,pfmev,rdmavg,bosc,th_kd,prange
      !temporary variables 
      integer:: i,k,nq,ia,ib,ind,nra,nrb,lra,lrb,jrel
      real(kd) ::pi,pk,rrel
      real(dp)::pii,pkk
      real(kd)::e_kd
      complex(kd),allocatable,dimension(:,:,:,:,:):: H0
      complex(kd),allocatable,dimension(:,:,:):: vint

      integer:: ipotin,icceff,inn,inp,ipp,ispin,iuc,iucmax,lbmax
      integer::lramin,mtot2,nrbmax,nna,nreffmax,sra
      double precision::cnp,cpp,cnn
      
      !quadrature 
      real(kd),allocatable::mesh_range(:),m_temp(:,:,:),mesh_we(:)
     & ,v_x_mesh(:)
      real(kd)::theta2,interr,regcut,regpow,temp_kd
      complex(kd),allocatable::u4(:,:,:)
      complex(kd):: wave,integ,reg,vrel
      integer,allocatable:: N_arr(:,:)
      integer::mesh_s,NN,totn,qn,n_leg
      real(kd),allocatable:: gau_leg_poi(:)
     & ,gau_leg_we(:)
      real(kd),allocatable:: roots(:,:,:)
      real(real128),allocatable:: mesh(:,:,:),step_arr(:,:,:)
      real(kd):: v_scale,v_theta

      ! eigen solver variables
      complex(kd),allocatable:: Heff_kd(:,:),eigr_kd(:),eigv_kd(:,:)
     &, Heff_kd2(:,:)
      complex(kind(0.d0)),allocatable::Heff(:,:),eigr(:),eigv(:,:)
      real:: vp
      real(kind(0.d0)):: vdp
      real(kd):: vkd
      complex(kd),allocatable,dimension(:,:,:,:):: sumvi
!----------------------------------------------------------------------------
      !common block for n3lo subroutine
      real*16:: XMEV,YMEV,V(6)
      LOGICAL HEFORM,SING,TRIP,COUP,ENDEP
      character*4 label
      integer::J,KREAD,KWRITE,KPUNCH,KDA(9),ipot
!----------------------------------------------------------------------------
!        THE USER NEEDS TO BE FAMILIAR WITH THE FOLLOWING
!        THREE COMMON BLOCKS:
!
      COMMON /CRDWRT_qp/ KREAD,KWRITE,KPUNCH,KDA
!
!        ARGUMENTS AND VALUES OF CDBONN:
      COMMON /CPOT_qp/   V,XMEV,YMEV
      COMMON /CSTATE_qp/ J,HEFORM,SING,TRIP,COUP,ENDEP,LABEL
      common /cnn_qp/ ipot
      interface
         subroutine mkl_zgeev(A,n,m,W,Z)
            integer, intent(in) :: n,m
            complex(kind(0.d0)),intent(in) :: A(n,n)
            complex(kind(0.d0)), intent(out) :: W(n),Z(n,n)
         end subroutine mkl_zgeev
      end interface
!----------------------------------------------------------------------------
!
!        ONCE FOR EVER, SET THE FOLLOWING 6 PARAMETERS AND FORGET THEM:
!        --------------------------------------------------------------
      KREAD=5
      KWRITE=6
      HEFORM=.FALSE.
      SING=.TRUE.
      TRIP=.TRUE.
      COUP=.TRUE.
!        IGNORE KPUNCH, KDA(9), ENDEP, and LABEL.
!        ----------------------------------------
      open(223,file="inter_2b_fit.out")
      call system_clock(count=time5, count_rate=ir)
!***********************
      !initialize variables
      ipot=ipotin
      prange=6.00
      Nf=9
      v_n=100
      v_scale=1._kd
      v_theta=0._kd
      e_kd=exp(1._kd)
      ipot=ipotin
      if (ipot==1) then
         ipp=1
         cpp=1.d0
         print *,' ipp=',ipp
         print *,' cpp=',cpp
      else
         ipp=0
      endif
      rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
      fmbyh2=rdmavg/hbc**2
      anu_kd=fmbyh2*hbo ! anu=(m*h_bar*omega)/ (h_bar*c)^2
      pfmev=hbc/dsqrt(2.d0)
      bsquare_kd=1.d0/anu_kd
      bosc=sqrt(real(bsquare_kd))
      hbc3=hbc*hbc*hbc/dsqrt(8.d0)
      th_kd=theta
      ! initialize gaus legendre points for n3lo_qp and store them in two arrays from n3lo_qp_module
      do i=25,125,25
      call gauleg2(i,-1._kd,1._kd,gaus_legendre_points(:i,i)
     &                      ,gaus_legendre_weights(:i,i))
      do j=1,i
         do lra=0,7
         call legendre_polynomial(legendre_polynom(j,i,lra)
     &  ,gaus_legendre_points(j,i),lra)
         enddo
      enddo
      enddo
     
!***********************
      !printing a description
      print*,"called inter_2b_fit"
      print*,"$$$$$ theta=",th_kd

!$OMP parallel
!$       n_threads=OMP_GET_MAX_THREADS()
!$OMP end parallel
      print*,"The number of threads being used with Open MP:"
     &                                            ,n_threads
      print *,'  iucdim=',iucdim

      if (icutnum==100) then
         write(2,*)
     &        'N3LO NN potential in momentum space'
      elseif (icutnum==102) then
         write(2,*)
     &        'NNLO(POUNDerS) NN potential in momentum space'
      elseif (icutnum==105) then
         write(2,*)
     &        'NNLOsat NN potential in momentum space'
      elseif (icutnum==106) then
         write(2,*)
     &        'N4LO500 NN potential in momentum space'
      elseif (icutnum==107) then
         write(2,*)
     &        'N3LO500new NN potential in momentum space'
      elseif (icutnum==108) then
         write(2,*)
     &      'N2LO500 NN potential in momentum space'
      elseif (icutnum==109) then
         write(2,*)
     &        'NLO500 NN potential in momentum space'
      elseif (icutnum==110) then
         write(2,*)
     &        'LO500 NN potential in momentum space'
      elseif (icutnum==115) then
         write(2,*)
     &        'chiral NN potential in momentum space from inputs'
         if(.not.coulomb_switch.and.ipp/=0) ipp=0
      endif

      if (ipp/=0) then
      endif
!*****************************
      theta2=theta
      !theta=0
      allocate(gau_leg_poi(100))
      allocate(gau_leg_we(100))
      if (allocated(roots)) deallocate(roots)
      allocate(roots(nrelm,nrelm,0:lrelm))
      do lra=0,lrelm
         call laguerre_rootfinder_kd(nrelm,lra+0.5_kd
     &                                 ,roots(1:nrelm,1:nrelm,lra))
      end do
      print*,"laguerre roots allocated"
      n_leg=50
      allocate(m_temp(0:nrelm,0:lrelm,100))
      m_temp=0
      NN=0
      allocate(mesh(0:5000,0:nrelm,0:lrelm),N_arr(0:nrelm,0:lrelm))
      allocate(step_arr(0:5000,0:nrelm,0:lrelm))
      call gauleg2(n_leg,0._kd,1._kd,gau_leg_poi,gau_leg_we)
         call cpu_time(time1)
         call system_clock(count=time3, count_rate=ir)
!$OMP parallel do 
!$OMP& private(nra,lra,totn,interr,sumv)
!$OMP& schedule(dynamic)
      do nra=0,nrelm
         do lra=0,lrelm
            totn=0
            interr=0
            wave=adapt_integ_leg(nra,nra,lra,0._kd,
     &   15._kd,-100._kd,n_leg,totn,interr)
            N_arr(nra,lra)=totn
            NN=max(NN,totn)
          if((mod(nra,20)==0 .or. nra==10) .and. lra<3) then
            print*,"nra,:",nra," integral: ", real(wave,8),
     &          " points ",totn,"error",real(interr,8)
          endif
         enddo
      enddo
!$OMP end parallel do
         call cpu_time(time2)
         call system_clock(count=time4, count_rate=ir)
      print*,'time', 
     &      real(time4 - time3,kind=8)/real(ir,kind=8),time2-time1
      print*,"NN= ",NN
      allocate(mesh_range(120))
      mesh_range=0
      mesh_s=0
      !create the boundaries of the integration ranges
      do nra=0,nrelm
         do lra=0,2
            do k=1,N_arr(nra,lra)/n_leg+1
         pi=m_temp(nra,lra,k)
         i=1
         do while( mesh_range(i)<pi .and. i<=mesh_s)
            i=i+1
         enddo
         if(i<=mesh_s .and. pi==mesh_range(i)) cycle
         do j=mesh_s+1,i+1,-1
            mesh_range(j)=mesh_range(j-1)
         enddo
         mesh_range(i)=pi
         mesh_s=mesh_s+1
            enddo
         enddo
      enddo
!      n_leg=100
!      deallocate(gau_leg_poi,gau_leg_we)
!      allocate(gau_leg_poi(n_leg),gau_leg_we(n_leg))
!      call gauleg2(n_leg,0._kd,1._kd,gau_leg_poi,gau_leg_we)
      !calculate the integration points
      print*,mesh_s
      do i=1,mesh_s
         print*,mesh_range(i)
      enddo
      allocate(v_x_mesh(n_leg*(mesh_s-1)),mesh_we(n_leg*(mesh_s-1)))
      do i=1,mesh_s-1
         temp_kd=mesh_range(i+1)-mesh_range(i)
         do k=1,n_leg
            v_x_mesh((i-1)*n_leg+k)=mesh_range(i)+gau_leg_poi(k)
     &             *temp_kd
            mesh_we((i-1)*n_leg+k)=temp_kd*gau_leg_we(k)
         enddo
      enddo
         qn=(mesh_s-1)*n_leg
         print*,"quadrature points:",qn

!
!      qn=2000
!      deallocate(gau_leg_poi,gau_leg_we,v_x_mesh,mesh_we)
!      allocate(gau_leg_poi(qn),v_x_mesh(qn),mesh_we(qn))
!      allocate(gau_leg_we(qn))
!      call gauleg2(qn,0._kd,15._kd,gau_leg_poi,gau_leg_we)
!      mesh_s=2
!      v_x_mesh(:)=gau_leg_poi
!      mesh_we=gau_leg_we
!
         theta=theta2
         allocate(u4(0:qn,0:nrelm,0:lrelm))
         print*,"bsquare", bsquare_kd
         do lra=0,lrelm
            do nra=0,nrelm
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,pi,wave)
!$OMP& SCHEDULE(DYNAMIC)
               do i=1,qn
!                  pp=gau_leg_poi(i)
                  pi=v_x_mesh(i)
                  !call waver(nr,lr,bsquare,pp,wave)
                  call wave_func_c(nra,lra,bsquare_kd,pi,wave
     &,roots(:nra,nra,lra),.false.)
                  u4(i,nra,lra)=wave*dble((-1)**nra)
     & *complex(cos(theta/2),sin(theta/2))
               end do
!$OMP end parallel do
            end do
         end do
!**************************
!**************************
         theta=theta2
         print*,"wave function created"
!*****************************
      print*,"hbc3",hbc3
      print*,"bsquare", bsquare_kd
!***********************
!calculating v_n3lo
      allocate(vint(0:qn,0:qn,6))
      allocate(sumvi(0:nrelm,0:2,0:qn,6))
      allocate(H0(0:nrelm,0:lrelm,0:nrelm,0:lrelm,6))
      do jrel=0,min(jrelm,mxeng+1)
         J=jrel
         print *,' doing J=',J
         call system_clock(count=time3, count_rate=ir)
         do i=1,qn
            pi=v_x_mesh(i)
            xmev=pi*pfmev
!**            do k=0,nstep
            do k=1,qn
!**               pk=pinf+k*pstep
               pk=v_x_mesh(k)
               ymev=pk*pfmev
               if (pk<=25.d0.and.pk>1.d-9.and.
     &             pi<=25.d0.and.pi>1.d-9     ) then
               if (icutnum==100) then
                  call n3lo_qp
               elseif (icutnum==102) then
                  call nnlo
               elseif (icutnum==105) then
                  call nnlo_sat
               elseif (icutnum==106) then
                  call n4lo500
               elseif (icutnum==107) then
                  call n3lo500new
               elseif (icutnum==108) then
                  call n2lo500
               elseif (icutnum==109) then
                  call nlo500
               elseif (icutnum==110) then
                  call lo500
               elseif (icutnum==115) then
                  call nnlo_input
               endif
               vint(i,k,:)=v(:)
!               vkd=-2.2_kd*1e-6_kd!*exp(
!     &         *exp(-(pi**2._kd+pk**2._kd)/(3.6_kd)**2)
!     &    *exp(-(pi/5._kd)**7._kd -(pk/5._kd)**7._kd) ! regulator
!               vdp=-2.2_dp*1e-6_dp!*exp(
!     &         *exp(-(pi**2._dp+pk**2._dp)/(3.6_dp)**2)
!               vp=-2.2*1e-6!*exp(
!     &         *exp(-(pi**2+pk**2)/(3.6)**2)
!               vint(i,k,:)=vkd
!     &  *exp(-(pi/3.6_kd)**7._kd -(pk/3.6_kd)**7._kd) ! regulator
!     &-3._kd*1e-6_kd*exp(
!     &         -(pi**1._kd+pk**1._kd)/2._kd)
!     &    *exp(-(pi/3.6_kd)**6._kd -(pk/3.6_kd)**6._kd) ! regulator
!     &-3._kd*1e-6_kd*exp(
!     &         -(pi**1._kd+pk**1._kd)/2._kd)
!     &    *exp(-(pi/3.6_kd)**6._kd -(pk/3.6_kd)**6._kd) ! regulator
!               pii=pi
!               pkk=pk
!               vint(i,k,:)=-20._dp*exp(
!     &         -(pii**2._dp+pkk**2._dp)/2._dp)/hbc3
!     &    *exp(-(pi/6._kd)**7._kd -(pk/6._kd)**7._kd) ! regulator
               else
               vint(i,k,:)=0.d0
               endif
            end do
         end do

         print *,' potential calculated'
         call system_clock(count=time4, count_rate=ir)
         print*,"time ",
     &      real(time4 - time3,kind=8)/real(ir,kind=8)
!*****************************************
       
       sumvi=complex(0._kd,0._kd)
       H0=complex(0._kd,0._kd)
       lramin=max(0,jrel-1)
!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(nra,lra,i,k,pk)
!$OMP& SCHEDULE(DYNAMIC)
         do nra=0,nrelm
            do lra=lramin,jrel+1
               if (2*nra+lra>2*nrelm) cycle
               do i=1,qn
                  do k=1,qn
                     pk=v_x_mesh(k)
                     sumvi(nra,lra-lramin,i,:)=
     &               sumvi(nra,lra-lramin,i,:)+mesh_we(k)*vint(i,k,:)
     &                 *u4(k,nra,lra)*pk
                  end do
               end do
            end do
         end do
!$OMP END PARALLEL DO
         print *,' initial integral calculated'

!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(nra,lra,nrb,lrb,sumv,i,pi)
!$OMP& SCHEDULE(DYNAMIC)
         do nra=0,nrelm
            do lra=lramin,jrel+1
               if (2*nra+lra>2*nrelm) cycle
               do nrb=nra,nrelm
                  do lrb=lramin,jrel+1
                     if (2*nrb+lrb>2*nrelm) cycle
                     if ((-1)**(lra+lrb)/=1) cycle

                     sumv=complex(0._kd,0._kd)
                     do i=1,qn
                        pi=v_x_mesh(i)
                        sumv(:)=sumv(:)
     &                    +sumvi(nrb,lrb-lramin,i,:)
     &                    *mesh_we(i)*u4(i,nra,lra)*pi
                     end do
                     sumv(:)=sumv(:)*real(hbc3,kd)
                     H0(nra,lra,nrb,lrb,:)=sumv(:)
                     ! fill the second half of H0
                     H0(nrb,lrb,nra,lra,:)=sumv(:)
                     H0(nrb,lrb,nra,lra,6)=sumv(5)
                     H0(nrb,lrb,nra,lra,5)=sumv(6)
                  end do
               end do
            end do
         end do
!$OMP END PARALLEL DO
         print *,' matrix elements calculated'
         do nra=0,nrelm
            do nrb=0,nrelm
               do lra=max(0,jrel-1),jrel+1
                  do lrb=max(0,lra-2),min(lra+2,jrel+1)
                     if((-1)**(lra+lrb)/=1) cycle
                     if(2*nra+lra>2*nrelm) cycle
                     if(2*nrb+lrb>2*nrelm) cycle

                    sumv(:)=H0(nra,lra,nrb,lrb,:)
                    ia=nra*2+(lra-jrel+1)/2
                    ib=nrb*2+(lrb-jrel+1)/2
                    !****************
                     if (jrel==0) then
                        if (lra==0) then
                           v_uc(nra,nrb,1)=v_uc(nra,nrb,1)
     &                                     +sumv(1)
                        elseif(lra==1) then
                           v_uc(nra,nrb,2)=v_uc(nra,nrb,2)
     &                                     +sumv(3)
                        endif
                     else
                        if (lra==jrel.and.lrb==jrel) then
                           v_uc(nra,nrb,2*jrel+1)=v_uc(nra,nrb,2*jrel+1
     &                                            )+sumv(1)
                           v_uc(nra,nrb,2*jrel+2)=v_uc(nra,nrb,2*jrel+2
     &                                           )+sumv(2)
                        elseif (lra==jrel+1.and.lrb==jrel+1) then
                           v_cc(ia,ib,jrel)=v_cc(ia,ib,jrel) 
     &                                      +sumv(3)
                        elseif (lra==jrel-1.and.lrb==jrel-1) then
                           v_cc(ia,ib,jrel)=v_cc(ia,ib,jrel)
     &                                      +sumv(4)
                        elseif (lra==jrel+1.and.lrb==jrel-1) then
                           v_cc(ia,ib,jrel)=v_cc(ia,ib,jrel)
     &                                      +sumv(5)
                        elseif (lra==jrel-1.and.lrb==jrel+1) then
                           v_cc(ia,ib,jrel)=v_cc(ia,ib,jrel)
     &                                      +sumv(6)
                        endif
                     endif
                    !****************

                  enddo
               enddo
            enddo
         enddo
      

      end do



!**********************************************
      close(223)
      return


!      stop
! end write
      close(223)

!************
      contains

      recursive function adapt_integ_leg(nra,nrb,lr,
     &    p_b,p_e,area,n,totn,interr) result(res)
         implicit none
         !integer,parameter:: kd=real128
         real(kd)::p_b,p_e,area,interr
         integer::n,nra,nrb,lr,totn

         real(kd):: res,tol,u1,u2
         real(kd)::p_m,area1,area2,step,pi,cf,tmp2,length,vrel
         complex(kd)::wave
         integer:: i,reste,tmp1
         area1=0._kd
         area2=0._kd
         tol=10._kd**(-25._kd)
         if(nra<70) tol=10._kd**(-26._kd)
         if(nra<40) tol=10._kd**(-28._kd)
         if(nra<10) tol=10._kd**(-30._kd)
         tol=10._kd**(-17._kd)
         p_m=(p_b+p_e)/2._kd
!******************************
         length=p_m-p_b
         do i=1,n
            pi=gau_leg_poi(i)*(length)+p_b
            call wave_func_c(nra,lr,bsquare_kd,pi,wave,roots(:nra,
     &               nra,lr),.false.)
            u1=wave*(-1)**nra
            area1=area1+u1*pi*gau_leg_we(i)*length!*vrel
!     &* -10._kd*1e-6_kd
     & *exp(-1._kd*(pi/3.6_kd)**4._kd)
         enddo
!******************************
         length=p_e-p_m
         do i=1,n
            pi=gau_leg_poi(i)*(length)+p_m
            call wave_func_c(nra,lr,bsquare_kd,pi,wave,roots(:nra,
     &               nra,lr),.false.)
            u1=wave*(-1)**nra
            area2=area2+u1*pi*gau_leg_we(i)*length!*vrel
!     &* -10._kd*1e-6_kd
     & *exp(-1._kd*(pi/3.6_kd)**4._kd)
         enddo
!******************************
!         if(abs(area1+area2-area)>tol*(p_e-p_b)/20._kd) then
         if(abs(area1+area2-area)>tol) then
            tmp2=0
            res=adapt_integ_leg(nra,nrb,lr,p_b,p_m,area1,n,totn,
     &             interr)
            res=res +adapt_integ_leg(nra,nrb,lr,p_m,p_e,area2,n,totn
     &             ,interr)
            interr=interr+tmp2
         else
            length=p_e-p_b
            if(totn==0) then
               m_temp(nra,lr,1)=p_b
            endif
            !print*,p_b,p_e
            m_temp(nra,lr,int(totn/n)+2)=p_e
            do i=1,n
               if(p_b/=0 .and. i==0) then
                  step_arr(totn+i,nra,lr)=(step_arr(totn+i,nra,lr)
     &                   +step)/2._kd
               else
                  mesh(totn+i,nra,lr)=p_b+gau_leg_poi(i)*length
                  step_arr(totn+i,nra,lr)=gau_leg_we(i)*length
               endif
            enddo
            totn=totn+n
            interr=interr+abs(area1+area2-area)
            res=area
         endif
      end function adapt_integ_leg
!**********************
      end subroutine inter_2b_quadrature_qp

!**********************
!**********************
      subroutine inter_2b_quadrature_qp_c(ipotin)
      use constants
      !$ use OMP_LIB
      use paramdef
      use harmosc
      use SRG_module
      use iso_fortran_env
      use nninteraction
      use pless_eft_param, only: coulomb_switch
      !use trep_integration
      use omp_lib
      use H_matrices
      use quadratures
      use HO_wf
      use n3lo_qp_c_module
      implicit none
      integer,parameter :: kd=real128
      integer,parameter :: dp=kind(0.d0)
!*      parameter (nstep=500,pc=30.d0,pinf=0.000000001)
      complex(kd) :: vnp(6),vpp(6),vnn(6),sumv(6)  
      !file variables
      character(len=5):: string1,string2,string3,string4
     & ,string5,string6
      character(len=100):: v_file_name
      logical:: filetest=.false.
      !timer variables
      real:: time1,time2
      integer:: ir,time3,time4,time5,time6
      !parameters 
      integer:: Nf,Nmax,Nd,v_n,n_threads
      real(kd):: fmbyh2,hbc3,pfmev,rdmavg,bosc,th_kd,prange
      !temporary variables 
      integer:: i,k,nq,ia,ib,ind,nra,nrb,lra,lrb,jrel
      real(kd) ::pi,pk,rrel
      real(dp)::pii,pkk
      real(kd)::e_kd
      complex(kd),allocatable,dimension(:,:,:,:,:):: H0
      complex(kd),allocatable,dimension(:,:,:):: vint

      integer:: ipotin,icceff,inn,inp,ipp,ispin,iuc,iucmax,lbmax
      integer::lramin,mtot2,nrbmax,nna,nreffmax,sra
      double precision::cnp,cpp,cnn
      
      !quadrature 
      real(kd),allocatable::mesh_range(:),m_temp(:,:,:),mesh_we(:)
     & ,v_x_mesh(:)
      real(kd)::theta2,interr,regcut,regpow,temp_kd
      complex(kd),allocatable::u4(:,:,:)
      complex(kd):: wave,integ,reg,vrel
      integer,allocatable:: N_arr(:,:)
      integer::mesh_s,NN,totn,qn,n_leg
      real(kd),allocatable:: gau_leg_poi(:)
     & ,gau_leg_we(:)
      real(kd),allocatable:: roots(:,:,:)
      real(real128),allocatable:: mesh(:,:,:),step_arr(:,:,:)
      real(kd):: v_scale,v_theta

      ! eigen solver variables
      complex(kd),allocatable:: Heff_kd(:,:),eigr_kd(:),eigv_kd(:,:)
     &, Heff_kd2(:,:)
      complex(kind(0.d0)),allocatable::Heff(:,:),eigr(:),eigv(:,:)

      complex(kd):: expo,expoo,expo1,reg1,reg2,reg_comp,connect,cutf
      complex(kd),allocatable,dimension(:,:,:,:):: sumvi
!----------------------------------------------------------------------------
      !common block for n3lo subroutine
      complex*16:: XMEV,YMEV,V(6)
      LOGICAL HEFORM,SING,TRIP,COUP,ENDEP
      character*4 label
      integer::J,KREAD,KWRITE,KPUNCH,KDA(9),ipot
!----------------------------------------------------------------------------
!        THE USER NEEDS TO BE FAMILIAR WITH THE FOLLOWING
!        THREE COMMON BLOCKS:
!
      COMMON /CRDWRT_qp/ KREAD,KWRITE,KPUNCH,KDA
!
!        ARGUMENTS AND VALUES OF CDBONN:
      COMMON /CPOT_qp/   V,XMEV,YMEV
      COMMON /CSTATE_qp/ J,HEFORM,SING,TRIP,COUP,ENDEP,LABEL
      common /cnn_qp/ ipot
      interface
         subroutine mkl_zgeev(A,n,m,W,Z)
            integer, intent(in) :: n,m
            complex(kind(0.d0)),intent(in) :: A(n,n)
            complex(kind(0.d0)), intent(out) :: W(n),Z(n,n)
         end subroutine mkl_zgeev
      end interface
!----------------------------------------------------------------------------
!
!        ONCE FOR EVER, SET THE FOLLOWING 6 PARAMETERS AND FORGET THEM:
!        --------------------------------------------------------------
      KREAD=5
      KWRITE=6
      HEFORM=.FALSE.
      SING=.TRUE.
      TRIP=.TRUE.
      COUP=.TRUE.
!        IGNORE KPUNCH, KDA(9), ENDEP, and LABEL.
!        ----------------------------------------
      open(223,file="inter_2b_fit.out")
      call system_clock(count=time5, count_rate=ir)
!***********************
      !initialize variables
      ipot=ipotin
      prange=6.00
      Nf=9
      v_n=100
      v_scale=1._kd
      v_theta=0._kd
      e_kd=exp(1._kd)
      ipot=ipotin
      if (ipot==1) then
         ipp=1
         cpp=1.d0
         print *,' ipp=',ipp
         print *,' cpp=',cpp
      else
         ipp=0
      endif
      rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
      fmbyh2=rdmavg/hbc**2
      anu_kd=fmbyh2*hbo ! anu=(m*h_bar*omega)/ (h_bar*c)^2
      pfmev=hbc/dsqrt(2.d0)
      bsquare_kd=1.d0/anu_kd
      bosc=sqrt(real(bsquare_kd))
      hbc3=hbc*hbc*hbc/dsqrt(8.d0)
      th_kd=theta
      ! initialize gaus legendre points for n3lo_qp and store them in two arrays from n3lo_qp_module
      do i=25,125,25
      call gauleg2(i,-1._kd,1._kd,gaus_legendre_points(:i,i)
     &                      ,gaus_legendre_weights(:i,i))
      do j=1,i
         do lra=0,7
         call legendre_polynomial(legendre_polynom(j,i,lra)
     &  ,gaus_legendre_points(j,i),lra)
         enddo
      enddo
      enddo
     
!***********************
      !printing a description
      print*,"called inter_2b_fit"
      print*,"$$$$$ theta=",th_kd
!$OMP parallel
!$       n_threads=OMP_GET_MAX_THREADS()
!$OMP end parallel
      print*,"The number of threads being used with Open MP:"
     &                                            ,n_threads
      print *,'  iucdim=',iucdim

      if (icutnum==100) then
         write(2,*)
     &        'N3LO NN potential in momentum space'
      elseif (icutnum==102) then
         write(2,*)
     &        'NNLO(POUNDerS) NN potential in momentum space'
      elseif (icutnum==105) then
         write(2,*)
     &        'NNLOsat NN potential in momentum space'
      elseif (icutnum==106) then
         write(2,*)
     &        'N4LO500 NN potential in momentum space'
      elseif (icutnum==107) then
         write(2,*)
     &        'N3LO500new NN potential in momentum space'
      elseif (icutnum==108) then
         write(2,*)
     &      'N2LO500 NN potential in momentum space'
      elseif (icutnum==109) then
         write(2,*)
     &        'NLO500 NN potential in momentum space'
      elseif (icutnum==110) then
         write(2,*)
     &        'LO500 NN potential in momentum space'
      elseif (icutnum==115) then
         write(2,*)
     &        'chiral NN potential in momentum space from inputs'
         if(.not.coulomb_switch.and.ipp/=0) ipp=0
      endif

      if (ipp/=0) then
      endif
!*****************************
      theta2=theta
      theta=0
      allocate(gau_leg_poi(100))
      allocate(gau_leg_we(100))
      if (allocated(roots)) deallocate(roots)
      allocate(roots(nrelm,nrelm,0:lrelm))
      do lra=0,lrelm
         call laguerre_rootfinder_kd(nrelm,lra+0.5_kd
     &                                 ,roots(1:nrelm,1:nrelm,lra))
      end do
      print*,"laguerre roots allocated"
      n_leg=50
      allocate(m_temp(0:nrelm,0:lrelm,100))
      m_temp=0
      NN=0
      allocate(mesh(0:5000,0:nrelm,0:lrelm),N_arr(0:nrelm,0:lrelm))
      allocate(step_arr(0:5000,0:nrelm,0:lrelm))
      call gauleg2(n_leg,0._kd,1._kd,gau_leg_poi,gau_leg_we)
         call cpu_time(time1)
         call system_clock(count=time3, count_rate=ir)
!$OMP parallel do 
!$OMP& private(nra,lra,totn,interr,sumv)
!$OMP& schedule(dynamic)
      do nra=0,nrelm
         do lra=0,lrelm
            totn=0
            interr=0
            wave=adapt_integ_leg(nra,nra,lra,0._kd,
     &   15._kd,-100._kd,n_leg,totn,interr)
            N_arr(nra,lra)=totn
            NN=max(NN,totn)
          if((mod(nra,20)==0 .or. nra==10) .and. lra<3) then
            print*,"nra,:",nra," integral: ", real(wave,8),
     &          " points ",totn,"error",real(interr,8)
          endif
         enddo
      enddo
!$OMP end parallel do
         call cpu_time(time2)
         call system_clock(count=time4, count_rate=ir)
      print*,'time', 
     &      real(time4 - time3,kind=8)/real(ir,kind=8),time2-time1
      print*,"NN= ",NN
      allocate(mesh_range(120))
      mesh_range=0
      mesh_s=0
      !create the boundaries of the integration ranges
      do nra=0,nrelm
         do lra=0,2
            do k=1,N_arr(nra,lra)/n_leg+1
         pi=m_temp(nra,lra,k)
         i=1
         do while( mesh_range(i)<pi .and. i<=mesh_s)
            i=i+1
         enddo
         if(i<=mesh_s .and. pi==mesh_range(i)) cycle
         do j=mesh_s+1,i+1,-1
            mesh_range(j)=mesh_range(j-1)
         enddo
         mesh_range(i)=pi
         mesh_s=mesh_s+1
            enddo
         enddo
      enddo
!      n_leg=100
!      deallocate(gau_leg_poi,gau_leg_we)
!      allocate(gau_leg_poi(n_leg),gau_leg_we(n_leg))
!      call gauleg2(n_leg,0._kd,1._kd,gau_leg_poi,gau_leg_we)
      !calculate the integration points
      print*,mesh_s
      do i=1,mesh_s
         print*,mesh_range(i)
      enddo
      allocate(v_x_mesh(n_leg*(mesh_s-1)),mesh_we(n_leg*(mesh_s-1)))
      do i=1,mesh_s-1
         temp_kd=mesh_range(i+1)-mesh_range(i)
         do k=1,n_leg
            v_x_mesh((i-1)*n_leg+k)=mesh_range(i)+gau_leg_poi(k)
     &             *temp_kd
            mesh_we((i-1)*n_leg+k)=temp_kd*gau_leg_we(k)
         enddo
      enddo
         qn=(mesh_s-1)*n_leg
         print*,"quadrature points:",qn

!
      qn=400
      deallocate(gau_leg_poi,gau_leg_we,v_x_mesh,mesh_we)
      allocate(gau_leg_poi(qn),v_x_mesh(qn),mesh_we(qn))
      allocate(gau_leg_we(qn))
      call gauleg2(qn,0._kd,10._kd,gau_leg_poi,gau_leg_we)
      mesh_s=2
      v_x_mesh(:)=gau_leg_poi
      mesh_we=gau_leg_we
!
         theta=0
         allocate(u4(0:qn,0:nrelm,0:lrelm))
         print*,"bsquare", bsquare_kd
         do lra=0,lrelm
            do nra=0,nrelm
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,pi,wave)
!$OMP& SCHEDULE(DYNAMIC)
               do i=1,qn
!                  pp=gau_leg_poi(i)
                  pi=v_x_mesh(i)
                  !call waver(nr,lr,bsquare,pp,wave)
                  call wave_func_c(nra,lra,bsquare_kd,pi,wave
     &,roots(:nra,nra,lra),.false.)
                  u4(i,nra,lra)=wave*dble((-1)**nra)
     & *complex(cos(theta/2),sin(theta/2))
               end do
!$OMP end parallel do
            end do
         end do
!**************************
!**************************
         theta=theta2
         print*,"wave function created"
!*****************************
      print*,"hbc3",hbc3
      print*,"bsquare", bsquare_kd
!***********************
!calculating v_n3lo
      allocate(vint(0:qn,0:qn,6))
      allocate(sumvi(0:nrelm,0:2,0:qn,6))
      allocate(H0(0:nrelm,0:lrelm,0:nrelm,0:lrelm,6))
      do jrel=0,min(jrelm,mxeng+1)
         J=jrel
         print *,' doing J=',J
         call system_clock(count=time3, count_rate=ir)
         do i=1,qn
            pi=v_x_mesh(i)
            xmev=pi*pfmev*complex(cos(th_kd),sin(-1._kd*th_kd))
!**            do k=0,nstep
            do k=1,qn
!**               pk=pinf+k*pstep
               pk=v_x_mesh(k)
               ymev=pk*pfmev*complex(cos(th_kd),sin(-1._kd*th_kd))
               if (pk<=25.d0.and.pk>1.d-9.and.
     &             pi<=25.d0.and.pi>1.d-9     ) then
               if (icutnum==100) then
                  call n3lo_qp_c
               elseif (icutnum==102) then
                  call nnlo
               elseif (icutnum==105) then
                  call nnlo_sat
               elseif (icutnum==106) then
                  call n4lo500
               elseif (icutnum==107) then
                  call n3lo500new
               elseif (icutnum==108) then
                  call n2lo500
               elseif (icutnum==109) then
                  call nlo500
               elseif (icutnum==110) then
                  call lo500
               elseif (icutnum==115) then
                  call nnlo_input
               endif
!               expo= (xmev/500._kd)**7._kd + (ymev/500._kd)**7._kd
!               cutf=3.6_kd*pfmev*complex(cos(th_kd),sin(-1._kd*th_kd))
!               if(abs(real(expo))>30) then
!                  expo=expo*30._kd/abs(real(expo))
!               endif
!               expo1=expo
!               if(real(expo)>0) expo1=-expo
!               reg1=(1._kd+exp(4*expo1))/(1._kd+exp(6*expo1))
!     & * exp(expo1)
!               reg1=1._kd/(1._kd+expo)
!               expo= ((xmev-cutf)/500._kd)**3._kd
!     & + ((ymev-cutf)/500._kd)**3._kd
!               if(abs(real(expo))>30) then
!                  expo=expo*30._kd/abs(real(expo))
!               endif
!               expo1=expo
!               if(real(expo)>0) expo1=-expo
!               reg2=(1._kd+exp(4*expo1))/(1._kd+exp(6*expo1))
!     & * exp(expo1)
!              
!               connect=exp(-(pi/4)**4._kd-(pk/4)**4._kd)
!               reg_comp=connect*reg1+(1._kd-connect)*reg2
!
!
               vint(i,k,:)=v(:)
!               vint(i,k,:)=-4._kd*1e-6_kd*exp(
!     &         -(xmev**2._kd+ymev**2._kd)/(500._kd)**2)
     &   *complex(cos(3._kd*th_kd),sin(-3._kd*th_kd))
!     &* reg1
!     &  *exp(-(pi/3.6_kd)**7._kd -(pk/3.6_kd)**7._kd) ! regulator
!     &-3._kd*1e-6_kd*exp(
!     &         -(pi**1._kd+pk**1._kd)/2._kd)
!     &    *exp(-(pi/3.6_kd)**6._kd -(pk/3.6_kd)**6._kd) ! regulator
!     &-3._kd*1e-6_kd*exp(
!     &         -(pi**1._kd+pk**1._kd)/2._kd)
!     &    *exp(-(pi/3.6_kd)**6._kd -(pk/3.6_kd)**6._kd) ! regulator
!               pii=pi
!               pkk=pk
!               vint(i,k,:)=-20._dp*exp(
!     &         -(pii**2._dp+pkk**2._dp)/2._dp)/hbc3
!     &    *exp(-(pi/6._kd)**7._kd -(pk/6._kd)**7._kd) ! regulator
               else
               vint(i,k,:)=0.d0
               endif
            end do
         end do

         print *,' potential calculated'
         call system_clock(count=time4, count_rate=ir)
         print*,"time ",
     &      real(time4 - time3,kind=8)/real(ir,kind=8)
!*****************************************
       
       sumvi=complex(0._kd,0._kd)
       H0=complex(0._kd,0._kd)
       lramin=max(0,jrel-1)
!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(nra,lra,i,k,pk)
!$OMP& SCHEDULE(DYNAMIC)
         do nra=0,nrelm
            do lra=lramin,jrel+1
               if (2*nra+lra>2*nrelm) cycle
               do i=1,qn
                  do k=1,qn
                     pk=v_x_mesh(k)
                     sumvi(nra,lra-lramin,i,:)=
     &               sumvi(nra,lra-lramin,i,:)+mesh_we(k)*vint(i,k,:)
     &                 *u4(k,nra,lra)*pk
                  end do
               end do
            end do
         end do
!$OMP END PARALLEL DO
         print *,' initial integral calculated'

!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(nra,lra,nrb,lrb,sumv,i,pi)
!$OMP& SCHEDULE(DYNAMIC)
         do nra=0,nrelm
            do lra=lramin,jrel+1
               if (2*nra+lra>2*nrelm) cycle
               do nrb=nra,nrelm
                  do lrb=lramin,jrel+1
                     if (2*nrb+lrb>2*nrelm) cycle
                     if ((-1)**(lra+lrb)/=1) cycle

                     sumv=complex(0._kd,0._kd)
                     do i=1,qn
                        pi=v_x_mesh(i)
                        sumv(:)=sumv(:)
     &                    +sumvi(nrb,lrb-lramin,i,:)
     &                    *mesh_we(i)*u4(i,nra,lra)*pi
                     end do
                     sumv(:)=sumv(:)*real(hbc3,kd)
                     sumv(:)=sumv(:)
                     H0(nra,lra,nrb,lrb,:)=sumv(:)
                     ! fill the second half of H0
                     H0(nrb,lrb,nra,lra,:)=sumv(:)
                     H0(nrb,lrb,nra,lra,6)=sumv(5)
                     H0(nrb,lrb,nra,lra,5)=sumv(6)
                  end do
               end do
            end do
         end do
!$OMP END PARALLEL DO
         print *,' matrix elements calculated'
         do nra=0,nrelm
            do nrb=0,nrelm
               do lra=max(0,jrel-1),jrel+1
                  do lrb=max(0,lra-2),min(lra+2,jrel+1)
                     if((-1)**(lra+lrb)/=1) cycle
                     if(2*nra+lra>2*nrelm) cycle
                     if(2*nrb+lrb>2*nrelm) cycle

                    sumv(:)=H0(nra,lra,nrb,lrb,:)
                    ia=nra*2+(lra-jrel+1)/2
                    ib=nrb*2+(lrb-jrel+1)/2
                    !****************
                     if (jrel==0) then
                        if (lra==0) then
                           v_uc(nra,nrb,1)=v_uc(nra,nrb,1)
     &                                     +sumv(1)
                        elseif(lra==1) then
                           v_uc(nra,nrb,2)=v_uc(nra,nrb,2)
     &                                     +sumv(3)
                        endif
                     else
                        if (lra==jrel.and.lrb==jrel) then
                           v_uc(nra,nrb,2*jrel+1)=v_uc(nra,nrb,2*jrel+1
     &                                            )+sumv(1)
                           v_uc(nra,nrb,2*jrel+2)=v_uc(nra,nrb,2*jrel+2
     &                                           )+sumv(2)
                        elseif (lra==jrel+1.and.lrb==jrel+1) then
                           v_cc(ia,ib,jrel)=v_cc(ia,ib,jrel) 
     &                                      +sumv(3)
                        elseif (lra==jrel-1.and.lrb==jrel-1) then
                           v_cc(ia,ib,jrel)=v_cc(ia,ib,jrel)
     &                                      +sumv(4)
                        elseif (lra==jrel+1.and.lrb==jrel-1) then
                           v_cc(ia,ib,jrel)=v_cc(ia,ib,jrel)
     &                                      +sumv(5)
                        elseif (lra==jrel-1.and.lrb==jrel+1) then
                           v_cc(ia,ib,jrel)=v_cc(ia,ib,jrel)
     &                                      +sumv(6)
                        endif
                     endif
                    !****************

                  enddo
               enddo
            enddo
         enddo
      

      end do



!**********************************************
      close(223)
      return


!      stop
! end write
      close(223)

!************
      contains

      recursive function adapt_integ_leg(nra,nrb,lr,
     &    p_b,p_e,area,n,totn,interr) result(res)
         implicit none
         !integer,parameter:: kd=real128
         real(kd)::p_b,p_e,area,interr
         integer::n,nra,nrb,lr,totn

         real(kd):: res,tol,u1,u2
         real(kd)::p_m,area1,area2,step,pi,cf,tmp2,length,vrel
         complex(kd)::wave
         integer:: i,reste,tmp1
         area1=0._kd
         area2=0._kd
         tol=10._kd**(-25._kd)
         if(nra<70) tol=10._kd**(-26._kd)
         if(nra<40) tol=10._kd**(-28._kd)
         if(nra<10) tol=10._kd**(-30._kd)
         tol=10._kd**(-30._kd)
         p_m=(p_b+p_e)/2._kd
!******************************
         length=p_m-p_b
         do i=1,n
            pi=gau_leg_poi(i)*(length)+p_b
            call wave_func_c(nra,lr,bsquare_kd,pi,wave,roots(:nra,
     &               nra,lr),.false.)
            u1=wave*(-1)**nra
            area1=area1+u1*pi*gau_leg_we(i)*length!*vrel
!     &* -10._kd*1e-6_kd
     & *exp(-1._kd*(pi/4._kd)**2._kd)
         enddo
!******************************
         length=p_e-p_m
         do i=1,n
            pi=gau_leg_poi(i)*(length)+p_m
            call wave_func_c(nra,lr,bsquare_kd,pi,wave,roots(:nra,
     &               nra,lr),.false.)
            u1=wave*(-1)**nra
            area2=area2+u1*pi*gau_leg_we(i)*length!*vrel
!     &* -10._kd*1e-6_kd
     & *exp(-1._kd*(pi/4._kd)**2._kd)
         enddo
!******************************
!         if(abs(area1+area2-area)>tol*(p_e-p_b)/20._kd) then
         if(abs(area1+area2-area)>tol) then
            tmp2=0
            res=adapt_integ_leg(nra,nrb,lr,p_b,p_m,area1,n,totn,
     &             interr)
            res=res +adapt_integ_leg(nra,nrb,lr,p_m,p_e,area2,n,totn
     &             ,interr)
            interr=interr+tmp2
         else
            length=p_e-p_b
            if(totn==0) then
               m_temp(nra,lr,1)=p_b
            endif
            !print*,p_b,p_e
            m_temp(nra,lr,int(totn/n)+2)=p_e
            do i=1,n
               if(p_b/=0 .and. i==0) then
                  step_arr(totn+i,nra,lr)=(step_arr(totn+i,nra,lr)
     &                   +step)/2._kd
               else
                  mesh(totn+i,nra,lr)=p_b+gau_leg_poi(i)*length
                  step_arr(totn+i,nra,lr)=gau_leg_we(i)*length
               endif
            enddo
            totn=totn+n
            interr=interr+abs(area1+area2-area)
            res=area
         endif
      end function adapt_integ_leg
!**********************
      end subroutine inter_2b_quadrature_qp_c

!**********************
      end module v2b_quadrature
