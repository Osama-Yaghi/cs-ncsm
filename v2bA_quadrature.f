      module v2bA_quadrature
      use iso_fortran_env
      contains
!**********************
!**********************
      subroutine inter_2b_quadrature_qp_c(icutnum)
      use constants
      !$ use OMP_LIB
      use paramdef
      use harmosc
      use SRG_module
      use iso_fortran_env
      use pless_eft_param, only: coulomb_switch
      !use trep_integration
      use omp_lib
      use pointers
      use ist
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
      integer:: Nf,Nmax,Nd,v_n,n_threads,icutnum
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
      if (nucleons==3.and.threeff) then
!         anu=hbo*(dble(neut3eff)*fmsn+dble(iprot3eff)*fmsp)
!     +        /dble(nucl3eff)/hbc**2
         MTot2=neut3eff-iprot3eff
         ipp=(nucl3eff-Mtot2)*(nucl3eff-Mtot2-2)
         inp=nucl3eff*(nucl3eff-2)+itot23eff*(itot23eff+2)
     &        -2*Mtot2**2
         inn=(nucl3eff+Mtot2)*(nucl3eff+Mtot2-2)
         cpp=dble(ipp)
     &        /dble(3*nucl3eff*(nucl3eff-2)+iTot23eff*(iTot23eff+2))
         cnp=dble(inp)
     &        /dble(3*nucl3eff*(nucl3eff-2)+iTot23eff*(iTot23eff+2))
         cnn=dble(inn)
     &        /dble(3*nucl3eff*(nucl3eff-2)+iTot23eff*(iTot23eff+2))

!**********
      else
!         anu=hbo*(dble(neutrons)*fmsn+dble(iprotons)*fmsp)
!     +        /dble(nucleons)/hbc**2
         MTot2=neutrons-iprotons
         ipp=(nucleons-Mtot2)*(nucleons-Mtot2-2)
         inp=nucleons*(nucleons-2)+itot2*(itot2+2)-2*Mtot2**2
         inn=(nucleons+Mtot2)*(nucleons+Mtot2-2)
         cpp=dble(ipp)
     &        /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
         cnp=dble(inp)
     &        /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
         cnn=dble(inn)
     &        /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
      endif
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
      !allocate(gau_leg_poi(100))
      !allocate(gau_leg_we(100))
      if (allocated(roots)) deallocate(roots)
      allocate(roots(nrelm,nrelm,0:lrelm))
      do lra=0,lrelm
         call laguerre_rootfinder_kd(nrelm,lra+0.5_kd
     &                                 ,roots(1:nrelm,1:nrelm,lra))
      end do
      print*,"laguerre roots allocated"
!**********************************************
!     adaptive integration section (ommitted)
      !      n_leg=50
!      allocate(m_temp(0:nrelm,0:lrelm,100))
!      m_temp=0
!      NN=0
!      allocate(mesh(0:5000,0:nrelm,0:lrelm),N_arr(0:nrelm,0:lrelm)) !      allocate(step_arr(0:5000,0:nrelm,0:lrelm))
!      call gauleg2(n_leg,0._kd,1._kd,gau_leg_poi,gau_leg_we)
!         call cpu_time(time1)
!         call system_clock(count=time3, count_rate=ir)
!!$OMP parallel do 
!!$OMP& private(nra,lra,totn,interr,sumv)
!!$OMP& schedule(dynamic)
!      do nra=0,nrelm
!         do lra=0,lrelm
!            totn=0
!            interr=0
!            wave=adapt_integ_leg(nra,nra,lra,0._kd,
!     &   15._kd,-100._kd,n_leg,totn,interr)
!            N_arr(nra,lra)=totn
!            NN=max(NN,totn)
!          if((mod(nra,20)==0 .or. nra==10) .and. lra<3) then
!            print*,"nra,:",nra," integral: ", real(wave,8),
!     &          " points ",totn,"error",real(interr,8)
!          endif
!         enddo
!      enddo
!!$OMP end parallel do
!         call cpu_time(time2)
!         call system_clock(count=time4, count_rate=ir)
!      print*,'time', 
!     &      real(time4 - time3,kind=8)/real(ir,kind=8),time2-time1
!      print*,"NN= ",NN
!      allocate(mesh_range(120))
!      mesh_range=0
!      mesh_s=0
!      !create the boundaries of the integration ranges
!      do nra=0,nrelm
!         do lra=0,2
!            do k=1,N_arr(nra,lra)/n_leg+1
!         pi=m_temp(nra,lra,k)
!         i=1
!         do while( mesh_range(i)<pi .and. i<=mesh_s)
!            i=i+1
!         enddo
!         if(i<=mesh_s .and. pi==mesh_range(i)) cycle
!         do j=mesh_s+1,i+1,-1
!            mesh_range(j)=mesh_range(j-1)
!         enddo
!         mesh_range(i)=pi
!         mesh_s=mesh_s+1
!            enddo
!         enddo
!      enddo
!!      n_leg=100
!!      deallocate(gau_leg_poi,gau_leg_we)
!!      allocate(gau_leg_poi(n_leg),gau_leg_we(n_leg))
!!      call gauleg2(n_leg,0._kd,1._kd,gau_leg_poi,gau_leg_we)
!      !calculate the integration points
!      print*,mesh_s
!      do i=1,mesh_s
!         print*,mesh_range(i)
!      enddo
!      allocate(v_x_mesh(n_leg*(mesh_s-1)),mesh_we(n_leg*(mesh_s-1)))
!      do i=1,mesh_s-1
!         temp_kd=mesh_range(i+1)-mesh_range(i)
!         do k=1,n_leg
!            v_x_mesh((i-1)*n_leg+k)=mesh_range(i)+gau_leg_poi(k)
!     &             *temp_kd
!            mesh_we((i-1)*n_leg+k)=temp_kd*gau_leg_we(k)
!         enddo
!      enddo
!         qn=(mesh_s-1)*n_leg
!         print*,"quadrature points:",qn
!
!******************************************
      qn=400
      !deallocate(gau_leg_poi,gau_leg_we,v_x_mesh,mesh_we)
      allocate(gau_leg_poi(qn),v_x_mesh(qn),mesh_we(qn))
      allocate(gau_leg_we(qn))
      call gauleg2(qn,0._kd,10._kd,gau_leg_poi,gau_leg_we)
      mesh_s=2
      v_x_mesh(:)=gau_leg_poi
      mesh_we=gau_leg_we
      print*,"Gauss-Legendre quadrature points:",qn
!
!******************************************
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
               !
               ipot=2
               if (icutnum==100) then
                  call n3lo_qp_c
               endif
               vnp(:)=v(:)
               !
               !
               ipot=1
               if (icutnum==100) then
                  call n3lo_qp_c
               endif
               vpp(:)=v(:)
               !
               !
               ipot=3
               if (icutnum==100) then
                  call n3lo_qp_c
               endif
               vnn(:)=v(:)
               !
               if (jrel==0) then
                  vint(i,k,:)=vnp(:)*cnp+vnn(:)*cnn+vpp(:)*cpp
               elseif (jrel==2*(jrel/2)) then
                  vint(i,k,1)=vnp(1)*cnp+vnn(1)*cnn+vpp(1)*cpp
                  vint(i,k,2)=vnp(2)
                  vint(i,k,3)=vnp(3)*cnp+vnn(3)*cnn+vpp(3)*cpp
                  vint(i,k,4)=vnp(4)*cnp+vnn(4)*cnn+vpp(4)*cpp
                  vint(i,k,5)=vnp(5)*cnp+vnn(5)*cnn+vpp(5)*cpp
                  vint(i,k,6)=vnp(6)*cnp+vnn(6)*cnn+vpp(6)*cpp
               else
                  vint(i,k,1)=vnp(1)
                  vint(i,k,2)=vnp(2)*cnp+vnn(2)*cnn+vpp(2)*cpp
                  vint(i,k,3)=vnp(3)
                  vint(i,k,4)=vnp(4)
                  vint(i,k,5)=vnp(5)
                  vint(i,k,6)=vnp(6)
               endif
!               vint(i,k,:)=-10._kd*1e-6_kd*exp(
!     &         -5._kd*(xmev**2._kd+ymev**2._kd)/(500._kd)**2)
               vint(i,k,:)=vint(i,k,:)
     &   *complex(cos(3._kd*th_kd),sin(-3._kd*th_kd))
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
                         v_uc(nra,nrb,1)=sumv(1)
                      elseif(lra==1) then
                         v_uc(nra,nrb,2)=sumv(3)
                      endif
                   else
                      if (lra==jrel.and.lrb==jrel) then
                         v_uc(nra,nrb,2*jrel+1)=sumv(1)
                         v_uc(nra,nrb,2*jrel+2)=sumv(2)
                      elseif (lra==jrel+1.and.lrb==jrel+1) then
                         v_cc(ia,ib,jrel)=sumv(3)
                      elseif (lra==jrel-1.and.lrb==jrel-1) then
                         v_cc(ia,ib,jrel)=sumv(4)
                      elseif (lra==jrel+1.and.lrb==jrel-1) then
                         v_cc(ia,ib,jrel)=sumv(5)
                      elseif (lra==jrel-1.and.lrb==jrel+1) then
                         v_cc(ia,ib,jrel)=sumv(6)
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
         tol=10._kd**(-25._kd)
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
      end subroutine inter_2b_quadrature_qp_c
!**********************
      subroutine inter_c_old
      use constants
      use paramdef
      use ist
      use pointers
      use nninteraction
      use harmosc
!      use potdef
      use gognypot
      use H_matrices
      use HO_wf
      use quadratures
      implicit double precision (a-h,o-z)
      integer,parameter :: kd=real128
      integer,parameter :: dp=kind(0.d0)
      parameter (intstep=2000)
      parameter (nstep=intstep,rc=1.d-6,rinf=22.d0)

       double precision,allocatable,dimension(:,:,:):: u
      complex(kd),allocatable,dimension(:,:):: vuncp,
     & vjm1,vcpl,vjp1
      complex(kd) vpot(2,2),vnnjp1,vnpjp1,vppjp1,vnnjm1,vnpjm1
     & ,vppjm1,vnncpl,vnpcpl,vppcpl,vnn,vnp,vpp
      complex(kd):: wave,r,xmult,over1,over2,fact,over3,over4
      real(kd),allocatable:: roots(:,:,:),rstep
     &,rcr,rinfr,rstepr
      integer:: ipot

      if (chirp) av18=.true.

      if (av18) then
!*         write(2,*)
!*     +        'Argonne V18 NN potential in coordinate space'
         if (chirp) then
            write(2, "(' local chiral NN potential: lpot=',i3)") lpot
         else
            write(2,"(' Argonne V18 NN potential: lpot=',i3)") lpot
         endif
         if (nucleons==3.and.threeff) then
!*            anu=hbo*(dble(neut3eff)*fmsn+dble(iprot3eff)*fmsp)
!*     +           /dble(nucl3eff)/hbc**2
            rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
            fmbyh2=rdmavg/hbc**2
            anu=fmbyh2*hbo
            MTot2=neut3eff-iprot3eff
            ipp=(nucl3eff-Mtot2)*(nucl3eff-Mtot2-2)
            inp=nucl3eff*(nucl3eff-2)+itot23eff*(itot23eff+2)
     &           -2*Mtot2**2
            inn=(nucl3eff+Mtot2)*(nucl3eff+Mtot2-2)
            cpp=dble(ipp)
     &        /dble(3*nucl3eff*(nucl3eff-2)+iTot23eff*(iTot23eff+2))
            cnp=dble(inp)
     &        /dble(3*nucl3eff*(nucl3eff-2)+iTot23eff*(iTot23eff+2))
            cnn=dble(inn)
     &        /dble(3*nucl3eff*(nucl3eff-2)+iTot23eff*(iTot23eff+2))
         else
!*            anu=hbo*(dble(neutrons)*fmsn+dble(iprotons)*fmsp)
!*     +           /dble(nucleons)/hbc**2
            MTot2=neutrons-iprotons
            ipp=(nucleons-Mtot2)*(nucleons-Mtot2-2)
            inp=nucleons*(nucleons-2)+itot2*(itot2+2)-2*Mtot2**2
            inn=(nucleons+Mtot2)*(nucleons+Mtot2-2)
            cpp=dble(ipp)
     &           /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
            cnp=dble(inp)
     &           /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
            cnn=dble(inn)
     &           /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
!*******************
!            cnp=1.d0
!            cnn=0.d0
!            cpp=0.d0
            rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
            fmbyh2=rdmavg/hbc**2
            anu=fmbyh2*hbo
!*******************
         endif
         print *,' ipp,inn,inp:',ipp,inn,inp
         print *,' cpp,cnn,cnp:',cpp,cnn,cnp
      elseif (MTV) then
         write(2,*)
!     +        'Malfliet-Tjon NN potential in coordinate space'
     &        'Minnesota NN potential in coordinate space'
!     +        'Gogny NN potential in coordinate space'
         rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
         fmbyh2=rdmavg/hbc**2
!         anu=fmbyh2*hbo
         anu=hbo/41.47d0
         if (nucleons==3.and.threeff) then
            MTot2=neut3eff-iprot3eff
            ipp=(nucl3eff-Mtot2)*(nucl3eff-Mtot2-2)
            inp=nucl3eff*(nucl3eff-2)+itot23eff*(itot23eff+2)
     &           -2*Mtot2**2
            inn=(nucl3eff+Mtot2)*(nucl3eff+Mtot2-2)
            cpp=dble(ipp)
     &        /dble(3*nucl3eff*(nucl3eff-2)+iTot23eff*(iTot23eff+2))
            cnp=dble(inp)
     &        /dble(3*nucl3eff*(nucl3eff-2)+iTot23eff*(iTot23eff+2))
            cnn=dble(inn)
     &        /dble(3*nucl3eff*(nucl3eff-2)+iTot23eff*(iTot23eff+2))
            rcoul=cpp*hbc/alpha
         else
            MTot2=neutrons-iprotons
            ipp=(nucleons-Mtot2)*(nucleons-Mtot2-2)
            inp=nucleons*(nucleons-2)+itot2*(itot2+2)-2*Mtot2**2
            inn=(nucleons+Mtot2)*(nucleons+Mtot2-2)
            cpp=dble(ipp)
     &              /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
            cnp=dble(inp)
     &              /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
            cnn=dble(inn)
     &              /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
            rcoul=cpp*hbc/alpha
         endif
         print *,' ipp,inn,inp:',ipp,inn,inp
         print *,' cpp,cnn,cnp:',cpp,cnn,cnp
      else
         write(2,*)
     &        'Argonne V8p NN potential in coordinate space'
         rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
         fmbyh2=rdmavg/hbc**2
         anu=fmbyh2*hbo
      endif

      write(2,1246) nrelm,lrelm,jrelm
 1246 format(' nrelm=',i4,'    lrelm=',i3,'    jrelm=',i3)


      bsquare=1.d0/anu
      write(2,1325) hbo,dsqrt(bsquare)
 1325 format(' hbar*Omega=',f8.4,'    b_HO=',f8.4)

      print *,' hbo3=',hbo3,'  intstep=',intstep,'  iucdim=',iucdim

      allocate(vuncp(0:intstep,iucdim))
      vuncp=0.d0
      RSTEP=(RINF-RC)/REAL(NSTEP)

      write(2,1457) rc,rinf,nstep,rstep
 1457 format(' rc=',e14.7,'     rinf=',f8.4,'     nstep=',i5,
     &                    '     rstep=',f8.4)
!*************************      
      if (allocated(roots)) deallocate(roots)
      allocate(roots(nrelm,nrelm,0:lrelm))
      do lra=0,lrelm
         call laguerre_rootfinder_kd(nrelm,lra+0.5_kd
     &                                 ,roots(1:nrelm,1:nrelm,lra))
      end do
      print*,"laguerre roots allocated"
!*************************      

      vnp=0.d0
      vpp=0.d0
      vnn=0.d0

      DO JPW=0,jrelm
         lpw=jpw
         do ispw=0,1
            if (ispw==1.and.jpw==0) lpw=1
            if ((-1)**(lpw+ispw)==1) then
               itpw=1
            else
               itpw=0
            endif
            iuc=2*jpw+ispw+1
            DO I=0,NSTEP
               R=(RC+I*RSTEP)
     &*complex(cos(1._kd*theta),sin(1._kd*theta))

               if (av18) then
                  print*,"not cs yet"
                  stop
                  if (inp/=0) then
                     CALL av18pw(lpot,lpw,ispw,jpw,itpw,1,-1,r,vpot)
                     vnp=vpot(1,1)
                  endif
                  if (ipp/=0.and.itpw==1) then
                     CALL av18pw(lpot,lpw,ispw,jpw,itpw,1,1,r,vpot)
                     vpp=vpot(1,1)
                  endif
                  if (inn/=0.and.itpw==1) then
                     CALL av18pw(lpot,lpw,ispw,jpw,itpw,-1,-1,r,vpot)
                     vnn=vpot(1,1)
                  endif
                  if (itpw==1) then
                     vuncp(i,iuc)=cnp*vnp+cpp*vpp+cnn*vnn
                  else
                     vuncp(i,iuc)=vnp
                  endif
               elseif (MTV) then
!                  VUNCP(I,iuc)=v(r)
                   vpp=v_c(r,1,ispw)
                   vnp=v_c(r,2,ispw)
                   vnn=v_c(r,3,ispw)
                  if (itpw==1) then
                     vuncp(i,iuc)=cnp*vnp+cpp*vpp+cnn*vnn
                  else
                     vuncp(i,iuc)=vnp
                  endif
                  if(lpw/=0) vuncp(i,iuc)=0
                  vuncp(i,iuc)=vuncp(i,iuc)+dble(itpw)*rcoul/r
!                  VUNCP(I,iuc)=vmn_c(r,ispw,itpw)+dble(itpw)*rcoul/r
!                  CALL gognypotpw(lpw,ispw,jpw,itpw,r,vpot)
!                  VUNCP(I,iuc)=VPOT(1,1)
               else
                  print*,"not cs yet"
                  stop
                  CALL av8ppw(lpw,ispw,jpw,itpw,r,vpot)
                  VUNCP(I,iuc)=VPOT(1,1)
               endif

            end do
         ENDDO
      end do
      print *,' uncoupled potential computed'

      allocate(vjm1(0:intstep,jrelm))
      allocate(vcpl(0:intstep,jrelm))
      allocate(vjp1(0:intstep,jrelm))
      vjm1=0.d0
      vcpl=0.d0
      vjp1=0.d0

      vnpjm1=0.d0
      vnpcpl=0.d0
      vnpjp1=0.d0
      vppjm1=0.d0
      vppcpl=0.d0
      vppjp1=0.d0
      vnnjm1=0.d0
      vnncpl=0.d0
      vnnjp1=0.d0

      do jpw=1,jrelm
         lpw=jpw-1
         ispw=1
         if ((-1)**(lpw+ispw)==1) then
            itpw=1
         else
            itpw=0
         endif
         DO I=0,NSTEP
            R=(RC+I*RSTEP)
     &*complex(cos(1._kd*theta),sin(1._kd*theta))
               if (av18) then
                  print*,"not cs yet"
                  stop
                  if (inp/=0) then
                     CALL av18pw(lpot,lpw,ispw,jpw,itpw,1,-1,r,vpot)
                     vnpjm1=vpot(1,1)
                     vnpcpl=vpot(1,2)
                     vnpjp1=vpot(2,2)
                  endif
                  if (ipp/=0.and.itpw==1) then
                     CALL av18pw(lpot,lpw,ispw,jpw,itpw,1,1,r,vpot)
                     vppjm1=vpot(1,1)
                     vppcpl=vpot(1,2)
                     vppjp1=vpot(2,2)
                  endif
                  if (inn/=0.and.itpw==1) then
                     CALL av18pw(lpot,lpw,ispw,jpw,itpw,-1,-1,r,vpot)
                     vnnjm1=vpot(1,1)
                     vnncpl=vpot(1,2)
                     vnnjp1=vpot(2,2)
                  endif
                  if (itpw==1) then
                     vjm1(i,jpw)=cnp*vnpjm1+cpp*vppjm1+cnn*vnnjm1
                     vcpl(i,jpw)=cnp*vnpcpl+cpp*vppcpl+cnn*vnncpl
                     vjp1(i,jpw)=cnp*vnpjp1+cpp*vppjp1+cnn*vnnjp1
                  else
                     vjm1(i,jpw)=vnpjm1
                     vcpl(i,jpw)=vnpcpl
                     vjp1(i,jpw)=vnpjp1
                  endif
               elseif (MTV) then
!                  vjm1(i,jpw)=v(r)
!                  vjp1(i,jpw)=v(r)
                  vppjm1=v_c(r,1,ispw)
                  vnpjm1=v_c(r,2,ispw)
                  vnnjm1=v_c(r,3,ispw)
                  vppjp1=v_c(r,1,ispw)
                  vnpjp1=v_c(r,2,ispw)
                  vnnjp1=v_c(r,3,ispw)
                  if (itpw==1) then
                     vjm1(i,jpw)=cnp*vnpjm1+cpp*vppjm1+cnn*vnnjm1
                     vjp1(i,jpw)=cnp*vnpjp1+cpp*vppjp1+cnn*vnnjp1
                  else
                     vjm1(i,jpw)=vnpjm1
                     vjp1(i,jpw)=vnpjp1
                  endif
                  vjp1(i,jpw)=0! take only s wave
                  if(jpw>1) then
                        vjm1(i,jpw)=0
                        vjp1(i,jpw)=0
                  endif
                  
                  vcpl(i,jpw)=0.d0
                  vjm1(i,jpw)=vjm1(i,jpw)+dble(itpw)*rcoul/r
                  vjp1(i,jpw)=vjp1(i,jpw)+dble(itpw)*rcoul/r
!                  vjm1(i,jpw)=vmn_c(r,ispw,itpw)+dble(itpw)*rcoul/r
!                  vjp1(i,jpw)=vmn_c(r,ispw,itpw)+dble(itpw)*rcoul/r
!                  CALL gognypotpw(lpw,ispw,jpw,itpw,r,vpot)
!                  vjm1(i,jpw)=VPOT(1,1)
!                  vcpl(i,jpw)=VPOT(1,2)
!                  vjp1(i,jpw)=VPOT(2,2)
               else
                  print*,"not cs yet"
                  stop
                  CALL av8ppw(lpw,ispw,jpw,itpw,r,vpot)
                  vjm1(i,jpw)=VPOT(1,1)
                  vcpl(i,jpw)=VPOT(1,2)
                  vjp1(i,jpw)=VPOT(2,2)
               endif
         ENDDO
      end do

      print *,' coupled potential computed'
      print *,' rc=',rc,' rinf=',rinf,' nstep=',nstep

      allocate(u(0:nstep,0:nrelm,0:lrelm))
      rcr=rc/dsqrt(2.d0)
      rinfr=rinf/dsqrt(2.d0)
      rstepr=(rinfr-rcr)/real(nstep)
*      do lr=0,lrelm
*         do nr=0,nrelm
*            if (2*nr+lr>2*nrelm) cycle
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,rr,wave)
!$OMP& SCHEDULE(DYNAMIC)
            do i=0,nstep
               rr=rcr+i*rstepr
!*               call waver(nr,lr,anu,rr,wave)
               call waver_gen(nrelm,lrelm,u(i,:,:),anu,rr)
!*               u(i,nr,lr)=wave
            end do
!$OMP END PARALLEL DO
*         end do
*      end do

         bsquare_kd=1.d0/anu
         anu_kd=anu
!         do lra=0,lrelm
!            do nra=0,nrelm
!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,pi,wave)
!!$OMP& SCHEDULE(DYNAMIC)
!               do i=0,nstep
!               rr=rcr+i*rstepr
!!                  pp=gau_leg_poi(i)
!                  !call waver(nr,lr,bsquare,pp,wave)
!                  call wave_func_c(nra,lra,bsquare_kd,rr,wave
!     &,roots(:max(nra,1),max(nra,1),lra),.false.)
!                  u(i,nra,lra)=wave*dble((-1)**nra)
!               end do
!!$OMP end parallel do
!            end do
!         end do
      print *,' u calculated'

      !allocate(vcc(0:iccdim,0:iccdim,jrelm))
      !allocate(vuc(0:nrelm,0:nrelm,iucdim))

      !vcc=0.d0
      !vuc=0.d0

      do nra=0,nrelm
         do lra=0,jrelm+1
            nna=nra+nra+lra
            if (nna>2*nrelm) cycle
            do nrb=nra,nrelm
               if (lra.eq.0) then
                   lrbmin=0
               elseif (lra.eq.1) then
                   lrbmin=1
               else
                   lrbmin=lra-2
               endif
               do lrb=lrbmin,min(lra+2,jrelm+1),2
                  if (2*nrb+lrb>2*nrelm) cycle
                  if (lrb.eq.lra) then

                      if (lra.eq.0) then

                          over1=0.d0
                          over2=0.d0
                          fact=4.d0
                          do i=1,nstep-1
                             xmult=u(i,nra,0)*u(i,nrb,0)*fact
                             over1=over1+xmult*vuncp(i,1)
                             over2=over2+xmult*vjm1(i,1)
                             fact=6.d0-fact
                          enddo
                          v_uc(nra,nrb,1)=
     &   (over1+u(0,nra,0)*u(0,nrb,0)*vuncp(0,1)
     & +u(nstep,nra,0)*u(nstep,nrb,0)*vuncp(nstep,1))*rstepr/3.d0
                          v_cc(2*nra,2*nrb,1)=
     &   (over2+u(0,nra,0)*u(0,nrb,0)*vjm1(0,1)
     & +u(nstep,nra,0)*u(nstep,nrb,0)*vjm1(nstep,1))*rstepr/3.d0
                          v_uc(nrb,nra,1)=v_uc(nra,nrb,1)
                          v_cc(2*nrb,2*nra,1)=v_cc(2*nra,2*nrb,1)

                      elseif (lra.eq.1) then

                          over1=0.d0
                          over2=0.d0
                          over3=0.d0
                          over4=0.d0
                          fact=4.d0
                          do i=1,nstep-1
                             xmult=u(i,nra,1)*u(i,nrb,1)*fact
                             over1=over1+xmult*vuncp(i,2)
                             over2=over2+xmult*vjm1(i,2)
                             over3=over3+xmult*vuncp(i,3)
                             over4=over4+xmult*vuncp(i,4)
                             fact=6.d0-fact
                          enddo

                          v_uc(nra,nrb,2)=
     &   (over1+u(0,nra,1)*u(0,nrb,1)*vuncp(0,2)
     & +u(nstep,nra,1)*u(nstep,nrb,1)*vuncp(nstep,2))*rstepr/3.d0
                          v_cc(2*nra,2*nrb,2)=
     &   (over2+u(0,nra,1)*u(0,nrb,1)*vjm1(0,2)
     & +u(nstep,nra,1)*u(nstep,nrb,1)*vjm1(nstep,2))*rstepr/3.d0
                          v_uc(nra,nrb,3)=
     &   (over3+u(0,nra,1)*u(0,nrb,1)*vuncp(0,3)
     & +u(nstep,nra,1)*u(nstep,nrb,1)*vuncp(nstep,3))*rstepr/3.d0
                          v_uc(nra,nrb,4)=
     &   (over4+u(0,nra,1)*u(0,nrb,1)*vuncp(0,4)
     & +u(nstep,nra,1)*u(nstep,nrb,1)*vuncp(nstep,4))*rstepr/3.d0
                          v_uc(nrb,nra,2)=v_uc(nra,nrb,2)
                          v_cc(2*nrb,2*nra,2)=v_cc(2*nra,2*nrb,2)
                          v_uc(nrb,nra,3)=v_uc(nra,nrb,3)
                          v_uc(nrb,nra,4)=v_uc(nra,nrb,4)

                      elseif (lra.eq.jrelm+1) then

                          over2=0.d0
                          fact=4.d0

                          do i=1,nstep-1
                             xmult=u(i,nra,lra)
     &                            *vjp1(i,jrelm)*u(i,nrb,lra)*fact
                             over2=over2+xmult
                             fact=6.d0-fact
                          enddo
                          v_cc(2*nra+1,2*nrb+1,jrelm)=
     &   (over2+u(0,nra,lra)*u(0,nrb,lra)*vjp1(0,jrelm)
     & +u(nstep,nra,lra)*u(nstep,nrb,lra)*vjp1(nstep,jrelm))
     &  *rstepr/3.d0
                          v_cc(2*nrb+1,2*nra+1,jrelm)=
     &                    v_cc(2*nra+1,2*nrb+1,jrelm)

                      elseif (lra.eq.jrelm) then

                          over1=0.d0
                          over2=0.d0
                          over4=0.d0
                          fact=4.d0
                          do i=1,nstep-1
                             xmult=u(i,nra,lra)*u(i,nrb,lra)*fact
                             over1=over1+xmult*vuncp(i,2*lra+1)
                             over2=over2+xmult*vuncp(i,2*lra+2)
                             over4=over4+xmult*vjp1(i,lra-1)
                             fact=6.d0-fact
                          enddo

                          v_uc(nra,nrb,2*lra+1)=
     &   (over1+u(0,nra,lra)*u(0,nrb,lra)*vuncp(0,2*lra+1)
     & +u(nstep,nra,lra)*u(nstep,nrb,lra)*vuncp(nstep,2*lra+1))
     &  *rstepr/3.d0
                          v_uc(nra,nrb,2*lra+2)=
     &   (over2+u(0,nra,lra)*u(0,nrb,lra)*vuncp(0,2*lra+2)
     & +u(nstep,nra,lra)*u(nstep,nrb,lra)*vuncp(nstep,2*lra+2))
     &  *rstepr/3.d0
                          v_cc(2*nra+1,2*nrb+1,lra-1)=
     &   (over4+u(0,nra,lra)*u(0,nrb,lra)*vjp1(0,lra-1)
     & +u(nstep,nra,lra)*u(nstep,nrb,lra)*vjp1(nstep,lra-1))
     &  *rstepr/3.d0
                          v_uc(nrb,nra,2*lra+1)=v_uc(nra,nrb,2*lra+1)
                          v_uc(nrb,nra,2*lra+2)=v_uc(nra,nrb,2*lra+2)
                          v_cc(2*nrb+1,2*nra+1,lra-1)=
     &                    v_cc(2*nra+1,2*nrb+1,lra-1)

                      else

                          over1=0.d0
                          over2=0.d0
                          over3=0.d0
                          over4=0.d0
                          fact=4.d0

                          do i=1,nstep-1
                             xmult=u(i,nra,lra)*u(i,nrb,lra)*fact
                             over1=over1+xmult*vuncp(i,2*lra+1)
                             over2=over2+xmult*vuncp(i,2*lra+2)
                             over3=over3+xmult*vjm1(i,lra+1)
                             over4=over4+xmult*vjp1(i,lra-1)
                             fact=6.d0-fact
                          enddo

                          v_uc(nra,nrb,2*lra+1)=
     &   (over1+u(0,nra,lra)*u(0,nrb,lra)*vuncp(0,2*lra+1)
     & +u(nstep,nra,lra)*u(nstep,nrb,lra)*vuncp(nstep,2*lra+1))
     &  *rstepr/3.d0
                          v_uc(nra,nrb,2*lra+2)=
     &   (over2+u(0,nra,lra)*u(0,nrb,lra)*vuncp(0,2*lra+2)
     & +u(nstep,nra,lra)*u(nstep,nrb,lra)*vuncp(nstep,2*lra+2))
     &  *rstepr/3.d0
                          v_cc(2*nra,2*nrb,lra+1)=
     &   (over3+u(0,nra,lra)*u(0,nrb,lra)*vjm1(0,lra+1)
     & +u(nstep,nra,lra)*u(nstep,nrb,lra)*vjm1(nstep,lra+1))
     &  *rstepr/3.d0
                          v_cc(2*nra+1,2*nrb+1,lra-1)=
     &   (over4+u(0,nra,lra)*u(0,nrb,lra)*vjp1(0,lra-1)
     & +u(nstep,nra,lra)*u(nstep,nrb,lra)*vjp1(nstep,lra-1))
     &  *rstepr/3.d0
                          v_uc(nrb,nra,2*lra+1)=v_uc(nra,nrb,2*lra+1)
                          v_uc(nrb,nra,2*lra+2)=v_uc(nra,nrb,2*lra+2)
                          v_cc(2*nrb,2*nra,lra+1)=
     &                    v_cc(2*nra,2*nrb,lra+1)
                          v_cc(2*nrb+1,2*nra+1,lra-1)=
     &                    v_cc(2*nra+1,2*nrb+1,lra-1)

                      endif

                   elseif(lrb.eq.lra-2) then
                      jrela=(lra+lrb)/2
                      ibc=2*nrb
                      iac=2*nra+1
                      over2=0.d0
                      fact=4.d0

                      do i=1,nstep-1
                         xmult=u(i,nra,lra)
     &                        *u(i,nrb,lrb)*fact
                         over2=over2+xmult*vcpl(i,jrela)
                         fact=6.d0-fact
                      enddo
                      v_cc(iac,ibc,jrela)=
     &   (over2+u(0,nra,lra)*u(0,nrb,lrb)*vcpl(0,jrela)
     & +u(nstep,nra,lra)*u(nstep,nrb,lrb)*vcpl(nstep,jrela))
     &  *rstepr/3.d0
                      v_cc(ibc,iac,jrela)=v_cc(iac,ibc,jrela)
                   elseif(lrb.eq.lra+2) then
                      jrela=(lra+lrb)/2
                      ibc=2*nrb+1
                      iac=2*nra
                      over2=0.d0
                      fact=4.d0

                      do i=1,nstep-1
                         xmult=u(i,nra,lra)
     &                        *u(i,nrb,lrb)*fact
                         over2=over2+xmult*vcpl(i,jrela)
                         fact=6.d0-fact
                      enddo
                      v_cc(iac,ibc,jrela)=
     &   (over2+u(0,nra,lra)*u(0,nrb,lrb)*vcpl(0,jrela)
     & +u(nstep,nra,lra)*u(nstep,nrb,lrb)*vcpl(nstep,jrela))
     &  *rstepr/3.d0
                      v_cc(ibc,iac,jrela)=v_cc(iac,ibc,jrela)
                   endif
               end do
            end do
         end do
      end do
      !v_cc=v_cc*complex(cos(1.355d0*theta),sin(-1.355d0*theta))
      !v_uc=v_uc*complex(cos(1.355d0*theta),sin(-1.355d0*theta))
      print *,' inter calculated'

      deallocate(u)
      deallocate(vuncp)
      deallocate(vjm1)
      deallocate(vjp1)
      deallocate(vcpl)

      if (chirp) av18=.false.
      contains

         complex(real128) function v_c(r,ipot,is)
         complex(real128), intent(in)::r
         complex(real128):: rmu1,rmu2,VMT1,VMT2
         integer,intent(in)::ipot,is
!            parameter (rmu1=3.11d0,rmu2=1.55d0)
!            parameter (VMT1=1458.05d0,VMT2=-578.09d0)

            rmu1=3.11d0
            rmu2=1.55d0
            VMT1=1438.720d0
         if(ipot==2) then
            if(is==0) then
               VMT2=-513.968d0
            else
               VMT2=-626.885d0
            endif
         else
            if(is==0) then
               VMT2=-509.4d0
               VMT2=-513.968d0
            else
               VMT2=-509.4d0
               VMT2=-626.885d0
            endif
         endif
            !rmu1=3.11d0
            !rmu2=1.55d0
            !VMT1=1458.05d0
            !VMT2=-626.885d0
               v_c=VMT1*exp(-rmu1*r)/(r)
     &          +VMT2*exp(-rmu2*r)/(r)
         end function v_c
         complex(real128) function vmn_c(r,is,it)
         complex(real128) r
         integer is,it
         complex(real128) rmu1,rmu2,rmu3,VMN1,VMN2,VMN3
         parameter (rmu1=1.487_kd,rmu2=0.639_kd,rmu3=0.465_kd)
         parameter (VMN1=200.0_kd,VMN2=-178.0_kd,VMN3=-91.85_kd)
         vmn_c=0.5_kd*VMN1*exp(-rmu1*r*r)
     &        *(1._kd-dble((-1)**(is+it)))
     &               +0.25_kd*VMN2*exp(-rmu2*r*r)
     &        *(1._kd+dble(-(-1)**is+(-1)**it-(-1)**(is+it)))
     &               +0.25_kd*VMN3*exp(-rmu3*r*r)
     &        *(1._kd+dble((-1)**is-(-1)**it-(-1)**(is+it)))
         end function vmn_c
      end subroutine inter_c_old
!**********************
!**********************
      subroutine inter_c
      use constants
      use paramdef
      use ist
      use pointers
      use nninteraction
      use harmosc
!      use potdef
      use gognypot
      use H_matrices
      use HO_wf
      use quadratures
      implicit double precision (a-h,o-z)
      integer,parameter :: kd=real128
      integer,parameter :: dp=kind(0.d0)
      parameter (intstep=500)
      parameter (nstep=intstep,rc=1.d-6,rinf=22.d0)

       double precision,allocatable,dimension(:,:,:):: u
      complex(kd),allocatable,dimension(:,:):: vuncp,
     & vjm1,vcpl,vjp1
      complex(kd) vpot(2,2),vnnjp1,vnpjp1,vppjp1,vnnjm1,vnpjm1
     & ,vppjm1,vnncpl,vnpcpl,vppcpl,vnn,vnp,vpp
      complex(kd):: wave,r,xmult,over1,over2,fact,over3,over4
      real(kd),allocatable:: roots(:,:,:),rstep
     &,rcr,rinfr,rstepr
      integer:: ipot,qn
      real(kd),allocatable:: gau_leg_poi(:)
     & ,gau_leg_we(:),mesh_we(:),v_x_mesh(:)

      if (chirp) av18=.true.

      if (av18) then
!*         write(2,*)
!*     +        'Argonne V18 NN potential in coordinate space'
         if (chirp) then
            write(2, "(' local chiral NN potential: lpot=',i3)") lpot
         else
            write(2,"(' Argonne V18 NN potential: lpot=',i3)") lpot
         endif
         if (nucleons==3.and.threeff) then
!*            anu=hbo*(dble(neut3eff)*fmsn+dble(iprot3eff)*fmsp)
!*     +           /dble(nucl3eff)/hbc**2
            rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
            fmbyh2=rdmavg/hbc**2
            anu=fmbyh2*hbo
            MTot2=neut3eff-iprot3eff
            ipp=(nucl3eff-Mtot2)*(nucl3eff-Mtot2-2)
            inp=nucl3eff*(nucl3eff-2)+itot23eff*(itot23eff+2)
     &           -2*Mtot2**2
            inn=(nucl3eff+Mtot2)*(nucl3eff+Mtot2-2)
            cpp=dble(ipp)
     &        /dble(3*nucl3eff*(nucl3eff-2)+iTot23eff*(iTot23eff+2))
            cnp=dble(inp)
     &        /dble(3*nucl3eff*(nucl3eff-2)+iTot23eff*(iTot23eff+2))
            cnn=dble(inn)
     &        /dble(3*nucl3eff*(nucl3eff-2)+iTot23eff*(iTot23eff+2))
         else
!*            anu=hbo*(dble(neutrons)*fmsn+dble(iprotons)*fmsp)
!*     +           /dble(nucleons)/hbc**2
            MTot2=neutrons-iprotons
            ipp=(nucleons-Mtot2)*(nucleons-Mtot2-2)
            inp=nucleons*(nucleons-2)+itot2*(itot2+2)-2*Mtot2**2
            inn=(nucleons+Mtot2)*(nucleons+Mtot2-2)
            cpp=dble(ipp)
     &           /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
            cnp=dble(inp)
     &           /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
            cnn=dble(inn)
     &           /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
!*******************
!            cnp=1.d0
!            cnn=0.d0
!            cpp=0.d0
            rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
            fmbyh2=rdmavg/hbc**2
            anu=fmbyh2*hbo
!*******************
         endif
         print *,' ipp,inn,inp:',ipp,inn,inp
         print *,' cpp,cnn,cnp:',cpp,cnn,cnp
      elseif (MTV) then
         write(2,*)
!     +        'Malfliet-Tjon NN potential in coordinate space'
     &        'Minnesota NN potential in coordinate space'
!     +        'Gogny NN potential in coordinate space'
         rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
         fmbyh2=rdmavg/hbc**2
!         anu=fmbyh2*hbo
         anu=hbo/41.47d0
         MTot2=neutrons-iprotons
         ipp=(nucleons-Mtot2)*(nucleons-Mtot2-2)
         inp=nucleons*(nucleons-2)+itot2*(itot2+2)-2*Mtot2**2
         inn=(nucleons+Mtot2)*(nucleons+Mtot2-2)
         cpp=dble(ipp)
     &           /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
         cnp=dble(inp)
     &           /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
         cnn=dble(inn)
     &           /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
         rcoul=cpp*hbc/alpha
         print *,' ipp,inn,inp:',ipp,inn,inp
         print *,' cpp,cnn,cnp:',cpp,cnn,cnp
      else
         write(2,*)
     &        'Argonne V8p NN potential in coordinate space'
         rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
         fmbyh2=rdmavg/hbc**2
         anu=fmbyh2*hbo
      endif

      write(2,1246) nrelm,lrelm,jrelm
 1246 format(' nrelm=',i4,'    lrelm=',i3,'    jrelm=',i3)


      bsquare=1.d0/anu
      write(2,1325) hbo,dsqrt(bsquare)
 1325 format(' hbar*Omega=',f8.4,'    b_HO=',f8.4)

      print *,' hbo3=',hbo3,'  intstep=',intstep,'  iucdim=',iucdim

      allocate(vuncp(intstep,iucdim))
      vuncp=0.d0
      RSTEP=(RINF-RC)/REAL(NSTEP)

      write(2,1457) rc,rinf,nstep,rstep
 1457 format(' rc=',e14.7,'     rinf=',f8.4,'     nstep=',i5,
     &                    '     rstep=',f8.4)
!*************************      
      if (allocated(roots)) deallocate(roots)
      allocate(roots(nrelm,nrelm,0:lrelm))
      do lra=0,lrelm
         call laguerre_rootfinder_kd(nrelm,lra+0.5_kd
     &                                 ,roots(1:nrelm,1:nrelm,lra))
      end do
      print*,"laguerre roots allocated"
!*************************      
      qn=intstep
      !deallocate(gau_leg_poi,gau_leg_we,v_x_mesh,mesh_we)
      allocate(gau_leg_poi(qn),v_x_mesh(qn),mesh_we(qn))
      allocate(gau_leg_we(qn))
      call gauleg2(qn,1e-8_kd,22._kd,gau_leg_poi,gau_leg_we)
      v_x_mesh(:)=gau_leg_poi
      mesh_we=gau_leg_we
!*************************      

      vnp=0.d0
      vpp=0.d0
      vnn=0.d0

      DO JPW=0,jrelm
         lpw=jpw
         do ispw=0,1
            if (ispw==1.and.jpw==0) lpw=1
            if ((-1)**(lpw+ispw)==1) then
               itpw=1
            else
               itpw=0
            endif
            iuc=2*jpw+ispw+1
            DO I=1,qn
               !R=RC+I*RSTEP!*complex(cos(theta),sin(theta))
               r=v_x_mesh(i)*complex(cos(theta),sin(theta))

               if (av18) then
                  print*,"not cs yet"
                  stop
                  if (inp/=0) then
                     CALL av18pw(lpot,lpw,ispw,jpw,itpw,1,-1,r,vpot)
                     vnp=vpot(1,1)
                  endif
                  if (ipp/=0.and.itpw==1) then
                     CALL av18pw(lpot,lpw,ispw,jpw,itpw,1,1,r,vpot)
                     vpp=vpot(1,1)
                  endif
                  if (inn/=0.and.itpw==1) then
                     CALL av18pw(lpot,lpw,ispw,jpw,itpw,-1,-1,r,vpot)
                     vnn=vpot(1,1)
                  endif
                  if (itpw==1) then
                     vuncp(i,iuc)=cnp*vnp+cpp*vpp+cnn*vnn
                  else
                     vuncp(i,iuc)=vnp
                  endif
               elseif (MTV) then
!                  VUNCP(I,iuc)=v(r)
                   vpp=v_c(r,1,ispw)
                   vnp=v_c(r,2,ispw)
                   vnn=v_c(r,3,ispw)
                  if (itpw==1) then
                     vuncp(i,iuc)=cnp*vnp+cpp*vpp+cnn*vnn
                  else
                     vuncp(i,iuc)=vnp
                  endif
                  if(lpw/=0) vuncp(i,iuc)=0
                  !vuncp(i,iuc)=vuncp(i,iuc)+dble(itpw)*rcoul/r
!                  VUNCP(I,iuc)=vmn_c(r,ispw,itpw)+dble(itpw)*rcoul/r
!                  CALL gognypotpw(lpw,ispw,jpw,itpw,r,vpot)
!                  VUNCP(I,iuc)=VPOT(1,1)
               else
                  print*,"not cs yet"
                  stop
                  CALL av8ppw(lpw,ispw,jpw,itpw,r,vpot)
                  VUNCP(I,iuc)=VPOT(1,1)
               endif

            end do
         ENDDO
      end do
      print *,' uncoupled potential computed'

      allocate(vjm1(intstep,jrelm))
      allocate(vcpl(intstep,jrelm))
      allocate(vjp1(intstep,jrelm))
      vjm1=0.d0
      vcpl=0.d0
      vjp1=0.d0

      vnpjm1=0.d0
      vnpcpl=0.d0
      vnpjp1=0.d0
      vppjm1=0.d0
      vppcpl=0.d0
      vppjp1=0.d0
      vnnjm1=0.d0
      vnncpl=0.d0
      vnnjp1=0.d0

      do jpw=1,jrelm
         lpw=jpw-1
         ispw=1
         if ((-1)**(lpw+ispw)==1) then
            itpw=1
         else
            itpw=0
         endif
         DO I=1,qn
            r=v_x_mesh(i)*complex(cos(theta),sin(theta))
               if (av18) then
                  print*,"not cs yet"
                  stop
                  if (inp/=0) then
                     CALL av18pw(lpot,lpw,ispw,jpw,itpw,1,-1,r,vpot)
                     vnpjm1=vpot(1,1)
                     vnpcpl=vpot(1,2)
                     vnpjp1=vpot(2,2)
                  endif
                  if (ipp/=0.and.itpw==1) then
                     CALL av18pw(lpot,lpw,ispw,jpw,itpw,1,1,r,vpot)
                     vppjm1=vpot(1,1)
                     vppcpl=vpot(1,2)
                     vppjp1=vpot(2,2)
                  endif
                  if (inn/=0.and.itpw==1) then
                     CALL av18pw(lpot,lpw,ispw,jpw,itpw,-1,-1,r,vpot)
                     vnnjm1=vpot(1,1)
                     vnncpl=vpot(1,2)
                     vnnjp1=vpot(2,2)
                  endif
                  if (itpw==1) then
                     vjm1(i,jpw)=cnp*vnpjm1+cpp*vppjm1+cnn*vnnjm1
                     vcpl(i,jpw)=cnp*vnpcpl+cpp*vppcpl+cnn*vnncpl
                     vjp1(i,jpw)=cnp*vnpjp1+cpp*vppjp1+cnn*vnnjp1
                  else
                     vjm1(i,jpw)=vnpjm1
                     vcpl(i,jpw)=vnpcpl
                     vjp1(i,jpw)=vnpjp1
                  endif
               elseif (MTV) then
!                  vjm1(i,jpw)=v(r)
!                  vjp1(i,jpw)=v(r)
                  vppjm1=v_c(r,1,ispw)
                  vnpjm1=v_c(r,2,ispw)
                  vnnjm1=v_c(r,3,ispw)
                  vppjp1=v_c(r,1,ispw)
                  vnpjp1=v_c(r,2,ispw)
                  vnnjp1=v_c(r,3,ispw)
                  if (itpw==1) then
                     vjm1(i,jpw)=cnp*vnpjm1+cpp*vppjm1+cnn*vnnjm1
                     vjp1(i,jpw)=cnp*vnpjp1+cpp*vppjp1+cnn*vnnjp1
                  else
                     vjm1(i,jpw)=vnpjm1
                     vjp1(i,jpw)=vnpjp1
                  endif
                  vjp1(i,jpw)=0! take only s wave
                  if(jpw>1) then
                        vjm1(i,jpw)=0
                        vjp1(i,jpw)=0
                  endif
                  
                  vcpl(i,jpw)=0.d0
                  !vjm1(i,jpw)=vjm1(i,jpw)+dble(itpw)*rcoul/r
!                  vjm1(i,jpw)=vmn_c(r,ispw,itpw)+dble(itpw)*rcoul/r
!                  vjp1(i,jpw)=vmn_c(r,ispw,itpw)+dble(itpw)*rcoul/r
!                  CALL gognypotpw(lpw,ispw,jpw,itpw,r,vpot)
!                  vjm1(i,jpw)=VPOT(1,1)
!                  vcpl(i,jpw)=VPOT(1,2)
!                  vjp1(i,jpw)=VPOT(2,2)
               else
                  print*,"not cs yet"
                  stop
                  CALL av8ppw(lpw,ispw,jpw,itpw,r,vpot)
                  vjm1(i,jpw)=VPOT(1,1)
                  vcpl(i,jpw)=VPOT(1,2)
                  vjp1(i,jpw)=VPOT(2,2)
               endif
         ENDDO
      end do

      print *,' coupled potential computed'
      print *,' rc=',rc,' rinf=',rinf,' nstep=',nstep

      allocate(u(nstep,0:nrelm,0:lrelm))
      rcr=rc/dsqrt(2.d0)
      rinfr=rinf/dsqrt(2.d0)
      rstepr=(rinfr-rcr)/real(nstep)
*      do lr=0,lrelm
*         do nr=0,nrelm
*            if (2*nr+lr>2*nrelm) cycle
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,rr,wave)
!$OMP& SCHEDULE(DYNAMIC)
         DO I=1,qn
            rr=v_x_mesh(i)/sqrt(2.d0)
!*               call waver(nr,lr,anu,rr,wave)
               call waver_gen(nrelm,lrelm,u(i,:,:),anu,rr)
!*               u(i,nr,lr)=wave
            end do
!$OMP END PARALLEL DO
*         end do
*      end do

         bsquare_kd=1.d0/anu
         anu_kd=anu
!         do lra=0,lrelm
!            do nra=0,nrelm
!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,pi,wave)
!!$OMP& SCHEDULE(DYNAMIC)
!               do i=0,nstep
!               rr=rcr+i*rstepr
!!                  pp=gau_leg_poi(i)
!                  !call waver(nr,lr,bsquare,pp,wave)
!                  call wave_func_c(nra,lra,bsquare_kd,rr,wave
!     &,roots(:max(nra,1),max(nra,1),lra),.false.)
!                  u(i,nra,lra)=wave*dble((-1)**nra)
!               end do
!!$OMP end parallel do
!            end do
!         end do
      print *,' u calculated'

      !allocate(vcc(0:iccdim,0:iccdim,jrelm))
      !allocate(vuc(0:nrelm,0:nrelm,iucdim))

      !vcc=0.d0
      !vuc=0.d0

      do nra=0,nrelm
         do lra=0,jrelm+1
            nna=nra+nra+lra
            if (nna>2*nrelm) cycle
            do nrb=nra,nrelm
               if (lra.eq.0) then
                   lrbmin=0
               elseif (lra.eq.1) then
                   lrbmin=1
               else
                   lrbmin=lra-2
               endif
               do lrb=lrbmin,min(lra+2,jrelm+1),2
                  if (2*nrb+lrb>2*nrelm) cycle
                  if (lrb.eq.lra) then

                      if (lra.eq.0) then

                          over1=0.d0
                          over2=0.d0
                     do i=1,qn
                        pi=v_x_mesh(i)
                        over1=over1+
     &                    +vuncp(i,1)
     &                    *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                        over2=over2+
     &                    +vjm1(i,1)
     &                    *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                     end do
                          v_uc(nra,nrb,1)=over1
                          v_cc(2*nra,2*nrb,1)=over2
                          v_uc(nrb,nra,1)=v_uc(nra,nrb,1)
                          v_cc(2*nrb,2*nra,1)=v_cc(2*nra,2*nrb,1)

                      elseif (lra.eq.1) then

                          over1=0.d0
                          over2=0.d0
                          over3=0.d0
                          over4=0.d0
                          do i=1,qn
                             pi=v_x_mesh(i)
                             over1=over1+
     &                         +vuncp(i,2)
     &                         *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                             over2=over2+
     &                         +vjm1(i,2)
     &                         *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                             over3=over3+
     &                         +vuncp(i,3)
     &                         *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                             over4=over4+
     &                         +vuncp(i,4)
     &                         *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                          end do

                          v_uc(nra,nrb,2)=over1
                          v_cc(2*nra,2*nrb,2)=over2
                          v_uc(nra,nrb,3)=over3
                          v_uc(nra,nrb,4)=over4
                          v_uc(nrb,nra,2)=v_uc(nra,nrb,2)
                          v_cc(2*nrb,2*nra,2)=v_cc(2*nra,2*nrb,2)
                          v_uc(nrb,nra,3)=v_uc(nra,nrb,3)
                          v_uc(nrb,nra,4)=v_uc(nra,nrb,4)

                      elseif (lra.eq.jrelm+1) then

                          over2=0.d0

                          do i=1,qn
                             pi=v_x_mesh(i)
                             over2=over2+
     &                         +vjp1(i,2)
     &                         *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                          end do
                          v_cc(2*nra+1,2*nrb+1,jrelm)=over2
                          v_cc(2*nrb+1,2*nra+1,jrelm)=
     &                    v_cc(2*nra+1,2*nrb+1,jrelm)

                      elseif (lra.eq.jrelm) then

                          over1=0.d0
                          over2=0.d0
                          over4=0.d0
                          do i=1,qn
                             pi=v_x_mesh(i)
                             over1=over1+
     &                         +vuncp(i,2*lra+1)
     &                         *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                             over2=over2+
     &                         +vuncp(i,2*lra+2)
     &                         *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                             over4=over4+
     &                         +vjp1(i,lra-1)
     &                         *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                          end do

                          v_uc(nra,nrb,2*lra+1)=over1
                          v_uc(nra,nrb,2*lra+2)=over2
                          v_cc(2*nra+1,2*nrb+1,lra-1)=over4
                          v_uc(nrb,nra,2*lra+1)=v_uc(nra,nrb,2*lra+1)
                          v_uc(nrb,nra,2*lra+2)=v_uc(nra,nrb,2*lra+2)
                          v_cc(2*nrb+1,2*nra+1,lra-1)=
     &                    v_cc(2*nra+1,2*nrb+1,lra-1)

                      else

                          over1=0.d0
                          over2=0.d0
                          over3=0.d0
                          over4=0.d0

                          do i=1,qn
                             pi=v_x_mesh(i)
                             over1=over1+
     &                         +vuncp(i,2*lra+1)
     &                         *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                             over2=over2+
     &                         +vuncp(i,2*lra+2)
     &                         *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                             over3=over3+
     &                         +vjm1(i,lra+1)
     &                         *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                             over4=over4+
     &                         +vjp1(i,lra-1)
     &                         *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                          end do
                          v_uc(nra,nrb,2*lra+1)=over1
                          v_uc(nra,nrb,2*lra+2)=over2
                          v_cc(2*nra,2*nrb,lra+1)=over3
                          v_cc(2*nra+1,2*nrb+1,lra-1)=over4
                          v_uc(nrb,nra,2*lra+1)=v_uc(nra,nrb,2*lra+1)
                          v_uc(nrb,nra,2*lra+2)=v_uc(nra,nrb,2*lra+2)
                          v_cc(2*nrb,2*nra,lra+1)=
     &                    v_cc(2*nra,2*nrb,lra+1)
                          v_cc(2*nrb+1,2*nra+1,lra-1)=
     &                    v_cc(2*nra+1,2*nrb+1,lra-1)

                      endif

                   elseif(lrb.eq.lra-2) then
                      jrela=(lra+lrb)/2
                      ibc=2*nrb
                      iac=2*nra+1
                      over2=0.d0

                      do i=1,qn
                         pi=v_x_mesh(i)
                         over2=over2+
     &                     +vcpl(i,jrela)
     &                     *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                      end do
                      v_cc(iac,ibc,jrela)=over2
                      v_cc(ibc,iac,jrela)=v_cc(iac,ibc,jrela)
                   elseif(lrb.eq.lra+2) then
                      jrela=(lra+lrb)/2
                      ibc=2*nrb+1
                      iac=2*nra
                      over2=0.d0
                      do i=1,qn
                         pi=v_x_mesh(i)
                         over2=over2+
     &                     +vcpl(i,jrela)
     &                     *mesh_we(i)*u(i,nra,lra)*u(i,nrb,lrb)
                      end do
                      v_cc(iac,ibc,jrela)=over2
                      v_cc(ibc,iac,jrela)=v_cc(iac,ibc,jrela)
                   endif
               end do
            end do
         end do
      end do

      !v_cc=v_cc*complex(cos(1.d0*theta),sin(-1.d0*theta))
      !v_uc=v_uc*complex(cos(1.d0*theta),sin(-1.d0*theta))
      print *,' inter calculated'

      deallocate(u)
      deallocate(vuncp)
      deallocate(vjm1)
      deallocate(vjp1)
      deallocate(vcpl)

      if (chirp) av18=.false.
      contains

         complex(real128) function v_c(r,ipot,it)
         complex(real128), intent(in)::r
         complex(real128):: rmu1,rmu2,VMT1,VMT2
         integer,intent(in)::ipot,it
!            parameter (rmu1=3.11d0,rmu2=1.55d0)
!            parameter (VMT1=1458.05d0,VMT2=-578.09d0)

            rmu1=3.11d0
            rmu2=1.55d0
            VMT1=1438.720d0
         if(ipot==2) then
            if(it==0) then
               VMT2=-513.968d0
            else
               VMT2=-626.885d0
            endif
         else
            if(it==0) then
               VMT2=-509.4d0
            else
               VMT2=-509.4d0
               VMT2=-626.885d0
            endif
         endif
            !rmu1=3.11d0
            !rmu2=1.55d0
            !VMT1=1458.05d0
            !VMT2=-626.885d0
               v_c=VMT1*exp(-rmu1*r)/(r)
     &          +VMT2*exp(-rmu2*r)/(r)
         end function v_c
         complex(real128) function vmn_c(r,is,it)
         complex(real128) r
         integer is,it
         complex(real128) rmu1,rmu2,rmu3,VMN1,VMN2,VMN3
         parameter (rmu1=1.487_kd,rmu2=0.639_kd,rmu3=0.465_kd)
         parameter (VMN1=200.0_kd,VMN2=-178.0_kd,VMN3=-91.85_kd)
         vmn_c=0.5_kd*VMN1*exp(-rmu1*r*r)
     &        *(1._kd-dble((-1)**(is+it)))
     &               +0.25_kd*VMN2*exp(-rmu2*r*r)
     &        *(1._kd+dble(-(-1)**is+(-1)**it-(-1)**(is+it)))
     &               +0.25_kd*VMN3*exp(-rmu3*r*r)
     &        *(1._kd+dble((-1)**is-(-1)**it-(-1)**(is+it)))
         end function vmn_c
      end subroutine inter_c
!**********************
      end module v2bA_quadrature
