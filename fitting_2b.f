! This module contains the subroutines that can be used to create the 2-body Hamiltonian
! matrix elements, using the fitting approach.
! where the interaction is approximated with a set of functions of the form:
!                 fk= c0 pi^bi pj^bj exp(-ci pi*pi/2) exp(-cj pj*pj/2)
! in particular:
!    subroutine generate_interaction_mesh: generates the n3lo interaction and stores
!              its values over a mesh in a file, that can be used by a python script
!              as the data to be fitted with the set of function.
!    subroutine inter_2b_fit: is used to calculate the matrix elements using the fit
!              parameters that it reads from files.

      module fitting_2b
      use iso_fortran_env

      integer::Nf,v_n=100
      real(real128)::prange=6
         contains

!**********************
      subroutine inter_2b_fit_quadrature(ipotin)
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
      use HO_wf
      use quadratures
      implicit none
      integer,parameter :: kd=real128
!*      parameter (nstep=500,pc=30.d0,pinf=0.000000001)
      complex(kd) :: vnp(6),vpp(6),vnn(6),sumv(6)  
      !file variables
      character(len=5):: string1,string2,string3,string4
     & ,string5,string6
      character(len=100):: v_file_name,string7
      logical:: filetest=.false.
      !timer variables
      real:: time1,time2
      integer:: ir,time3,time4,time5,time6
      !parameters 
      integer:: Nmax,Nd,n_threads,J
      real(kd):: fmbyh2,hbc3,pfmev,rdmavg,bosc,th_kd
      !temporary variables 
      integer:: i,k,ij,i_f,ia,ib,ind,nra,nrb,lra,lrb,jrel
      real(kd) ::pi,pk,rrel
      real(kd)::c,c0,c2,x0,x2,p0,ci,ck,ni,nk,e_kd,Hmax
      complex(kd),allocatable::H_part1(:,:),H_part2(:,:)
      complex(kd),allocatable,dimension(:,:,:,:,:):: H0

      integer:: ipotin,icceff,inn,inp,ipp,ispin,iuc,iucmax,lbmax
      integer::lbmin,lrbm,lrbmin,mtot2,nrbmax,nna,nreffmax,sra
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
!----------------------------------------------------------------------------
      !common block for n3lo subroutine
      real(kind(0.d0))::ipot
      common /cnn/ ipot
      interface
         subroutine mkl_zgeev(A,n,m,W,Z)
            integer, intent(in) :: n,m
            complex(kind(0.d0)),intent(in) :: A(n,n)
            complex(kind(0.d0)), intent(out) :: W(n),Z(n,n)
         end subroutine mkl_zgeev
      end interface
!----------------------------------------------------------------------------
      open(223,file="inter_2b_fit.out")
      call system_clock(count=time5, count_rate=ir)
!***********************
      !initialize variables
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
           string7="n3lo"
      elseif (icutnum==102) then
         write(2,*)
     &        'NNLO(POUNDerS) NN potential in momentum space'
      elseif (icutnum==105) then
         write(2,*)
     &        'NNLOsat NN potential in momentum space'
           string7="nnlo_sat"
      elseif (icutnum==106) then
         write(2,*)
     &        'N4LO500 NN potential in momentum space'
           string7="n4lo500"
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
     &,roots(:max(nra,1),max(nra,1),lra),.false.)
                  u4(i,nra,lra)=wave*dble((-1)**nra)
     & *complex(cos(theta/2),sin(theta/2))
               end do
!$OMP end parallel do
            end do
         end do
         theta=theta2
         print*,"wave function created"
!*****************************
      print*,"hbc3",hbc3
      print*,"bsquare", bsquare_kd
!***********************
!******************************
      allocate(H_part1(0:nrelm,0:lrelm))
      allocate(H_part2(0:nrelm,0:lrelm))
      if(allocated(H0))deallocate(H0)
      allocate(H0(0:nrelm,0:lrelm,0:nrelm,0:lrelm,6))
      vnp=0.d0
      vpp=0.d0
      vnn=0.d0
      !loop over J
      do jrel=0,min(jrelm,mxeng+1)
      !***************
       print *,' doing J=',jrel
       write(223,*),"***********************"
       write(223,*),' J=',jrel
       H0=complex(0._kd,0._kd)
!*****************************
       !loop over interaction channels
       do ij=1,6
       !*****************
         !read from file
        print*,"**************************"
        print*,"starting channel ",ij
        write(223,*),"starting channel ",ij
        write(string2,'(i5)') int (jrel)
        write(string1,'(i5)') int (ij)
        
        write(string3,'(i5)')int(Nf)
        write(string4,'(i5)')int(v_n)
        write(string6,'(F5.2)') prange 
        select case(int(ipot))
         case(1) 
           string5="pp"
         case(2)
           string5="np"
         case(3)
           string5="nn"
         case default
           print*,"illegal value for ipot"
           stop
        end select
        string1=adjustl(string1)
        string2=adjustl(string2)
        string3=adjustl(string3)
        string4=adjustl(string4)
        string5=adjustl(string5)
        string6=adjustl(string6)
        string7=adjustl(string7)
        !v_file_name="2b_fitting/fit_parameters/v_fit3_params_"//
        v_file_name="2b_fitting/lmfit_parameters/v_lmfit_params_v5_"//
     &trim(string5)//"_"//trim(string7)//
     &"_v_n"//trim(string4)//"_range"//trim(string6)//"_J"//
     &trim(string2)//"_ch"//
     &trim(string1)//"_Nf"//trim(string3)//'.fit'
!     &trim(string5)//
!     &"_n3lo_v_n"//trim(string4)//"_J"//trim(string2)//"_ch"//
!     &trim(string1)//"_Nf"//trim(string3)//'.fit'
        inquire(file=v_file_name,exist=filetest)
        if(.not. filetest) then
           print*,"Error,potential fit parameters not found"
           print*,"file:",v_file_name
           stop
        endif
        print*,"reading the parameters from:",v_file_name
        open(222,file=v_file_name)
        read(222,*) Nf
        !loop over fitting functions
        do i_f=1,Nf
        !******************
         read(222,*) c,x0,x2,c0,c2
            c=c*10._kd**-6
            write(223,"('function ',I3,4X,'parameters',4X,ES14.3
     &                        ,4F14.4)"),i_f,c,x0,x2,c0,c2
         p0=c
         ci=2._kd/(c0*c0)
         ck=2._kd/(c2*c2)
         ni=ci*x0*x0
         nk=ck*x2*x2
         if(x0/=0) then
            p0=p0*((x0*x0/e_kd)**(-ni/2._kd))
            !p0=p0*((ni/ci/e_kd)**(-ni/2._kd))
         endif
         if(x2/=0) then
            p0=p0*((x2*x2/e_kd)**(-nk/2._kd))
         endif
            write(223,"('function ',I3,4X,'parameters',4X,ES14.3
     &                        ,4F14.4)"),i_f,p0,ni,nk,ci,ck
         Hmax=0
!*****************************
         lbmin=jrel
         lbmax=jrel
         H_part1=complex(0._kd,0._kd)
         H_part2=complex(0._kd,0._kd)
         if(ij>2) then
            lbmin=max(0,jrel-1)
            lbmax=min(lrelm,jrel+1)
         endif
         call cpu_time(time1)
         call system_clock(count=time3, count_rate=ir)
!$OMP parallel do collapse(2) 
!$OMP& private(vrel,integ,reg,regcut,regpow,i,nra,lra,pi)
         do nra=0,nrelm
            do lra=lbmin,lbmax
               H_part1(nra,lra)=0._kd
               H_part2(nra,lra)=0._kd
!**************************
!                call HO_integral_arb(0,nra,0,lra,1/ci*bsquare_kd
!     & ,th_kd,real(lra+1+ni,kd),temp,temp3)
!               sumv(1)=complex(get_float128(temp3(1)),
!     &                get_float128(temp3(2)))
!               sumv(1)=sumv(1)*sqrt(2._kd)
!     & *((bsquare_kd/ci)
!     & **(lra/2.0_kd))
!               sumv(1)=sumv(1)*((0.5/ci)**(1.5_kd))
!     &*complex(cos(th_kd),sin(th_kd))**(lra+1.5_kd)
!               H_part1(nra,lra)=sumv(1)*(ci**(-ni/2._kd+0.75_kd))
!**************************
!     &*complex(cos(0.5_kd*th_kd),sin(-0.5_kd*th_kd))
               integ=0
               regcut=20.0_kd
               regpow=7
               do i=1,qn!0,N_arr(nra,lra)
                  !pi=mesh(i,nra,lra)
                  pi=v_x_mesh(i)
                  vrel=sqrt(abs(p0))*(pi*complex(cos(th_kd)
     &,sin(-th_kd)))**ni!pp/step/ck
     &   * exp(pi**2*complex(cos(2._kd*th_kd),sin(-2._kd*th_kd)) 
     &*-0.5_kd*ci)

               reg=3._kd/(1._kd+ exp((pi*complex(cos(th_kd)
     &,sin(-th_kd)) /regcut)**regpow*-1._kd)
     &                     + exp((pi*complex(cos(th_kd),sin(-th_kd))
     &  /regcut)**regpow))
                  vrel=vrel*reg
                  !vrel=vrel*exp(sin(-th_kd)*pi)
                  vrel=vrel!*v_scale*complex(cos(v_theta),sin(v_theta))
                  integ=integ+vrel*pi
     &                   *u4(i,nra,lra)*mesh_we(i)
               enddo
               integ=integ*complex(cos(1.5_kd*th_kd)
     &                               ,sin(-1.5_kd*th_kd))
               H_part1(nra,lra)=integ
!!**************************
!                call HO_integral_arb(0,nra,0,lra,1/ck*bsquare_kd
!     & ,th_kd,real(lra+1+nk,kd),temp,temp3)
!               sumv(1)=complex(get_float128(temp3(1)),
!     &                get_float128(temp3(2)))
!               sumv(1)=sumv(1)*sqrt(2._kd)
!     & *((bsquare_kd/ck)
!     & **(lra/2.0_kd))
!               sumv(1)=sumv(1)*((0.5/ck)**(1.5_kd))
!     &*complex(cos(th_kd),sin(th_kd))**(lra+1.5_kd)
!               H_part2(nra,lra)=sumv(1)*(ck**(-nk/2._kd+0.75_kd))
!                   if(abs(H_part1(nra,lra))>Hmax) then
!                      Hmax=abs(H_part1(nra,lra))
!                   endif
               integ=0
               do i=1,qn!0,N_arr(nra,lra)
                  !pi=mesh(i,nra,lra)
                  pi=v_x_mesh(i)
                  vrel=sqrt(abs(p0))*(pi*complex(cos(th_kd)
     &,sin(-th_kd)))**nk!pp/step/ck
     &   * exp(pi**2*complex(cos(2._kd*th_kd),sin(-2._kd*th_kd)) 
     &*-0.5_kd*ck)
               reg=3._kd/(1._kd+ exp((pi*complex(cos(th_kd)
     &,sin(-th_kd)) /regcut)**regpow*-1._kd)
     &                     + exp((pi*complex(cos(th_kd),sin(-th_kd))
     &  /regcut)**regpow))
                  vrel=vrel*reg
                  integ=integ+vrel*pi
     &                   *u4(i,nra,lra)*mesh_we(i)
               enddo
               integ=integ*complex(cos(1.5_kd*th_kd)
     &                               ,sin(-1.5_kd*th_kd))
               H_part2(nra,lra)=integ
            enddo
         enddo
!$OMP end parallel do
         call system_clock(count=time4, count_rate=ir)
         call cpu_time(time2)
         Hmax=0
         lrbm=max(0,jrel-1)
!***********************************************
         call cpu_time(time1)
         call system_clock(count=time3, count_rate=ir)
!$OMP PARALLEL DO DEFAULT(SHARED) collapse(2) num_threads(1)
!$OMP& PRIVATE(nra,lra,nrb,lrb,sumv)
!$OMP& SCHEDULE(dynamic)
         do nra=0,nrelm
            do nrb=0,nrelm
               do lra=max(0,jrel-1),min(jrel+1,jrelm+1)
                  do lrb=max(0,jrel-1),min(jrel+1,jrelm+1)

         !*********************************
                     if ((-1)**(lra+lrb)/=1) cycle
               lrbmin=max(0,jrel-1)
               sumv=cmplx(0.d0,0.d0)
         !*********************************
         ! coupled integral
!                sumv(ij)=H_part1(nra,lra)*H_part2(nrb,lrb)*c1
                if(ij<5) then               
                sumv(ij)=(H_part1(nra,lra)*H_part2(nrb,lrb)+
     &    H_part1(nrb,lrb)*H_part2(nra,lra))*0.5_kd
                else
                sumv(ij)=H_part1(nra,lra)*H_part2(nrb,lrb)
                endif
                if(p0<0) sumv(ij)=sumv(ij)*-1._kd
         !*********************************
                if(abs(sumv(ij))*hbc3>Hmax) then
                   Hmax=abs(sumv(ij))*hbc3
                endif
                     H0(nra,lra,nrb,lrb,ij)=H0(nra,lra,nrb,lrb,ij)
     &                           +sumv(ij)*real(hbc3,kd)!+vrel
         !*********************************
                  enddo
               enddo
               !print*,"nra,nrb ",nra,nrb
            enddo
         enddo
!$OMP END PARALLEL DO
         call system_clock(count=time4, count_rate=ir)
         call cpu_time(time2)
!***********************************************
         call cpu_time(time1)
         enddo
         enddo

         if(Jrel>0) then
            do nra=0,nrelm
               do nrb=0,nrelm
                  lra=jrel-1
                  lrb=jrel+1
                   if ((-1)**(lra+lrb)/=1) cycle
                   sumv(1)=(H0(nra,lra,nrb,lrb,6)+H0(nrb,lrb,nra
     &                                            ,lra,5))*0.5_kd
                   H0(nra,lra,nrb,lrb,6)=sumv(1)
                   H0(nrb,lrb,nra,lra,5)=sumv(1)
               enddo
            enddo
         endif
         print*,"H calculated"

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
      print *,' inter calculated'
      deallocate(H_part1,H_part2)
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
            call wave_func_c(nra,lr,bsquare_kd,pi,wave
     & ,roots(:max(nra,1),
     &               max(nra,1),lr),.false.)
            u1=wave*(-1)**nra
            area1=area1+u1*pi*gau_leg_we(i)*length!*vrel
     & *exp(-(pi/7)**6)
         enddo
!******************************
         length=p_e-p_m
         do i=1,n
            pi=gau_leg_poi(i)*(length)+p_m
            call wave_func_c(nra,lr,bsquare_kd,pi,wave
     &,roots(:max(nra,1),
     &               max(nra,1),lr),.false.)
            u1=wave*(-1)**nra
            area2=area2+u1*pi*gau_leg_we(i)*length!*vrel
     & *exp(-(pi/7)**6)
         enddo
!******************************
         if(abs(area1+area2-area)>tol*(p_e-p_b)/20._kd) then
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
      end subroutine inter_2b_fit_quadrature

!**********************
      subroutine add_residual_interaction(ipotin)
      
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
      integer::Nmax,Nd,n_threads
      real(kd):: fmbyh2,hbc3,pfmev,rdmavg,bosc,th_kd
      !temporary variables 
      integer:: i,k,ij,i_f,ia,ib,ind,nra,nrb,lra,lrb,jrel
      real(kd) ::pi,pk,rrel
      real(kd)::c,c0,c2,x0,x2,p0,ci,ck,ni,nk,e_kd
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
      complex(kd):: wave,integ,reg,vrel,vrel1,vrel2
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
      th_kd=0
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
      allocate(m_temp(0:nrelm,0:lrelm,80))
      m_temp=0
      NN=0
      allocate(mesh(0:5000,0:nrelm,0:lrelm),N_arr(0:nrelm,0:lrelm))
      allocate(step_arr(0:5000,0:nrelm,0:lrelm))
      call gauleg2(50,0._kd,1._kd,gau_leg_poi,gau_leg_we)
         call cpu_time(time1)
         call system_clock(count=time3, count_rate=ir)
!$OMP parallel do num_threads(1)
!$OMP& private(nra,lra,totn,interr,sumv)
!$OMP& schedule(dynamic)
      do nra=0,nrelm
         do lra=0,lrelm
            totn=0
            interr=0
            wave=adapt_integ_leg(nra,nra,lra,0._kd,
     &   15._kd,complex(-100._kd,0._kd),50,totn,interr)
            N_arr(nra,lra)=totn
            NN=max(NN,totn)
          if((mod(nra,20)>0 .or. nra==10) .and. lra<3) then
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
         do lra=0,lrelm
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
               !vint(i,k,:)=-100*exp(-pi-pk)/hbc3
               else
               vint(i,k,:)=0.d0
               endif
            end do
         end do

         print *,' potential calculated'
!*****************************************
      !Subtract the fitted form of the interaction
       print *,' doing J=',jrel
       write(223,*),"***********************"
       write(223,*),' J=',jrel
!*****************************
       !loop over interaction channels
       do ij=1,6
       !*****************
         !read from file
        print*,"**************************"
        print*,"starting channel ",ij
        write(223,*),"starting channel ",ij
        write(string2,'(i5)') int (jrel)
        write(string1,'(i5)') int (ij)
        
        write(string3,'(i5)')int(Nf)
        write(string4,'(i5)')int(v_n)
        write(string6,'(F5.2)') prange 
        select case(int(ipot))
         case(1) 
           string5="pp"
         case(2)
           string5="np"
         case(3)
           string5="nn"
         case default
           print*,"illegal value for ipot"
           stop
        end select
        string1=adjustl(string1)
        string2=adjustl(string2)
        string3=adjustl(string3)
        string4=adjustl(string4)
        string5=adjustl(string5)
        string6=adjustl(string6)
        !v_file_name="2b_fitting/fit_parameters/v_fit3_params_"//
        v_file_name="2b_fitting/lmfit_parameters/v_lmfit_params_"//
     &trim(string5)//
     &"_n3lo_v_n"//trim(string4)//"_range"//trim(string6)//"_J"//
     &trim(string2)//"_ch"//
     &trim(string1)//"_Nf"//trim(string3)//'.fit'
!     &trim(string5)//
!     &"_n3lo_v_n"//trim(string4)//"_J"//trim(string2)//"_ch"//
!     &trim(string1)//"_Nf"//trim(string3)//'.fit'
        inquire(file=v_file_name,exist=filetest)
        if(.not. filetest) then
           print*,"Error,potential fit parameters not found"
           print*,"file:",v_file_name
           stop
        endif
        print*,"reading the parameters from:",v_file_name
        open(222,file=v_file_name)
        read(222,*) Nf
        !loop over fitting functions
        do i_f=1,Nf
        !******************
         read(222,*) c,x0,x2,c0,c2
            c=c*10._kd**-6
            write(223,"('function ',I3,4X,'parameters',4X,ES14.3
     &                        ,4F14.4)"),i_f,c,x0,x2,c0,c2
         c0=max(c0,x0/2.4)
         c2=max(c2,x2/2.4)
         p0=c
         ci=2._kd/(c0*c0)
         ck=2._kd/(c2*c2)
         ni=ci*x0*x0
         nk=ck*x2*x2
         if(x0/=0) then
            p0=p0*((x0*x0/e_kd)**(-ni/2._kd))
            !p0=p0*((ni/ci/e_kd)**(-ni/2._kd))
         endif
         if(x2/=0) then
            p0=p0*((x2*x2/e_kd)**(-nk/2._kd))
         endif
            write(223,"('function ',I3,4X,'parameters',4X,ES14.3
     &                        ,4F14.4)"),i_f,p0,ni,nk,ci,ck
         
         regcut=6.7_kd
         regpow=7
         do i=1,qn
            pi=v_x_mesh(i)
            vrel1=p0*(pi*complex(cos(th_kd)
     &,sin(-th_kd)))**ni
     &   * exp(pi**2*complex(cos(2._kd*th_kd),sin(-2._kd*th_kd)) 
     &*-0.5_kd*ci)
           
!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(k,pk,vrel,vrel2,reg)
!$OMP& SCHEDULE(static)
            do k=1,qn

               pk=v_x_mesh(k)

               vrel2=(pk*complex(cos(th_kd)
     &,sin(-th_kd)))**nk
     &   * exp(pk**2*complex(cos(2._kd*th_kd),sin(-2._kd*th_kd))
     &*-0.5_kd*ck)
               vrel=vrel1*vrel2
               !regulator
               reg=3._kd/(1._kd+ exp((pi*complex(cos(th_kd)
     &,sin(-th_kd)) /regcut)**regpow*-1._kd)
     &                     + exp((pi*complex(cos(th_kd),sin(-th_kd))
     &  /regcut)**regpow))
     &            *3._kd/(1._kd+ exp((pk*complex(cos(th_kd)
     &,sin(-th_kd)) /regcut)**regpow*-1._kd)
     &                     + exp((pk*complex(cos(th_kd),sin(-th_kd))
     &  /regcut)**regpow))
               vrel=vrel*reg
               !
               if(ij<5) then 
               vint(i,k,ij)=vint(i,k,ij)-vrel/2
               vint(k,i,ij)=vint(k,i,ij)-vrel/2
               else
               vint(i,k,ij)=vint(i,k,ij)-vrel/2
               vint(k,i,11-ij)=vint(k,i,11-ij)-vrel/2
               endif

            end do
!$OMP END PARALLEL DO
         end do

         enddo!i_f
         enddo!ij
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
         real(kd)::p_b,p_e,interr
         integer::n,nra,nrb,lr,totn
         complex(kd)::area


         real(kd):: tol
         real(kd)::p_m,step,pi,cf,tmp2,length,vrel
         complex(kd)::wave,u1,u2,res,area1,area2
         integer:: i,reste,tmp1
         area1=0._kd
         area2=0._kd
         tol=10._kd**(-25._kd)
         if(nra<70) tol=10._kd**(-26._kd)
         if(nra<40) tol=10._kd**(-28._kd)
         if(nra<10) tol=10._kd**(-30._kd)
         tol=10._kd**(-13._kd)
         p_m=(p_b+p_e)/2._kd
!******************************
         if((p_e-p_b)<10._kd**-28) then
         print*,nra,lr,p_e,p_b,p_e-p_b
            area1=0
            area2=0
            length=p_e-p_b
            do i=1,n
               pi=gau_leg_poi(i)*(length)+p_b
               call wave_func_c(nra,lr,bsquare_kd,pi,wave,roots(:nra,
     &               nra,lr),.false.)
               u1=wave*(-1)**nra
               print*,u1,(u1-area2)/(pi-area1)!*vrel
               area1=pi
               area2=u1
!     & *exp(-(pi/5._kd)**7._kd)
!     & *u1
            enddo
            stop
         endif
!*****************************
         length=p_m-p_b
         do i=1,n
            pi=gau_leg_poi(i)*(length)+p_b
            call wave_func_c(nra,lr,bsquare_kd,pi,wave,roots(:nra,
     &               nra,lr),.false.)
            u1=wave*(-1)**nra
            area1=area1+u1*pi*gau_leg_we(i)*length!*vrel
     & *exp(-(pi/5._kd)**7._kd)
!     & *u1
         enddo
!******************************
         length=p_e-p_m
         do i=1,n
            pi=gau_leg_poi(i)*(length)+p_m
            call wave_func_c(nra,lr,bsquare_kd,pi,wave,roots(:nra,
     &               nra,lr),.false.)
            u1=wave*(-1)**nra
            area2=area2+u1*pi*gau_leg_we(i)*length!*vrel
     & *exp(-(pi/5._kd)**7._kd)
!     & *u1
         enddo
!******************************
         if(abs(area1+area2-area)>tol*(p_e-p_b)/20._kd
     &.and.abs(area1+area2-area)>tol) then
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



      end subroutine add_residual_interaction



!**********************
!*******************************
      ! the interaction is scaled by 10^6 to make the fit more accurate
      ! the amplitude retrieved by the fitting procedure must be scaled by 10^_6 to compensate.
      subroutine generate_interaction_mesh(v_n,j_rel_max,p_range,
     & icutnum,ipotin)
         use constants
         use iso_fortran_env
         implicit none
         integer,intent(in):: v_n,j_rel_max,ipotin,icutnum
         integer,parameter :: dp=kind(0.d0)

         !*********
         real(dp):: pfmev,pi,pk,v_x_mesh(v_n),v_mesh(v_n,v_n,6)
         real(dp):: p_range
         integer:: i,k,jrel
         real:: time1,time2
         character(len=10)::string1,string2,string3,iso_name
     & ,string4
         character(len=100):: int_name,file_name
         !*********
         !common block for n3lo subroutine
         real(dp):: KREAD,KWRITE,KPUNCH,XMEV,YMEV
     &,ipot,V(6),KDA(9)
         integer::J
         character*4 label
         LOGICAL HEFORM,SING,TRIP,COUP,ENDEP
         COMMON /CRDWRT/ KREAD,KWRITE,KPUNCH,KDA
!        ARGUMENTS AND VALUES OF CDBONN:
         COMMON /CPOT/   V,XMEV,YMEV
         COMMON /CSTATE/ J,HEFORM,SING,TRIP,COUP,ENDEP,LABEL
         common /cnn/ ipot
         KREAD=5
         KWRITE=6
         HEFORM=.FALSE.
         SING=.TRUE.
         TRIP=.TRUE.
         COUP=.TRUE.
         !*********

         ipot=ipotin
         call cpu_time(time1)
         pfmev=hbc/dsqrt(2.d0)
         print*,"v_n, jrelm,p_range",v_n, j_rel_max,p_range
      !*****************
      !determining the interaction requested
         if (icutnum==100) then
            int_name='n3lo'
         elseif (icutnum==102) then
            int_name="nnlo"
         elseif (icutnum==105) then
            int_name="nnlo_sat"
         elseif (icutnum==106) then
            int_name="n4lo500"
         elseif (icutnum==107) then
            int_name="n3lo500new"
         elseif (icutnum==108) then
            int_name="n2lo500"
         elseif (icutnum==109) then
            int_name="nlo500"
         elseif (icutnum==110) then
            int_name="lo500"
         elseif (icutnum==115) then
            int_name="nnlo_input"
         endif
         print*,"using interaction: ",int_name
         do i=1,v_n
            v_x_mesh(i)=i*(p_range/v_n)
         enddo
         if(ipotin==2) iso_name="np"
         if(ipotin==1) iso_name="pp"
         if(ipotin==3) iso_name="nn"
      !*****************
         do jrel=0,j_rel_max
            J=jrel
            print *,' doing J=',int(J)
            !*****************
            do i=1,v_n
               pi=v_x_mesh(i)
               xmev=pi*pfmev
               do k=1,v_n
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
                  v_mesh(i,k,:)=v(:)
                  else
                  v_mesh(i,k,:)=0
                  endif
               end do
            end do
            !*****************
            write(string3,'(i3)'),int(v_n)
            string3=adjustl(string3)
            write(string4,'(f5.2)'),p_range
            string4=adjustl(string4)
            open(222,file="2b_fitting/fit_input/fit-mesh_v_n"//
     &trim(string3)//"_prange"//trim(string4)
     &//".dat")
            do i=1,v_n
               write(222,"(F14.7)"),v_x_mesh(i)
            enddo
            close(222)
            do k=1,6
            write(string1,'(i2)'),int(k)
            write(string2,'(i2)'),int(jrel)
            string1=adjustl(string1)
            string2=adjustl(string2)
            file_name="2b_fitting/fit_input/interaction_mesh_"//
     &trim(iso_name)
     &//"_"//trim(int_name)//"_v_n"//trim(string3)//"_range"//
     &trim(string4)//"_J"//trim(string2)
     &//"_ch"//trim(string1)//".dat"
            print*,"writing to file: ",file_name
            open(222,file=file_name)
            do i=1,v_n
             write(222,"(10000F14.7)"),(10.d0**6*v_mesh(i,j,k),j=1,v_n)
            enddo
            close(222)
            enddo
         enddo
      !*****************
         call cpu_time(time2)
         print *,'2-body interaction mesh generated',time2-time1
   
      end subroutine generate_interaction_mesh
      end module fitting_2b
