      module SRG_module
      use iso_fortran_env
      logical :: srg=.false.,NN_SRG_for3b_extract=.false.,
     &     NNN_SRG=.false.,NNNN_SRG=.false.
      real(kind(0.d0)) :: lambda,delta_E=-1.d0
!!     changed (Conor-Nuiok) - start
      integer :: dimcut,delta_N1=0,Ngen=-1,dimgen=0
      real(kind(0.d0)),allocatable :: H(:,:),G(:,:),DYDS1(:,:),
     &     Hgen(:,:),ETA1(:,:),DYDS2(:,:),ETA2(:,:),DYDS(:,:),
     &     ETA(:,:),H_s(:),TkinSave(:,:)
      complex(kind(0.d0)),allocatable:: HH(:,:),GG(:,:), DDYDS1(:,:),
     &     HHGEN(:,:),EETA1(:,:),DDYDS2(:,:),EETA2(:,:),DDYDS(:,:),
     &     EETA(:,:), HH_s(:),TTkinSave(:,:),ckin(:,:)
      integer :: ena,enb,Ndiag
      logical :: bswitch=.false.
!!     changed (Conor-Nuiok) - finish
      real(kind(0.d0)) :: ATOL=1.d-6,RTOL=1.d-7
      complex(real128),allocatable :: HH_kd(:,:),GG_kd(:,:)
     & ,DDYDS1_kd(:,:),HHGEN_kd(:,:),EETA1_kd(:,:),DDYDS2_kd(:,:)
     &,EETA2_kd(:,:),DDYDS_kd(:,:),EETA_kd(:,:),HH_s_kd(:)
     &,TTkinSave_kd(:,:)
      complex(kind(0.d0)),allocatable:: GG2(:,:)
      real(kind(0.d0)):: srg_f=1 ! proportion of T in the generator
      integer:: block_size=60
      contains
!********************************
!!!!!! by E. Jurgenson; modified
!*******************************
      !***********************************
      subroutine csrg_evolution(dim,H_evolve,Tkin,lambda_in)
      use paramdef, only: theta
      use constants
      use dvodef90              ! an ode solver
      implicit none
      integer,intent(IN) :: dim
      integer,parameter:: dp=kind(1.d0)
      complex(dp),intent(IN) :: Tkin(dim,dim)
      complex(dp),intent(INOUT) :: H_evolve(dim,dim)
      real(kind(0.d0)),intent(in),optional:: lambda_in
      integer :: i,j,k
      double precision :: rdmavg,lam_scale,lambda_sc, RSTATS(22)
      complex(dp),allocatable:: Zwork(:)
      double precision, allocatable:: Rwork(:)
      integer, allocatable:: Iwork(:)
      integer:: Lzw,Lrw,Liw,mf
!! ODE vars
      type(VODE_OPTS) :: OPTIONS
      real(dp) :: T, TOUT
      integer :: ISTATE, ITASK, NEQ, ISTATS(31)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!! set some ODE constants for the solver routines
      itask = 1                 ! or maybe 4, see dvodef90.f90 comments for details
      istate = 1                ! initial call, hereafter, this  should return as 2
      OPTIONS = SET_OPTS(METHOD_FLAG=10,ABSERR=ATOL,RELERR=RTOL,
     &     MXSTEP=50000)
      dimcut=dim
      NEQ=dim*(dim+1)/2
      rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
      lam_scale = hbc/sqrt(rdmavg) !hbar_sq_over_Mn ! make sure the srg units are right
      allocate(HH(dim,dim),GG(dim,dim),DDYDS1(dim,dim),EETA1(dim,dim),
     &     DDYDS2(dim,dim),EETA2(dim,dim),DDYDS(dim,dim),EETA(dim,dim))
      allocate(HH_s(NEQ))
      allocate(GG2(dim,dim))

!     changed (Conor) - start
      GG=Tkin
      allocate(ckin(dim,dim))
      Ckin(:,:)=Tkin(:,:)*complex(cos(2.d0*theta),sin(-2.d0*theta))
!     G=Hgen
!     changed (Conor) - finish
      do i=1,dim
         do j=i,dim
            k=i+j*(j-1)/2
            HH_s(k)=H_evolve(j,i)
         end do
      end do
      T=0
      lambda_sc = lambda*lam_scale
      if(present(lambda_in))then
         if(lambda_in>0.1) lambda_sc = lambda_in*lam_scale
      endif!lambda_sc=
      TOUT = 1/(lambda_sc*lambda_sc*lambda_sc*lambda_sc) !TOUT = 1/(lambda**4) 
      write(*,*) 's_in: ',T,'   s_out: ',TOUT,'   lam out: ',lambda
     &                                        ,' srg_f: ',srg_f
! new osama
      mf=10 ! depends on whether the problem is 'stiff', see documentation of ZVODE
      ! the follwing lengths should depend on MF, see documentation of ZVODE
      Lzw=15*NEQ 
      Lrw=20+NEQ
      Liw= 30
      allocate(Zwork(Lzw),Rwork(Lrw),Iwork(Liw))
      Iwork=0
      Rwork=0
      Zwork=0
      Iwork(6)=150000 ! the maximum number of steps allowed in the zvode routine
      ! incase of an error "maximum number of steps taken before reaching TOUT",
      ! either increase Iwork(6) or increase lambda
! end
      print *,' ZVODE_F90 called; ISTATE,matrix dim=',ISTATE,dim
      if(dimcut<305)then
      ATOL=1e-4
      RTOL=1e-5
      !ATOL=1e-6
      !RTOL=1e-7
      print *,' Atol, Rtol =',Atol,Rtol
!      call ZVODE(cSRG_RHS,NEQ,HH_s,T,TOUT,1,RTOL,ATOL,ITASK,ISTATE
!     & ,1,Zwork,Lzw,Rwork,Lrw,Iwork,Liw,srg_jac,mf )
       call ZVODE(cSRG_RHS_modified,NEQ,HH_s,T,TOUT,1,RTOL,ATOL
     & ,ITASK,ISTATE,1,Zwork,Lzw,Rwork,Lrw,Iwork,Liw,srg_jac,mf )
      else
      ATOL=1e-4
      RTOL=1e-5
      !ATOL=1e-6
      !RTOL=1e-7
      print *,' Atol, Rtol =',Atol,Rtol
!     block_size=max(int(dimcut/50),20)
      print*,"block size: ",block_size
      call ZVODE(cSRG_RHS_block,NEQ,HH_s,T,TOUT,1,RTOL,ATOL,ITASK,ISTATE
     & ,1,Zwork,Lzw,Rwork,Lrw,Iwork,Liw,srg_jac,mf )
      endif
      WRITE(6,60) Iwork(11),Iwork(12),Iwork(13),Iwork(20), 
     &     Iwork(21),Iwork(22),Iwork(23)
 60   FORMAT(/,'  No. steps =',I6,'   No. f-s =',I6, 
     &     '  No. J-s =',I4,'   No. LU-s =',I4,/,
     &     '  No. nonlinear iterations =',I4,/,        
     &     '  No. nonlinear convergence failures =',I4,/, 
     &     '  No. error test failures =',I4,/)
!      print*,' last step size',Rwork(11),Rwork(12)
      if (ISTATE/=2) then
         print *,'***error in DVODE_F90'
         stop
      endif
      do i=1,dim
         do j=i,dim
            k=i+j*(j-1)/2
            H_evolve(j,i)=HH_s(k)
            H_evolve(i,j)=H_evolve(j,i)
         end do
      end do
      deallocate(HH_s)
      deallocate(HH,GG,DDYDS1,EETA1,DDYDS2,EETA2,DDYDS,EETA)
      deallocate(ckin)
      deallocate(GG2)
      end subroutine csrg_evolution


      subroutine cSRG_RHS (NEQ, TIN, Y, YDOT)  

! counter variables
      integer :: i,j,k,NEQ,time3,time4,ir,i1,i2
! these are the inputs and outputs
      real(kind(1.d0)) :: TIN,sum1,sum2
      complex(kind(1.d0)) :: Y(NEQ), YDOT(NEQ)
!write(*,*) 'unfolding Y and RPAR'  
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j !dimcut
            k = i+j*(j-1)/2     ! Petr's half of the matrix
!            HH(i,j)=Y(k)         ! Petr's half of the matrix
!            HH(j,i)=Y(k)         ! get the other half too
            HH(i,j)=Y(k)!-ckin(i,j)       
            HH(j,i)=HH(i,j)
         end do
      end do
!$OMP end PARALLEL DO
! do the first matrix multiplication
      EETA1=0.d0
!*********************
!      if(dimcut<100000) then
!         call system_clock(count=time3, count_rate=ir)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) Collapse(2)
!$OMP& schedule(static)
         do  j = 1,dimcut
            do i = 1,dimcut
               !EETA1(i,j) = SUM(GG(:,i)*HH(:,j))
               EETA1(i,j) = SUM(GG(:,i)*HH(:,j))
            end do
         end do
!$OMP end PARALLEL DO
!         call system_clock(count=time4, count_rate=ir)
!         print *,' multiplication'
!     &   ,real(time4 - time3,kind=8)/real(ir,kind=8)
!*********************
!      else
!         !call system_clock(count=time3, count_rate=ir)
!         call zsymm('L','U',dimcut,dimcut,complex(1.d0,0.d0),GG,dimcut
!     &                     ,HH,dimcut,complex(0.d0,0.d0),EETA1,dimcut)
!         !call system_clock(count=time4, count_rate=ir)
!!               print *,'zsymm'
!!     &   ,real(time4 - time3,kind=8)/real(ir,kind=8)
!      endif
!*********************
      EETA2 = TRANSPOSE(EETA1)
      EETA = EETA1 - EETA2
      DDYDS1=0.d0
      sum1=0
      !do i=1,dimcut
      !   do j=1,dimcut
      !      sum1=max(sum1,abs(EETA(i,j)))
      !   enddo
      !enddo
!*********************
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j !dimcut
            k = i+j*(j-1)/2     ! Petr's half of the matrix
            HH(i,j)=Y(k)         ! Petr's half of the matrix
            HH(j,i)=Y(k)         ! get the other half too
         end do
      end do
!$OMP end PARALLEL DO
!*********************
!      if(dimcut<1000000) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) Collapse(2)
!$OMP& schedule(static)
         do  j = 1,dimcut
            do i = 1,dimcut
               DDYDS1(i,j) = -SUM(EETA(:,i)*HH(:,j))
            end do
         end do
!$OMP end PARALLEL DO
!*********************
!      else
!         call zsymm('L','U',dimcut,dimcut,complex(1.d0,0.d0),HH,dimcut
!     &                   ,EETA,dimcut,complex(0.d0,0.d0),DDYDS1,dimcut)
!      endif
      DDYDS2 = TRANSPOSE(DDYDS1)
      DDYDS = DDYDS1 + DDYDS2      ! because \eta = -\eta^\dagger --> H*eta = - (eta*H)^\dagger
      sum2=0
      !do i=1,dimcut
      !   do j=1,dimcut
      !      sum2=max(sum2,abs(DDYDS(i,j)))
      !   enddo
      !enddo
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j  !dimcut
            k = i+j*(j-1)/2
            YDOT(k) = DDYDS(i,j)
         end do
      end do
!$OMP end PARALLEL DO
      !sum1=0
      !sum2=0
      !do i=1,NEQ
      !   sum1=max(sum1,abs(Y(i)))
      !   sum2=max(sum2,abs(YDOT(i)))
      !enddo

      !print*,"iteration, T=",TIN,sum1,sum2
      end subroutine cSRG_RHS

      subroutine cSRG_RHS_block (NEQ, TIN, Y, YDOT)  

! counter variables
      integer :: i,j,k,NEQ,time3,time4,ir,i1,i2,ii,jj,kk
! these are the inputs and outputs
      real(kind(1.d0)) :: TIN,sum1,sum2
      complex(kind(1.d0)) :: Y(NEQ), YDOT(NEQ)
!write(*,*) 'unfolding Y and RPAR'
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j !dimcut
            k = i+j*(j-1)/2     ! Petr's half of the matrix
!            HH(i,j)=Y(k)         ! Petr's half of the matrix
!            HH(j,i)=Y(k)         ! get the other half too
            HH(i,j)=Y(k)!-ckin(i,j)
            HH(j,i)=HH(i,j)
         end do
      end do
!$OMP end PARALLEL DO
! do the first matrix multiplication
      EETA1=0.d0
!*********************
!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) Collapse(2)
!!$OMP& schedule(static)
!      do  j = 1,dimcut
!         do i = 1,dimcut
!            !EETA1(i,j) = SUM(GG(:,i)*HH(:,j))
!            EETA1(i,j) = SUM(GG(:,i)*HH(:,j))
!         end do
!      end do
!!$OMP end PARALLEL DO
      EETA1=0
      !
        GG2(:,:)=srg_f*GG(:,:)
        do j=1,dimcut
           do i=j,j
              if(abs(GG(j,i))>1e-4) then
           GG2(j,i)=GG2(j,i)+(1-srg_f)*conjg(HH(j,i))
              endif
           enddo
        enddo
      !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,ii,jj,kk) Collapse(2)
!$OMP& schedule(static)
      do  j = 1,dimcut,block_size
         do i = 1,dimcut,block_size
            do k=1,dimcut,block_size
               do jj=j,min(j+block_size-1,dimcut)
                  do ii=i,min(i+block_size-1,dimcut)
                     kk=min(k+block_size-1,dimcut)
                     EETA1(ii,jj) = EETA1(ii,jj)
     &                  +SUM(GG2(k:kk,ii)*HH(k:kk,jj))
                  enddo
               enddo
            enddo
         end do
      end do
!$OMP end PARALLEL DO
!*********************
      EETA2 = TRANSPOSE(EETA1)
      EETA = EETA1 - EETA2
      DDYDS1=0.d0
      !do i=1,dimcut
      !   do j=1,dimcut
      !      sum1=max(sum1,abs(EETA(i,j)))
      !   enddo
      !enddo
!*********************
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j !dimcut
            k = i+j*(j-1)/2     ! Petr's half of the matrix
            HH(i,j)=Y(k)         ! Petr's half of the matrix
            HH(j,i)=Y(k)         ! get the other half too
         end do
      end do
!$OMP end PARALLEL DO
!*********************
      DDYDS1=0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,ii,jj,kk) Collapse(2)
!$OMP& schedule(static)
      do  j = 1,dimcut,block_size
         do i = 1,dimcut,block_size
            do k=1,dimcut,block_size
               do jj=j,min(j+block_size-1,dimcut)
                  do ii=i,min(i+block_size-1,dimcut)
                     kk=min(k+block_size-1,dimcut)
                     DDYDS1(ii,jj) = DDYDS1(ii,jj)
     &                  -SUM(EETA(k:kk,ii)*HH(k:kk,jj))
                  enddo
               enddo
            enddo
         end do
      end do
!$OMP end PARALLEL DO
!*********************
      DDYDS2 = TRANSPOSE(DDYDS1)
      DDYDS = DDYDS1 + DDYDS2      ! because \eta = -\eta^\dagger --> H*eta = - (eta*H)^\dagger
      sum2=0
      !do i=1,dimcut
      !   do j=1,dimcut
      !      sum2=max(sum2,abs(DDYDS(i,j)))
      !   enddo
      !enddo
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j  !dimcut
            k = i+j*(j-1)/2
            YDOT(k) = DDYDS(i,j)
         end do
      end do
!$OMP end PARALLEL DO
      !sum1=0
      !sum2=0
      !do i=1,NEQ
      !   sum1=max(sum1,abs(Y(i)))
      !   sum2=max(sum2,abs(YDOT(i)))
      !enddo

      !print*,"iteration, T=",TIN,sum1,sum2
      end subroutine cSRG_RHS_block
!**********************************
      subroutine cSRG_RHS_modified (NEQ, TIN, Y, YDOT)  

! counter variables
      integer :: i,j,k,NEQ,time3,time4,ir,i1,i2
! these are the inputs and outputs
      real(kind(1.d0)) :: TIN,sum1,sum2
      complex(kind(1.d0)) :: Y(NEQ), YDOT(NEQ)
!write(*,*) 'unfolding Y and RPAR'  
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j !dimcut
            k = i+j*(j-1)/2     ! Petr's half of the matrix
!            HH(i,j)=Y(k)         ! Petr's half of the matrix
!            HH(j,i)=Y(k)         ! get the other half too
            HH(i,j)=Y(k)!-ckin(i,j)       
            HH(j,i)=HH(i,j)
         end do
      end do
!$OMP end PARALLEL DO
! do the first matrix multiplication
      EETA1=0.d0
!*********************
!      if(dimcut<100000) then
!         call system_clock(count=time3, count_rate=ir)
        GG2(:,:)=srg_f*GG(:,:)
        do j=1,dimcut
           do i=1,dimcut
              if(abs(GG(j,i))>1e-4) then
              !if(abs(i-j)<3) then
           GG2(j,i)=srg_f*GG(j,i)+(1.d0-srg_f)*conjg(HH(j,i))
              endif
           enddo
        enddo
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) Collapse(2)
!$OMP& schedule(static)
         do  j = 1,dimcut
            do i = 1,dimcut
               !EETA1(i,j) = SUM(GG(:,i)*HH(:,j))
               EETA1(i,j) = SUM(GG2(:,i)*HH(:,j))
            end do
         end do
!$OMP end PARALLEL DO
!         call system_clock(count=time4, count_rate=ir)
!         print *,' multiplication'
!     &   ,real(time4 - time3,kind=8)/real(ir,kind=8)
!*********************
!      else
!         !call system_clock(count=time3, count_rate=ir)
!         call zsymm('L','U',dimcut,dimcut,complex(1.d0,0.d0),GG,dimcut
!     &                     ,HH,dimcut,complex(0.d0,0.d0),EETA1,dimcut)
!         !call system_clock(count=time4, count_rate=ir)
!!               print *,'zsymm'
!!     &   ,real(time4 - time3,kind=8)/real(ir,kind=8)
!      endif
!*********************
      EETA2 = TRANSPOSE(EETA1)
      EETA = EETA1 - EETA2
      DDYDS1=0.d0
      sum1=0
      !do i=1,dimcut
      !   do j=1,dimcut
      !      sum1=max(sum1,abs(EETA(i,j)))
      !   enddo
      !enddo
!*********************
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j !dimcut
            k = i+j*(j-1)/2     ! Petr's half of the matrix
            HH(i,j)=Y(k)         ! Petr's half of the matrix
            HH(j,i)=Y(k)         ! get the other half too
         end do
      end do
!$OMP end PARALLEL DO
!*********************
!      if(dimcut<1000000) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) Collapse(2)
!$OMP& schedule(static)
         do  j = 1,dimcut
            do i = 1,dimcut
               DDYDS1(i,j) = -SUM(EETA(:,i)*HH(:,j))
            end do
         end do
!$OMP end PARALLEL DO
!*********************
!      else
!         call zsymm('L','U',dimcut,dimcut,complex(1.d0,0.d0),HH,dimcut
!     &                   ,EETA,dimcut,complex(0.d0,0.d0),DDYDS1,dimcut)
!      endif
      DDYDS2 = TRANSPOSE(DDYDS1)
      DDYDS = DDYDS1 + DDYDS2      ! because \eta = -\eta^\dagger --> H*eta = - (eta*H)^\dagger
      sum2=0
      !do i=1,dimcut
      !   do j=1,dimcut
      !      sum2=max(sum2,abs(DDYDS(i,j)))
      !   enddo
      !enddo
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j  !dimcut
            k = i+j*(j-1)/2
            YDOT(k) = DDYDS(i,j)
         end do
      end do
!$OMP end PARALLEL DO
      !sum1=0
      !sum2=0
      !do i=1,NEQ
      !   sum1=max(sum1,abs(Y(i)))
      !   sum2=max(sum2,abs(YDOT(i)))
      !enddo

      !print*,"iteration, T=",TIN,sum1,sum2
      end subroutine cSRG_RHS_modified
!**********************************
      subroutine cSRG_RHS_modified2 (NEQ, TIN, Y, YDOT)  

! counter variables
      integer :: i,j,k,NEQ,time3,time4,ir,i1,i2
! these are the inputs and outputs
      real(kind(1.d0)) :: TIN,sum1,sum2
      complex(kind(1.d0)) :: Y(NEQ), YDOT(NEQ),phase1,phase2
!write(*,*) 'unfolding Y and RPAR'  
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j !dimcut
            k = i+j*(j-1)/2     ! Petr's half of the matrix
!            HH(i,j)=Y(k)         ! Petr's half of the matrix
!            HH(j,i)=Y(k)         ! get the other half too
            HH(i,j)=Y(k)!-ckin(i,j)       
            HH(j,i)=HH(i,j)
         end do
      end do
!$OMP end PARALLEL DO
! do the first matrix multiplication
      EETA1=0.d0
!*********************
!      if(dimcut<100000) then
!         call system_clock(count=time3, count_rate=ir)
        GG2(:,:)=srg_f*GG(:,:)
        do j=1,dimcut
           do i=1,dimcut
              if(abs(GG(j,i))>1e-4) then
           GG2(j,i)=srg_f*GG(j,i)+(1.d0-srg_f)*conjg(HH(j,i))
              endif
           enddo
        enddo
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) Collapse(2)
!$OMP& schedule(static)
         do  j = 1,dimcut
            do i = 1,dimcut
               !EETA1(i,j) = SUM(GG(:,i)*HH(:,j))
               EETA1(i,j) = SUM(GG2(:,i)*HH(:,j))
            end do
         end do
!$OMP end PARALLEL DO
!         call system_clock(count=time4, count_rate=ir)
!         print *,' multiplication'
!     &   ,real(time4 - time3,kind=8)/real(ir,kind=8)
!*********************
!      else
!         !call system_clock(count=time3, count_rate=ir)
!         call zsymm('L','U',dimcut,dimcut,complex(1.d0,0.d0),GG,dimcut
!     &                     ,HH,dimcut,complex(0.d0,0.d0),EETA1,dimcut)
!         !call system_clock(count=time4, count_rate=ir)
!!               print *,'zsymm'
!!     &   ,real(time4 - time3,kind=8)/real(ir,kind=8)
!      endif
!*********************
      EETA2 = TRANSPOSE(EETA1)
      EETA = EETA1 - EETA2
      DDYDS1=0.d0
      sum1=0
      !do i=1,dimcut
      !   do j=1,dimcut
      !      sum1=max(sum1,abs(EETA(i,j)))
      !   enddo
      !enddo
!*********************
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j !dimcut
            k = i+j*(j-1)/2     ! Petr's half of the matrix
            HH(i,j)=Y(k)         ! Petr's half of the matrix
            HH(j,i)=Y(k)         ! get the other half too
         end do
      end do
!$OMP end PARALLEL DO
!*********************
!      if(dimcut<1000000) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) Collapse(2)
!$OMP& schedule(static)
         do  j = 1,dimcut
            do i = 1,dimcut
               DDYDS1(i,j) = -SUM(EETA(:,i)*HH(:,j))
            end do
         end do
!$OMP end PARALLEL DO
!*********************
!      else
!         call zsymm('L','U',dimcut,dimcut,complex(1.d0,0.d0),HH,dimcut
!     &                   ,EETA,dimcut,complex(0.d0,0.d0),DDYDS1,dimcut)
!      endif
      DDYDS2 = TRANSPOSE(DDYDS1)
      DDYDS = DDYDS1 + DDYDS2      ! because \eta = -\eta^\dagger --> H*eta = - (eta*H)^\dagger
      sum2=0
      !do i=1,dimcut
      !   do j=1,dimcut
      !      sum2=max(sum2,abs(DDYDS(i,j)))
      !   enddo
      !enddo
      GG2(:,:)=complex(0.d0,0.d0)
      do j=1,dimcut
         do i=1,dimcut
            if(abs(i-j)<4) cycle
            GG2(i,j)=DDYDS(i,j)
     &       *conjg(HH(i,j))*abs(i-j)!/(abs(HH(i,j))+1e-16)
         enddo
      enddo
      phase1=SUM(GG2(:,:))
      phase2=-1.d0* conjg(phase1)/(abs(phase1)+1e-16)
      DDYDS(:,:)=DDYDS(:,:)*phase2
      print*,TIN,phase1,phase2
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j  !dimcut
            k = i+j*(j-1)/2
            YDOT(k) = DDYDS(i,j)
         end do
      end do
!$OMP end PARALLEL DO
      !sum1=0
      !sum2=0
      !do i=1,NEQ
      !   sum1=max(sum1,abs(Y(i)))
      !   sum2=max(sum2,abs(YDOT(i)))
      !enddo

      !print*,"iteration, T=",TIN,sum1,sum2
      end subroutine cSRG_RHS_modified2
!*******************************
      !***********************************
      subroutine csrg_evolution_dynamic(dim_i,dim_f,H_evolve,Tkin
     &                                    ,n_stepsi,lambda_in)
      use paramdef, only: theta
      use constants
      use dvodef90              ! an ode solver
      implicit none
      integer,intent(IN) :: dim_i,dim_f,n_stepsi
      integer:: n_steps
      integer,parameter:: dp=kind(1.d0)
      complex(dp),intent(IN) :: Tkin(dim_i,dim_i)
      complex(dp),intent(INOUT) :: H_evolve(dim_i,dim_i)
      real(kind(0.d0)),intent(in),optional:: lambda_in
      integer :: i,j,k
      double precision :: rdmavg,lam_scale,lambda_sc, RSTATS(22)
!! ODE vars
      type(VODE_OPTS) :: OPTIONS
      real(dp) :: T, TOUT,TOUT_j,lambda_j,T_step
      complex(dp),allocatable:: H_temp(:,:),Tkin_temp(:,:)
      integer ::dim_j 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!! set some ODE constants for the solver routines
      n_steps=6
      rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
      lam_scale = hbc/sqrt(rdmavg) !hbar_sq_over_Mn ! make sure the srg units are right
      if(present(lambda_in))then
         if(lambda_in>0.1) lambda=lambda_in
      endif!lambda_sc=
      lambda_sc = lambda*lam_scale
      TOUT = 1/(lambda_sc*lambda_sc*lambda_sc*lambda_sc) !TOUT = 1/(lambda**4) 
      T=0
      write(*,*) 's_in: ',T,'   s_out: ',TOUT,'   lam out: ',lambda

      print*,"n_steps:",n_steps,"  dim_i:",dim_i,"  dim_f:",dim_f
      T_step=(TOUT-T)/n_steps
      dim_j=dim_i
      TOUT_j=0
      do i=1,n_steps
        if(i==1) then
           T_step=min(TOUT,1.8E-5)
        elseif(i==2) then
           T_step=(TOUT-TOUT_j)/(n_steps-1)
        endif
        TOUT_j=TOUT_j+T_step
        lambda_j=T_step**(-0.25d0)
        lambda_j=lambda_j/lam_scale
        allocate(H_temp(dim_j,dim_j),Tkin_temp(dim_j,dim_j))
        H_temp(:dim_j,:dim_j)=H_evolve(:dim_j,:dim_j)
        Tkin_temp(:dim_j,:dim_j)=Tkin(:dim_j,:dim_j)
        print*,"lambda_j:",lambda_j
        call csrg_evolution(dim_j,H_temp,Tkin_temp,lambda_j)
        H_evolve(:dim_j,:dim_j)=H_temp(:dim_j,:dim_j)
        deallocate(Tkin_temp,H_temp)
        dim_j= evolved_dim(i)
      enddo
      contains
        integer function evolved_dim(i)
           implicit none
           integer,intent(in):: i
           !evolved_dim=(dim_f-dim_i)/(TOUT-T)*T_step*i+dim_i
           evolved_dim=int(dim_f+(dim_i-dim_f)*exp(
     &                       -(1.2d0*TOUT_j/TOUT)**4.d0))
        end function evolved_dim


      end subroutine csrg_evolution_dynamic

!*******************************
      subroutine cSRG_RHS_test (NEQ, TIN, Y, YDOT)  

! counter variables
      integer :: i,j,k,NEQ,time3,time4,ir,i1,i2
! these are the inputs and outputs
      real(kind(1.d0)) :: TIN,sum1,sum2
      complex(kind(1.d0)) :: Y(NEQ), YDOT(NEQ)
!write(*,*) 'unfolding Y and RPAR' 
      !GG(:,:)=0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j !dimcut
            k = i+j*(j-1)/2     ! Petr's half of the matrix
            HH(i,j)=Y(k)!-ckin(i,j)       
            HH(j,i)=HH(i,j)
            if(j<0) then
               if(i==j) then
                  GG(i,j)=HH(i,j)
               else
                  GG(i,j)=0
               endif
            else
               if(GG(i,j)/=0) then
                  GG(i,j)=HH(i,j)
               else
                  GG(i,j)=0
               endif
            endif 
            GG(j,i)=GG(i,j)
         end do
      end do
!$OMP end PARALLEL DO
! do the first matrix multiplication
      EETA1=0.d0
!*********************
!         call system_clock(count=time3, count_rate=ir)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) Collapse(2)
!$OMP& schedule(static)
         do  j = 1,dimcut
            do i = 1,dimcut
               !EETA1(i,j) = SUM(GG(:,i)*HH(:,j))
               !EETA1(i,j) = SUM(GG(:,i)*HH(:,j))
               EETA1(i,j) = GG(i,i)*HH(i,j)
            end do
         end do
!$OMP end PARALLEL DO
!         call system_clock(count=time4, count_rate=ir)
!         print *,' multiplication'
!     &   ,real(time4 - time3,kind=8)/real(ir,kind=8)
!*********************
!*********************
      EETA2 = TRANSPOSE(EETA1)
      EETA = EETA1 - EETA2
      DDYDS1=0.d0
      sum1=0
      !do i=1,dimcut
      !   do j=1,dimcut
      !      sum1=max(sum1,abs(EETA(i,j)))
      !   enddo
      !enddo
!*********************
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j !dimcut
            k = i+j*(j-1)/2     ! Petr's half of the matrix
            HH(i,j)=Y(k)         ! Petr's half of the matrix
            HH(j,i)=Y(k)         ! get the other half too
         end do
      end do
!$OMP end PARALLEL DO
!*********************
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) Collapse(2)
!$OMP& schedule(static)
         do  j = 1,dimcut
            do i = 1,dimcut
               DDYDS1(i,j) = -SUM(EETA(:,i)*HH(:,j))
            end do
         end do
!$OMP end PARALLEL DO
!*********************
      DDYDS2 = TRANSPOSE(DDYDS1)
      DDYDS = DDYDS1 + DDYDS2      ! because \eta = -\eta^\dagger --> H*eta = - (eta*H)^\dagger
      sum2=0
      !do i=1,dimcut
      !   do j=1,dimcut
      !      sum2=max(sum2,abs(DDYDS(i,j)))
      !   enddo
      !enddo
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j  !dimcut
            k = i+j*(j-1)/2
            YDOT(k) = DDYDS(i,j)
         end do
      end do
!$OMP end PARALLEL DO
      !sum1=0
      !sum2=0
      !do i=1,NEQ
      !   sum1=max(sum1,abs(Y(i)))
      !   sum2=max(sum2,abs(YDOT(i)))
      !enddo

      !print*,"iteration, T=",TIN,sum1,sum2
      end subroutine cSRG_RHS_test
!*******************************
      subroutine csrg_evolution_kd(dim,H_evolve,Tkin,lambda_in)
      use iso_fortran_env
      use constants
      use dvodef90             ! an ode solver
      implicit none
      integer,intent(IN) :: dim
      integer,parameter:: kd=real128
      complex(kd),intent(IN) :: Tkin(dim,dim)
      complex(kd),intent(INOUT) :: H_evolve(dim,dim)
      real(kd),intent(in),optional:: lambda_in
      integer :: i,j,k
      real(kd) :: rdmavg,lam_scale,lambda_sc, RSTATS(22)
      complex(kd),allocatable:: Zwork(:)
      real(kd), allocatable:: Rwork(:)
      integer, allocatable:: Iwork(:)
      integer:: Lzw,Lrw,Liw,mf
!! ODE vars
      type(VODE_OPTS) :: OPTIONS
      real(kd) :: T, TOUT
      real(real128),parameter :: ATOL_kd=1.d-6,RTOL_kd=1.d-7
      integer :: ISTATE, ITASK, NEQ, ISTATS(31)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!! set some ODE constants for the solver routines
      itask = 1                 ! or maybe 4, see dvodef90.f90 comments for details
      istate = 1                ! initial call, hereafter, this  should return as 2
      OPTIONS = SET_OPTS(METHOD_FLAG=10,ABSERR=ATOL,RELERR=RTOL,
     &     MXSTEP=50000)
      dimcut=dim
      NEQ=dim*(dim+1)/2
      rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn)
      lam_scale = hbc/sqrt(rdmavg) !hbar_sq_over_Mn ! make sure the srg units are right
      allocate(HH_kd(dim,dim),GG_kd(dim,dim),DDYDS1_kd(dim,dim)
     &,EETA1_kd(dim,dim),DDYDS2_kd(dim,dim),EETA2_kd(dim,dim)
     & ,DDYDS_kd(dim,dim),EETA_kd(dim,dim))
      allocate(HH_s_kd(NEQ))
!     changed (Conor) - start
      GG_kd=Tkin
!     G=Hgen
!     changed (Conor) - finish
      do i=1,dim
         do j=i,dim
            k=i+j*(j-1)/2
            HH_s_kd(k)=H_evolve(j,i)
         end do
      end do
      T=0
      if(present(lambda_in))then
         if(lambda_in>0.1) lambda=lambda_in
      endif!lambda_sc=
      lambda_sc = lambda*lam_scale
      TOUT = 1/(lambda_sc*lambda_sc*lambda_sc*lambda_sc) !TOUT = 1/(lambda**4) 
      write(*,*) 's_in: ',T,'   s_out: ',TOUT,'   lam out: ',lambda
! new osama
      mf=10 ! depends on whether the problem is 'stiff', see documentation of ZVODE
      ! the follwing lengths should depend on MF, see documentation of ZVODE
      Lzw=15*NEQ 
      Lrw=20+NEQ
      Liw= 30
      allocate(Zwork(Lzw),Rwork(Lrw),Iwork(Liw))
      Iwork=0
      Rwork=0
      Zwork=0
      Iwork(6)=160000 ! the maximum number of steps allowed in the zvode routine
      call CVODE(cSRG_RHS_kd,NEQ,HH_s_kd,T,TOUT,1,RTOL_kd,ATOL_kd,ITASK
     & ,ISTATE,1,Zwork,Lzw,Rwork,Lrw,Iwork,Liw,srg_jac,mf )
      print *,' CVODE called; ISTATE,matrix dim=',ISTATE,dim
      WRITE(6,60) Iwork(11),Iwork(12),Iwork(13),Iwork(20), 
     &     Iwork(21),Iwork(22),Iwork(23)
 60   FORMAT(/,'  No. steps =',I5,'   No. f-s =',I5, 
     &     '  No. J-s =',I4,'   No. LU-s =',I4,/,
     &     '  No. nonlinear iterations =',I4,/,        
     &     '  No. nonlinear convergence failures =',I4,/, 
     &     '  No. error test failures =',I4,/)
      if (ISTATE/=2) then
         print *,'***error in DVODE_F90'
         stop
      endif
      do i=1,dim
         do j=i,dim
            k=i+j*(j-1)/2
            H_evolve(j,i)=HH_s_kd(k)
            H_evolve(i,j)=H_evolve(j,i)
         end do
      end do
      deallocate(HH_s_kd)
      deallocate(HH_kd,GG_kd,DDYDS1_kd,EETA1_kd,DDYDS2_kd,EETA2_kd
     &  ,DDYDS_kd,EETA_kd)
      end subroutine csrg_evolution_kd


      subroutine cSRG_RHS_kd (NEQ, TIN, Y, YDOT)  
      use iso_fortran_env

! counter variables
      integer :: i,j,k,NEQ,time3,time4,ir

! these are the inputs and outputs
      real(real128) :: TIN
      complex(real128) :: Y(NEQ), YDOT(NEQ)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j !dimcut
            k = i+j*(j-1)/2     ! Petr's half of the matrix
            HH_kd(i,j)=Y(k)         ! Petr's half of the matrix
            HH_kd(j,i)=Y(k)         ! get the other half too
         end do
      end do
!$OMP end PARALLEL DO
! do the first matrix multiplication
      EETA1_kd=0.d0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) Collapse(2)
!$OMP& schedule(static)
      do  j = 1,dimcut
         do i = 1,dimcut
            EETA1_kd(i,j) = SUM(GG_kd(:,i)*HH_kd(:,j))
         end do
      end do
!$OMP end PARALLEL DO
! other part of commutator is the transpose
      EETA2_kd = TRANSPOSE(EETA1_kd)
      EETA_kd = EETA1_kd - EETA2_kd
      DDYDS1_kd=0.d0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) Collapse(2)
!$OMP& schedule(static)
      do  j = 1,dimcut
         do i = 1,dimcut
            DDYDS1_kd(i,j) = -SUM(EETA_kd(:,i)*HH_kd(:,j))
         end do
      end do
!$OMP end PARALLEL DO
      DDYDS2_kd = TRANSPOSE(DDYDS1_kd)
      DDYDS_kd = DDYDS1_kd + DDYDS2_kd      ! because \eta = -\eta^\dagger --> H*eta = - (eta*H)^\dagger
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) 
      do j = 1,dimcut
         do i = 1,j  !dimcut
            k = i+j*(j-1)/2
            YDOT(k) = DDYDS_kd(i,j)
         end do
      end do
!$OMP end PARALLEL DO

      end subroutine cSRG_RHS_kd
!*******************************
      subroutine SRG_JAC

      end subroutine SRG_JAC

      end module SRG_module
