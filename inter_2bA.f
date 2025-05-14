!*******************************
      module inter_2bA
         use H_matrices
      contains
!**********************
      subroutine two_body_matrix(icutnum_in,v2b_method)
      use paramdef
      use harmosc
      use H_matrices
      use fitting_2bA
!      use CS_matrix_2bA
      use nninteraction
      use v2bA_quadrature
      implicit none
      integer:: icutnum_in,v2b_method
      integer::ia,ib
      print*,"iccdim=",iccdim,"iucdim",iucdim
      if (allocated(v_cc)) deallocate(v_cc)
      allocate(v_cc(0:iccdim,0:iccdim,jrelm))
      allocate(vcc(0:iccdim,0:iccdim,jrelm))
      if (allocated(v_uc)) deallocate(v_uc)
      allocate(v_uc(0:nrelm,0:nrelm,iucdim))
      allocate(vuc(0:nrelm,0:nrelm,iucdim))
      v_cc=0
      v_uc=0
      call mkl_set_num_threads(240)
      if(icutnum_in==100.or. (icutnum_in>=105 
     &  .and.icutnum_in<=110)) then
       !********
         if(v2b_method==2) then
      !call inter_2bA_fit(icutnum_in)
         call inter_2bA_fit_quadrature(
     &                           icutnum_in)
         elseif (v2b_method==3) then
!            call inter_CS_multistep(icutnum_in)
             print*,"This routine is not integrated in the code anymore"
             stop
         elseif (v2b_method==4) then
            call inter_2b_quadrature_qp_c(icutnum_in)
         endif
       !********
         else
         call inter_c_old
      endif
      if(MTV.neqv..true.) then! inter_c adds coul contribution on its own
         call add_coulomb()
      endif
      vcc(:,:,:)=v_cc(:,:,:)
      vuc(:,:,:)=v_uc(:,:,:)
      end subroutine two_body_matrix
      

      subroutine v2b_eff
      use H_matrices
      use paramdef
      use ist
      use nninteraction
      use effinteraction
      use harmosc
      use SRG_module, only: Ngen,srg
     & ,NN_SRG_for3b_extract 
      use pointers, only: jtot2,itot2
      implicit none
      logical uncpl
      character(len=3):: string1,string2
      character(len=20):: fname
      integer:: i1,i2,ia1,ic2,ii,is1,is2,it1,it2,iuc,j1,j2,l1,l2,n1
     & ,icceff,n2,nnl,nreffmax
      integer:: iucmax
      nreffmax=mxeng/2
      icceff=2*nreffmax+1
      nreffmax=mxeng/2
      icceff=2*nreffmax+1
      iucmax=2*max(jrelm,mxeng/2)+2
      if (allocated(inter_cc)) deallocate(inter_cc)
      if (allocated(inter_uc)) deallocate(inter_uc)
      if(multi) then
         allocate(inter_uc(0:nreffmax,0:nreffmax,iucmax,0:mxeng/2))
         allocate(inter_cc(0:icceff,0:icceff,jrelm,0:mxeng/2))
      else
         allocate(inter_uc(0:nreffmax,0:nreffmax,iucmax,0:0))
         allocate(inter_cc(0:icceff,0:icceff,jrelm,0:0))
      endif
      inter_cc=0.d0
      inter_uc=0.d0
      if (nucleons==3.and.threeff) then
!c         hbo3=hbo/dble(nucl3eff) ! for SRG
         hbo3=hbo/dble(nucleons)  ! for SRG
      else
         hbo3=hbo/dble(nucleons)
      endif
      
      !call eff_cc_inter_2b_test()
      call eff_cc_inter_2b()
      print *,' effcoupl called'
      call eff_uc_inter_2b()
      print *,' effuncpl called'
      if (allocated(v2bc)) deallocate(v2bc)
      if (allocated(t2bc)) deallocate(t2bc)
      if (allocated(Pv2bPc)) deallocate(Pv2bPc)
      if (multi) then
         allocate(v2bc(ist2dim*(ist2dim+1)/2,0:mxeng/2))
      else
         allocate(v2bc(ist2dim*(ist2dim+1)/2,0:0))
         allocate(t2bc(ist2dim*(ist2dim+1)/2))
         allocate(Pv2bPc(ist2dim*(ist2dim+1)/2))
      endif   
      v2bc=cmplx(0.d0,0.d0)
      t2bc=cmplx(0.d0,0.d0)
      Pv2bPc=cmplx(0.d0,0.d0)
      do i1=1,ist2dim
         n1=ist2(1,i1)
         l1=ist2(2,i1)
         is1=ist2(3,i1)
         j1=ist2(4,i1)
         it1=ist2(5,i1)      
         nnl=2*n1+l1
         if (j1==l1.or.j1==0) then
            uncpl=.true.
         else
            uncpl=.false.
         endif
         do i2=i1,ist2dim
            n2=ist2(1,i2)
            l2=ist2(2,i2)
            is2=ist2(3,i2)
            j2=ist2(4,i2)
            it2=ist2(5,i2)  
            if (is2/=is1) cycle
            if (it2/=it1) cycle
            if (j2/=j1) cycle
            ii=i1+i2*(i2-1)/2
            if (j1<=jrelm) then
               if (uncpl) then
                  iuc=2*j1+1+is1 
                  v2bc(ii,:)=inter_uc(n1,n2,iuc,:)
                  if (2*n2+l2<=Ngen.and.2*n1+l1<=Ngen) Pv2bPc(ii)=
     &                 v_uc(n1,n2,iuc)
               else
                  if (l1==j1+1) then 
                     ia1=2*n1+1
                  else
                     ia1=2*n1
                  endif
                  if (l2==j1+1) then 
                     ic2=2*n2+1
                  else
                     ic2=2*n2
                  endif
                  v2bc(ii,:)=inter_cc(ia1,ic2,j1,:)
                  if (2*n2+l2<=Ngen.and.2*n1+l1<=Ngen) Pv2bPc(ii)=
     &                 v_cc(ia1,ic2,j1)
               endif
            else
               if (l1==l2) then
                  if (n1==n2) then
                     if (noho) then
                        v2bc(ii,:)=(nnl+1.5d0)*hbo3
     & *complex(cos(2.d0*theta),sin(-2.d0*theta))
                     else
                        v2bc(ii,:)=-(nnl+1.5d0)*hbo3
     & *complex(cos(2.d0*theta),sin(-2.d0*theta))
                     endif   
                  elseif (n1==n2+1) then
                     v2bc(ii,:)=
     &                    sqrt((n2+1.d0)*(n2+l2+1.5d0))*hbo3
     & *complex(cos(2.d0*theta),sin(-2.d0*theta))
                  elseif (n1==n2-1) then
                     v2bc(ii,:)=
     &                    sqrt((n1+1.d0)*(n1+l1+1.5d0))*hbo3
     & *complex(cos(2.d0*theta),sin(-2.d0*theta))
                  endif
               endif   
            endif
            if (l1==l2) then
               if (n1==n2) then
                  t2bc(ii)=(nnl+1.5d0)*hbo3
     & *complex(cos(2.d0*theta),sin(-2.d0*theta))
               elseif (n1==n2+1) then
                  t2bc(ii)=
     &                 sqrt((n2+1.d0)*(n2+l2+1.5d0))*hbo3
     & *complex(cos(2.d0*theta),sin(-2.d0*theta))
               elseif (n1==n2-1) then
                  t2bc(ii)=
     &                 sqrt((n1+1.d0)*(n1+l1+1.5d0))*hbo3
     & *complex(cos(2.d0*theta),sin(-2.d0*theta))
               endif
            endif     
         end do   
      end do
      deallocate(inter_uc,inter_cc)
      deallocate(v_cc,vcc)
      deallocate(v_uc,vuc)
!*****************************
!      open(225,file="v2b.test")
!      write(225,"(1000ES14.4)"),(v2bc(ii,0),ii=1
!     &                   ,ist2dim+ist2dim*(ist2dim-1)/2)
!      close(225)
!*****************************
      if(NN_SRG_for3b_extract) then
         fname='v2b_j'
         write(string1,'(i3)'),jtot2
         write(string2,'(i2)'),itot2
         string1=adjustl(string1)
         string2=adjustl(string2)
         open(14156,file=trim(fname)//trim(string1)//'_T'//
     &    trim(string2)//'.out',status='unknown')
         write(14156,*),'jtot2=',jtot2,'itot2=',itot2
         write(14156,*),'nrelm=',nrelm,'srg=',srg
         write(14156,*),'leesuzuki=',leesuz
         write(14156,*),'i1 ,  i2  ,  v2b(i1+i2(i2-1)/2,0:6)'
         do i1=1,ist2dim
            do i2=1,ist2dim
               ii=i1+i2*(i2-1)/2
               write(14156,'(2i5,2E12.5)'),i1,i2,v2bc(ii,0)
            enddo
         enddo
         close(14156)
      endif
      end subroutine
!**********************
      subroutine eff_cc_inter_2b()
      use paramdef
      use nninteraction
      use harmosc
      use SRG_module!, only: srg,Hgen,Ngen,dimgen,TkinSave,
!     &    Ndiag,delta_E,ena,enb,bswitch,HHgen,TTkinSave
      use HOtrap, only: HO_trap,V_HO
      use H_matrices
      implicit none
      !**********
      integer,parameter:: kd=real128
      integer:: lx1,lx2,jrel,jfmin,jfmax,jspin,ispin,
     & NE_max,num1_states,num2_states,ipmax,ipot
      real(kd):: theta_kd
      !lx1: orbital momentum of bra
      !lx2: orbital momentum of ket
      !jrel: total angular momentum of bra and ket
      !jspin: spin (S)= 1 (spin triplet) or 0 (spin singlet)
      !ispin: isospin (T)= 1 or 0
      !num1_states: the number of states with L=J-1 and E<Emax
      !num2_states: the number of states with L=J+1 and E<Emax
      !NE_max: the highest allowed state energetically 
      !ipmax: the number of states with L=J+-1 and E<Emax
      !**********
      complex(kd),allocatable:: Heff(:,:),H_cc(:,:),eigr(:),eigv(:,:)
     & ,Tkin(:,:),TkinA(:,:),G_cc(:,:)
      !states space in Heff: 
      !  (n,l): (0,lx1) (0,lx2) (1,lx) (1,lx2) ...
      complex(kd):: rrel
      integer:: ia,ib,nra,nrb,lra,lrb
      integer:: max_ncut,ncut,ipmax_eff,Nmx_eff
      integer:: time3,time4,ir
      logical:: convergence
      !temporary untile writing csrg_kd
      complex(kind(0.d0)),allocatable:: Heff_db(:,:),G_db(:,:)
      integer:: rep
      !
!      select case(ipot)
!      case(2)
!         open(44,file='V_pn.mat',status='unknown')
!      case(1)
!         open(44,file='V_pp.mat',status='unknown')
!      case(3)
!         open(44,file='V_nn.mat',status='unknown')
!      case default
!         open(44,file='V_NN.mat',status='unknown')
!      end select
      print*,"************************"
      print*,"eff_cc_inter_2b started"
      jfmin=1
      jfmax=min(jrelm,mxeng+1)
      theta_kd=real(theta,kd)
      print*,"theta=",theta_kd
      convergence=multi
      open(223,file="inter_2b_fit.out",Access='append')
      do jrel=jfmin,jfmax
         lx1=jrel-1
         lx2=jrel+1
         ! in coupled channels S=1 :
         jspin=1
         ispin=0
         !maintain antisymmetry
         if( (-1)**(lx1+jspin+ispin) /= (-1)) ispin=1
         if ((ipot==1.or.ipot==3).and.ispin==0) cycle
         !calculating the dimensions of the matrix starting from the energy cutoff: NE_max
         ! and NE=2n+l
         if (srg.and.NN_SRG_for3b_extract) then
            NE_max=mxeng
         else
            NE_max=2*nrelm
         endif
         num1_states=(NE_max-lx1)/2+1
         if (lx1>NE_max) num1_states=0
         num2_states=(NE_max-lx2)/2+1
         if (lx2>NE_max) num2_states=0
         ipmax=num1_states+num2_states
         !
         num1_states=(Ngen-lx1)/2+1
         if (lx1>Ngen) num1_states=0
         num2_states=(Ngen-lx2)/2+1
         if (lx2>Ngen) num2_states=0
         dimgen=num1_states+num2_states
         !
         print *,'*****'
         print *,' Ngen, jpw, dimgen=',Ngen,jrel,dimgen
         print *,' Ndiag, delta_E, ipmiax=',Ndiag,delta_E,ipmax
         if(allocated(H_cc)) deallocate(H_cc,Tkin,TkinA,G_cc)
        

         allocate(H_cc(ipmax,ipmax),Tkin(ipmax,ipmax)
     &          ,G_cc(ipmax,ipmax),TkinA(ipmax,ipmax))
         G_cc=0
         write(223,*),"************************"
         write(223,*),' S=',jspin,' j=',jrel,' T=',ispin
         write(223,*),' Dimension=', ipmax
         !**********
         !filling Heff and Tkin
         do ia=1,ipmax
            nra=(ia-1)/2
            if(mod(ia,2)==1) then
               lra=lx1
            else 
               lra=lx2
            endif
            do ib=1,ipmax
               nrb=(ib-1)/2
               if(mod(ib,2)==1) then
                  lrb=lx1
               else 
                  lrb=lx2
               endif
               !***************
               !calculate the kinetic matrix elements
               if (lra==lrb) then
                  if (nra==nrb) then
                     rrel=(2*nra+lra+1.5_kd)
                  elseif (nrb==nra+1) then
                     rrel=sqrt((nra+1_kd)*(nra+lra+1.5d0))
                  elseif (nra==nrb+1) then
                     rrel=sqrt((nrb+1_kd)*(nrb+lrb+1.5d0))
                  else
                     rrel=0._kd
                  endif
               else
                  rrel=complex(0._kd,0._kd)
               endif
               rrel=rrel
               !***************
               Tkin(ia,ib)= rrel*hbo/2!**2*0.1
               TkinA(ia,ib)= rrel*hbo3!**2*0.1
               if(srg)then

                  H_cc(ia,ib)=v_cc(ia-1,ib-1,jrel)+Tkin(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
               else
                  H_cc(ia,ib)=v_cc(ia-1,ib-1,jrel)+TkinA(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
               endif
               if(ia==ib.or..true.) G_cc(ia,ib)=Tkin(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
               H_cc(ib,ia)=H_cc(ia,ib)
            enddo
         enddo
!         if(jrel==2) then
!         print*,ipmax
!         open(226,file="Htest2.test")
!         do ia=1,min(140,ipmax)
!            write(226,"(500ES14.4)"),(H_cc(ia,ib),ib=1
!     &             ,min(140,ipmax))
!         enddo
!         close(226)
           !TEST 
!            open(224,file="2b_spect.dat",Access='append')
!            write(224,*),"******************************"
!            write(224,*),"Two-body hamiltonian exact eigenvalues"
!            write(224,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
!     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
!     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
!     & -2*theta))),ia=1,min(ipmax,40))
!            close(224)
!            stop
!         endif
         !**********
         ! print the exact eigen values before SRG
            allocate(eigr(ipmax),eigv(ipmax,ipmax))
            print*,"before cgeev"

         call system_clock(count=time3, count_rate=ir)
            call mkl_Cgeev(H_cc,ipmax,eigr,eigv)
         call system_clock(count=time4, count_rate=ir)
         print *,'cgeev finished',
     &   real(time4 - time3,kind=8)/real(ir,kind=8)
         write(223,*),"Two-body hamiltonian exact eigenvalues"
         write(223,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax,80))
         !**********
         deallocate(eigr,eigv)
         !**********
         !**********
         !performing SRG
         if(srg) then
            print *,'srg evolution is started'
            if(.false.) then
               call system_clock(count=time3, count_rate=ir)
                  call csrg_evolution_kd(ipmax,H_cc,Tkin)
               call system_clock(count=time4, count_rate=ir)
               print *,'srg evolution is finished',
     &         real(time4 - time3,kind=8)/real(ir,kind=8)
            else
               allocate(Heff_db(ipmax,ipmax),G_db(ipmax,ipmax))
               Heff_db(:,:)=H_cc(:,:)
     &*complex(cos(2._kd*theta_kd),+sin(2._kd*theta_kd))
               G_db(:,:)=G_cc(:,:)
     &*complex(cos(2._kd*theta_kd),+sin(2._kd*theta_kd))
               do ia=1,ipmax
                  do ib=ia+1,ipmax
               if(abs(Heff_db(ia,ib)-Heff_db(ib,ia))>10.d0**-13)then
                     print*,'non-sym',ia,ib,abs(Heff_db(ia,ib)
     &  -Heff_db(ib,ia))
                   endif
                  enddo
               enddo
               call system_clock(count=time3, count_rate=ir)
                  call csrg_evolution(ipmax,Heff_db,G_db)
               call system_clock(count=time4, count_rate=ir)
               print *,'srg finished',
     &         real(time4 - time3,kind=8)/real(ir,kind=8)
               H_cc(:,:)=Heff_db(:,:)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
               deallocate(Heff_db,G_db)
               !********************
            endif
         endif
         ! print the exact eigen values before SRG
         allocate(eigr(ipmax),eigv(ipmax,ipmax))
         call system_clock(count=time3, count_rate=ir)
            call mkl_Cgeev(H_cc,ipmax,eigr,eigv)
         call system_clock(count=time4, count_rate=ir)
         print *,'cgeev finished',
     &      real(time4 - time3,kind=8)/real(ir,kind=8)
         write(223,*),"Two-body hamiltonian exact eigenvalues
     &                                                  after srg"
         write(223,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax,80))
         deallocate(eigr,eigv)
         !**********
         !print the spectrum and the convergence trend if requested
         write(223,*),"model-space diagonalization"
         max_ncut=0
         if(convergence) max_ncut=mxeng-2
         do Nmx_eff=mxeng-max_ncut,mxeng,2
         !do Nmx_eff=2,60,2
            ncut=mxeng-Nmx_eff
            !ncut=NE_max-Nmx_eff
            !**********
            !calculating the dimension of the effective Hamiltonian
            num1_states=(Nmx_eff-lx1)/2+1
            if (lx1>Nmx_eff) num1_states=0
            num2_states=(Nmx_eff-lx2)/2+1
            if (lx2>Nmx_eff) num2_states=0
            ipmax_eff=num1_states+num2_states

            allocate( Heff(ipmax_eff,ipmax_eff))
            Heff(1:ipmax_eff,1:ipmax_eff)=H_cc(1:ipmax_eff
     &                                       ,1:ipmax_eff)
!            Heff(1:ipmax_eff,1:ipmax_eff)=Tkin(1:ipmax_eff
!     &                                       ,1:ipmax_eff)
!            Heff(1:min(ipmax_eff,14),1:min(ipmax_eff,14))=
!     & H_cc(1:min(ipmax_eff,14),1:min(ipmax_eff,14))
            allocate(eigr(ipmax_eff),eigv(ipmax_eff,ipmax_eff))
            call mkl_Cgeev(Heff,ipmax_eff,eigr,eigv)
            write(223,*),"Two-body effective hamiltonian eigenvalues
     &                        for N_2=",Nmx_eff,"dimension",ipmax_eff
            write(223,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax_eff,40))
            deallocate(eigr,eigv)

            do ia=1,ipmax_eff
               do ib=ia,ipmax_eff
                  if(srg)then

                     inter_cc(ia-1,ib-1,jrel,ncut/2)=H_cc(ia,ib)
     &                              +TkinA(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
     &-Tkin(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
                  else
                     inter_cc(ia-1,ib-1,jrel,ncut/2)=H_cc(ia,ib)
                  endif
                  inter_cc(ib-1,ia-1,jrel,ncut/2)=inter_cc(ia-1
     &                                          ,ib-1,jrel,ncut/2)
               enddo
            enddo


            deallocate(Heff)
         enddo
         !**********
      enddo
      close(223)
      end subroutine eff_cc_inter_2b
!*******************************
!**********************
      subroutine eff_uc_inter_2b()
      use paramdef
      use nninteraction
      use harmosc
      use SRG_module!, only: srg,Hgen,Ngen,dimgen,TkinSave,
!     &    Ndiag,delta_E,ena,enb,bswitch,HHgen,TTkinSave
      use HOtrap, only: HO_trap,V_HO
      use H_matrices
      implicit none
      !**********
      integer,parameter:: kd=real128
      integer::jrel,jfmin,jfmax,jspin,ispin,
     & NE_max,ipmax,ipot
      real(kd):: theta_kd
      !lx1: orbital momentum of bra
      !lx2: orbital momentum of ket
      !jrel: total angular momentum of bra and ket
      !jspin: spin (S)= 1 (spin triplet) or 0 (spin singlet)
      !ispin: isospin (T)= 1 or 0
      !num1_states: the number of states with L=J-1 and E<Emax
      !num2_states: the number of states with L=J+1 and E<Emax
      !NE_max: the highest allowed state energetically 
      !ipmax: the number of states with L=J+-1 and E<Emax
      !**********
      complex(kd),allocatable:: Heff(:,:),H_uc(:,:),eigr(:),eigv(:,:)
     & ,Tkin(:,:),TkinA(:,:),G_cc(:,:)
      !states space in Heff: 
      !  (n,l): (0,lx1) (0,lx2) (1,lx) (1,lx2) ...
      complex(kd):: rrel
      integer:: ia,ib,nra,nrb,lra,lrb,iuc,lrel
      integer:: max_ncut,ncut,ipmax_eff,Nmx_eff
      integer:: time3,time4,ir
      logical:: convergence
      !temporary untile writing csrg_kd
      complex(kind(0.d0)),allocatable:: Heff_db(:,:),G_db(:,:)
      !
!      select case(ipot)
!      case(2)
!         open(44,file='V_pn.mat',status='unknown')
!      case(1)
!         open(44,file='V_pp.mat',status='unknown')
!      case(3)
!         open(44,file='V_nn.mat',status='unknown')
!      case default
!         open(44,file='V_NN.mat',status='unknown')
!      end select
      print*,"************************"
      print*,"eff_uc_inter_2b started"
      jfmin=0
      jfmax=min(jrelm,mxeng+1)
      theta_kd=real(theta,kd)
      convergence=multi
      open(223,file="inter_2b_fit.out",Access='append')
      do jrel=jfmin,jfmax
      do jspin=0,1
         ! in uncoupled channels J=l1=l2 :
         lrel=jrel 
         if(jrel==0 .and. jspin==1) lrel=1
         ispin=0
         !maintain antisymmetry
         if( (-1)**(lrel+jspin+ispin) /= (-1)) ispin=1
         if ((ipot==1.or.ipot==3).and.ispin==0) cycle
         !calculating the dimensions of the matrix starting from the energy cutoff: NE_max
         ! and NE=2n+l
         if (srg.and.NN_SRG_for3b_extract) then
            NE_max=mxeng
         else
            NE_max=2*nrelm
         endif
         if(lrel>2*nrelm) then
            ipmax=0
         else
            ipmax=(NE_max-lrel)/2+1
         endif
         !
         if(lrel>Ngen) then
            dimgen=0
         else
            dimgen=(Ngen-lrel)/2+1
         endif
         !

         print *,'*****'
         print *,' Ngen, jpw, dimgen=',Ngen,jrel,dimgen
         print *,' Ndiag, delta_E, ipmiax=',Ndiag,delta_E,ipmax
         if(allocated(H_uc)) deallocate(H_uc,Tkin,TkinA,G_cc)
         
         allocate(H_uc(ipmax,ipmax),Tkin(ipmax,ipmax)
     &         ,G_cc(ipmax,ipmax),TkinA(ipmax,ipmax))
         G_cc=0 
         write(223,*),"************************"
         write(223,*),' S=',jspin,' j=',jrel,' T=',ispin,' l=',lrel
         write(223,*),' Dimension=', ipmax
         !**********
         !filling Heff and Tkin
         iuc=2*jrel+jspin+1
         do ia=1,ipmax
            nra=ia-1
            do ib=1,ipmax
               nrb=ib-1
               !***************
               !calculate the kinetic matrix elements
               if (nra==nrb) then
                  rrel=(2*nra+lrel+1.5_kd)
               elseif (nrb==nra+1) then
                  rrel=sqrt((nra+1_kd)*(nra+lrel+1.5_kd))
               elseif (nra==nrb+1) then
                  rrel=sqrt((nrb+1_kd)*(nrb+lrel+1.5_kd))
               else
                  rrel=0._kd
               endif
               rrel=rrel
               !***************
               Tkin(ia,ib)=rrel*hbo/2
               TkinA(ia,ib)=rrel*hbo3
               if(srg) then
               H_uc(ia,ib)=v_uc(nra,nrb,iuc)+Tkin(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
               else
               H_uc(ia,ib)=v_uc(nra,nrb,iuc)+TkinA(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
               endif
               if(ia==ib.or..true.) G_cc(ia,ib)= Tkin(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
            enddo
         enddo
         !**********
         ! print the exact eigen values before SRG
         allocate(eigr(ipmax),eigv(ipmax,ipmax))
         call mkl_Cgeev(H_uc,ipmax,eigr,eigv)
         write(223,*),"Two-body hamiltonian exact eigenvalues"
         write(223,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax,80))
         deallocate(eigr,eigv)
         !**********
         !performing SRG
         allocate(Heff_db(ipmax,ipmax),G_db(ipmax,ipmax))
         if(srg) then
            Heff_db(:,:)=H_uc(:,:)
     &*complex(cos(2._kd*theta_kd),+sin(2._kd*theta_kd))
            G_db(:,:)=G_cc(:,:)
     &*complex(cos(2._kd*theta_kd),+sin(2._kd*theta_kd))
            !call csrg_evolution(ipmax,H_uc,Tkin)
         call system_clock(count=time3, count_rate=ir)
            call csrg_evolution(ipmax,Heff_db,G_db)
         call system_clock(count=time4, count_rate=ir)
         print *,'srg evolution is finished',
     &      real(time4 - time3,kind=8)/real(ir,kind=8)
            H_uc(:,:)=Heff_db(:,:)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
         endif
         deallocate(Heff_db,G_db)
         allocate(eigr(ipmax),eigv(ipmax,ipmax))
         call system_clock(count=time3, count_rate=ir)
            call mkl_Cgeev(H_uc,ipmax,eigr,eigv)
         call system_clock(count=time4, count_rate=ir)
         print *,'cgeev finished',
     &      real(time4 - time3,kind=8)/real(ir,kind=8)
         write(223,*),"Two-body hamiltonian exact eigenvalues
     &                                                  after srg"
         write(223,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax,80))
         deallocate(eigr,eigv)
         !**********
         !print the spectrum and the convergence trend if requested
         write(223,*),"model-space diagonalization"
         max_ncut=0
         if(convergence) max_ncut=mxeng-2
         do Nmx_eff=mxeng-max_ncut,mxeng,2
            ncut=mxeng-Nmx_eff
            !**********
            !calculating the dimension of the effective Hamiltonian
            ipmax_eff=(Nmx_eff-lrel)/2+1
            if (lrel>Nmx_eff) ipmax_eff=0
            print*,"ipmax_eff",ipmax_eff

            allocate( Heff(ipmax_eff,ipmax_eff))
            Heff(1:ipmax_eff,1:ipmax_eff)=H_uc(1:ipmax_eff
     &                                       ,1:ipmax_eff)
            allocate(eigr(ipmax_eff),eigv(ipmax_eff,ipmax_eff))
            call mkl_Cgeev(Heff,ipmax_eff,eigr,eigv)
            write(223,*),"Two-body effective hamiltonian eigenvalues
     &                        for N_2=",Nmx_eff,"dimension",ipmax_eff
            write(223,*),"Model-space diagonalization"
            write(223,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax_eff,40))
            deallocate(eigr,eigv)

            do ia=1,ipmax_eff
               do ib=ia,ipmax_eff
                  if(srg) then
                     inter_uc(ia-1,ib-1,iuc,ncut/2)=H_uc(ia,ib)
     &                                 +TkinA(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
     & -Tkin(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
                  else
                     inter_uc(ia-1,ib-1,iuc,ncut/2)=H_uc(ia,ib)
                  endif
                  inter_uc(ib-1,ia-1,iuc,ncut/2)=inter_uc(ia-1
     &                                          ,ib-1,iuc,ncut/2)
               enddo
            enddo


            deallocate(Heff)
         enddo
         !**********
      enddo
      enddo
      close(223)
      end subroutine eff_uc_inter_2b
!*******************************
      subroutine eff_cc_inter_2b_test()
      use paramdef
      use nninteraction
      use harmosc
      use SRG_module!, only: srg,Hgen,Ngen,dimgen,TkinSave,
!     &    Ndiag,delta_E,ena,enb,bswitch,HHgen,TTkinSave
      use HOtrap, only: HO_trap,V_HO
      use H_matrices
      implicit none
      !**********
      integer,parameter:: kd=real128
      integer:: lx1,lx2,jrel,jfmin,jfmax,jspin,ispin,
     & NE_max,num1_states,num2_states,ipmax,ipot
      real(kd):: theta_kd
      !lx1: orbital momentum of bra
      !lx2: orbital momentum of ket
      !jrel: total angular momentum of bra and ket
      !jspin: spin (S)= 1 (spin triplet) or 0 (spin singlet)
      !ispin: isospin (T)= 1 or 0
      !num1_states: the number of states with L=J-1 and E<Emax
      !num2_states: the number of states with L=J+1 and E<Emax
      !NE_max: the highest allowed state energetically 
      !ipmax: the number of states with L=J+-1 and E<Emax
      !**********
      complex(kd),allocatable:: Heff(:,:),H_cc(:,:),eigr(:),eigv(:,:)
     & ,Tkin(:,:),TkinA(:,:),G_cc(:,:),U(:,:)
      !states space in Heff: 
      !  (n,l): (0,lx1) (0,lx2) (1,lx) (1,lx2) ...
      complex(kd):: rrel
      integer:: ia,ib,nra,nrb,lra,lrb,i
      integer:: max_ncut,ncut,ipmax_eff,Nmx_eff
      integer:: time3,time4,ir
      logical:: convergence
      !temporary untile writing csrg_kd
      complex(kind(0.d0)),allocatable:: Heff_db(:,:),G_db(:,:)
      integer:: rep
      !
!      select case(ipot)
!      case(2)
!         open(44,file='V_pn.mat',status='unknown')
!      case(1)
!         open(44,file='V_pp.mat',status='unknown')
!      case(3)
!         open(44,file='V_nn.mat',status='unknown')
!      case default
!         open(44,file='V_NN.mat',status='unknown')
!      end select
      print*,"************************"
      print*,"eff_cc_inter_2b started"
      jfmin=1
      jfmax=min(jrelm,mxeng+1)
      theta_kd=real(theta,kd)
      print*,"theta=",theta_kd
      convergence=multi
      open(223,file="inter_2b_fit.out",Access='append')
      do jrel=jfmin,jfmax
         lx1=jrel-1
         lx2=jrel+1
         ! in coupled channels S=1 :
         jspin=1
         ispin=0
         !maintain antisymmetry
         if( (-1)**(lx1+jspin+ispin) /= (-1)) ispin=1
         if ((ipot==1.or.ipot==3).and.ispin==0) cycle
         !calculating the dimensions of the matrix starting from the energy cutoff: NE_max
         ! and NE=2n+l
         if (srg.and.NN_SRG_for3b_extract) then
            NE_max=mxeng
         else
            NE_max=2*nrelm
         endif
         num1_states=(NE_max-lx1)/2+1
         if (lx1>NE_max) num1_states=0
         num2_states=(NE_max-lx2)/2+1
         if (lx2>NE_max) num2_states=0
         ipmax=num1_states+num2_states
         !
         num1_states=(Ngen-lx1)/2+1
         if (lx1>Ngen) num1_states=0
         num2_states=(Ngen-lx2)/2+1
         if (lx2>Ngen) num2_states=0
         dimgen=num1_states+num2_states
         !
         print *,'*****'
         print *,' Ngen, jpw, dimgen=',Ngen,jrel,dimgen
         print *,' Ndiag, delta_E, ipmiax=',Ndiag,delta_E,ipmax
         if(allocated(H_cc)) deallocate(H_cc,Tkin,TkinA,G_cc)
        

         allocate(H_cc(ipmax,ipmax),Tkin(ipmax,ipmax)
     &          ,G_cc(ipmax,ipmax),TkinA(ipmax,ipmax))
         G_cc=0
         write(223,*),"************************"
         write(223,*),' S=',jspin,' j=',jrel,' T=',ispin
         write(223,*),' Dimension=', ipmax
         !**********
         !filling Heff and Tkin
         do ia=1,ipmax
            nra=(ia-1)/2
            if(mod(ia,2)==1) then
               lra=lx1
            else 
               lra=lx2
            endif
            do ib=1,ipmax
               nrb=(ib-1)/2
               if(mod(ib,2)==1) then
                  lrb=lx1
               else 
                  lrb=lx2
               endif
               !***************
               !calculate the kinetic matrix elements
               if (lra==lrb) then
                  if (nra==nrb) then
                     rrel=(2*nra+lra+1.5_kd)
                  elseif (nrb==nra+1) then
                     rrel=sqrt((nra+1_kd)*(nra+lra+1.5d0))
                  elseif (nra==nrb+1) then
                     rrel=sqrt((nrb+1_kd)*(nrb+lrb+1.5d0))
                  else
                     rrel=0._kd
                  endif
               else
                  rrel=complex(0._kd,0._kd)
               endif
               rrel=rrel
               !***************
               Tkin(ia,ib)= rrel*hbo/2!**2*0.1
               TkinA(ia,ib)= rrel*hbo3!**2*0.1
               if(srg)then

                  H_cc(ia,ib)=v_cc(ia-1,ib-1,jrel)+Tkin(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
               else
                  H_cc(ia,ib)=v_cc(ia-1,ib-1,jrel)+TkinA(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
               endif
               if(ia==ib.or..true.) G_cc(ia,ib)=Tkin(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
               H_cc(ib,ia)=H_cc(ia,ib)
            enddo
         enddo
         if(jrel==2) then
         print*,ipmax
         open(226,file="Htest2.test")
         do ia=1,min(140,ipmax)
            write(226,"(500ES14.4)"),(H_cc(ia,ib),ib=1
     &             ,min(140,ipmax))
         enddo
         close(226)
           !TEST 
!            open(224,file="2b_spect.dat",Access='append')
!            write(224,*),"******************************"
!            write(224,*),"Two-body hamiltonian exact eigenvalues"
!            write(224,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
!     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
!     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
!     & -2*theta))),ia=1,min(ipmax,40))
!            close(224)
!            stop
         endif
         !**********
         ! print the exact eigen values before SRG
            allocate(eigr(ipmax),eigv(ipmax,ipmax))
            print*,"before cgeev"

         call system_clock(count=time3, count_rate=ir)
            call mkl_Cgeev(H_cc,ipmax,eigr,eigv)
         call system_clock(count=time4, count_rate=ir)
         print *,'cgeev finished',
     &   real(time4 - time3,kind=8)/real(ir,kind=8)
         write(223,*),"Two-body hamiltonian exact eigenvalues"
         write(223,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax,80))
         !transform to momentum basis
         deallocate(eigr,eigv)
            allocate(eigr(ipmax),eigv(ipmax,ipmax))
         call system_clock(count=time3, count_rate=ir)
            call mkl_Cgeev(Tkin,ipmax,eigr,eigv)
         call system_clock(count=time4, count_rate=ir)
         print *,'cgeev finished'
         allocate (U(ipmax,ipmax))
         !transform H
         do ia=1,ipmax
            do ib=1,ipmax
               U(ia,ib)=0
               do i=1,ipmax
                  U(ia,ib)=U(ia,ib)+H_cc(ia,i)*eigv(i,ib)
               enddo
            enddo
         enddo
         do ia=1,ipmax
            do ib=1,ipmax
               !H_cc(ia,ib)=0
               do i=1,ipmax
                  !H_cc(ia,ib)=H_cc(ia,ib)+eigv(i,ia)*U(i,ib)
               enddo
            enddo
         enddo
         !transform T
         do ia=1,ipmax
            do ib=1,ipmax
               U(ia,ib)=0
               do i=1,ipmax
                  U(ia,ib)=U(ia,ib)+Tkin(ia,i)*eigv(i,ib)
               enddo
            enddo
         enddo
         do ia=1,ipmax
            do ib=1,ipmax
               !Tkin(ia,ib)=0
               do i=1,ipmax
                  !Tkin(ia,ib)=Tkin(ia,ib)+eigv(i,ia)*U(i,ib)
               enddo
            enddo
         enddo
         do ia=1,ipmax
            do ib=1,ipmax
               U(ia,ib)=0
               do i=1,ipmax
                  U(ia,ib)=U(ia,ib)+eigv(ia,i)*eigv(ib,i)
               enddo
            enddo
         enddo
         write(223,*),"T eigen values"
         write(223,"(100ES14.4)"),(eigr(ib),ib=1,50)
         call system_clock(count=time3, count_rate=ir)
            call mkl_Cgeev(Tkin,ipmax,eigr,eigv)
         call system_clock(count=time4, count_rate=ir)
         print *,'cgeev finished'
         write(223,*),"T eigen values"
         write(223,"(100ES14.4)"),(eigr(ib),ib=1,50)

         write(223,*),"H after transformation"
         do ia=1,50
            write(223,"(100ES14.4)"),(H_cc(ia,ib),ib=1,50)
         enddo
         write(223,*),"T after transformation"
         do ia=1,50
            write(223,"(100ES14.4)"),(Tkin(ia,ib),ib=1,50)
         enddo
         write(223,*),"Norm of the transformation"
         do ia=1,50
            write(223,"(100ES14.4)"),(U(ia,ib),ib=1,50)
         enddo
         deallocate(U)
         !!
         !**********
         !deallocate(eigr,eigv)
         !**********
         !**********
         !performing SRG
         if(srg) then
            print *,'srg evolution is started'
            if(.false.) then
               call system_clock(count=time3, count_rate=ir)
                  call csrg_evolution_kd(ipmax,H_cc,Tkin)
               call system_clock(count=time4, count_rate=ir)
               print *,'srg evolution is finished',
     &         real(time4 - time3,kind=8)/real(ir,kind=8)
            else
               allocate(Heff_db(ipmax,ipmax),G_db(ipmax,ipmax))
               Heff_db(:,:)=H_cc(:,:)
     &*complex(cos(2._kd*theta_kd),+sin(2._kd*theta_kd))
               G_db(:,:)=G_cc(:,:)
     &*complex(cos(2._kd*theta_kd),+sin(2._kd*theta_kd))
               G_db(:,:)=Tkin(:,:)
               do ia=1,ipmax
                  do ib=ia+1,ipmax
               if(abs(Heff_db(ia,ib)-Heff_db(ib,ia))>10.d0**-13)then
                     print*,'non-sym',ia,ib,abs(Heff_db(ia,ib)
     &  -Heff_db(ib,ia))
                   endif
                  enddo
               enddo
               call system_clock(count=time3, count_rate=ir)
                  call csrg_evolution(ipmax,Heff_db,G_db)
               call system_clock(count=time4, count_rate=ir)
               print *,'srg finished',
     &         real(time4 - time3,kind=8)/real(ir,kind=8)
               H_cc(:,:)=Heff_db(:,:)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
               deallocate(Heff_db,G_db)
               !********************
         !transform to momentum basis
         allocate (U(ipmax,ipmax))
         !transform H
         do ia=1,ipmax
            do ib=1,ipmax
               U(ia,ib)=0
               do i=1,ipmax
                  U(ia,ib)=U(ia,ib)+H_cc(ia,i)*eigv(ib,i)
               enddo
            enddo
         enddo
         do ia=1,ipmax
            do ib=1,ipmax
               !H_cc(ia,ib)=0
               do i=1,ipmax
                  !H_cc(ia,ib)=H_cc(ia,ib)+eigv(ia,i)*U(i,ib)
               enddo
            enddo
         enddo
         deallocate(U)
         !!
         !**********
         deallocate(eigr,eigv)
            endif
         endif
         ! print the exact eigen values before SRG
         if(nucleons==2.or..true.) then
            allocate(eigr(ipmax),eigv(ipmax,ipmax))
         call system_clock(count=time3, count_rate=ir)
            call mkl_Cgeev(H_cc,ipmax,eigr,eigv)
         call system_clock(count=time4, count_rate=ir)
         print *,'cgeev finished',
     &      real(time4 - time3,kind=8)/real(ir,kind=8)
            write(223,*),"Two-body hamiltonian exact eigenvalues
     &                                                  after srg"
            write(223,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax,40))
            deallocate(eigr,eigv)
         endif
         !**********
         !print the spectrum and the convergence trend if requested
         write(223,*),"model-space diagonalization"
         max_ncut=0
         if(convergence) max_ncut=mxeng-2
         do Nmx_eff=mxeng-max_ncut,mxeng,2
         !do Nmx_eff=2,60,2
            ncut=mxeng-Nmx_eff
            !ncut=NE_max-Nmx_eff
            !**********
            !calculating the dimension of the effective Hamiltonian
            num1_states=(Nmx_eff-lx1)/2+1
            if (lx1>Nmx_eff) num1_states=0
            num2_states=(Nmx_eff-lx2)/2+1
            if (lx2>Nmx_eff) num2_states=0
            ipmax_eff=num1_states+num2_states

            allocate( Heff(ipmax_eff,ipmax_eff))
            Heff(1:ipmax_eff,1:ipmax_eff)=H_cc(1:ipmax_eff
     &                                       ,1:ipmax_eff)
!            Heff(1:ipmax_eff,1:ipmax_eff)=Tkin(1:ipmax_eff
!     &                                       ,1:ipmax_eff)
!            Heff(1:min(ipmax_eff,14),1:min(ipmax_eff,14))=
!     & H_cc(1:min(ipmax_eff,14),1:min(ipmax_eff,14))
            allocate(eigr(ipmax_eff),eigv(ipmax_eff,ipmax_eff))
            call mkl_Cgeev(Heff,ipmax_eff,eigr,eigv)
            write(223,*),"Two-body effective hamiltonian eigenvalues
     &                        for N_2=",Nmx_eff,"dimension",ipmax_eff
            write(223,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax_eff,40))
            deallocate(eigr,eigv)

            do ia=1,ipmax_eff
               do ib=ia,ipmax_eff
                  if(srg)then

                     inter_cc(ia-1,ib-1,jrel,ncut/2)=H_cc(ia,ib)
     &                              +TkinA(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
     &-Tkin(ia,ib)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
                  else
                     inter_cc(ia-1,ib-1,jrel,ncut/2)=H_cc(ia,ib)
                  endif
                  inter_cc(ib-1,ia-1,jrel,ncut/2)=inter_cc(ia-1
     &                                          ,ib-1,jrel,ncut/2)
               enddo
            enddo


            deallocate(Heff)
         enddo
         !**********
      enddo
      close(223)
      end subroutine eff_cc_inter_2b_test

!*******************************
!*******************************
      subroutine add_coulomb()
      use iso_fortran_env
      use paramdef
      use harmosc
      use nninteraction
      use quadratures
      use ist, only: mscheme,threeff
      use pointers, only: itot2
      use HO_wf
      implicit none
      integer,parameter :: kd=real128
      real(kd),allocatable:: roots(:,:,:)
      double precision::cnp,cpp,cnn
      integer:: nra,lra,nrb,jrel,ispin,nrbmax,ipp,iuc,Mtot2
      real(kd):: th_kd
      integer:: ia,ib
      ipp=0
      if (nucleons==3.and.threeff) then
         MTot2=neut3eff-iprot3eff
         ipp=(nucl3eff-Mtot2)*(nucl3eff-Mtot2-2)
         cpp=dble(ipp)
     &        /dble(3*nucl3eff*(nucl3eff-2)+iTot23eff*(iTot23eff+2))
      else
         MTot2=neutrons-iprotons
         ipp=(nucleons-Mtot2)*(nucleons-Mtot2-2)
         cpp=dble(ipp)
     &        /dble(3*nucleons*(nucleons-2)+iTot2*(iTot2+2))
      endif
      if(ipp==0) return
      th_kd=theta
      if (allocated(roots)) deallocate(roots)
      allocate(roots(nrelm,nrelm,0:lrelm))
      do lra=0,lrelm
         call laguerre_rootfinder_kd(nrelm,real(lra+0.5_kd,kd)
     &                         ,roots(1:nrelm,1:nrelm,lra))
      end do
      if(allocated(vcoul)) deallocate(vcoul)
!      !test
!         do nra=1,3
!            write(*,"(100F14.4)"),(roots(lra,nra,0),lra=1,nra)
!            write(*,"(100F14.4)"),(roots(lra,nra,1),lra=1,nra)
!         enddo
!         print*,"**************"
!      !
      call coulomb_c(nrelm,lrelm,roots)
      if (mscheme) vcoul=0.d0
!***********************************
      do nra=0,nrelm
         do lra=0,jrelm+1
            if (2*nra+lra>2*nrelm) cycle
            if (ipp==0) then
               nrbmax=min(nra+1,nrelm)
            else
               nrbmax=nrelm
               ispin=0
               if ((-1)**(lra+1)==1) ispin=1
            endif
            do nrb=nra,nrbmax
               if (2*nrb+lra>2*nrelm) cycle

               if (ipp/=0) then
                  do jrel=iabs(lra-ispin),min(lra+ispin,jrelm)
                     if (jrel==lra.or.jrel==0) then
                        iuc=2*jrel+ispin+1
                        v_uc(nra,nrb,iuc)=
     &                  v_uc(nra,nrb,iuc)+cpp*vcoul(nra,nrb,lra)
     &                                *complex(cos(th_kd),-sin(th_kd))
                       v_uc(nrb,nra,iuc)=v_uc(nra,nrb,iuc)
                     else
                        if (lra==jrel+1) then
                           v_cc(2*nra+1,2*nrb+1,jrel)=
     &         v_cc(2*nra+1,2*nrb+1,jrel)+cpp*vcoul(nra,nrb,lra)
     &                                *complex(cos(th_kd),-sin(th_kd))
                           v_cc(2*nrb+1,2*nra+1,jrel)=
     &                            v_cc(2*nra+1,2*nrb+1,jrel)
                        elseif (lra==jrel-1) then
                           v_cc(2*nra,2*nrb,jrel)=
     &          v_cc(2*nra  ,2*nrb  ,jrel)+cpp*vcoul(nra,nrb,lra)
     &                                 *complex(cos(th_kd),-sin(th_kd))
                   v_cc(2*nrb,2*nra,jrel)=
     &                            v_cc(2*nra,2*nrb,jrel)
                        endif
                     endif
                  end do
               endif

            end do
         end do
      end do

      end subroutine add_coulomb
!*******************************
      subroutine coulomb_c(nrelm,lrelm,roots)
      use constants
      use harmosc
      use nninteraction
      use paramdef,only : theta
      use iso_fortran_env
      use HO_wf
      use quadratures
      implicit double precision (a-h,o-z)
      integer,parameter :: kd=real128
      real(kd) :: roots(nrelm,nrelm,0:lrelm),rr,over
      parameter (nstep=800)
      real(kd),allocatable,dimension(:,:,:):: uc
      real(kd),allocatable :: gaupoi_pp(:),gauwe_pp(:)

      rcoul=hbc/alpha

      if (allocated(gaupoi_pp)) deallocate(gaupoi_pp)
      allocate(gaupoi_pp(nstep))
      if (allocated(gauwe_pp)) deallocate(gauwe_pp)
      allocate(gauwe_pp(nstep))
      call gen_laguerre_quad_qp(nstep,gaupoi_pp,gauwe_pp)
      allocate(uc(0:nstep,0:nrelm,0:lrelm))
      uc=0
      do lr=0,lrelm
         do nr=0,nrelm
            if (2*nr+lr>2*nrelm) cycle
            do i=1,nstep
               rr=gaupoi_pp(i)
               call wf_coulomb(nr,lr,anu_kd,rr,uc(i,nr,lr),
     &                         roots(:max(nr,1),max(nr,1),lr))
!               call wave_func_r(nr,lr,bsquare_kd,rr,uc(i,nr,lr),
!     &                 roots(:max(nr,1),max(nr,1),lr),.false.)
            end do
         end do
      end do
!      !test
!      write(*,"(100F14.4)"),(gaupoi_pp(lr),lr=1,10)
!         print*,"**************"
!      !
!      !test
!         do nr=1,3
!            write(*,"(100F14.4)"),(roots(lr,nr,0),lr=1,nr)
!            write(*,"(100F14.4)"),(roots(lr,nr,1),lr=1,nr)
!         enddo
!         print*,"**************"
!      !
!      !test
!      do lr=0,3
!         do nr=0,3
!            print*,uc(1,nr,lr),uc(10,nr,lr),uc(40,nr,lr),uc(100,nr,lr)
!         enddo
!      enddo
!      !
      uc=uc/sqrt(2._kd)
      uc=uc*sqrt(2._kd) ! change of variable
      print *,' u-coulomb calculated'
      if (allocated(vcoul)) deallocate(vcoul)
      allocate(vcoul(0:nrelm,0:nrelm,0:lrelm))
      vcoul=cmplx(0.d0,0.d0)
      do lr=0,lrelm
         do nra=0,nrelm
            if (2*nra+lr>2*nrelm) cycle
            do nrb=nra,nrelm
               if (2*nrb+lr>2*nrelm) cycle
               over=0._kd
               do i=1,nstep
                  rr=gaupoi_pp(i)*2 ! change of variable
                  if(.not.sqrt(rr/anu_kd)>1.d-6) cycle
                  over=over
     &            +uc(i,nra,lr)*uc(i,nrb,lr)*gauwe_pp(i)*rcoul/rr
               enddo
               vcoul(nra,nrb,lr)=over
               vcoul(nrb,nra,lr)=vcoul(nra,nrb,lr)
            end do
         end do
      end do
      vcoul=vcoul/sqrt(2.d0)

      print *, ' Coulomb calculated'

      deallocate(uc)

      end subroutine coulomb_c


      end module inter_2bA
