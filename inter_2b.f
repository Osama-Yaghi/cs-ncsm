!*******************************
      module inter_2b
         use H_matrices
      contains
!**********************
      subroutine two_body_matrix(ipot)
      use paramdef
      use harmosc
      use H_matrices
      use fitting_2b
      use v2b_quadrature
      implicit none
      integer:: nreffmax,icceff,iucmax,ipot
      !call inter_2b_fit(ipot)
      !*************
      if (allocated(v_cc)) deallocate(v_cc)
      if (allocated(v_uc)) deallocate(v_uc)
      allocate(v_cc(0:iccdim,0:iccdim,jrelm))
      allocate(v_uc(0:nrelm,0:nrelm,iucdim))
      v_cc=0
      v_uc=0
      call inter_2b_fit_quadrature(ipot)
      !call add_residual_interaction(ipot)
      !call inter_2b_quadrature_qp(ipot)
      !call inter_2b_quadrature_qp_c(ipot)

      if(ipot==1) then
         call add_coulomb()
         print*,"coulomb added"
      endif

      !*************
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
      call eff_cc_inter_2b(ipot)
      call eff_uc_inter_2b(ipot)
      deallocate(v_cc,v_uc)
      call write_interaction_tofile(ipot)
      !*************
      deallocate(inter_cc,inter_uc)
      end subroutine two_body_matrix
!**********************
      subroutine eff_cc_inter_2b(ipot)
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
     & NE_maxi,NE_maxf,num1_states,num2_states,ipmax,ipot
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
     & ,Tkin(:,:)
      !states space in Heff: 
      !  (n,l): (0,lx1) (0,lx2) (1,lx) (1,lx2) ...
      complex(kd):: rrel
      integer:: ia,ib,nra,nrb,lra,lrb,i
      integer:: max_ncut,ncut,ipmax_eff,Nmx_eff,ipmax_e
      integer:: time3,time4,ir
      logical:: convergence
      !temporary untile writing csrg_kd
      complex(kind(0.d0)),allocatable:: Heff_db(:,:),Tkin_db(:,:)
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
      print*,multi
      jfmin=1
      jfmax=min0(jrelm,mxeng+1)
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
         !calculating the dimensions of the matrix starting from the energy cutoff: NE_maxi and NE_maxf
         ! and NE=2n+l
         NE_maxi=2*nrelm
         num1_states=(NE_maxi-lx1)/2+1
         if (lx1>NE_maxi) num1_states=0
         num2_states=(NE_maxi-lx2)/2+1
         if (lx2>NE_maxi) num2_states=0
         ipmax=num1_states+num2_states
         !
         NE_maxf=mxeng
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
         if(allocated(H_cc)) deallocate(H_cc,Tkin)
        

         allocate(H_cc(ipmax,ipmax),Tkin(ipmax,ipmax))
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
               !*************** !calculate the kinetic matrix elements
               if (lra==lrb) then
                  if (nra==nrb) then
                     rrel=(2*nra+lra+1.5_kd)*hbo/real(nucleons,kd)
                  elseif (nrb==nra+1) then
                     rrel=sqrt((nra+1_kd)*(nra+lra+1.5d0))
     &                    *hbo/real(nucleons,kd)
                  elseif (nra==nrb+1) then
                     rrel=sqrt((nrb+1_kd)*(nrb+lrb+1.5d0))
     &                    *hbo/real(nucleons,kd)
                  else
                     rrel=0._kd
                  endif
               else
                  rrel=complex(0._kd,0._kd)
               endif
               rrel=rrel!*complex(cos(2._kd*theta_kd),
!     &                              -sin(2._kd*theta_kd))
               !***************
            H_cc(ia,ib)= v_cc(ia-1,ib-1,jrel)+rrel
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
            Tkin(ia,ib)=rrel!**2*0.1
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
               H_cc(ib,ia)=H_cc(ia,ib)
            enddo
         enddo
         !**********
      !************************************
       if(jrel==1) then
        open(224,file="H.test")
         do ia=1,ipmax
            write(224,"(1000ES14.4)"),(real(H_cc(ib,ia)),ib=1,ipmax)
         enddo
        close(224)
       endif
      !************************************
         ! print the exact eigen values before SRG
         if(nucleons==2) then
            allocate(eigr(ipmax),eigv(ipmax,ipmax))
         call system_clock(count=time3, count_rate=ir)
            call mkl_Cgeev(H_cc,ipmax,eigr,eigv)
         call system_clock(count=time4, count_rate=ir)
         print *,'cgeev finished',
     &      real(time4 - time3,kind=8)/real(ir,kind=8)
            write(223,*),"Two-body hamiltonian exact eigenvalues"
            write(223,*),"Spectrum for Nmax=",ipmax-1
            write(223,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax,40))
            write(223,*),"End of spectrum"

           !TEST 
            open(224,file="2b_spect.dat",Access='append')
            write(224,*),"******************************"
            write(224,*),"Two-body hamiltonian exact eigenvalues"
            write(224,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax,40))
            close(224)
      !************************************
       if(jrel==1) then
        open(224,file="H.test")
         do ia=1,ipmax
            write(224,"(1000ES14.4)"),(real(H_cc(ib,ia)),ib=1,ipmax)
         enddo
        close(224)
       endif
      !************************************
            deallocate(eigr,eigv)
         endif
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
               allocate(Heff_db(ipmax,ipmax),Tkin_db(ipmax,ipmax))
               Heff_db(:,:)=H_cc(:,:)
!     &*complex(cos(2._kd*theta_kd),+sin(2._kd*theta_kd))
               Tkin_db(:,:)=Tkin(:,:)
     &*complex(cos(4._kd*theta_kd),+sin(4._kd*theta_kd))
               do ia=1,ipmax
                  do ib=ia+1,ipmax
               if(abs(Heff_db(ia,ib)-Heff_db(ib,ia))>10.d0**-13)then
                     print*,'non-sym',ia,ib,abs(Heff_db(ia,ib)
     &  -Heff_db(ib,ia))
                   endif
                  enddo
               enddo
               call system_clock(count=time3, count_rate=ir)
                  call csrg_evolution(ipmax,Heff_db,Tkin_db)
               call system_clock(count=time4, count_rate=ir)
               print *,'srg finished',
     &         real(time4 - time3,kind=8)/real(ir,kind=8)
               do ia=1,ipmax
                  do ib=ia+1,ipmax
               if(abs(Heff_db(ia,ib)-Heff_db(ib,ia))>10.d0**-13)then
                     print*,'non-sym',ia,ib,abs(Heff_db(ia,ib)
     &  -Heff_db(ib,ia))
                   endif
                  enddo
               enddo
               H_cc(:,:)=Heff_db(:,:)
!     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
               deallocate(Heff_db,Tkin_db)
            endif
         endif
         !**********
!         do ia=1,min(40,ipmax)
!            write(223,"(200ES14.4)"),(H_cc(ia,ib),ib=1
!     &             ,min(40,ipmax))
!         enddo
!         open(225,file="Hfit.test")
!         do ia=1,ipmax
!            write(225,"(1000F14.4)"), ( H_cc(ia,ib),ib=1,ipmax)
!         enddo
!         close(225)
         !**********
         ! print the exact eigen values before SRG
         if(nucleons==2) then
            allocate(eigr(ipmax),eigv(ipmax,ipmax))
         call system_clock(count=time3, count_rate=ir)
            call mkl_Cgeev(H_cc,ipmax,eigr,eigv)
         call system_clock(count=time4, count_rate=ir)
         print *,'cgeev finished',
     &      real(time4 - time3,kind=8)/real(ir,kind=8)
            write(223,*),"Two-body hamiltonian  eigenvalues after srg"
            write(223,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax,40))
      !************************************
       if(jrel==1) then
        open(224,file="eigv_aftersrg.test")
         do ia=1,5
            write(224,"(2ES13.3)"),(eigv(ib,ia),ib=1,ipmax)
            write(224,*),"*******************"
         enddo
        close(224)
       endif
      !************************************
            deallocate(eigr,eigv)
         endif
         !**********
         !print the spectrum and the convergence trend if requested
         max_ncut=0
         if(convergence) max_ncut=NE_maxf-6
         do Nmx_eff=NE_maxf-max_ncut,NE_maxf,2
            ncut=NE_maxf-Nmx_eff
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
            write(223,*),"Spectrum for Nmax=",ipmax_eff-1
            write(223,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax_eff,40))
            write(223,*),"End of spectrum"
            deallocate(eigr,eigv)

            do ia=1,ipmax_eff
               do ib=ia,ipmax_eff
                  inter_cc(ia-1,ib-1,jrel,ncut/2)=H_cc(ia,ib)
     &                                                -Tkin(ia,ib)
                  inter_cc(ib-1,ia-1,jrel,ncut/2)=inter_cc(ia-1
     &                                          ,ib-1,jrel,ncut/2)
               enddo
            enddo


            deallocate(Heff)
         enddo
         !**********
      
      enddo ! end of Jrel loop
      close(223)
      end subroutine eff_cc_inter_2b
!*******************************
!**********************
      subroutine eff_uc_inter_2b(ipot)
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
     & NE_maxi,NE_maxf,ipmax,ipot
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
     & ,Tkin(:,:)
      !states space in Heff: 
      !  (n,l): (0,lx1) (0,lx2) (1,lx) (1,lx2) ...
      complex(kd):: rrel
      integer:: ia,ib,nra,nrb,lra,lrb,iuc,lrel
      integer:: max_ncut,ncut,ipmax_eff,Nmx_eff
      integer:: time3,time4,ir
      logical:: convergence
      !temporary untile writing csrg_kd
      complex(kind(0.d0)),allocatable:: Heff_db(:,:),Tkin_db(:,:)
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
      print*,multi
      jfmin=0
      jfmax=min0(jrelm,mxeng+1)
      theta_kd=real(theta,kd)
      convergence=multi
      open(223,file="inter_2b_fit.out",Access='append')
      do jrel=jfmin,jfmax
      do jspin=0,1
         ! in uncoupled channels J=l1=l2 :
         lrel=jrel 
         ispin=0
         !maintain antisymmetry
         if( (-1)**(lrel+jspin+ispin) /= (-1)) ispin=1
         if ((ipot==1.or.ipot==3).and.ispin==0) cycle
         !calculating the dimensions of the matrix starting from the energy cutoff: NE_max
         ! and NE=2n+l
         NE_maxi=2*nrelm
         if(lrel>2*nrelm) then
            ipmax=0
         else
            ipmax=(NE_maxi-lrel)/2+1
         endif
         !
         NE_maxf=mxeng
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
         if(allocated(H_uc)) deallocate(H_uc,Tkin)
         
         allocate(H_uc(ipmax,ipmax),Tkin(ipmax,ipmax))
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
                  rrel=(2*nra+lrel+1.5_kd)*hbo/real(nucleons,kd)
               elseif (nrb==nra+1) then
                  rrel=sqrt((nra+1_kd)*(nra+lrel+1.5_kd))
     &                 *hbo/real(nucleons,kd)
               elseif (nra==nrb+1) then
                  rrel=sqrt((nrb+1_kd)*(nrb+lrel+1.5_kd))
     &                 *hbo/real(nucleons,kd)
               else
                  rrel=0._kd
               endif
               rrel=rrel*complex(cos(2._kd*theta_kd),
     &                              -sin(2._kd*theta_kd))
               !***************
               H_uc(ia,ib)=v_uc(nra,nrb,iuc)+rrel
               Tkin(ia,ib)=rrel
            enddo
         enddo
         !**********
         ! print the exact eigen values before SRG
         if(nucleons==2) then
            allocate(eigr(ipmax),eigv(ipmax,ipmax))
            call mkl_Cgeev(H_uc,ipmax,eigr,eigv)
            write(223,*),"Two-body hamiltonian exact eigenvalues"
            write(223,"(3ES16.5,3X,3X,3X,I3,'%')") ,(eigr(ia),atan
     &(aimag(eigr(ia))/real(eigr(ia))),abs(100-int(
     & 100*atan(aimag(eigr(ia))/real(eigr(ia)))/(
     & -2*theta))),ia=1,min(ipmax,40))
            deallocate(eigr,eigv)
         endif
         !**********
         !performing SRG
         allocate(Heff_db(ipmax,ipmax),Tkin_db(ipmax,ipmax))
         if(srg) then
            Heff_db(:,:)=H_uc(:,:)
     &*complex(cos(2._kd*theta_kd),+sin(2._kd*theta_kd))
            Tkin_db(:,:)=Tkin(:,:)
     &*complex(cos(2._kd*theta_kd),+sin(2._kd*theta_kd))
            !call csrg_evolution(ipmax,H_uc,Tkin)
         call system_clock(count=time3, count_rate=ir)
            call csrg_evolution(ipmax,Heff_db,Tkin_db)
         call system_clock(count=time4, count_rate=ir)
         print *,'srg evolution is finished',
     &      real(time4 - time3,kind=8)/real(ir,kind=8)
            H_uc(:,:)=Heff_db(:,:)
     &*complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
         endif
         deallocate(Heff_db,Tkin_db)
         !**********
         !print the spectrum and the convergence trend if requested
         max_ncut=0
         if(convergence) max_ncut=NE_maxf-6
         do Nmx_eff=NE_maxf-max_ncut,NE_maxf,2
            ncut=NE_maxf-Nmx_eff
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
                  inter_uc(ia-1,ib-1,iuc,ncut/2)=H_uc(ia,ib)
     &                                                -Tkin(ia,ib)
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
      subroutine write_interaction_tofile(ipot)
      use paramdef
      use nninteraction
      use harmosc
      use HOtrap, only: HO_trap,V_HO
      use H_matrices
      use quadratures
      use HO_wf
      implicit none
      integer,parameter:: kd=real128
      integer,parameter:: dp=kind(0.d0)
      integer:: tz,ipot,ipmax,jrel,jfmin,jfmax,jspin,ispin
      integer:: i,ia,ib,nra,lra,nrb,lrb,lx1,lx2,lrel,iuc
     & ,num1_states,num2_states,NE_max,max_ncut
      complex(dp),allocatable:: Tkin
      real(kd)::theta_kd
      logical:: convergence


      real(kd),allocatable:: roots(:,:,:)
      complex(kd):: zero_cmplx
      zero_cmplx=complex(0.d0,0.d0)
      jfmin=1
      jfmax=min0(jrelm,mxeng+1)
      max_ncut=0
      if(multi) max_ncut=2*nrelm-6
      theta_kd=real(theta,kd)
      !*******************
      ! generate coulomb matrix: vcoul
      if (allocated(roots)) deallocate(roots)
      allocate(roots(nrelm,nrelm,0:lrelm))
      do lra=0,lrelm
         call laguerre_rootfinder_kd(nrelm,real(lra+0.5_kd,kd)
     &                         ,roots(1:nrelm,1:nrelm,lra))
      end do
      if (ipot==1) then
         if(allocated(vcoul)) deallocate(vcoul)
         call coulomb_c(mxeng/2,min(mxeng,lmax),roots)
      endif
      !*******************
      do jrel=jfmin,min(mxeng+1,lmax-1)
      
      !*****************
         lx1=jrel-1
         lx2=jrel+1
         jspin=1
         ispin=0

         if (ipot==1) then
            tz=1
         elseif(ipot==3) then
            tz=-1
         else
            tz=0
         endif
         if( (-1)**(lx1+jspin+ispin) /= (-1)) ispin=1
         if ((ipot==1.or.ipot==3).and.ispin==0) cycle
         !
         NE_max=mxeng
         num1_states=(NE_max-lx1)/2+1
         if (lx1>NE_max) num1_states=0
         num2_states=(NE_max-lx2)/2+1
         if (lx2>NE_max) num2_states=0
         ipmax=num1_states+num2_states
         !
      !*****************
         do ia=1,ipmax
            nra=(ia-1)/2
            if (mod(ia,2)==0) then
               lra=lx2
            else
               lra=lx1
            endif
            do ib=ia,ipmax
               nrb=(ib-1)/2
               if (mod(ib,2)==0) then
                  lrb=lx2
               else
                  lrb=lx1
               endif
               !*************** !calculate the kinetic matrix elements
               if (lra==lrb) then
                  if (nra==nrb) then
                     Tkin=(2*nra+lra+1.5_kd)*hbo/real(nucleons,kd)
                  elseif (nrb==nra+1) then
                     Tkin=sqrt((nra+1_kd)*(nra+lra+1.5d0))
     &                    *hbo/real(nucleons,kd)
                  elseif (nra==nrb+1) then
                     Tkin=sqrt((nrb+1_kd)*(nrb+lrb+1.5d0))
     &                    *hbo/real(nucleons,kd)
                  else
                     Tkin=0._kd
                  endif
               else
                  Tkin=complex(0._kd,0._kd)
               endif
               Tkin=Tkin
     & *complex(cos(2._kd*theta_kd),-sin(2._kd*theta_kd))
               !***************
               if(jrel<=jfmax) then
                  write(3) jspin,jrel,ispin,
     & nra,lra,nrb,lrb,(cmplx(inter_cc(ia-1,ib-1,jrel,i),    
     & KIND=dp),i=0,max_ncut/2)
                  write(2,6767) jspin,jrel,ispin,
     & nra,lra,nrb,lrb,(cmplx(inter_cc(ia-1,ib-1,jrel,i),
     & KIND=dp),i=0,max_ncut/2)
 6767             format(7i4,30f15.8)
                  write(37) jspin,jrel,ispin,nra,lra,nrb,lrb,
     &         cmplx(inter_cc(ia-1,ib-1,jrel,0)+Tkin,KIND=dp)
                  write(377,6767) nra,lra,nrb,lrb,jspin,jrel,-tz,
     &         cmplx(inter_cc(ia-1,ib-1,jrel,0)+Tkin,KIND=dp)
               else
                     !noho is assumed True and 0 is printed
                  write(3) jspin,jrel,ispin,
     &              nra,lra,nrb,lrb,(zero_cmplx,i=0,max_ncut/2)
                  write(2,6767) jspin,jrel,ispin,
     &              nra,lra,nrb,lrb,(zero_cmplx,i=0,max_ncut/2)
                  if(ipot==1.and.lra==lrb) then
                     write(37) jspin,jrel,ispin,nra,lra,nrb,lrb
     &                        ,cmplx(vcoul(nra,nrb,lra),KIND=dp)
                  endif
               endif
               !***************
            end do
         end do
      !*****************
      enddo
      ! the uncpled part
      jfmin=0
      jfmax=min0(jrelm,mxeng+1)
      do jrel=jfmin,min(mxeng,lmax)
      do jspin=0,1 
      !*****************
         lrel=jrel
         lra=lrel
         lrb=lrel
         ispin=0
         if( (-1)**(lrel+jspin+ispin) /= (-1)) ispin=1
         if ((ipot==1.or.ipot==3).and.ispin==0) cycle
         NE_max=mxeng
         !
         if(lrel>2*nrelm) then
            ipmax=0
         else
            ipmax=(NE_max-lrel)/2+1
         endif
         !
         iuc=2*jrel+jspin+1
      !*****************
         do ia=1,ipmax
            nra=ia-1
            do ib=1,ipmax
               nrb=ib-1
               !***************
               !calculate the kinetic matrix elements
               if (nra==nrb) then
                  Tkin=(2*nra+lrel+1.5_kd)*hbo/real(nucleons,kd)
               elseif (nrb==nra+1) then
                  Tkin=sqrt((nra+1_kd)*(nra+lrel+1.5_kd))
     &                 *hbo/real(nucleons,kd)
               elseif (nra==nrb+1) then
                  Tkin=sqrt((nrb+1_kd)*(nrb+lrel+1.5_kd))
     &                 *hbo/real(nucleons,kd)
               else
                  Tkin=0._kd
               endif
               Tkin=Tkin*complex(cos(2._kd*theta_kd),
     &                              -sin(2._kd*theta_kd))
               !***************
               if(jrel<=jfmax) then
                        write(3)      jspin,jrel,ispin,
     & nra,lra,nrb,lrb,(cmplx(inter_uc(nra,nrb,iuc,i)
     &,KIND=dp),i=0,max_ncut/2)
                        write(2,6767) jspin,jrel,ispin,
     & nra,lra,nrb,lrb,(cmplx(inter_uc(nra,nrb,iuc,i),
     & KIND=dp),i=0,max_ncut/2)
                        write(37) jspin,jrel,ispin,nra,lra,nrb,lrb,
     & cmplx(inter_uc(nra,nrb,iuc,0)+Tkin,KIND=dp)
                        write(377,6767) nra,lra,nrb,lrb,jspin,jrel,-tz,
     & cmplx(inter_uc(nra,nrb,iuc,0)+Tkin,KIND=dp)
               else

                  write(3) jspin,jrel,ispin,
     & nra,lra,nrb,lrb,(zero_cmplx,i=0,max_ncut/2)
                  write(2,6767) jspin,jrel,ispin,
     & nra,lra,nrb,lrb,(zero_cmplx,i=0,max_ncut/2)
!     +           ,-urel(0)/hbo3,-urel(0)/hbo3

                  if (ipot==1.and.lra==lrb) then
                     write(37) jspin,jrel,ispin,nra,lra,nrb,lrb,
     &               cmplx(vcoul(nra,nrb,lra),KIND=dp)
                  endif
               endif
               !***************
            enddo
         enddo
      !*****************
      enddo
      enddo
      if(allocated(vcoul)) deallocate(vcoul)
      endsubroutine write_interaction_tofile
!*******************************
!*******************************
      subroutine add_coulomb()
      use iso_fortran_env
      use paramdef
      use harmosc
      use nninteraction
      use quadratures
      use HO_wf
      implicit none
      integer,parameter :: kd=real128
      real(kd),allocatable:: roots(:,:,:)
      double precision::cnp,cpp,cnn
      integer:: nra,lra,nrb,jrel,ispin,nrbmax,ipp,iuc
      real(kd):: th_kd
      cpp=1.d0
      ipp=1
      th_kd=theta
      if (allocated(roots)) deallocate(roots)
      allocate(roots(nrelm,nrelm,0:lrelm))
      do lra=0,lrelm
         call laguerre_rootfinder_kd(nrelm,real(lra+0.5_kd,kd)
     &                         ,roots(1:nrelm,1:nrelm,lra))
      end do
      if (ipp/=0) then
         if(allocated(vcoul)) deallocate(vcoul)
      endif
      print*,"nrelm,lrelm",nrelm,lrelm
      !*********************
      !*********************
      call coulomb_c(nrelm,lrelm,roots)
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
      if(allocated(vcoul)) deallocate(vcoul)
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
      !*********************
      !*********************
      allocate(uc(0:nstep,0:nrelm,0:lrelm))
      do lr=0,lrelm
         do nr=0,nrelm
            if (2*nr+lr>2*nrelm) cycle
            do i=1,nstep
               rr=gaupoi_pp(i)
               call wf_coulomb(nr,lr,anu_kd,rr,uc(i,nr,lr),
     &                               roots(:nr,nr,lr))
            end do
         end do
      end do
      uc=uc/sqrt(2._kd)

      print *,' u-coulomb calculated'
      !*********************
      !*********************

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
                  rr=gaupoi_pp(i)
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
      !*********************
      !*********************

      print *, ' Coulomb calculated'

      deallocate(uc)
      end subroutine coulomb_c

      end module inter_2b
