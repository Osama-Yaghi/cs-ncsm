!*************** n3lo *********************************************
       module n3lo_qp_c_module
       use iso_fortran_env
       real(real128):: gaus_legendre_points(1:200,1:200),
     & gaus_legendre_weights(1:200,1:200)
     &,legendre_polynom(1:200,1:200,0:8)
       real(real64):: reg_n=10
       contains
         subroutine regulator_poly(px,py,reg,lambda,n)
            use iso_fortran_env
            implicit none
            integer,parameter:: qp=real64
            complex(qp),intent(in)::px,py
            complex(qp),intent(out):: reg
            complex(qp)::lambda,n
            complex(qp)::power_factor,factor_p,factor_n
            complex(qp)::m,cutoff
            m=n
            if(real(m)>5.9) then
            m=6._qp
            cutoff=lambda*1.25_qp
            power_factor=0.5_qp*( (px/cutoff)**m+(py/cutoff)**m)
     &    + 0.5_qp*( (px/cutoff)**(m-1)+(py/cutoff)**(m-1))
            elseif(real(m)>4.9) then
            m=5._qp
            cutoff=lambda*1.2_qp
            power_factor= (px/cutoff)**m+(py/cutoff)**m
            else
            m=4.7_qp
            cutoff=lambda*1.2_qp
            cutoff=lambda*1.229_qp
            power_factor= (px/cutoff)**m+(py/cutoff)**m
            endif
            reg=(1._qp)/(1._qp+ power_factor)**3._qp
         end subroutine regulator_poly
         subroutine regulator_poly_2(px,py,reg,lambda,n)
            use iso_fortran_env
            implicit none
            integer,parameter:: qp=real64
            complex(qp),intent(in)::px,py
            complex(qp),intent(out):: reg
            complex(qp)::lambda,n
            complex(qp)::power_factor,factor_p,factor_n
            complex(qp)::m,cutoff
            m=n
            if(real(m)>2.9) then
            m=3._qp
            cutoff=lambda*1.56_qp
            power_factor=0.5_qp*( (px/cutoff)**m+(py/cutoff)**m)
     &    + 0.5_qp*( (px/cutoff)**(m-0.5)+(py/cutoff)**(m-0.5))
            elseif(real(m)>2.4) then
            m=2.5_qp
            cutoff=lambda*1.44_qp
            power_factor= (px/cutoff)**m+(py/cutoff)**m
            else
            m=2.35_qp
            cutoff=lambda*1.44_qp
            cutoff=lambda*1.5104_qp
            power_factor= (px/cutoff)**m+(py/cutoff)**m
            endif
            reg=(1._qp)/(1._qp+ power_factor)**3._qp
         end subroutine regulator_poly_2
         subroutine regulator_c(px,py,reg,lambda,n)
            use iso_fortran_env
            implicit none
            integer,parameter:: qp=real64
            complex(qp),intent(in)::px,py
            complex(qp),intent(out):: reg
            complex(qp)::lambda,n
            complex(qp)::power_factor,factor_p,factor_n
            complex(qp)::m
            m=n
            if(real(m)>4) m=4._qp
            power_factor= (px/lambda)**m+(py/lambda)**m
            if (abs(real(power_factor))>80) then
               power_factor=80._qp/abs(real(power_factor))
     &                           *power_factor
            endif
            if(real(power_factor)>0) then
               factor_p=power_factor
               factor_n=-power_factor
            else
               factor_p=-power_factor
               factor_n=power_factor
            endif
            reg=(1._qp+exp(10*factor_n))/(1._qp+ exp(12*factor_n))
     &               * exp(factor_n)
         end subroutine regulator_c
         subroutine regulator_c2(exp_factor,reg,n)
            use iso_fortran_env
            implicit none
            integer,parameter:: qp=real64
            complex(qp),intent(in)::exp_factor
            complex(qp),intent(out):: reg
            complex(qp)::lambda,n
            complex(qp)::power_factor,factor_p,factor_n
            complex(qp)::m
            m=n
            if(real(m)>2.5) m=2.5_qp
            power_factor= exp_factor**m
            if (abs(real(power_factor))>80) then
               power_factor=80._qp/abs(real(power_factor))
     &                           *power_factor
            endif
            if(real(power_factor)>0) then
               factor_p=power_factor
               factor_n=-power_factor
            else
               factor_p=-power_factor
               factor_n=power_factor
            endif
            reg=(1._qp+exp(10*factor_n))/(1._qp+ exp(12*factor_n))
     &               * exp(factor_n)
         end subroutine regulator_c2
         subroutine regulator_c3(px,py,reg,lambda,n)
            use iso_fortran_env
            implicit none
            integer,parameter:: qp=real64
            complex(qp),intent(in)::px,py
            complex(qp),intent(out):: reg
            complex(qp)::lambda,n
            complex(qp)::power_factor,factor_p,factor_n
            complex(qp)::m
            m=n
            if(real(m)>2.5) m=2.5_qp
            power_factor= (px/lambda)**m+(py/lambda)**m
            if (abs(real(power_factor))>80) then
               power_factor=80._qp/abs(real(power_factor))
     &                           *power_factor
            endif
            if(real(power_factor)>0) then
               factor_p=power_factor
               factor_n=-power_factor
            else
               factor_p=-power_factor
               factor_n=power_factor
            endif
            reg=(1._qp+exp(10*factor_n))/(1._qp+ exp(12*factor_n))
     &               * exp(factor_n)
         end subroutine regulator_c3
         subroutine n3lo_qp_c
!
!******************************************************************
!
!        FEBRUARY 2003
!        (FORTRAN improved 6/30/05)
!
!******************************************************************
!
!        This code computes the
!
!        Charge-Dependent Chiral NN Potential at Order Four (N3LO)
!        ---------------------------------------------------------
!        in momentum space.
!
!        this package is self-contained and includes
!        all subroutines needed.
!        only `n3lo' needs to be called by the user.
!        all codes are consistently in double precision.
!        when working on an UNIX/LINUX system, it is recommended
!        to compile this code with the  -static  option.
!        more information on the code is given below.
!
!*******************************************************************
!
!        authors:     D. R. Entem and R. Machleidt
!                     department of physics
!                     university of idaho
!                     moscow, idaho 83844
!                     u. s. a.
!                     e-mails: dentem@uidaho.edu
!                              machleid@uidaho.edu
!
!        publication: PRC 68, 041001 (2003).
!
!*******************************************************************
!
!
!
!
      use iso_fortran_env
      implicit complex*16 (a-h,o-z)
      integer,parameter:: qp=real64
!
!
      common /crdwrt_qp/ kread,kwrite,kpunch,kda(9)
!
!        arguments and values of this subroutine
!
      common /cpot_qp/   v,xmev,ymev
      common /cstate_qp/ j,heform,sing,trip,coup,endep,label
      common /cnn_qp/ inn
!
!
!        this has been the end of the common-blocks containing
!        the arguments and values of this subroutine.
!
!        specifications for these common blocks
!
      complex*16:: v(6),xmev,ymev
      logical heform,sing,trip,coup,endep
      character*4 label
!
!
!*********************************************************************
!        THE ABOVE FOUR COMMON BLOCKS IS ALL THE USER NEEDS
!        TO BE FAMILIAR WITH.
!*********************************************************************
!
!        here are now some explanations of what those common blocks contain:
!        -------------------------------------------------------------------
!
!        xmev and ymev are the final and initial relative momenta,
!        respectively, in units of mev/c.
!        v is the potential in units of mev**(-2).
!        concerning units, factors of pi, etc.,
!        cf. with the partial-wave Lippmann-Schwinger equation, Eq. (A25),
!        and with the phase shift relation, Eq. (A33), given in Appendix A
!        of the article: R. Machleidt, PRC 63, 024001 (2001).
!
!        the partial-wave Lippmann-Schwinger equation for the
!        K-matrix reads:
!
!        K(q',q) = V(q',q) + M P \int dk k^2 V(q',k) K(k,q)/(q^2-k^2)
!
!        with M the nucleon mass in MeV and P denoting the principal value;
!        V(q',q) as provided by this code in common block /cpot/;
!        all momenta in MeV.
!
!        the phase-shift relation is:
!
!        tan \delta_L = -(pi/2) M q K_L(q,q)
!
!        with M and q in units of MeV, K_L in MeV**(-2) like V.
!
!
!        if heform=.true., v contains the 6 matrix elements
!        associated with one j in the helicity formalism
!        in the following order:
!        0v, 1v, 12v, 34v, 55v, 66v
!        (for notation see Appendix A of above article).
!
!        if heform=.false., v contains the 6 matrix elements
!        associated with one j in the lsj formalism
!        in the following order:
!        0v(singlet), 1v(uncoupled triplet), v++, v--, v+-, v-+ (coupled)
!        (for notation, see explanations given in the above article
!        below Eq. (A31)).
!
!        j is the total angular momentum. there is essentially no upper
!        limit for j.
!        sing, trip, and coup should in general be .true..
!        endep and label can be ignored.
!        it is customary, to set kread=5 and kwrite=6;
!        ignore kpunch and kda(9).
!
!        the meaning of the parameter inn in the common block
!
!                  common /cnn/ inn
!        is
!                  inn=1  means pp potential,
!                  inn=2  means np potential, and
!                  inn=3  means nn potential.
!
!        the user needs to include this common block in his/her code,
!        and specify which potential he/she wants to use.
!
!
!        THIS IS ESSENTIALLY ALL THE USER NEEDS TO KNOW.
!
!        if you have further questions, do not hesitate to contact one
!        of the authors (see e-mail addresses above).
!
!**********************************************************************
!
!
!        common block for all chi-subroutines
!
      common /cchi_qp/ vj(32,270),c(20,270),fff,ff,f(52),aa(200)
     &,ai(19,30),
     &                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     &                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(200),
     &wt(200),
     &                ic(20,270),ift(3),mint(3),maxt(3),nt,
     &                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     &                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     &                indc(2,270),indpar(3),indxy
!
!         specifications for this common block
!
      logical indc,indxy,indpar
!
      common /comlsj_qp/ clsj(15,50),cutlsj(15,50),indlsj
      logical indlsj
!
      common /crrr_qp/ rrr
!
!
!        further specifications
!
      dimension vl(4),adminv(4,4),ldminv(4),mdminv(4)
      dimension vv0(6),vv2(6),vv4(6)
      character*4 nucnuc(3)
      character*4 mesong(40)
      logical index
      logical indmg(40)
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     &            '1-  ','1-t ','1-tt','1-st','1-ss',
     &            'c   ','ss  ','ls  ','sq  ','sk  ',
     &            'sl  ',
     &         24*'    '/
      data index/.false./
      data indmg/40*.false./
      data jj/-1/
!      data pi/3.141592653589793_qp/
      data pi/3.141592653589793238462643383279502884197_qp/
      data innn/-1/
      data nucnuc/'N3pp','N3np','N3nn'/
      save
!
!
!
!
      if (index) go to 10
      index=.true.
!
!
!        call subroutine chipar once and only once
!
!
      call chipar_qp
!     -----------
!
!
!        if you want the potential to be zero for very large momenta,
!        choose rrr=1000.
!        if you want no technical problems in the calculation of the deuteron
!        wave functions, choose rrr=80.
!
      rrr=80_qp
!      rrr=1000._qp
!
   10 continue
!
!
!
!
      if (inn.lt.1.or.inn.gt.3) then
!        choose the np potential as the default:
      inn=2
      endif
      if (j.lt.0) then
      write (kwrite,19002)
19002 format (////' error in n3lo: total angular momentum j',
     &' is negative.'/' execution terminated.'////)
      stop
      endif
!
!
!
!
!        set the inn dependent parameters
!
      if (inn.eq.innn) go to 30
      innn=inn
      inter=inn
      label=nucnuc(inter)
!
!**      go to (21,22,23), inter
!**   21 write (kwrite,10001)
!**10001 format (' The pp potential is used.')
!**      go to 24
!**   22 write (kwrite,10002)
!**10002 format (' The np potential is used.')
!**      go to 24
!**   23 write (kwrite,10003)
!**10003 format (' The nn potential is used.')
!**   24 write (kwrite,10004)
!**10004 format (' -------------------------'/)
!
!
      iftgo=ift(inter)+1
      dwn=1._qp/wnn(inter)
!
!
!        prepare constant over-all factor
!
      fac=pi/(2._qp*pi)**3._qp*dwn*dwn
!     ---------------------------
!
!
!
      iman=imaa(inter)
      imen=imea(inter)
!
      imanm1=iman-1
!
      iman1=imanm1+1
      iman2=imanm1+2
      iman3=imanm1+3
      iman4=imanm1+4
      iman5=imanm1+5
      iman6=imanm1+6
      iman7=imanm1+7
      iman8=imanm1+8
      iman9=imanm1+9
      imen24=imen-24
      imen23=imen-23
      imen22=imen-22
      imen21=imen-21
      imen15=imen-15
      imen14=imen-14
!
!
      ez1=ezz1(inter)
      ez2=ezz2(inter)
!
!
!
   30 if (j.eq.jj) go to 50
      jj=j
      if (j.eq.0) go to 50
      aj=real(j,qp)
      aj1=real(j+1,qp)
      a2j1=real(2*j+1,qp)
      aaj6=sqrt(aj*aj1)
!
!        coefficient matrix for the translations into lsj formalism
!
      adminv(1,1)=aj1
      adminv(1,2)=aj
      adminv(1,3)=-aaj6
      adminv(1,4)=-aaj6
      adminv(2,1)=aj
      adminv(2,2)=aj1
      adminv(2,3)=aaj6
      adminv(2,4)=aaj6
      adminv(3,1)=aaj6
      adminv(3,2)=-aaj6
      adminv(3,3)=aj1
      adminv(3,4)=-aj
      adminv(4,1)=aaj6
      adminv(4,2)=-aaj6
      adminv(4,3)=-aj
      adminv(4,4)=aj1
!
!       inversion
!
      call dminv_qp (adminv,4,deter,ldminv,mdminv)
!
!
!
!
!        prepare expressions depending on x and y
!        ----------------------------------------
!        ----------------------------------------
!
!
!
!
   50 xa=xmev*dwn
      ya=ymev*dwn
      indxy=.false.
      x=xa
      xx=x*x
      y=ya
      yy=y*y
      xy2=x*y*2._qp
      xxpyy=xx+yy
      ex=sqrt(1._qp+xx)
      ey=sqrt(1._qp+yy)
      eem12=(ex*ey-1._qp)*2._qp
!
!
      xy=xy2*0.5_qp
      ee=ex*ey
      ree=sqrt(ee)
      eem1=ee-1._qp
      eme=ex-ey
      emeh=eme*0.5_qp
      emehq=emeh*emeh
      eep1=ee+1._qp
       epe=ex+ey
      xxyy=xx*yy
!
!
      xxpyyh=xxpyy*0.5_qp
      xy3=xy*3._qp
      xy4=xy*4._qp
!
!
!
!
      do 63 iv=1,6
      vv0(iv)=0._qp
      vv2(iv)=0._qp
      vv4(iv)=0._qp
   63 v(iv)=0._qp
      do 65 il=iman,imen
      do 65 iv=1,32
   65 vj(iv,il)=0._qp
!
!
!
!
!        prepare over-all factor
!
!
      go to (70,71,72,71,72,75,76),iftgo
!
!        no additional factor
!
   70 fff=fac
      go to 80
!
!        minimal relativity
!
   71 fff=fac/ree
      go to 80
!
!        factor m/e*m/e
!
   72 fff=fac/ee
      go to 80
!
!        sharp cutoff
!
   75 if (real(xmev).gt.real(ez1).or.real(ymev).gt.real(ez1)) then
      return
      else
      fff=fac
      end if
      go to 80
!
!        exponential cutoff
!
   76 expo=(xmev/ez1)**(2._qp*ez2)+(ymev/ez1)**(2._qp*ez2)
      if (real(expo).gt.real(rrr)) then
      expo=rrr
      end if
      expexp=exp(-expo)
      call regulator_poly(xmev,ymev,expexp,ez1,2._qp*ez2) ! call regulator function that stores the factor in expexp
      fff=fac*expexp
!
!
   80 continue
!
!
!
!
!        contributions
!        -------------
!        -------------
!
!
!
!
      do 5995 img=1,mge
      mg=mggo(img,inter)
      if (mg.gt.16) go to 9000
      if (mg.eq.0) go to 8000
      me=mgg(mg,inter)
      go to (9000,9000,9000,9000,9000,9000,9000,9000,9000,9000,
     &       1100,1200,1300,1400,1500,1600),mg
!
!
!
!
!        c   , central force
!        -------------------
!
!
!
!
 1100 mc=1
!
      ff=1._qp
      f(1)=2._qp
      f(2)=0._qp
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-f(1)
      f(8)=f(7)
!
       !test
      call chistr_qp(1,1,me)
      go to 5995
!
!
!
!
!        ss  , spin-spin force
!        ---------------------
!
!
!
!
 1200 mc=1
!
      ff=1._qp
      f(1)=-6._qp
      f(2)=0._qp
      f(3)=2._qp
      f(4)=0._qp
      f(5)=0._qp
      f(6)=f(3)
      f(7)=-f(3)
      f(8)=f(7)
!
      call chistr_qp(1,1,me)
      go to 5995
!
!
!
!
!        ls  , spin-orbit force
!        ----------------------
!
!
!
!
 1300 mc=1
!
      ff=1._qp
      f(1)=0._qp
      f(2)=0._qp
      f(3)=0._qp
      f(4)=-xy2
      f(5)=-xy2
      f(6)=0._qp
      f(7)=0._qp
      f(8)=0._qp
      f(9)=0._qp
      f(10)=+xy2
      f(11)=-xy2
!
      call chistr_qp(2,1,me)
      go to 5995
!
!
!
!
!        sq  , sq tensor force (where q denotes the momentum transfer)
!        ---------------------
!
!
!
!
 1400 mc=1
!
      ff=1._qp
      f(1)=-xxpyy*2.0_qp
      f(2)=xy*4._qp
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=(xx-yy)*2.0_qp
      f(8)=-f(7)
!
      call chistr_qp(1,1,me)
      go to 5995
!
!
!
!
!        sk  , sk tensor force (where k denotes the average momentum)
!        ---------------------
!
!
!
!
 1500 mc=1
!
      ff=0.25_qp
      f(1)=-xxpyy*2.0_qp
      f(2)=-xy*4._qp
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=(xx-yy)*2.0_qp
      f(8)=-f(7)
!
      call chistr_qp(1,1,me)
      go to 5995
!
!
!
!
!        sl  , "quadratic spin-orbit force"
!               or sigma-l operator
!        ----------------------------------
!
!
!
!
 1600 mc=1
!
      ff=1._qp
      f(1)=-xxyy*2._qp
      f(2)=0._qp
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=-f(1)
      f(7)=f(1)
      f(8)=f(7)
      f(9)=f(6)*2._qp
!
      call chistr_qp(4,1,me)
      go to 5995
!
!
!
!
!
!        this has been the end of the contributions of mesons
!        ----------------------------------------------------
!
!
!
!
!        errors and warnings
!        -------------------
!
!
!
!
 9000 if (indmg(mg)) go to 5995
!**** write (kwrite,19000) mesong(mg)
19000 format(1h ////' warning in chinn: contribution ',a4,' does not exi
     &st in this program.'/' contribution ignored. execution continued.'
     &////)
      indmg(mg)=.true.
!
!
!
!
 5995 continue

!
!
!
!
!        add up contributions
!        --------------------
!
!
!
!
 8000 continue
!
!
!        charge-dependent OPE contribution
!        ---------------------------------
!
      if (mod(j,2).eq.1) go to 8020
!
!        j even
!
      v(1)=-vj(1,iman1)+2._qp*vj(1,iman5)
      v(1)=v(1)-vj(1,iman2)+2._qp*vj(1,iman6)
      v(1)=v(1)-vj(1,iman3)+2._qp*vj(1,iman7)
      v(1)=v(1)-vj(1,iman4)+2._qp*vj(1,iman8)
!
      v(2)=-vj(2,iman1)-2._qp*vj(2,iman5)
      v(2)=v(2)-vj(2,iman2)-2._qp*vj(2,iman6)
      v(2)=v(2)-vj(2,iman3)-2._qp*vj(2,iman7)
      v(2)=v(2)-vj(2,iman4)-2._qp*vj(2,iman8)
!
      do 8015 iv=3,6
      v(iv)=-vj(iv,iman1)+2._qp*vj(iv,iman5)
      v(iv)=v(iv)-vj(iv,iman2)+2._qp*vj(iv,iman6)
      v(iv)=v(iv)-vj(iv,iman3)+2._qp*vj(iv,iman7)
      v(iv)=v(iv)-vj(iv,iman4)+2._qp*vj(iv,iman8)
 8015 continue
      go to 8030
!
!        j odd
!
 8020 continue
      v(1)=-vj(1,iman1)-2._qp*vj(1,iman5)
      v(1)=v(1)-vj(1,iman2)-2._qp*vj(1,iman6)
      v(1)=v(1)-vj(1,iman3)-2._qp*vj(1,iman7)
      v(1)=v(1)-vj(1,iman4)-2._qp*vj(1,iman8)
!
      v(2)=-vj(2,iman1)+2._qp*vj(2,iman5)
      v(2)=v(2)-vj(2,iman2)+2._qp*vj(2,iman6)
      v(2)=v(2)-vj(2,iman3)+2._qp*vj(2,iman7)
      v(2)=v(2)-vj(2,iman4)+2._qp*vj(2,iman8)
!
      do 8025 iv=3,6
      v(iv)=-vj(iv,iman1)-2._qp*vj(iv,iman5)
      v(iv)=v(iv)-vj(iv,iman2)-2._qp*vj(iv,iman6)
      v(iv)=v(iv)-vj(iv,iman3)-2._qp*vj(iv,iman7)
      v(iv)=v(iv)-vj(iv,iman4)-2._qp*vj(iv,iman8)
 8025 continue
!
!
 8030 continue
!
!
      if (iman9.gt.imen) go to 8500
!
!
      if (.not.indlsj) then
      do 8105 il=iman9,imen
      do 8105 iv=1,6
 8105 v(iv)=v(iv)+vj(iv,il)
      else
!
!
!        there are contact terms
!        -----------------------
!
      if (iman9.gt.imen24) go to 8200
!
!        the non-contact terms
!
      do 8155 il=iman9,imen24
      do 8155 iv=1,6
 8155 v(iv)=v(iv)+vj(iv,il)
!
!        contact contributions
!        ---------------------
!
 8200 continue
!
!        Q^0 contacts
      do 8205 il=imen23,imen22
      do 8205 iv=1,6
 8205 vv0(iv)=vv0(iv)+vj(iv,il)
!
!        Q^2 contacts
      do 8215 il=imen21,imen15
      do 8215 iv=1,6
 8215 vv2(iv)=vv2(iv)+vj(iv,il)
!
!        Q^4 contacts
      do 8225 il=imen14,imen
      do 8225 iv=1,6
 8225 vv4(iv)=vv4(iv)+vj(iv,il)
!
!
!        ------------------------------------------------------
!        NOTE: partial-wave potentials that add-up to zero need
!        to be cutoff, because they diverge for large momenta.
!        ------------------------------------------------------
!
!        use 3d3 cutoff as default for all j.gt.3 partial waves
!
      if (j.gt.3) then
      if (cutlsj(1,15).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(2,15))**(2._qp*cutlsj(1,15))
     &    +(ymev/cutlsj(2,15))**(2._qp*cutlsj(1,15))
      if (real(expo).gt.real(rrr)) expo=rrr
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      expexp=exp(-expo)
      call regulator_poly(xmev,ymev,expexp,cutlsj(2,15)
     &                     ,2._qp*cutlsj(1,15)) ! call regulator function that stores the factor in expexp
      end if
!
      do 8275 iv=1,6
      vv0(iv)=vv0(iv)*expexp
      vv2(iv)=vv2(iv)*expexp
 8275 vv4(iv)=vv4(iv)*expexp
      go to 8400
      end if
!
!
!        look into individual partial waves and
!        multiply with partial-wave dependent cutoffs
!        --------------------------------------------
!
      j1=j+1
      go to (8310,8320,8330,8340),j1
!
!
!        j=0
!        ---
!        ---
!
 8310 continue
!
!        1s0
!        ---
!        Q^0 term
!
      if (cutlsj(1,1).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(2,1))**(2._qp*cutlsj(1,1))
     &    +(ymev/cutlsj(2,1))**(2._qp*cutlsj(1,1))
      if (real(expo).gt.real(rrr)) expo=rrr
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      expexp=exp(-expo)
      call regulator_poly(xmev,ymev,expexp,cutlsj(2,1)
     &                     ,2._qp*cutlsj(1,1)) ! call regulator function that stores the factor in expexp
      end if
      vv0(1)=vv0(1)*expexp
!
!        Q^2 terms
!
      if (cutlsj(3,1).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(4,1))**(2._qp*cutlsj(3,1))
     &    +(ymev/cutlsj(4,1))**(2._qp*cutlsj(3,1))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(4,1)
     &                     ,2._qp*cutlsj(3,1)) ! call regulator function that stores the factor in expexp
      end if
      vv2(1)=vv2(1)*expexp
!
!        Q^4 terms
!
      if (cutlsj(5,1).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(6,1))**(2._qp*cutlsj(5,1))
     &    +(ymev/cutlsj(6,1))**(2._qp*cutlsj(5,1))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(6,1)
     &                     ,2._qp*cutlsj(5,1)) ! call regulator function that stores the factor in expexp
      end if
      vv4(1)=vv4(1)*expexp
!
!        3p0
!        ---
!        Q^2 term
!
      if (cutlsj(1,2).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(2,2))**(2._qp*cutlsj(1,2))
     &    +(ymev/cutlsj(2,2))**(2._qp*cutlsj(1,2))
      if (real(expo).gt.real(rrr)) expo=rrr
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      expexp=exp(-expo)
      call regulator_poly(xmev,ymev,expexp,cutlsj(2,2)
     &                     ,2._qp*cutlsj(1,2)) ! call regulator function that stores the factor in expexp
      end if
      vv2(3)=vv2(3)*expexp
      vv0(3)=vv0(3)*expexp
!
!        Q^4 term
!
      if (cutlsj(3,2).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(4,2))**(2._qp*cutlsj(3,2))
     &    +(ymev/cutlsj(4,2))**(2._qp*cutlsj(3,2))
      if (real(expo).gt.real(rrr)) expo=rrr
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      expexp=exp(-expo)
      call regulator_poly(xmev,ymev,expexp,cutlsj(4,2)
     &                     ,2._qp*cutlsj(3,2)) ! call regulator function that stores the factor in expexp
      end if
      vv4(3)=vv4(3)*expexp
!
      go to 8400
!
!
!        j=1
!        ---
!        ---
!
 8320 continue
!
         !test
         !vv0(:)=-1e-6_qp
         !vv2(:)=-1e-6_qp
         !vv4(:)=-1e-6_qp
!        1p1
!        ---
!        Q^2 term
!      

      if (cutlsj(1,3).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(2,3))**(2._qp*cutlsj(1,3))
     &    +(ymev/cutlsj(2,3))**(2._qp*cutlsj(1,3))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(2,3)
     &                     ,2._qp*cutlsj(1,3)) ! call regulator function that stores the factor in expexp
      end if
      vv2(1)=vv2(1)*expexp
      vv0(1)=vv0(1)*expexp
!
!        Q^4 term
!
      if (cutlsj(3,3).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(4,3))**(2._qp*cutlsj(3,3))
     &    +(ymev/cutlsj(4,3))**(2._qp*cutlsj(3,3))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(4,3)
     &                     ,2._qp*cutlsj(3,3)) ! call regulator function that stores the factor in expexp
      end if
      vv4(1)=vv4(1)*expexp
!
!        3p1
!        ---
!        Q^2 term
!
      if (cutlsj(1,4).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(2,4))**(2._qp*cutlsj(1,4))
     &    +(ymev/cutlsj(2,4))**(2._qp*cutlsj(1,4))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(2,4)
     &                     ,2._qp*cutlsj(1,4)) ! call regulator function that stores the factor in expexp
      end if
      vv2(2)=vv2(2)*expexp
      vv0(2)=vv0(2)*expexp
!
!        Q^4 term
!
      if (cutlsj(3,4).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(4,4))**(2._qp*cutlsj(3,4))
     &    +(ymev/cutlsj(4,4))**(2._qp*cutlsj(3,4))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(4,4)
     &                     ,2._qp*cutlsj(3,4)) ! call regulator function that stores the factor in expexp
      end if
      vv4(2)=vv4(2)*expexp
!
!        3s1
!        ---
!        Q^0 term
!
      if (cutlsj(1,5).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(2,5))**(2._qp*cutlsj(1,5))
     &    +(ymev/cutlsj(2,5))**(2._qp*cutlsj(1,5))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(2,5)
     &                     ,2._qp*cutlsj(1,5)) ! call regulator function that stores the factor in expexp
      end if
      vv0(4)=vv0(4)*expexp
!
!        Q^2 terms
!
      if (cutlsj(3,5).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(4,5))**(2._qp*cutlsj(3,5))
     &    +(ymev/cutlsj(4,5))**(2._qp*cutlsj(3,5))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(4,5)
     &                     ,2._qp*cutlsj(3,5)) ! call regulator function that stores the factor in expexp
      end if
      vv2(4)=vv2(4)*expexp
!
!        Q^4 terms
!
      if (cutlsj(5,5).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(6,5))**(2._qp*cutlsj(5,5))
     &    +(ymev/cutlsj(6,5))**(2._qp*cutlsj(5,5))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(6,5)
     &                     ,2._qp*cutlsj(5,5)) ! call regulator function that stores the factor in expexp
      end if
      vv4(4)=vv4(4)*expexp
!
!        3d1
!        ---
!        Q^4 term
!
      if (cutlsj(1,6).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(2,6))**(2._qp*cutlsj(1,6))
     &    +(ymev/cutlsj(2,6))**(2._qp*cutlsj(1,6))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(2,6)
     &                     ,2._qp*cutlsj(1,6)) ! call regulator function that stores the factor in expexp
      end if
      vv4(3)=vv4(3)*expexp
      vv2(3)=vv2(3)*expexp
      vv0(3)=vv0(3)*expexp
!
!        3s/d1
!        -----
!        Q^2 term
!
      if (cutlsj(1,7).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(2,7))**(2._qp*cutlsj(1,7))
     &    +(ymev/cutlsj(2,7))**(2._qp*cutlsj(1,7))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(2,7)
     &                     ,2._qp*cutlsj(1,7)) ! call regulator function that stores the factor in expexp
      end if
      vv2(5)=vv2(5)*expexp
      vv2(6)=vv2(6)*expexp
      vv0(5)=vv0(5)*expexp
      vv0(6)=vv0(6)*expexp
!
!        Q^4 term
!
      if (cutlsj(3,7).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(4,7))**(2._qp*cutlsj(3,7))
     &    +(ymev/cutlsj(4,7))**(2._qp*cutlsj(3,7))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(4,7)
     &                     ,2._qp*cutlsj(3,7)) ! call regulator function that stores the factor in expexp
      end if
      vv4(5)=vv4(5)*expexp
      vv4(6)=vv4(6)*expexp
!
      go to 8400
!
!
!        j=2
!        ---
!        ---
!
 8330 continue
!
!        1d2
!        ---
!        Q^4 term
!
      if (cutlsj(1,8).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(2,8))**(2._qp*cutlsj(1,8))
     &    +(ymev/cutlsj(2,8))**(2._qp*cutlsj(1,8))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(2,8)
     &                     ,2._qp*cutlsj(1,8)) ! call regulator function that stores the factor in expexp
      end if
      vv4(1)=vv4(1)*expexp
      vv2(1)=vv2(1)*expexp
      vv0(1)=vv0(1)*expexp
!
!        3d2
!        ---
!        Q^4 term
!
      if (cutlsj(1,9).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(2,9))**(2._qp*cutlsj(1,9))
     &    +(ymev/cutlsj(2,9))**(2._qp*cutlsj(1,9))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(2,9)
     &                     ,2._qp*cutlsj(1,9)) ! call regulator function that stores the factor in expexp
      end if
      vv4(2)=vv4(2)*expexp
      vv2(2)=vv2(2)*expexp
      vv0(2)=vv0(2)*expexp
!
!        3p2
!        ---
!        Q^2 term
!
      if (cutlsj(1,10).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(2,10))**(2._qp*cutlsj(1,10))
     &    +(ymev/cutlsj(2,10))**(2._qp*cutlsj(1,10))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(2,10)
     &                     ,2._qp*cutlsj(1,10)) ! call regulator function that stores the factor in expexp
      end if
      vv2(4)=vv2(4)*expexp
      vv0(4)=vv0(4)*expexp
      vv2(3)=vv2(3)*expexp
      vv0(3)=vv0(3)*expexp
!
!        Q^4 terms
!
      if (cutlsj(3,10).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(4,10))**(2._qp*cutlsj(3,10))
     &    +(ymev/cutlsj(4,10))**(2._qp*cutlsj(3,10))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(4,10)
     &                     ,2._qp*cutlsj(3,10)) ! call regulator function that stores the factor in expexp
      end if
      vv4(4)=vv4(4)*expexp
      vv4(3)=vv4(3)*expexp
!
!        3p/f2
!        -----
!        Q^4 term
!
      if (cutlsj(1,12).eq.0._qp) then
      expexp=1._qp
      else
      expo=(xmev/cutlsj(2,12))**(2._qp*cutlsj(1,12))
     &    +(ymev/cutlsj(2,12))**(2._qp*cutlsj(1,12))
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
      call regulator_poly(xmev,ymev,expexp,cutlsj(2,12)
     &                     ,2._qp*cutlsj(1,12)) ! call regulator function that stores the factor in expexp
      end if
      vv4(5)=vv4(5)*expexp
      vv4(6)=vv4(6)*expexp
      vv2(5)=vv2(5)*expexp
      vv2(6)=vv2(6)*expexp
      vv0(5)=vv0(5)*expexp
      vv0(6)=vv0(6)*expexp
!
      go to 8400
!
!
!        j=3
!        ---
!        ---
!
 8340 continue
!
!        3d3
!        ---
!        special cutoff: for the exponent, the parameter
!        from the list is not used. Instead it is used
!        what you see below.
!        Q^4 term
!
      if (cutlsj(1,15).eq.0._qp) then
      expexp=1._qp
      else
      !test
      expo1=(xmev/cutlsj(2,15))**(2._qp*2.0_qp)
     &    +(ymev/cutlsj(2,15))**(2._qp*2.0_qp)
!      expo1=(xmev/cutlsj(2,15))**(2._qp*1.0_qp)
!     &    +(ymev/cutlsj(2,15))**(2._qp*1.0_qp)
      if (real(expo1).gt.real(rrr)) expo1=rrr
      expexp1=exp(-expo1)
!      expexp1=(exp(-(reg_n-1._qp)*expo1)
!     &+exp(+(reg_n-1._qp)*expo1))
!     &/(exp(-(reg_n)*expo1)+exp(+(reg_n)*expo1))
      expo2=(xmev/cutlsj(2,15))**(2._qp*3.0_qp)
     &    +(ymev/cutlsj(2,15))**(2._qp*3.0_qp)
      if (real(expo2).gt.real(rrr)) expo2=rrr
      expexp2=exp(-expo2)
!      expexp2=(exp(-(reg_n-1._qp)*expo2)
!     &+exp(+(reg_n-1._qp)*expo2))
!     &/(exp(-(reg_n)*expo2)+exp(+(reg_n)*expo2))
      call regulator_poly(xmev,ymev,expexp1,cutlsj(2,15)
     &                     ,2._qp*complex(2._qp,0._qp)) ! call regulator function that stores the factor in expexp
      call regulator_poly(xmev,ymev,expexp2,cutlsj(2,15)
     &                     ,2._qp*complex(3._qp,0._qp)) ! call regulator function that stores the factor in expexp
      expexp=0.5_qp*(expexp1+expexp2)
      end if
!
!        use 3d3 cutoff for all j.eq.3 partial waves
!
      do 8345 iv=1,6
      vv0(iv)=vv0(iv)*expexp
      vv2(iv)=vv2(iv)*expexp
 8345 vv4(iv)=vv4(iv)*expexp
!
!
!
!
!
!
!        final add up
!        ------------
!
 8400 do 8405 iv=1,6
      !test
      !v(iv)=0._qp
 8405 v(iv)=v(iv)+vv0(iv)+vv2(iv)+vv4(iv)
!

      end if
!
!
!
!
 8500 if (j.eq.0.or..not.heform) go to 8900
!
!
!         translation into (combinations of) helicity states
!
!
      do 8505 i=1,4
 8505 vl(i)=v(i+2)
!
      do 8520 ii=1,4
      iii=ii+2
      v(iii)=0._qp
!
      do 8515 i=1,4
 8515 v(iii)=v(iii)+adminv(ii,i)*vl(i)
 8520 v(iii)=v(iii)*a2j1
!
!
!
!
 8900 continue

      return
      end
      subroutine chipar_qp
!
!        chipar reads, writes, and stores the parameter for all
!        chi-subroutines.
!
!
      use iso_fortran_env
      implicit complex*16 (a-h,o-z)
      integer,parameter:: qp=real64
!
!
      common /crdwrt_qp/ kread,kwrite,kpunch,kda(9)
!
      common /cstate_qp/ j,heform,sing,trip,coup,endep,label
      common /cnn_qp/ inn
      logical heform,sing,trip,coup,endep
!
!
!        common block for all chi-subroutines
!
      common /cchi_qp/ vj(32,270),c(20,270),fff,ff,f(52),aa(200)
     &,ai(19,30),
     &                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     &                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(200)
     &,wt(200),
     &                ic(20,270),ift(3),mint(3),maxt(3),nt,
     &                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     &                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     &                indc(2,270),indpar(3),indxy
!
!         specifications for this common block
!
      logical indc,indxy,indpar
!
      common /compar_qp/ cb1a(3),cb2a(3),cb3a(3),cb4a(3),
     &                cd12a(3),cd3a(3),cd5a(3),cd145a(3)
!
      common /comlsj_qp/ clsj(15,50),cutlsj(15,50),indlsj
      logical indlsj
!
!
!        further specifications
!
      dimension cc(5),cca(5)
      dimension clec(15,50)
      dimension a(1024),b(32)
      dimension ttab(5,131),tab(5,131)
      dimension topepp(5,2),topenp(5,2),topenn(5,2)
      dimension t1s0pp(5),t1s0np(5),t1s0nn(5)
      real*16 eps
      character*4 name(3)
      character*4 ntab(3,131)
      integer imga(3)
      character*4 cut,cuta,fun,lsj,lec,end
      character*4 mesong(40)
      logical index
      logical zerocp
      logical indlec
      logical indca,indlca
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     &            '1-  ','1-t ','1-tt','1-st','1-ss',
     &            'c   ','ss  ','ls  ','sq  ','sk  ',
     &            'sl  ',
     &         24*'    '/
      data index/.false./
      data zerocp/.true./
!      data pi/3.141592653589793_qp/
      data pi/3.141592653589793238462643383279502884197_qp/
      data eps/1.e-32/
      data cut/'cut '/,cuta/'cuta'/
      data fun/'fun '/,lsj/'lsj '/,lec/'lec '/,end/'end '/
!
!
!
!
!        parameter tables
!        ----------------
!        ----------------
!
!
!        identification table
!        --------------------
!
!
      data ntab/
     & 'cuta','ll  ','  ',
     & 'sq  ',' ope','p ',
     & '    ','    ','  ',
     & '    ','    ','  ',
     & '    ','    ','  ',
     & 'sq  ',' ope','p ',
     & 'sq  ',' pi-','g ',
     & 'fun ','    ','  ',
     & '    ','    ','  ',
     & '    ','    ','  ',
     & 'cuta','ll  ','  ',
     & 'sq  ',' tpn','1 ',
     & 'fun ','    ','  ',
     & 'ss  ',' tpn','1 ',
     & 'fun ','    ','  ',
     & 'fun ','    ','  ',
     & 'c   ',' tpn','2 ',
     & 'fun ','    ','  ',
     & 'sq  ',' tpn','2 ',
     & 'fun ','    ','  ',
     & 'ss  ',' tpn','2 ',
     & 'fun ','    ','  ',
     & 'fun ','    ','  ',
     & 'ls  ',' tpn','2 ',
     & 'fun ','    ','  ',
     & 'c   ',' tpn','3 ',
     & 'fun ','    ','  ',
     & 'c   ',' tpn','3 ',
     & 'fun ','    ','  ',
     & 'ls  ',' tpn','3 ',
     & 'fun ','    ','  ',
     & 'c   ',' tpn','32',
     & 'fun ','    ','  ',
     & 'ss  ',' tpn','32',
     & 'fun ','    ','  ',
     & 'fun ','    ','  ',
     & 'sq  ',' tpn','32',
     & 'fun ','    ','  ',
     & 'c   ',' tpn','3m',
     & 'fun ','    ','  ',
     & 'ss  ',' tpn','3m',
     & 'fun ','    ','  ',
     & 'fun ','    ','  ',
     & 'sq  ',' tpn','3m',
     & 'fun ','    ','  ',
     & 'ls  ',' tpn','3m',
     & 'fun ','    ','  ',
     & 'sl  ',' tpn','3m',
     & 'fun ','    ','  ',
     & 'c   ',' tpn','2c',
     & 'fun ','    ','  ',
     & 'sq  ',' tpn','2c',
     & 'fun ','    ','  ',
     & 'ss  ',' tpn','2c',
     & 'fun ','    ','  ',
     & 'fun ','    ','  ',
     & 'c   ',' tpn','1 ',
     & 'fun ','    ','  ',
     & 'c   ',' tpn','2 ',
     & 'fun ','    ','  ',
     & 'sq  ',' tpn','2 ',
     & 'fun ','    ','  ',
     & 'ss  ',' tpn','2 ',
     & 'fun ','    ','  ',
     & 'fun ','    ','  ',
     & 'ls  ',' tpn','2 ',
     & 'fun ','    ','  ',
     & 'sq  ',' tpn','3 ',
     & 'fun ','    ','  ',
     & 'ss  ',' tpn','3 ',
     & 'fun ','    ','  ',
     & 'fun ','    ','  ',
     & 'c   ',' tpn','3 ',
     & 'fun ','    ','  ',
     & 'fun ','    ','  ',
     & 'ss  ',' tpn','3 ',
     & 'fun ','    ','  ',
     & 'fun ','    ','  ',
     & 'sq  ',' tpn','3 ',
     & 'fun ','    ','  ',
     & 'ls  ',' tpn','3 ',
     & 'fun ','    ','  ',
     & 'ss  ',' tpn','32',
     & 'fun ','    ','  ',
     & 'fun ','    ','  ',
     & 'sq  ',' tpn','32',
     & 'fun ','    ','  ',
     & 'c   ',' tpn','32',
     & 'fun ','    ','  ',
     & 'c   ',' tpn','3m',
     & 'fun ','    ','  ',
     & 'sq  ',' tpn','3m',
     & 'fun ','    ','  ',
     & 'ss  ',' tpn','3m',
     & 'fun ','    ','  ',
     & 'fun ','    ','  ',
     & 'ls  ',' tpn','3m',
     & 'fun ','    ','  ',
     & 'c   ',' tpn','2c',
     & 'fun ','    ','  ',
     & 'sq  ',' tpn','2c',
     & 'fun ','    ','  ',
     & 'ss  ',' tpn','2c',
     & 'fun ','    ','  ',
     & 'fun ','    ','  ',
     & 'cuta','ll  ','  ',
     & 'lsj ',' 1S0','  ',
     & 'lsj ',' 1S0','  ',
     & 'lsj ',' 1S0','  ',
     & 'lsj ',' 1S0','  ',
     & 'lsj ',' 3P0','  ',
     & 'lsj ',' 3P0','  ',
     & 'lsj ',' 1P1','  ',
     & 'lsj ',' 1P1','  ',
     & 'lsj ',' 3P1','  ',
     & 'lsj ',' 3P1','  ',
     & 'lsj ',' 3S1','  ',
     & 'lsj ',' 3S1','  ',
     & 'lsj ',' 3S1','  ',
     & 'lsj ',' 3S1','  ',
     & 'lsj ',' 3D1','  ',
     & 'lsj ',' 3S-','D1',
     & 'lsj ',' 3S-','D1',
     & 'lsj ',' 3S-','D1',
     & 'lsj ',' 1D2','  ',
     & 'lsj ',' 3D2','  ',
     & 'lsj ',' 3P2','  ',
     & 'lsj ',' 3P2','  ',
     & 'lsj ',' 3P-','F2',
     & 'lsj ',' 3D3','  ',
     & 'end ','para','m.'/
!
!
!        parameters
!        ----------
!
!
      data tab/
     &    6.000000_qp,   0.0_qp,    4.0000_qp,  500.0_qp,   0.0_qp,
     &   -1.290000_qp,  92.4_qp,  134.9766_qp,    0.0_qp,   0.0_qp,
     &    0.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    0.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    0.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.290000_qp,  92.4_qp,  139.5702_qp,    0.0_qp,   0.0_qp,
     &   -0.062170_qp,  92.4_qp,  139.5702_qp,    0.0_qp,   0.0_qp,
     &   10.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    0.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    0.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    6.000000_qp,   0.0_qp,    2.0000_qp,  500.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   14.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   11.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   14.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   15.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   17.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   11.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   17.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   19.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   30.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   39.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   38.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   40.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   42.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   11.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   42.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   48.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   50.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   11.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   50.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   53.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   54.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.500000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   55.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.500000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   56.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.500000_qp,   0.0_qp,  138.0390_qp,    0.0_qp,  -1.0_qp,
     &   11.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   56.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   13.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   16.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   18.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   11.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   18.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   20.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   31.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   11.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   31.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   36.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   11.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   37.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   11.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   37.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    4.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   36.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   41.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   11.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   41.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   43.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   49.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   51.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   51.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   11.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   52.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   55.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   56.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.000000_qp,   0.0_qp,  138.0390_qp,    1.0_qp,  -1.0_qp,
     &   11.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &   56.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp,
     &    0.000000_qp,   0.0_qp,    2.0000_qp,  500.0_qp,   0.0_qp,
     &   -0.147167_qp,   3.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    2.380000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &   -2.545000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &  -16.000000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    1.487000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    0.245000_qp,   3.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    0.656000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    5.250000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &   -0.630000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    2.350000_qp,   4.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &   -0.118972496_qp,3.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    0.760000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    7.000000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    6.550000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &   -2.800000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    0.826000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    2.250000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    6.610000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.770000_qp,   4.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &   -1.460000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &   -0.538000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    2.295000_qp,   2.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &   -0.465000_qp,   4.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    5.660000_qp,   2.5_qp,  500.0000_qp,    0.0_qp,   0.0_qp,
     &    0.000000_qp,   0.0_qp,    0.0000_qp,    0.0_qp,   0.0_qp/
!
!
      data topepp/
     &   -1.290000_qp,  92.4_qp,  134.9766_qp,    0.0_qp,   0.0_qp,
     &    0.000000_qp,   0.0_qp,  139.5702_qp,    0.0_qp,   0.0_qp/
!
!
      data topenp/
     &   -1.290000_qp,  92.4_qp,  139.5702_qp,    0.0_qp,   0.0_qp,
     &   -0.062170_qp,  92.4_qp,  139.5702_qp,    0.0_qp,   0.0_qp/
!
!
      data topenn/
     &   -1.290000_qp,  92.4_qp,  134.9766_qp,    0.0_qp,   0.0_qp,
     &    0.000000_qp,   0.0_qp,  139.5702_qp,    0.0_qp,   0.0_qp/
!
!
      data t1s0pp/
     &   -0.145286_qp,   3.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp/
!
!
      data t1s0np/
     &   -0.147167_qp,   3.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp/
!
!
      data t1s0nn/
     &   -0.146285_qp,   3.0_qp,  500.0000_qp,    0.0_qp,   0.0_qp/
!
!
!        this has been the end of tables
!
!
      save
!
!
!
!
!
10004 format (1h ,2a4,a2,f12.6,f10.6,1x,f10.5,2(f7.1,3x))
10005 format (1h ,2a4,a2,f4.1,1x,2f10.5,f13.5,f10.5)
10007 format (1h ,2a4,a2,3i3)
10008 format (1h ,63(1h-))
10010 format (1h ,2a4,a2,i3,2f10.2)
10011 format (//' n3lo: Charge-Dependent Chiral NN Potential ',
     & ' at Order Four (N3LO), with complex scaling')
10020 format (1h ,2a4,a2,f12.6,4f9.2)
10021 format (1h ,2a4,a2,5f10.6)
!
!
!
!
      if (index) go to 50
      index=.true.
!
      x=-1._qp
      y=-1._qp
!
!
!
!
!        maxima of certain indices related to the dimension as follows:
!        dimension c(mme,imee),ic(mice,imee),indc(mindce,imee),
!                  mgg(mge,3),mggo(mge,3),mesong(mge),vj(32,imee),
!                  ima(mee,mge,3)
!
      mge=40
      mee=30
      mme=20
      mice=20
      mindce=2
      imb=1
      ime=0
      imee=270
!        mme always ge mice, mindce
!
!        set all parameters and indices to zero or .false.
!
      do 1 int=1,3
      imga(int)=0
      indpar(int)=.false.
      do 1 mgx=1,mge
      mgg(mgx,int)=0
    1 mggo(mgx,int)=0
!
!
      do 2 il=1,imee
      do 2 mm=1,mme
      if (mm.le.mindce) indc(mm,il)=.false.
      if (mm.le.mice) ic(mm,il)=0
    2 c(mm,il)=0._qp
      endep=.false.
!
!
      pi2=pi*pi
      pi4=pi2*pi2
!
!
!
!
!        start
!        -----
!        -----
!
!
!
!        write title
!
   50 continue
      write (kwrite,10011)
      write (kwrite,10008)
      write (kwrite,10008)
!
!
!        store systematically the parameter sets
!        for pp, np, nn.
!
!
      do 8999 inter=1,3
!
!
      indca=.false.
      indlca=.false.
      indlsj=.false.
      indlec=.false.
      ilsj=0
      ilec=0
      do 55 ii=1,50
      do 55 i=1,15
      clsj(i,ii)=0._qp
      cutlsj(i,ii)=0._qp
   55 clec(i,ii)=0._qp
!
!
!        fix index-parameters concerning the factor
!        and the cutoff for the potential as a whole
!
      ift(inter)=1
      ezz1(inter)=0._qp
      ezz2(inter)=0._qp
!
!**** write (kwrite,10010) name,ift(inter),ezz1(inter),ezz2(inter)
      iftyp=ift(inter)
      if (iftyp.lt.0.or.iftyp.gt.6) go to 9003
!
!
!        fix parameters for numerical integration
!
      mint(inter)=4
      maxt(inter)=48
!
!**** write (kwrite,10007) name,mint(inter),maxt(inter)
!
!        nucleon mass
!
      go to (51,52,53), inter
!        mass used for pp
   51 wn=938.272_qp
      go to 54
!        mass used for np
   52 wn=938.9182_qp
      go to 54
!        mass used for nn
   53 wn=939.5653_qp
   54 continue
!**** write (kwrite,10004) name,wn
      wnq=wn*wn
      dwn=1._qp/wn
      dwnq=dwn*dwn
      wnn(inter)=wn
!
!
!        ga and fpi
!
      ga=1.29_qp
      fpi=92.4_qp
!
!**** write (kwrite,10004) name,ga,fpi
      ga2=ga*ga
      ga4=ga2*ga2
      ga6=ga4*ga2
      fpi=fpi*dwn
      fpi2=fpi*fpi
      fpi4=fpi2*fpi2
      fpi6=fpi4*fpi2
      gaa(inter)=ga2
      fpia(inter)=fpi2
!
!        fix the LECs of the pi-N Lagrangian
!
!        the c_i LECs
      cc(1)=-0.81_qp
      cc(2)=2.8_qp
      cc(3)=-3.2_qp
      cc(4)=5.4_qp
!**** write (kwrite,10021) name,cc
      cb1a(inter)=cc(1)*wn*1e-3_qp
      cb2a(inter)=cc(2)*wn*1e-3_qp
      cb3a(inter)=cc(3)*wn*1e-3_qp
      cb4a(inter)=cc(4)*wn*1e-3_qp
!
!        the d_i LECs
      cc(1)=3.06_qp
      cc(2)=-3.27_qp
      cc(3)=0.45_qp
      cc(4)=-5.65_qp
!**** write (kwrite,10021) name,cc
      cd12a(inter)=cc(1)*wnq*1e-6_qp
      cd3a(inter)=cc(2)*wnq*1e-6_qp
      cd5a(inter)=cc(3)*wnq*1e-6_qp
      cd145a(inter)=cc(4)*wnq*1e-6_qp
!
      cb1=cb1a(inter)
      cb2=cb2a(inter)
      cb3=cb3a(inter)
      cb4=cb4a(inter)
      cd12=cd12a(inter)
      cd3=cd3a(inter)
      cd5=cd5a(inter)
      cd145=cd145a(inter)
!
!
!
!        prepare table
!
      do 56 ll=1,131
      do 56 i=1,5
   56 ttab(i,ll)=tab(i,ll)
!
!
!        charge-dependent modifications for pp
!
      if (inter.eq.1) then
      do 57 i=1,5
      ttab(i,6)=topepp(i,1)
      ttab(i,7)=topepp(i,2)
   57 ttab(i,107)=t1s0pp(i)
      end if
!
!
!        charge-dependent modifications for np
!
      if (inter.eq.2) then
      do 58 i=1,5
      ttab(i,6)=topenp(i,1)
      ttab(i,7)=topenp(i,2)
   58 ttab(i,107)=t1s0np(i)
      end if
!
!
!        charge-dependent modifications for nn
!
      if (inter.eq.3) then
      do 59 i=1,5
      ttab(i,6)=topenn(i,1)
      ttab(i,7)=topenn(i,2)
   59 ttab(i,107)=t1s0nn(i)
      end if
!
!
!
!
!        get parameters from tables, line by line
!        ----------------------------------------
!        ----------------------------------------
!
!
!
      line=0
!
   61 line=line+1
      do i=1,5
      if (i.le.3) then
      name(i)=ntab(i,line)
      end if
      cc(i)=ttab(i,line)
      end do
!
!        check if end of input
!
      if (name(1).eq.end) go to 7000
!
!        check if lsj or lec
!
      if (name(1).eq.lsj) go to 6000
      if (name(1).eq.lec) go to 6500
!
!        check if data-card just read contains cut-off or
!        function parameters
!
      if (name(1).eq.cut.or.name(1).eq.fun) go to 70
!
      if (name(1).eq.cuta) then
!**** write (kwrite,10005) name,cc
      indca=.true.
      do i=1,5
      cca(i)=cc(i)
      end do
      go to 61
      end if

!
!
!
!
!        write parameters which are no cut-off or function parameters
!        ------------------------------------------------------------
!
!
!
!
!**** write (kwrite,10004) name,cc
!
!        check if coupling constant is zero
!
!****    do not use zerocp anymore
!****    because the first eight input lines are always pions.
!****    these lines must never be skipped even when g_pi zero.
!**** if (cc(1).ne.0.d0) go to 62
!**** zerocp=.true.
!**** go to 61
!
   62 zerocp=.false.
!
!        find out number of contribution mg
!
      do 63 mg=1,mge
      if (name(1).eq.mesong(mg)) go to 64
   63 continue
      go to 9000
!
!
!
!
!        store parameters which are no cut-off or function parameters
!        ------------------------------------------------------------
!
!
!
!
   64 ime=ime+1
      if (ime.gt.imee) go to 9011
      mgg(mg,inter)=mgg(mg,inter)+1
      m=mgg(mg,inter)
      if (m.gt.mee) go to 9001
      ima(m,mg,inter)=ime
      if (m.ne.1) go to 65
      imga(inter)=imga(inter)+1
      mggo(imga(inter),inter)=mg
   65 continue
!
!
      c(1,ime)=cc(1)
!
!
      if (mg.le.10) then
      c(1,ime)=c(1,ime)*4._qp*pi
      end if

      if (mg.le.3.and.cc(2).ne.0._qp) then
      c(1,ime)=(cc(1)/cc(2)*wn)**2
      if (real(cc(1)).lt.0._qp) c(1,ime)=-c(1,ime)
      end if
!
!
      if (mg.ge.6.and.mg.le.10) then
!        store coupling constant f*g
      c(3,ime)=cc(2)*c(1,ime)
!        store coupling constant f**2
      c(2,ime)=cc(2)*c(3,ime)
      if (mg.eq.10)
     &  c(1,ime)=c(1,ime)+c(3,ime)*2._qp+c(2,ime)
      end if
!
!
      if (mg.ge.11.and.cc(2).ne.0._qp) then
      c(1,ime)=(cc(1)/(2._qp*cc(2))*wn)**2
      if (real(cc(1)).lt.0._qp) c(1,ime)=-c(1,ime)
      end if
!
!
!        store meson mass square in units of nucleon mass square
      c(4,ime)=cc(3)*cc(3)*dwnq
!
!        test iso-spin
      icc=cc(4)
      if (icc.ne.0.and.icc.ne.1) go to 9004
!         store isospin as logical constant
      if (icc.eq.1) indc(1,ime)=.true.
!        store and test iprsp
      icc=cc(5)
      ic(1,ime)=icc
      if (iabs(ic(1,ime)).gt.1) go to 9005
!
!        index values for further storing
      mi=4
      mm=5
!
!
!        check if there is a `cutall' cutoff
!
      if (indca) then
      name(1)=cut
      do i=1,5
      cc(i)=cca(i)
      end do
      go to 72
      else
      go to 61
      end if
!
!
!
!
!        write cut-off or function parameters
!        ------------------------------------
!
!
!
!
   70 continue
!**** write (kwrite,10005) name,cc
!
      if (zerocp) go to 61
!
   72 continue
!
!
!
!
!        store parameters
!        ----------------
!
!
!
      ityp=cc(1)
!
      if (ityp.eq.0) go to 5995
      if (ityp.lt.1.or.ityp.gt.56) go to 9002
!
      im=ime
!
!        store typ of cut-off or function
      ic(mi,im)=ityp
!
      if (ityp.le.10) then
!        store and test typ of propagator of cut-off
      ic(mi+1,im)=cc(2)
      if (ic(mi+1,im).lt.0.or.ic(mi+1,im).gt.1) go to 9006
      end if
!
      go to (100,100,300,9002,500,600,9002,9002,9002,1000,
     & 1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,
     & 2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,
     & 3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,
     & 4100,4200,4300,9002,9002,9002,9002,4800,4900,5000,
     & 5100,5200,5300,5400,5500,5600),ityp
!
!
!
!
!        cut-off of dipole type
!        **********************
!
!
!        store and test exponent of cut-off
  100 ic(mi+2,im)=cc(3)
      if (ic(mi+2,im).lt.0) go to 9009
      if (ic(mi+2,im).gt.0) go to 101
!        exponent is zero, omit cut-off
      ic(mi,im)=0
      ic(mi+1,im)=0
      go to 5995
!        store cut-off mass for denominator
  101 c(mm+1,im)=cc(4)*cc(4)*dwnq
!        store numerator of cut-off
      c(mm,im)=c(mm+1,im)
      if (ityp.eq.2)     c(mm,im)=c(mm,im)-c(4,im)
      mi=mi+3
      mm=mm+2
      go to 5995
!
!
!
!
!        exponential form factor of momentum transfer
!        ********************************************
!
!
!        check exponent
  300 if (real(cc(3)).lt.0._qp) go to 9009
      if (real(cc(3)).gt.0._qp) go to 301
!        exponent is zero, omit cutoff
      ic (mi,im)=0
      ic (mi+1,im)=0
      go to 5995
!        store exponent
  301 c(mm+1,im)=cc(3)
!        compute constant factor for argument of exponential function
      c(mm,im)=wnq/(cc(4)*cc(4))
      mi=mi+2
      mm=mm+2
      go to 5995
!
!
!
!
!        sharp cutoff in x and y
!        ***********************
!
!
  500 c(mm,im)=cc(4)*dwn
      mi=mi+2
      mm=mm+1
      go to 5995
!
!
!
!
!        exponential form factor of xx and yy
!        ************************************
!
!
!        check exponent
  600 if (real(cc(3)).lt.0._qp) go to 9009
      if (real(cc(3)).gt.0._qp) go to 601
!        exponent is zero, omit cutoff
      ic (mi,im)=0
      ic (mi+1,im)=0
      go to 5995
!        store exponent
  601 c(mm+1,im)=cc(3)
!        compute constant factor for argument of exponential function
      c(mm,im)=wnq/(cc(4)*cc(4))
      mi=mi+2
      mm=mm+2
      go to 5995
!
!
!
!
!        pi-gamma potential
!        ******************
!
!
 1000 c(mm,im)=cc(3)
      mi=mi+2
      mm=mm+1
      go to 5995
!
!
!
!
!        function q^2 (momentum-transfer squared)
!        ************
!
!
 1100 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        function k^2 (average-momentum squared)
!        ************
!
!
 1200 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        function 1 for tpn1 (=NLO)
!        **************************
!
!
 1300 c(mm,im)=-1._qp/(384._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        function 2 for tpn1
!        *******************
!
!
 1400 c(mm,im)=-3._qp*ga4/(64._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn2 (=N^2LO), function 1
!        *************************
!
!
 1500 c(mm,im)=-3._qp*ga2/(16._qp*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn2, function 2
!        ****************
!
!
 1600 c(mm,im)=-ga2/(128._qp*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn2, function 3
!        ****************
!
!
 1700 c(mm,im)=9._qp*ga4/(512._qp*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn2, function 4
!        ****************
!
!
 1800 c(mm,im)=-ga2/(32._qp*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!
!        tpn2, function 5
!        ****************
!
!
 1900 c(mm,im)=6._qp*ga4/(64._qp*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!        tpn2, function 6
!        ****************
!
!
 2000 c(mm,im)=2._qp*ga2*(1._qp-ga2)/(64._qp*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        function q^4 (momentum-transfer to the power of 4)
!        ************
!
 2100 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!        function k^4 (average-momentum to the power of 4)
!        ************
!
 2200 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!        function +q^2*k^2
!        *****************
!
 2300 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!        function  (\vec q x \vec k)^2
!        *****************************
!
 2400 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!        function xy
!        ***********
!
 2500 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!        function xx+yy
!        **************
!
 2600 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!        function xx*xx+yy*yy
!        ********************
!
 2700 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!        function xx
!        ***********
!
 2800 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!        function yy
!        ***********
!
 2900 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn3 (= N^3LO with one loop), function 1
!        ****************************************
!
!
 3000 c(mm,im)=3._qp/(16._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn3, function 2
!        ****************
!
!
 3100 continue
      c(mm,im)=cb4*cb4/(96._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!        function 1._qp
!        *************
!
 3200 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!        function 1-q^2/8-k^2/2
!        *************************
!
 3300 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!        function 1-q^2/8
!        ****************
!
 3400 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!        function 1+k^2/2
!        ****************
!
 3500 continue
      c(mm,im)=cc(2)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn3, function 3
!        ****************
!
!
 3600 c(mm,im)=-cb4/(192._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn3, function 4
!        ****************
!
!
 3700 c(mm,im)=cb4/(192._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn3, function 5
!        ****************
!
!
 3800 continue
      c(mm,im)=cb2*ga2/(8._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn3, function 6
!        ****************
!
!
 3900 c(mm,im)=-ga2/(32._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn32 (=N^3LO with 2 loops), function 1
!        ***************************************
!
!
 4000 c(mm,im)=3._qp*ga4/(1024._qp*pi2*fpi6)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn32, function 2
!        *****************
!
!
 4100 c(mm,im)=-ga4/(2048._qp*pi2*fpi6)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn32, function 3
!        *****************
!
!
 4200 c(mm,im)=ga2*cd145/(32._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn32, function 4
!        *****************
!
!
 4300 c(mm,im)=1._qp/(18432._qp*pi4*fpi6)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn3m (= N^3LO, 1/M^2 terms), function 1
!        ****************************************
!
!
 4800 c(mm,im)=-ga4/(32._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn3m, function 2
!        *****************
!
!
 4900 c(mm,im)=-1._qp/(768._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn3m, function 3
!        *****************
!
!
 5000 c(mm,im)=ga4/(32._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn3m, function 4
!        *****************
!
!
 5100 c(mm,im)=1._qp/(1536._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn3m, function 5
!        *****************
!
!
 5200 c(mm,im)=1._qp/(256._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn3m, function 6
!        *****************
!
!
 5300 c(mm,im)=ga4/(4._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn3m, function 7
!        *****************
!
!
 5400 c(mm,im)=ga4/(32._qp*pi2*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn2c (correction for our it 2pi), function 1
!        *********************************************
!
!
 5500 c(mm,im)=ga4/(128._qp*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        tpn2c (correction for our it 2pi), function 2
!        *********************************************
!
!
 5600 c(mm,im)=-ga4/(256._qp*pi*fpi4)
      mi=mi+1
      mm=mm+1
      go to 5995
!
!
!
!
!        end cut-offs and functions
!        **************************
!
!        test dimensions
 5995 if (mi.gt.mice.or.mm.gt.mme) go to 9010
!
      if (indlca) go to 7800
!
      go to 61
!
!
!
!
!        partial wave LEC's
!        ------------------
!
!
 6000 continue
!**** write (kwrite,10020) name,cc
      indlsj=.true.
      ilsj=ilsj+1
      if (ilsj.le.4) iilsj=1
      if (ilsj.ge.5.and.ilsj.le.6) iilsj=2
      if (ilsj.ge.7.and.ilsj.le.8) iilsj=3
      if (ilsj.ge.9.and.ilsj.le.10) iilsj=4
      if (ilsj.ge.11.and.ilsj.le.14) iilsj=5
      if (ilsj.ge.15.and.ilsj.le.15) iilsj=6
      if (ilsj.ge.16.and.ilsj.le.18) iilsj=7
      if (ilsj.ge.19.and.ilsj.le.19) iilsj=8
      if (ilsj.ge.20.and.ilsj.le.20) iilsj=9
      if (ilsj.ge.21.and.ilsj.le.22) iilsj=10
      if (ilsj.ge.23.and.ilsj.le.23) iilsj=12
      if (ilsj.ge.24.and.ilsj.le.24) iilsj=15
      if (ilsj.eq.1) iord=0
      if (ilsj.eq.5) iord=0
      if (ilsj.eq.7) iord=0
      if (ilsj.eq.9) iord=0
      if (ilsj.eq.11) iord=0
      if (ilsj.eq.15) iord=0
      if (ilsj.eq.16) iord=0
      if (ilsj.eq.19) iord=0
      if (ilsj.eq.20) iord=0
      if (ilsj.eq.21) iord=0
      if (ilsj.eq.23) iord=0
      if (ilsj.eq.24) iord=0
      iord=iord+1
      !test
!      if(abs(cc(2)-2._qp)<0.1) cc(2)=2._qp
!      if(abs(cc(2)-2.5_qp)<0.1) cc(2)=4._qp
!      if(abs(cc(2)-3._qp)<0.1) cc(2)=2._qp
!      if(abs(cc(2)-4._qp)<0.1) cc(2)=4._qp
!      cc(2)=1._qp
      clsj(iord,iilsj)=cc(1)
      cutlsj(2*iord-1,iilsj)=cc(2)
      cutlsj(2*iord,iilsj)=cc(3)
      go to 61
!
!
!
!
!        lec LEC's
!        ---------
!
!
 6500 continue
!**** write (kwrite,10020) name,cc
      indlec=.true.
      ilec=ilec+1
      go to (6510,6510,6522,6540,6542,6544),ilec
 6510 do 6515 i=1,5
 6515 clec(ilec,i)=cc(i)
      go to 61
 6522 do 6523 i=1,2
 6523 clec(2,i+5)=cc(i)
      go to 61
 6540 do 6541 i=1,5
 6541 clec(3,i)=cc(i)
      go to 61
 6542 do 6543 i=1,5
 6543 clec(3,i+5)=cc(i)
      go to 61
 6544 do 6545 i=1,5
 6545 clec(3,i+10)=cc(i)
      go to 61
!
!
!
!
!        conclusions
!        -----------
!        -----------
!
!
!        write end
 7000 continue
!**** write (kwrite,10004) name
!**** write (kwrite,10008)
!**** write (kwrite,10008)
!
      if (indlsj) go to 7100
      if (indlec) go to 7500
      go to 8995
!
!
!        determine the low-energy constants (clec)
!        from the partial wave constants (clsj)
!
!
!        LEC's for Q^0 (LO)
!        ------------------
!
 7100 clec(1,1)=(clsj(1,1)+3._qp*clsj(1,5))*0.25_qp/(4._qp*pi)
      clec(1,2)=(clsj(1,5)  -   clsj(1,1))*0.25_qp/(4._qp*pi)
!
!
!        LEC's for Q^2 (NLO)
!        -------------------
!
!        vector b
!
!        1s0
      b(1)=clsj(2,1)
!        3p0
      b(2)=clsj(1,2)
!        1p1
      b(3)=clsj(1,3)
!        3p1
      b(4)=clsj(1,4)
!        3s1
      b(5)=clsj(2,5)
!        3s-d1
      b(6)=clsj(1,7)
!        3p2
      b(7)=clsj(1,10)
!
!
      do 7205 i=1,7
 7205 b(i)=b(i)/(4._qp*pi)
!
!
!        matrix a for the C parameters
!
!        1. column
      a(1)=1._qp
      a(2)=-2._qp/3._qp
      a(3)=a(2)
      a(4)=a(2)
      a(5)=1._qp
      a(6)=0._qp
      a(7)=a(2)
!
!        2. column
      a(8)=0.25_qp
      a(9)=1._qp/6._qp
      a(10)=a(9)
      a(11)=a(9)
      a(12)=0.25_qp
      a(13)=0._qp
      a(14)=a(9)
!
!        3. column
      a(15)=-3._qp
      a(16)=-2._qp/3._qp
      a(17)=2._qp
      a(18)=a(16)
      a(19)=1._qp
      a(20)=0._qp
      a(21)=a(16)
!
!        4. column
      a(22)=-0.75_qp
      a(23)=1._qp/6._qp
      a(24)=-0.5_qp
      a(25)=a(23)
      a(26)=0.25_qp
      a(27)=0._qp
      a(28)=a(23)
!
!        5. column
      a(29)=0._qp
      a(30)=-2._qp/3._qp
      a(31)=0._qp
      a(32)=-1._qp/3._qp
      a(33)=0._qp
      a(34)=0._qp
      a(35)=1._qp/3._qp
!
!        6. column
      a(36)=-1._qp
      a(37)=2._qp
      a(38)=2._qp/3._qp
      a(39)=-4._qp/3._qp
      a(40)=1._qp/3._qp
      a(41)=-2._qp*sqrt(2._qp)/3._qp
      a(42)=0._qp
!
!        7. column
      a(43)=-0.25_qp
      a(44)=-0.5_qp
      a(45)=-1._qp/6._qp
      a(46)=1._qp/3._qp
      a(47)=1._qp/12._qp
      a(48)=-sqrt(2._qp)/6._qp
      a(49)=0._qp
!
!
!
!
      call dgelg_qp (b,a,7,1,eps,ier)
!
      if (ier.ne.0) write (kwrite,19500) ier
19500 format (///' warning in chipar. the error index of dgelg is',
     & ' ier =',i5/' for the calculation of the C parameters.'///)
!

      do 7255 i=1,7
 7255 clec(2,i)=b(i)
!
!
!        LEC's for Q^4 (N^3LO)
!        -------------------
!
!        vector b
!
!        1s0
      b(1)=clsj(3,1)
!        1s0
      b(2)=clsj(4,1)
!        3p0
      b(3)=clsj(2,2)
!        1p1
      b(4)=clsj(2,3)
!        3p1
      b(5)=clsj(2,4)
!        3s1
      b(6)=clsj(3,5)
!        3s1
      b(7)=clsj(4,5)
!        3d1
      b(8)=clsj(1,6)
!        3s-d1
      b(9)=clsj(2,7)
!        3s-d1
      b(10)=clsj(3,7)
!        1d2
      b(11)=clsj(1,8)
!        3d2
      b(12)=clsj(1,9)
!        3p2
      b(13)=clsj(2,10)
!        3p-f2
      b(14)=clsj(1,12)
!        3d3
      b(15)=clsj(1,15)
!
!
      do 7305 i=1,15
 7305 b(i)=b(i)/(4._qp*pi)
!
!
!        matrix a for the D parameters
!
!        1. column
      a(1)=1._qp
      a(2)=10._qp/3._qp
      a(3)=-4._qp/3._qp
      a(4)=a(3)
      a(5)=a(3)
      a(6)=1._qp
      a(7)=a(2)
      a(8)=8._qp/15._qp
      a(9)=0._qp
      a(10)=0._qp
      a(11)=a(8)
      a(12)=a(8)
      a(13)=a(3)
      a(14)=0._qp
      a(15)=a(8)
!
!        2. column
      a(16)=1._qp/16._qp
      a(17)=5._qp/24._qp
      a(18)=1._qp/12._qp
      a(19)=a(18)
      a(20)=a(18)
      a(21)=a(16)
      a(22)=a(17)
      a(23)=1._qp/30._qp
      a(24)=0._qp
      a(25)=0._qp
      a(26)=a(23)
      a(27)=a(23)
      a(28)=a(18)
      a(29)=0._qp
      a(30)=a(23)
!
!        3. column
      a(31)=1._qp/4._qp
      a(32)=1._qp/6._qp
      a(33)=0._qp
      a(34)=0._qp
      a(35)=0._qp
      a(36)=a(31)
      a(37)=a(32)
      a(38)=-2._qp/15._qp
      a(39)=0._qp
      a(40)=0._qp
      a(41)=a(38)
      a(42)=a(38)
      a(43)=0._qp
      a(44)=0._qp
      a(45)=a(38)
!
!        4. column
      a(46)=0._qp
      a(47)=2._qp/3._qp
      a(48)=0._qp
      a(49)=0._qp
      a(50)=0._qp
      a(51)=0._qp
      a(52)=a(47)
      a(53)=-2._qp/15._qp
      a(54)=0._qp
      a(55)=0._qp
      a(56)=a(53)
      a(57)=a(53)
      a(58)=0._qp
      a(59)=0._qp
      a(60)=a(53)
!
!        5. column
      a(61)=-3._qp
      a(62)=-10._qp
      a(63)=-4._qp/3._qp
      a(64)=4._qp
      a(65)=a(63)
      a(66)=1._qp
      a(67)=10._qp/3._qp
      a(68)=8._qp/15._qp
      a(69)=0._qp
      a(70)=0._qp
      a(71)=-8._qp/5._qp
      a(72)=a(68)
      a(73)=a(63)
      a(74)=0._qp
      a(75)=a(68)
!
!        6. column
      a(76)=-3._qp/16._qp
      a(77)=-5._qp/8._qp
      a(78)=1._qp/12._qp
      a(79)=-1._qp/4._qp
      a(80)=a(78)
      a(81)=1._qp/16._qp
      a(82)=5._qp/24._qp
      a(83)=1._qp/30._qp
      a(84)=0._qp
      a(85)=0._qp
      a(86)=-1._qp/10._qp
      a(87)=a(83)
      a(88)=a(78)
      a(89)=0._qp
      a(90)=a(83)
!
!        7. column
      a(91)=-3._qp/4._qp
      a(92)=-1._qp/2._qp
      a(93)=0._qp
      a(94)=0._qp
      a(95)=0._qp
      a(96)=1._qp/4._qp
      a(97)=1._qp/6._qp
      a(98)=-2._qp/15._qp
      a(99)=0._qp
      a(100)=0._qp
      a(101)=2._qp/5._qp
      a(102)=a(98)
      a(103)=0._qp
      a(104)=0._qp
      a(105)=a(98)
!
!        8. column
      a(106)=0._qp
      a(107)=-2._qp
      a(108)=0._qp
      a(109)=0._qp
      a(110)=0._qp
      a(111)=0._qp
      a(112)=2._qp/3._qp
      a(113)=-2._qp/15._qp
      a(114)=0._qp
      a(115)=0._qp
      a(116)=2._qp/5._qp
      a(117)=a(113)
      a(118)=0._qp
      a(119)=0._qp
      a(120)=a(113)
!
!        9. column
      a(121)=0._qp
      a(122)=0._qp
      a(123)=-2._qp/3._qp
      a(124)=0._qp
      a(125)=-1._qp/3._qp
      a(126)=0._qp
      a(127)=0._qp
      a(128)=2._qp/5._qp
      a(129)=0._qp
      a(130)=0._qp
      a(131)=0._qp
      a(132)=2._qp/15._qp
      a(133)=1._qp/3._qp
      a(134)=0._qp
      a(135)=-4._qp/15._qp
!
!        10. column
      a(136)=0._qp
      a(137)=0._qp
      a(138)=-1._qp/6._qp
      a(139)=0._qp
      a(140)=-1._qp/12._qp
      a(141)=0._qp
      a(142)=0._qp
      a(143)=-1._qp/10._qp
      a(144)=0._qp
      a(145)=0._qp
      a(146)=0._qp
      a(147)=-1._qp/30._qp
      a(148)=1._qp/12._qp
      a(149)=0._qp
      a(150)=1._qp/15._qp
!
!        11. column
      a(151)=-1._qp
      a(152)=-10._qp/3._qp
      a(153)=8._qp/3._qp
      a(154)=4._qp/3._qp
      a(155)=-2._qp
      a(156)=1._qp/3._qp
      a(157)=10._qp/9._qp
      a(158)=-4._qp/9._qp
      a(159)=-2._qp*sqrt(2._qp)/3._qp
      a(160)=-14._qp*sqrt(2._qp)/9._qp
      a(161)=-8._qp/15._qp
      a(162)=4._qp/5._qp
      a(163)=-2._qp/15._qp
      a(164)=4._qp*sqrt(6._qp)/15._qp
      a(165)=0._qp
!
!        12. column
      a(166)=-1._qp/4._qp
      a(167)=-1._qp/6._qp
      a(168)=1._qp/3._qp
      a(169)=0._qp
      a(170)=-1._qp/6._qp
      a(171)=1._qp/12._qp
      a(172)=1._qp/18._qp
      a(173)=1._qp/9._qp
      a(174)=-sqrt(2._qp)/6._qp
      a(175)=sqrt(2._qp)/18._qp
      a(176)=2._qp/15._qp
      a(177)=-1._qp/5._qp
      a(178)=1._qp/30._qp
      a(179)=-sqrt(6._qp)/15._qp
      a(180)=0._qp
!
!        13. column
      a(181)=-1._qp/4._qp
      a(182)=-1._qp/6._qp
      a(183)=-1._qp/3._qp
      a(184)=0._qp
      a(185)=1._qp/6._qp
      a(186)=1._qp/12._qp
      a(187)=1._qp/18._qp
      a(188)=1._qp/9._qp
      a(189)=-sqrt(2._qp)/6._qp
      a(190)=sqrt(2._qp)/18._qp
      a(191)=2._qp/15._qp
      a(192)=-1._qp/5._qp
      a(193)=-1._qp/30._qp
      a(194)=sqrt(6._qp)/15._qp
      a(195)=0._qp
!
!        14. column
      a(196)=-1._qp/16._qp
      a(197)=-5._qp/24._qp
      a(198)=-1._qp/6._qp
      a(199)=-1._qp/12._qp
      a(200)=1._qp/8._qp
      a(201)=1._qp/48._qp
      a(202)=5._qp/72._qp
      a(203)=-1._qp/36._qp
      a(204)=-sqrt(2._qp)/24._qp
      a(205)=-7._qp*sqrt(2._qp)/72._qp
      a(206)=-1._qp/30._qp
      a(207)=1._qp/20._qp
      a(208)=1._qp/120._qp
      a(209)=-sqrt(6._qp)/60._qp
      a(210)=0._qp
!
!        15. column
      a(211)=0._qp
      a(212)=-2._qp/3._qp
      a(213)=0._qp
      a(214)=0._qp
      a(215)=0._qp
      a(216)=0._qp
      a(217)=2._qp/9._qp
      a(218)=-16._qp/45._qp
      a(219)=0._qp
      a(220)=2._qp*sqrt(2._qp)/9._qp
      a(221)=2._qp/15._qp
      a(222)=4._qp/15._qp
      a(223)=0._qp
      a(224)=0._qp
      a(225)=-2._qp/15._qp
!
!
!
!
      call dgelg_qp (b,a,15,1,eps,ier)
!
      if (ier.ne.0) write (kwrite,19501) ier
19501 format (///' warning in chipar. the error index of dgelg is',
     & ' ier =',i5/' for the calculation of the D parameters.'///)
!
!
      do 7355 i=1,15
 7355 clec(3,i)=b(i)
!
!
!
!        write LEC's
!        -----------
!
!
 7500 continue
!**** write (kwrite,10100)
10100 format (//' Low energy parameters (LEC):'/
     &          ' ----------------------------'/)
!
!        Q^0 (LO)
!**** write (kwrite,10101) (clec(1,i),i=1,2)
10101 format ('lec  CS,CT',2f10.6)
!
!        Q^2 (NLO)
!**** write (kwrite,10102) (clec(2,i),i=1,7)
10102 format ('lec  C_i  ',5f10.6)
!
!        Q^4 (N^3LO)
!**** write (kwrite,10103) (clec(3,i),i=1,15)
10103 format ('lec  D_i  ',5f10.6)
!
!
!
!
!        store LEC's appropriately
!        -------------------------
!
!
      iorder=0
 7600 iorder=iorder+1
!
!
      mg=10
      iterm=0
 7700 iterm=iterm+1
!
!
      if (iorder.eq.1.and.iterm.gt.2) go to 7600
      if (iorder.eq.2.and.iterm.gt.7) go to 7600
!
!
      mg=mg+1
!
      if (iorder.eq.2) then
      if (iterm.eq.2) mg=mg-1
      if (iterm.eq.4) mg=mg-1
      end if
!
      if (iorder.eq.3) then
      if (iterm.eq.2) mg=mg-1
      if (iterm.eq.3) mg=mg-1
      if (iterm.eq.4) mg=mg-1
      if (iterm.eq.6) mg=mg-1
      if (iterm.eq.7) mg=mg-1
      if (iterm.eq.8) mg=mg-1
      if (iterm.eq.10) mg=mg-1
      if (iterm.eq.12) mg=mg-1
      if (iterm.eq.14) mg=mg-1
      end if
!
!
      ime=ime+1
      if (ime.gt.imee) go to 9011
      mgg(mg,inter)=mgg(mg,inter)+1
      m=mgg(mg,inter)
      if (m.gt.mee) go to 9001
      ima(m,mg,inter)=ime
      if (m.eq.1) then
      imga(inter)=imga(inter)+1
      mggo(imga(inter),inter)=mg
      end if
!
!
      c(1,ime)=clec(iorder,iterm)*wnq*1e-2_qp
      ic(1,ime)=-1
!
!
      mi=4
      mm=5
!
!
      if (indca) then
      indlca=.true.
      name(1)=cut
      do i=1,5
      cc(i)=cca(i)
      end do
      go to 72
      end if
!
!
 7800 indlca=.false.
!
!
      if (iorder.eq.2) then
      c(1,ime)=c(1,ime)*wnq*1e-6_qp
      if (iterm.le.4) then
      imod=mod(iterm,2)
      if (imod.eq.0) imod=2
      ic(mi,ime)=10+imod
      end if
      end if
!
!
      if (iorder.eq.3) then
      c(1,ime)=c(1,ime)*(wnq*1e-6_qp)**2
      if (iterm.le.8) then
      imod=mod(iterm,4)
      if (imod.eq.0) imod=4
      ic(mi,ime)=20+imod
      end if
      if (iterm.ge.9.and.iterm.le.14) then
      imod=mod(iterm,2)
      if (imod.eq.0) imod=2
      ic(mi,ime)=10+imod
      end if
      end if
!
!
 7900 if (iterm.lt.15) go to 7700
      if (iorder.lt.3) go to 7600
!
!
 8995 imaa(inter)=imb
      imea(inter)=ime
      imb=ime+1
!
 8999 continue
!        this has been the end of the inter loop
!
      return
!
!
!
!        errors
!        ------
!        ------
!
!
!
!
 9000 write (kwrite,19000) name(1)
19000 format (1h ////' error in chipar:  contribution  ',a4,'   does not
     & exist in this program.'/' execution terminated.'////)
      go to 9999
!
!
 9001 write (kwrite,19001)
19001 format (1h ////' error in chipar:too many contributions within a g
     &roup with respect to'/' the given dimensions. execution terminated
     &.'////)
      go to 9999
!
!
 9002 write (kwrite,19002) cc(1)
19002 format (1h ////' error in chipar: cut/fun typ',f10.4,'  does not e
     &xist in this program.'/' execution terminated.'////)
      go to 9999
!
!
 9003 write (kwrite,19003) iftyp
19003 format (1h ////' error in chipar: factor typ has the non-permissib
     &le value',i4,' .'/' execution terminated.'////)
      go to 9999
!
!
 9004 write (kwrite,19004) cc(4)
19004 format (1h ////' error in chipar: isospin has the non-permissible
     &value',f10.4,'  .'/' execution terminated.'////)
      go to 9999
!
!
 9005 write (kwrite,19005) cc(5)
19005 format (1h ////' error in chipar: iprop/spe has the non-permissibl
     &e value',f10.4,'  .'/' execution terminated.'////)
      go to 9999
!
!
 9006 write (kwrite,19006) cc(2)
19006 format (1h ////' error in chipar: the index for the propagator of
     &the cut-off has the'/' non-permissible value',f10.4,'  . execution
     & terminated.'////)
      go to 9999
!
!
 9009 write (kwrite,19009)
19009 format (1h ////' error in chipar: the exponent of the cut-off is l
     &ess than zero.'/' execution terminated.'////)
      go to 9999
!
!
 9010 write (kwrite,19010)
19010 format (1h ////' error in chipar: too many cut/fun parameters with
     & respect to the given'/' dimensions. execution terminated.'////)
      go to 9999
!
!
 9011 write (kwrite,19011)
19011 format (1h ////' error in chipar:  too many contr. with respect to
     & the dimensions given'/' to this program. execution terminated.'
     &////)
      go to 9999
!
!
 9999 stop
      end
      subroutine chistr_qp (icase,max,mex)
!
!        chistr computes the structure of one-boson-exchanges
!
!
      use iso_fortran_env
      implicit complex*16 (a-h,o-z)
      integer,parameter:: qp=real64
!
!
!        common blocks
!
      common /crdwrt_qp/ kread,kwrite,kpunch,kda(9)
!
      common /cstate_qp/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
      common /cnn_qp/ inn
!
!
!        common block for all chi-subroutines
!
      common /cchi_qp/ vj(32,270),c(20,270),fff,ff,f(52),aa(200)
     &,ai(19,30),
     &                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     &                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(200)
     &,wt(200),
     &                ic(20,270),ift(3),mint(3),maxt(3),nt,
     &                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     &                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     &                indc(2,270),indpar(3),indxy
!
!         specifications for this common block
!
      logical indc,indxy,indpar
!
!     further specifications
!
      dimension vv(32)
      dimension tt(2,3)
      logical index
      logical indiso
      data jj/-1/
      data index/.false./
      save
!
!
!
!
      if (index) go to 50
      index=.true.
!
!
      do 1 ii=1,3
      tt(1,ii)=1._qp
    1 tt(2,ii)=-3._qp
!
!
!
!
!
   50 do 1095 m=max,mex
      im=ima(m,mg,inter)
!
!
      if (mc.ne.1) go to 60
!
!
!
!
!        call integrals
!        --------------
!
!
!
!
      call chiai_qp
!
!
!
!
   60 if (mc.lt.1) mc=1
!
      if (c(mc,im).eq.0._qp) go to 1095
!
!
!
!
!        nn-nn helicity amplitudes /combinations/
!        ----------------------------------------
!
!
!
!
!        basic structure (a factor of 2 is included in v5 and v6)
!
!
      ive=6
!
      vv(1)=f(1)*ai(1,m)+f(2)*ai(2,m)
      vv(2)=f(3)*ai(1,m)+f(4)*ai(3,m)
      vv(3)=f(5)*ai(1,m)+f(6)*ai(2,m)
      vv(4)=f(4)*ai(1,m)+f(3)*ai(3,m)
      vv(5)=f(7)*ai(4,m)
      vv(6)=f(8)*ai(4,m)
!
!
      go to (1000,120,130,140),icase
!
!
!        additional terms required for the tensor coupling
!        of the rho-meson or for certain operators,
!        like, the spin-orbit operator (`ls  ')
!
!
  120 vv(1)=vv(1)+f(9)*ai(5,m)
      vv(2)=vv(2)+f(10)*ai(2,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(10)*ai(5,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(10)*ai(6,m)
         e1=f(11)*ai(7,m)
      vv(5)=vv(5)+e1
      vv(6)=vv(6)+e1
      go to 1000
!
!
!        additional terms in case of 2+ mesons
!        not needed here
!
!
  130 continue
      go to 1000
!
!
!        additional terms needed for the sigma-l operator (`sl  ')
!
!
  140 vv(1)=vv(1)+f(6)*ai(5,m)
      vv(2)=vv(2)+f(1)*ai(5,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(1)*ai(11,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(1)*ai(12,m)
      vv(5)=vv(5)+f(6)*ai(13,m)
      vv(6)=vv(6)+f(6)*ai(13,m)
!
!
!
!
 1000 continue
!
!
!
!
!        set certain cases to zero
!
      if (j.ne.0) go to 1021
      vv(2)=0._qp
      vv(4)=0._qp
      vv(5)=0._qp
      vv(6)=0._qp
!
 1021 mmod=mod(j,2)
      if (.not.sing.or.(mmod.eq.1.and.inn.ne.2)) vv(1)=0._qp
      if (.not.trip.or.(mmod.eq.0.and.inn.ne.2)) vv(2)=0._qp
      if (coup.and.(mmod.eq.0.or.inn.eq.2)) go to 1030
      do 1025 iv=3,6
 1025 vv(iv)=0._qp
!
 1030 continue
!
!
!
!
!        transformation into lsj-formalism
!
      if (j.eq.jj) go to 1035
      jj=j
      aj=real(j,qp)
      aj1=real(j+1,qp)
      d2j1=1._qp/real(2*j+1,qp)
      arjj1=sqrt(aj*aj1)
!
 1035 v3=vv(3)
      v4=vv(4)
      v5=vv(5)
      v6=vv(6)
      v34=-arjj1*(v3-v4)
      v56=arjj1*(v5+v6)
      vv(3)=d2j1*(aj1*v3+aj*v4-v56)
      vv(4)=d2j1*(aj*v3+aj1*v4+v56)
      vv(5)=d2j1*(v34-aj1*v5+aj*v6)
      vv(6)=d2j1*(v34+aj*v5-aj1*v6)
!
!
!        possible different sign depending on the convention used
      vv(5)=-vv(5)
      vv(6)=-vv(6)
!
!
!
!
!        multiply with factors
!        ---------------------
!
!
!
!
 1040 is=mod(j,2)+1
      it=mod(is,2)+1
      indiso=indc(1,im)
      cmc=c(mc,im)
      fc=fff*ff*cmc
      do 1045 iv=1,ive
!
!        multiply with coupling-constant and factors fff and ff
!
      vv(iv)=vv(iv)*fc
!
!        multiply with isospin factor
!
      if (.not.indiso) go to 1045
      if (iv.eq.2) go to 1043
      vv(iv)=vv(iv)*tt(is,inter)
      go to 1045
 1043 vv(iv)=vv(iv)*tt(it,inter)
      !test
      !vv(iv)=1._qp
!     
!     add up in case of several couplings for one meson-exchange
!     and store
 1045 vj(iv,im)=vj(iv,im)+vv(iv)
!
!
 1095 continue
!
!
      return
      end
      subroutine chiai_qp
!
!        chiai integrates over theta
!
!
      use iso_fortran_env
      implicit complex*16 (a-h,o-z)
      integer,parameter:: qp=real64
!
      common /cpot_qp/   v,xmev,ymev
      common /cstate_qp/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
!
!
!        common block for all chi-subroutines
!
      common /cchi_qp/ vj(32,270),c(20,270),fff,ff,f(52),aa(200)
     &,ai(19,30),
     &                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     &                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(200)
     &,wt(200),
     &                ic(20,270),ift(3),mint(3),maxt(3),nt,
     &                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     &                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     &                indc(2,270),indpar(3),indxy
!
!         specifications for this common block
!
      logical indc,indxy,indpar
      complex*16:: v(6),xmev,ymev
!
!
!        further specifications
      dimension gi(7)
!
      dimension pj(7,200)
      complex*16 axy2,aomq,am
      logical indj
      data ige/7/
      data nnt/-1/,iinter/-1/,jj/-1/
      save
!
!
!
!
      if (inter.eq.iinter) go to 60
      iinter=inter
      min=mint(inter)
      max=maxt(inter)
!
      igeint=7
!
      wn=wnn(inter)
      dwn=1._qp/wn
      wnq=wn*wn
!
!
!
!
   60 if (j.eq.jj) go to 70
      jj=j
      indj=.false.
!
!
      aj=real(j,qp)
      aj1=real(j+1,qp)
      dj1=1._qp/aj1
      ajdj1=aj*dj1
      aaj=sqrt(ajdj1)
!
!
      aj2=real(j+2,qp)
      ajm1=real(j-1,qp)
!
!
      ajj1=aj*aj1
      ajj2=ajm1*aj2
      ajjb=aj*ajm1
!
      aajj=0._qp
      if (j.gt.1)
     &aajj=aj/sqrt(ajj1*ajj2)
!
      aaj1=aajj*ajm1
      aaj2=aajj*aj1
      aaj3=aajj*2._qp
!
      if (j.gt.1) go to 62
      aajj=0._qp
      go to 63
   62 aajj=1._qp/(aj1*sqrt(ajj2))
!
   63 aaj4=aajj*ajjb
      aaj5=aajj*aj1*2._qp
      aaj6=aajj*(ajj1+2._qp)
      aaj7=aajj*ajj2
!
!
!
!
!        find out appropriate number of gauss-points
!        -------------------------------------------
!
!
   70 c4=c(4,im)
      if (c4.eq.0._qp) then
      c4=(138._qp*dwn)**2
      end if
      iprsp=ic(1,im)
!
!
!        compute am
!
      axy2=xy2
      if (iprsp.ne.1) go to 91
      aomq=eem12+c4
      go to 92
   91 aomq=xxpyy+c4
!
   92 am=axy2/aomq
!
!
!        compute number of gausspoints (nt)
!
!
      if (real(am).gt.0.999) go to 94
!
!
      if (real(am).gt.0.85) am=am**(-log(1._qp-am)-0.9)
!
!
      nt=float(min)/(1._qp-am)+0.9
!
!
      if (nt.gt.max) nt=max
      go to 95
!
!
   94 nt=max
!
!
   95 nt=nt+j
!
!        compute nt, which is suitable for gset
!
      if (nt.le.16) go to 98
      if (nt.gt.24) go to 96
      nt=4*(nt/4)
      go to 98
   96 if (nt.gt.48) go to 97
      nt=8*(nt/8)
      go to 98
   97 nt=16*(nt/16)
      if (nt.gt.96) nt=96
!
   98 if (nt.eq.nnt.and.indj) go to 100
!
!
!
!
!        call gauss-points
!        -----------------
!
!
!
!
      nt=25*(int(nt/25)+1)
      if(nt>100) nt=100
      !call gset_qp (-1._qp,1._qp,nt,ct,wt)
      ct(1:nt)=gaus_legendre_points(1:nt,nt)
      wt(1:nt)=gaus_legendre_weights(1:nt,nt)
      nnt=nt
!
!
!
!
!        call legendre-polynoms if necessary
!        -----------------------------------
!
!
!
!
      indxy=.false.
      indj=.true.
!      do 99 i=1,nt
!      t=ct(i)
!      call legp_qp (pj(1,i),pj(3,i),t,j)
!      pj(2,i)=pj(1,i)*t
!      pj(4,i)=pj(2,i)*t
!      pj(6,i)=pj(4,i)*t
!      pj(5,i)=pj(3,i)*t
!   99 pj(7,i)=pj(5,i)*t
      pj(1,:nt)=legendre_polynom(:nt,nt,j)
      if(j==0) then
         pj(3,:nt)=0._qp
      else
         pj(3,:nt)=legendre_polynom(:nt,nt,j-1)
      endif
      pj(2,:nt)=pj(1,:nt)*ct(:nt)
      pj(4,:nt)=pj(2,:nt)*ct(:nt)
      pj(6,:nt)=pj(4,:nt)*ct(:nt)
      pj(5,:nt)=pj(3,:nt)*ct(:nt)
      pj(7,:nt)=pj(5,:nt)*ct(:nt)
!
!
!
!
!        call integrand
!        --------------
!
!
!
!
  100 call chiaa_qp
!
!
!
!
!        prepare for integration
!
!
!
!
      do 2001 ig=1,igeint
 2001 gi(ig)=0._qp
!
!
!
!
!        integration-loop of theta
!        -------------------------
!
!
!
!
      do 2005 i=1,nt
      do 2005 ig=1,igeint
 2005 gi(ig)=gi(ig)+pj(ig,i)*aa(i)
!         do ig=1,igeint
!            gi(ig)=SUM(pj(ig,:)*aa(:))
!         enddo
!
!
!
      if (j.ne.0) go to 2010
      gi(3)=0._qp
      gi(5)=0._qp
      gi(7)=0._qp
!
!
!
!
!        combinations of integrals
!        -------------------------
!
!
!
!
 2010 ai(1,m)=gi(1)
!
      ai(2,m)=gi(2)
      ai(3,m)= ajdj1*gi(2)+dj1*gi(3)
      gi23m  =gi(2)-gi(3)
      ai(4,m)=aaj*gi23m
!
!
      ai(5,m)=gi(4)
      ai(6,m)= ajdj1*gi(4)+dj1*gi(5)
      gi45m  =gi(4)-gi(5)
      ai(7,m)=aaj*gi45m
!
!
      ai( 8,m)= aaj1*gi(4)-aaj2*gi(1)+aaj3*gi(5)
      aai1    = aaj4*gi(4)+aaj5*gi(1)-aaj6*gi(5)
      aai2    = aaj7*gi23m
      ai( 9,m)= aai2+aai1
      ai(10,m)= aai2-aai1
!
!
      ai(11,m)=gi(6)
      ai(12,m)=ajdj1*gi(6)+dj1*gi(7)
      ai(13,m)=aaj*(gi(6)-gi(7))
      !test
      !ai(:,:)=1._qp
!
!
      return
      end
      subroutine chiaa_qp
!
!        chiaa computes propagators, cutoffs, and functions
!
!
      use iso_fortran_env
      implicit complex*16 (a-h,o-z)
      integer,parameter:: qp=real64
!
      common /crdwrt_qp/ kread,kwrite,kpunch,kda(9)
!
      common /cstate_qp/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
!
!
!        common block for all chi-subroutines
!
      common /cchi_qp/ vj(32,270),c(20,270),fff,ff,f(52),aa(200)
     & ,ai(19,30),
     &                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     &                gaa(3),fpia(3),ezz1(3),ezz2(3),ct(200)
     &,wt(200),
     &                ic(20,270),ift(3),mint(3),maxt(3),nt,
     &                mge,mgg(40,3),mggo(40,3),ima(30,40,3),
     &                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     &                indc(2,270),indpar(3),indxy
!
!         specifications for this common block
!
      logical indc,indxy,indpar
!
      common /compar_qp/ cb1a(3),cb2a(3),cb3a(3),cb4a(3),
     &                cd12a(3),cd3a(3),cd5a(3),cd145a(3)
!
!
      common /crrr_qp/ rrr
!
!
!
!        further specifications
      dimension deltaq(200,7)
      dimension ell(200),cpa(200),cpaa(200)
      logical indla
      data iinter/-1/
      data cc4/-1._qp/
!      data pi/3.141592653589793_qp/
      data pi/3.141592653589793238462643383279502884197_qp/
      save
!
!
!
!
      if (inter.eq.iinter) go to 10
      iinter=inter
      ga2=gaa(inter)
      ga4=ga2*ga2
      fpi2=fpia(inter)
!
      cb1=cb1a(inter)
      cb2=cb2a(inter)
      cb3=cb3a(inter)
      cb4=cb4a(inter)
      cd12=cd12a(inter)
      cd3=cd3a(inter)
      cd5=cd5a(inter)
      cd145=cd145a(inter)
!
      pi2=pi*pi
   10 continue
!
!
!
!
!        delta square
!        ------------
!
!
!
!
      if (indxy) go to 50
      indxy=.true.
      indla=.false.
      do 15 i=1,nt
      xy2t=xy2*ct(i)
!
!
!        function  -q^2 (- momentum-transfer-squared)
!        --------------
!
!        retardation ignored
!
      deltaq(i,1)=xy2t-xxpyy
!
!        retardation incorporated
!
      deltaq(i,2)=xy2t-eem12
!
!
!        function  +k^2 (average-momentum squared)
!        --------------
!
      deltaq(i,3)=(xy2t+xxpyy)*0.25_qp
!
!        function  q^4 (momentum-transfer to the power of 4)
!        -------------
!
      deltaq(i,4)=deltaq(i,1)*deltaq(i,1)
!
!        function  k^4 (average-momentum to the power of 4)
!        -------------
!
      deltaq(i,5)=deltaq(i,3)*deltaq(i,3)
!
!        function  +q^2*k^2
!        -----------------
!
      deltaq(i,6)=-deltaq(i,1)*deltaq(i,3)
!
!        function  (\vec q x \vec k)^2
!        -----------------------------
!
      deltaq(i,7)=xx*yy*(1._qp-ct(i)*ct(i))
!
   15 continue
      go to 50
!
!
!
!     calculate ell, cpa, and cpaa
!
   20 indla=.true.
      cc4=c4
      do 25 i=1,nt
      akk=-deltaq(i,1)
      ak=sqrt(akk)
      radi=4._qp*c4+akk
      root=sqrt(radi)
      deno=2._qp*sqrt(c4)
      ell(i)=root*log((root+ak)/deno)/ak
      cpa(i)=atan(ak/deno)/(2._qp*ak)
      cpaa(i)=(2._qp*c4+akk)*cpa(i)
   25 continue
      go to 6000
!
!
!
!
!        propagator
!        ----------
!        ----------
!
!
!
!
   50 c4=c(4,im)
      iprsp=ic(1,im)
      if (iprsp.lt.0) go to 60
      iret=iprsp+1
!
!         propagator for the nn case
      do 55 i=1,nt
   55 aa(i)=wt(i)/(c4-deltaq(i,iret))
      go to 80
!
!
!        "no propagator"
!
   60 do 65 i=1,nt
   65 aa(i)=wt(i)
!
!
   80 continue
!
!
!
!
!
!        cut-offs and functions
!        ----------------------
!        ----------------------
!
!
!
!
      mi=4
      mm=5
!
!
 5999 ityp=ic(mi,im)
      if (ityp.eq.0) go to 8000
      if (ityp.le.10) then
      iprspc=ic(mi+1,im)
      iret=iprspc+1
      end if
 6000 go to (100,100,300,9002,500,600,9002,9002,9002,1000,
     & 1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,
     & 2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,
     & 3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,
     & 4100,4200,4300,9002,9002,9002,9002,4800,4900,5000,
     & 5100,5200,5300,5400,5500,5600),ityp
!
!
!
!
!        cut-off of dipole type
!        **********************
!
!
  100 c5=c(mm,im)
      c6=c(mm+1,im)
      nexp=ic(mi+2,im)
!
      do 105 i=1,nt
!
      aaa=c5/(c6-deltaq(i,iret))
!     -------------------------
!
      do 105 ii=1,nexp
  105 aa(i)=aa(i)*aaa
!
!
      mi=mi+3
      mm=mm+2
      go to 5999
!
!
!
!
!        exponential form factor of momentum transfer
!        ********************************************
!
!
  300 c5=c(mm,im)
      c6=c(mm+1,im)
      !test
      !c6=2._qp
      do 305 i=1,nt
!
      expo=(c5*abs(deltaq(i,iret)))**c6
!     ----------------------------
!
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!     
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
!      call regulator_c2(c5*abs(deltaq(i,iret)),expexp
!     &,c6) 
      call regulator_poly_2(deltaq(i,iret)
     &,complex(0._qp,0._qp),expexp,1._qp/c5,c6)
! call regulator function that stores the factor in expexp
      aa(i)=aa(i)*expexp
!     ----------------------
!
  305 continue
      mi=mi+2
      mm=mm+2
      go to 5999
!
!
!
!
!        sharp cutoff of x and y
!        ***********************
!
!
  500 c5=c(mm,im)
!
      if (real(x).gt.real(c5).or.real(y).gt.real(c5)) then
!     ----------------------------
      do 505 i=1,nt
  505 aa(i)=0._qp
      end if
!
      mi=mi+2
      mm=mm+1
      go to 5999
!
!
!
!
!        exponential form factor of xx and yy
!        ************************************
!
!
  600 c5=c(mm,im)
      c6=c(mm+1,im)
!
      !test
      !c6=2._qp
      expo=(c5*xx)**c6+(c5*yy)**c6
!     ----------------------------
      if (real(expo).gt.real(rrr)) expo=rrr
      expexp=exp(-expo)
!      expexp=(exp(-(reg_n-1._qp)*expo)+exp(+(reg_n-1._qp)*expo))
!     &/(exp(-(reg_n)*expo)+exp(+(reg_n)*expo))
!      call regulator_c3(xx,yy,expexp,1._qp/c5
!     &                     ,c6) ! call regulator function that stores the factor in expexp
      call regulator_poly_2(xx,yy,expexp,1._qp/c5
     &                     ,c6) ! call regulator function that stores the factor in expexp
!     ------------------
!
      do 605 i=1,nt
  605 aa(i)=aa(i)*expexp
      mi=mi+2
      mm=mm+2
      go to 5999
!
!
!
!
!
!        pi-gamma potential
!        ******************
!
!
 1000 c5=c(mm,im)
      do 1055 i=1,nt
      betaq=-deltaq(i,1)/c4
      betaq1=betaq+1._qp
      aaa=-(1._qp-betaq)**2/(2._qp*betaq*betaq)*log(betaq1)
     &    +betaq1/(2._qp*betaq)
     &    -2._qp*c5
 1055 aa(i)=aa(i)*aaa
      mi=mi+2
      mm=mm+1
      go to 5999
!
!
!
!
!        function +q^2 (momentum-transfer squared)
!        *************
!
!
 1100 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      do 1105 i=1,nt
 1105 aa(i)=-aa(i)*deltaq(i,1)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        function k^2 (average-momentum squared)
!        ************
!
!
 1200 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      do 1205 i=1,nt
 1205 aa(i)=aa(i)*deltaq(i,3)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        function 1 for tpn1
!        *******************
!
!
 1300 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1305 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      brak=4._qp*c4*(5._qp*ga4-4._qp*ga2 -1._qp)
     &    +akk*(23._qp*ga4-10._qp*ga2-1._qp)
     &    +48._qp*ga4*c4*c4/radi
 1305 aa(i)=aa(i)*c5*ell(i)*brak
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        function 2 for tpn1
!        *******************
!
!
 1400 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1405 i=1,nt
 1405 aa(i)=aa(i)*c5*ell(i)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn2, function 1
!        ****************
!
!
 1500 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1505 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      term1=-ga2*c4**(2.5_qp)/(16._qp*radi)
      term2=(2._qp*c4*(2._qp*cb1-cb3)-akk*(cb3+3._qp/16._qp*ga2))
     &     *cpaa(i)
 1505 aa(i)=aa(i)*c5*(term1+term2)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn2, function 2
!        ****************
!
!
 1600 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1605 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      term1=-3._qp*ga2*c4**(2.5_qp)/radi
      term2=(4._qp*c4+2._qp*akk-ga2*(4._qp*c4+3._qp*akk))
     &     *cpaa(i)
 1605 aa(i)=aa(i)*c5*(term1+term2)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn2, function 3
!        ****************
!
!
 1700 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1705 i=1,nt
 1705 aa(i)=aa(i)*c5*cpaa(i)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn2, function 4
!        ****************
!
!
 1800 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1805 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      term1=(cb4+0.25_qp)*radi
      term2=-ga2/8._qp*(10._qp*c4+3._qp*akk)
 1805 aa(i)=aa(i)*c5*(term1+term2)*cpa(i)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn2, function 5
!        ****************
!
!
 1900 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 1905 i=1,nt
 1905 aa(i)=aa(i)*c5*cpaa(i)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn2, function 6
!        ****************
!
!
 2000 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 2005 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
 2005 aa(i)=aa(i)*c5*radi*cpa(i)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        function q^4 (momentum-transfer to the power of 4)
!        ************
!
!
 2100 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      do 2105 i=1,nt
 2105 aa(i)=aa(i)*deltaq(i,4)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        function k^4 (average-momentum to the power of 4)
!        ************
!
!
 2200 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      do 2205 i=1,nt
 2205 aa(i)=aa(i)*deltaq(i,5)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        function +q^2*k^2
!        *****************
!
!
 2300 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      do 2305 i=1,nt
 2305 aa(i)=aa(i)*deltaq(i,6)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        function  (\vec q x \vec k)^2
!        *****************************
!
!
 2400 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      do 2405 i=1,nt
 2405 aa(i)=aa(i)*deltaq(i,7)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!        function xy
!        ***********
!
 2500 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      aaxy=xy2*0.5_qp*c5
      do 2505 i=1,nt
 2505 aa(i)=aa(i)*aaxy
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!        function xx+yy
!        **************
!
 2600 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      do 2605 i=1,nt
 2605 aa(i)=aa(i)*xxpyy*c5
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!        function xx*xx+yy*yy
!        ********************
!
 2700 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      aaxy=(xx*xx+yy*yy)*c5
      do 2705 i=1,nt
 2705 aa(i)=aa(i)*aaxy
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!        function xx
!        ***********
!
 2800 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      do 2805 i=1,nt
 2805 aa(i)=aa(i)*xx*c5
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!        function yy
!        ***********
!
 2900 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      do 2905 i=1,nt
 2905 aa(i)=aa(i)*yy*c5
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn3, function 1
!        ****************
!
!
 3000 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 3005 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      brak=(cb2*radi/6._qp+cb3*(2._qp*c4+akk)-4._qp*cb1*c4)**2
     &    +(cb2*radi)**2/45._qp
 3005 aa(i)=aa(i)*c5*ell(i)*brak
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn3, function 2
!        ****************
!
!
 3100 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 3105 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
 3105 aa(i)=aa(i)*c5*ell(i)*radi
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!        function 1._qp
!        *************
!
 3200 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      do 3205 i=1,nt
 3205 aa(i)=aa(i)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!        function 1-q^2/8-k^/2
!        *********************
!
 3300 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      do 3305 i=1,nt
 3305 aa(i)=aa(i)*(1._qp+deltaq(i,1)/8._qp-deltaq(i,3)/2._qp)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!        function 1-q^2/8
!        ****************
!
 3400 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      do 3405 i=1,nt
 3405 aa(i)=aa(i)*(1._qp+deltaq(i,1)/8._qp)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!        function 1+k^/2
!        ***************
!
 3500 c5=c(mm,im)
      if (c5.eq.0._qp) c5=1._qp
      do 3505 i=1,nt
 3505 aa(i)=aa(i)*(1._qp+deltaq(i,3)/2._qp)*c5
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn3, function 3
!        ****************
!
!
 3600 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 3605 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      brak=radi+ga2*(8._qp*c4+5._qp*akk)
 3605 aa(i)=aa(i)*c5*ell(i)*brak
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn3, function 4
!        ****************
!
!
 3700 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 3705 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      brak=radi-ga2*(16._qp*c4+7._qp*akk)
 3705 aa(i)=aa(i)*c5*ell(i)*brak
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn3, function 5
!        ****************
!
!
 3800 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 3805 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
 3805 aa(i)=aa(i)*c5*ell(i)*radi
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn3, function 6
!        ****************
!
!
 3900 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 3905 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      brak=(cb2-6._qp*cb3)*akk*akk
     &     +4._qp*(6._qp*cb1+cb2-3._qp*cb3)*akk*c4
     &     +6._qp*(cb2-2._qp*cb3)*c4*c4
     &     +24._qp*(2._qp*cb1+cb3)*c4*c4*c4/radi
 3905 aa(i)=aa(i)*c5*ell(i)*brak
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn32, function 1
!        *****************
!
!
 4000 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 4005 i=1,nt
      akk=-deltaq(i,1)
      wpi=sqrt(c4)
      brak=(c4+2._qp*akk)*(2._qp*wpi+cpaa(i))
     &     +4._qp*ga2*wpi*(2._qp*c4+akk)
 4005 aa(i)=aa(i)*c5*brak*cpaa(i)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn32, function 2
!        *****************
!
!
 4100 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 4105 i=1,nt
      akk=-deltaq(i,1)
      wpi=sqrt(c4)
      radi=4._qp*c4+akk
      brak=radi*cpa(i)+2._qp*wpi
     &     +4._qp*ga2*wpi
 4105 aa(i)=aa(i)*c5*brak*radi*cpa(i)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn32, function 3
!        *****************
!
!
 4200 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 4205 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
 4205 aa(i)=aa(i)*c5*ell(i)*radi
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn32, function 4
!        *****************
!
!
 4300 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 4305 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      brak1=4._qp*c4*(2._qp*ga2+1._qp)+akk*(5._qp*ga2+1._qp)
      brak2=akk*(5._qp+13._qp*ga2)/3._qp+8._qp*c4*(1._qp+2._qp*ga2)
      brak=brak1*(brak1*ell(i)-brak2)
!
      brak3=384._qp*pi2*fpi2*(cd12*(2._qp*c4+akk)+4._qp*c4*cd5)
      brak4=2._qp*ga2*(2._qp*c4+akk)-3._qp/5._qp*(ga2-1._qp)*radi
      brak5=192._qp*pi2*fpi2*radi*cd3
      brakbrak=brak1*brak3+brak4*brak5
!
 4305 aa(i)=aa(i)*c5*(brak+brakbrak)*ell(i)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn3m (= N^3LO, 1/M^2 terms), function 1
!        ****************************************
!
!
 4800 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      c44=c4*c4
      c46=c44*c4
      c48=c46*c4
      do 4805 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      radiq=radi*radi
      brak1=c46/(2._qp*radi)
      brak2=(2._qp*c48/radiq+8._qp*c46/radi-akk*akk-2._qp*c44)*ell(i)
 4805 aa(i)=aa(i)*c5*(brak1+brak2)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn3m, function 2
!        *****************
!
!
 4900 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      c44=c4*c4
      c46=c44*c4
      c48=c46*c4
      yyoff=0.5_qp*(xx+yy)
      do 4905 i=1,nt
      akk=-deltaq(i,1)
      akkq=akk*akk
      radi=4._qp*c4+akk
      radiq=radi*radi
      brak1=(radi*(akk-4._qp*yyoff)
     &      +8._qp*ga2*(11._qp/4._qp*akkq+5._qp*c4*akk+3._qp*c44
     &       -6._qp*c46/radi-yyoff*(8._qp*c4+5._qp*akk))
     &      +4._qp*ga4*(yyoff*(20._qp*c4+7._qp*akk-16._qp*c44/radi)
     &       +16._qp*c48/radiq+12._qp*c46/radi-27._qp/4._qp*akkq
     &       -11._qp*c4*akk-6._qp*c44))*ell(i)
      brak2=16._qp*ga4*c46/radi
 4905 aa(i)=aa(i)*c5*(brak1+brak2)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn3m, function 3
!        *****************
!
!
 5000 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      c44=c4*c4
      yyoff=0.5_qp*(xx+yy)
      do 5005 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      brak=yyoff+3._qp/8._qp*akk+c44/radi
 5005 aa(i)=aa(i)*c5*brak*ell(i)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn3m, function 4
!        *****************
!
!
 5100 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      c44=c4*c4
      do 5105 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      brak=radi-32._qp*ga2*(c4+7._qp/16._qp*akk)
     &     +4._qp*ga4*(7._qp*c4+17._qp/4._qp*akk+4._qp*c44/radi)
 5105 aa(i)=aa(i)*c5*brak*ell(i)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn3m, function 5
!        *****************
!
!
 5200 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      c44=c4*c4
      do 5205 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      brak=-radi+16._qp*ga2*(c4+3._qp/8._qp*akk)
     &     +4._qp/3._qp*ga4*(-9._qp*c4-11._qp/4._qp*akk+4._qp*c44/radi)
 5205 aa(i)=aa(i)*c5*brak*ell(i)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn3m, function 6
!        *****************
!
!
 5300 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      c44=c4*c4
      do 5305 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      brak=11._qp/32._qp*akk+c44/radi
 5305 aa(i)=aa(i)*c5*brak*ell(i)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn3m, function 7
!        *****************
!
!
 5400 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 5405 i=1,nt
 5405 aa(i)=aa(i)*c5*ell(i)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn2c (correction for our it 2pi), function 1
!        *********************************************
!
!
 5500 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 5505 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      term1=sqrt(c4)*radi
      term2=(2._qp*c4+akk)*cpaa(i)
 5505 aa(i)=aa(i)*c5*(term1+term2)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
!        tpn2c (correction for our it 2pi), function 2
!        *********************************************
!
!
 5600 if (.not.indla.or.cc4.ne.c4) go to 20
      c5=c(mm,im)
      do 5605 i=1,nt
      akk=-deltaq(i,1)
      radi=4._qp*c4+akk
      term1=sqrt(c4)
      term2=radi*cpa(i)
 5605 aa(i)=aa(i)*c5*(term1+term2)
      mi=mi+1
      mm=mm+1
      go to 5999
!
!
!
!
 9002 write (kwrite,19002) ityp
19002 format (1h ////' error in chiaa:  cut/fun typ',i10  ,'  does not e
     &xist in this program.'/' execution terminated.'////)
      stop
!
!
!
!
!
 8000 return
      end
      subroutine legp_qp (pj,pjm1,x,j)
!
!
!        subroutine legp   computes the legendre polynominals
!
      use iso_fortran_env
      integer,parameter:: qp=real64
      complex*16 pj,pjm1,x,a,b
!
!
!
!        compute legendre polynom for j equals zero
!
!
      if (j.gt.0) go to 1
      pj=1._qp
      pjm1=0._qp
      if (j.lt.0) pj=0._qp
      return
!
!
!
!        compute legendre polynoms for j equals one
!
!
!
    1 pj=x
      pjm1=1._qp
      if (j.eq.1) return
!
!
!
!        compute legendre polynom for j greater or equal two
!
!
!
      do 2 i=2,j
      a=x*pj
      b=a-pjm1
      pjm1=pj
    2 pj=-b/real(i,qp)+b+a
!
!
      return
      end
      subroutine gset_qp(ax,bx,n,z,w)
!
!
!        this code has been obtained from the CERN computer library
!        in the year of the lord 1972.
!
!
      use iso_fortran_env
      implicit real*16 (a-h,o-z)
      integer,parameter:: qp=real64
!
      real*16::ax,bx
!     n-point gauss zeros and weights for the interval (ax,bx) are
!           stored in  arrays z and w respectively.
!
      dimension     a(273),x(273),ktab(96)
      dimension z(*),w(*)
!
!-----table of initial subscripts for n=2(1)16(4)96
      data ktab(2)/1/
      data ktab(3)/2/
      data ktab(4)/4/
      data ktab(5)/6/
      data ktab(6)/9/
      data ktab(7)/12/
      data ktab(8)/16/
      data ktab(9)/20/
      data ktab(10)/25/
      data ktab(11)/30/
      data ktab(12)/36/
      data ktab(13)/42/
      data ktab(14)/49/
      data ktab(15)/56/
      data ktab(16)/64/
      data ktab(20)/72/
      data ktab(24)/82/
      data ktab(28)/82/
      data ktab(32)/94/
      data ktab(36)/94/
      data ktab(40)/110/
      data ktab(44)/110/
      data ktab(48)/130/
      data ktab(52)/130/
      data ktab(56)/130/
      data ktab(60)/130/
      data ktab(64)/154/
      data ktab(68)/154/
      data ktab(72)/154/
      data ktab(76)/154/
      data ktab(80)/186/
      data ktab(84)/186/
      data ktab(88)/186/
      data ktab(92)/186/
      data ktab(96)/226/
!
!-----table of abscissae (x) and weights (a) for interval (-1,+1).
!
!**** n=2
      data x(1)/0.577350269189626_qp/, a(1)/1.000000000000000_qp/
!**** n=3
      data x(2)/0.774596669241483_qp/, a(2)/0.555555555555556_qp/
      data x(3)/0.000000000000000_qp/, a(3)/0.888888888888889_qp/
!**** n=4
      data x(4)/0.861136311594053_qp/, a(4)/0.347854845137454_qp/
      data x(5)/0.339981043584856_qp/, a(5)/0.652145154862546_qp/
!**** n=5
      data x(6)/0.906179845938664_qp/, a(6)/0.236926885056189_qp/
      data x(7)/0.538469310105683_qp/, a(7)/0.478628670499366_qp/
      data x(8)/0.000000000000000_qp/, a(8)/0.568888888888889_qp/
!**** n=6
      data x(9)/0.932469514203152_qp/, a(9)/0.171324492379170_qp/
      data x(10)/0.661209386466265_qp/, a(10)/0.360761573048139_qp/
      data x(11)/0.238619186083197_qp/, a(11)/0.467913934572691_qp/
!**** n=7
      data x(12)/0.949107912342759_qp/, a(12)/0.129484966168870_qp/
      data x(13)/0.741531185599394_qp/, a(13)/0.279705391489277_qp/
      data x(14)/0.405845151377397_qp/, a(14)/0.381830050505119_qp/
      data x(15)/0.000000000000000_qp/, a(15)/0.417959183673469_qp/
!**** n=8
      data x(16)/0.960289856497536_qp/, a(16)/0.101228536290376_qp/
      data x(17)/0.796666477413627_qp/, a(17)/0.222381034453374_qp/
      data x(18)/0.525532409916329_qp/, a(18)/0.313706645877887_qp/
      data x(19)/0.183434642495650_qp/, a(19)/0.362683783378362_qp/
!**** n=9
      data x(20)/0.968160239507626_qp/, a(20)/0.081274388361574_qp/
      data x(21)/0.836031107326636_qp/, a(21)/0.180648160694857_qp/
      data x(22)/0.613371432700590_qp/, a(22)/0.260610696402935_qp/
      data x(23)/0.324253423403809_qp/, a(23)/0.312347077040003_qp/
      data x(24)/0.000000000000000_qp/, a(24)/0.330239355001260_qp/
!**** n=10
      data x(25)/0.973906528517172_qp/, a(25)/0.066671344308688_qp/
      data x(26)/0.865063366688985_qp/, a(26)/0.149451349150581_qp/
      data x(27)/0.679409568299024_qp/, a(27)/0.219086362515982_qp/
      data x(28)/0.433395394129247_qp/, a(28)/0.269266719309996_qp/
      data x(29)/0.148874338981631_qp/, a(29)/0.295524224714753_qp/
!**** n=11
      data x(30)/0.978228658146057_qp/, a(30)/0.055668567116174_qp/
      data x(31)/0.887062599768095_qp/, a(31)/0.125580369464905_qp/
      data x(32)/0.730152005574049_qp/, a(32)/0.186290210927734_qp/
      data x(33)/0.519096129206812_qp/, a(33)/0.233193764591990_qp/
      data x(34)/0.269543155952345_qp/, a(34)/0.262804544510247_qp/
      data x(35)/0.000000000000000_qp/, a(35)/0.272925086777901_qp/
!**** n=12
      data x(36)/0.981560634246719_qp/, a(36)/0.047175336386512_qp/
      data x(37)/0.904117256370475_qp/, a(37)/0.106939325995318_qp/
      data x(38)/0.769902674194305_qp/, a(38)/0.160078328543346_qp/
      data x(39)/0.587317954286617_qp/, a(39)/0.203167426723066_qp/
      data x(40)/0.367831498998180_qp/, a(40)/0.233492536538355_qp/
      data x(41)/0.125233408511469_qp/, a(41)/0.249147045813403_qp/
!**** n=13
      data x(42)/0.984183054718588_qp/, a(42)/0.040484004765316_qp/
      data x(43)/0.917598399222978_qp/, a(43)/0.092121499837728_qp/
      data x(44)/0.801578090733310_qp/, a(44)/0.138873510219787_qp/
      data x(45)/0.642349339440340_qp/, a(45)/0.178145980761946_qp/
      data x(46)/0.448492751036447_qp/, a(46)/0.207816047536889_qp/
      data x(47)/0.230458315955135_qp/, a(47)/0.226283180262897_qp/
      data x(48)/0.000000000000000_qp/, a(48)/0.232551553230874_qp/
!**** n=14
      data x(49)/0.986283808696812_qp/, a(49)/0.035119460331752_qp/
      data x(50)/0.928434883663574_qp/, a(50)/0.080158087159760_qp/
      data x(51)/0.827201315069765_qp/, a(51)/0.121518570687903_qp/
      data x(52)/0.687292904811685_qp/, a(52)/0.157203167158194_qp/
      data x(53)/0.515248636358154_qp/, a(53)/0.185538397477938_qp/
      data x(54)/0.319112368927890_qp/, a(54)/0.205198463721296_qp/
      data x(55)/0.108054948707344_qp/, a(55)/0.215263853463158_qp/
!**** n=15
      data x(56)/0.987992518020485_qp/, a(56)/0.030753241996117_qp/
      data x(57)/0.937273392400706_qp/, a(57)/0.070366047488108_qp/
      data x(58)/0.848206583410427_qp/, a(58)/0.107159220467172_qp/
      data x(59)/0.724417731360170_qp/, a(59)/0.139570677926154_qp/
      data x(60)/0.570972172608539_qp/, a(60)/0.166269205816994_qp/
      data x(61)/0.394151347077563_qp/, a(61)/0.186161000015562_qp/
      data x(62)/0.201194093997435_qp/, a(62)/0.198431485327111_qp/
      data x(63)/0.000000000000000_qp/, a(63)/0.202578241925561_qp/
!**** n=16
      data x(64)/0.989400934991650_qp/, a(64)/0.027152459411754_qp/
      data x(65)/0.944575023073233_qp/, a(65)/0.062253523938648_qp/
      data x(66)/0.865631202387832_qp/, a(66)/0.095158511682493_qp/
      data x(67)/0.755404408355003_qp/, a(67)/0.124628971255534_qp/
      data x(68)/0.617876244402644_qp/, a(68)/0.149595988816577_qp/
      data x(69)/0.458016777657227_qp/, a(69)/0.169156519395003_qp/
      data x(70)/0.281603550779259_qp/, a(70)/0.182603415044924_qp/
      data x(71)/0.095012509837637_qp/, a(71)/0.189450610455069_qp/
!**** n=20
      data x(72)/0.993128599185094_qp/, a(72)/0.017614007139152_qp/
      data x(73)/0.963971927277913_qp/, a(73)/0.040601429800386_qp/
      data x(74)/0.912234428251325_qp/, a(74)/0.062672048334109_qp/
      data x(75)/0.839116971822218_qp/, a(75)/0.083276741576704_qp/
      data x(76)/0.746331906460150_qp/, a(76)/0.101930119817240_qp/
      data x(77)/0.636053680726515_qp/, a(77)/0.118194531961518_qp/
      data x(78)/0.510867001950827_qp/, a(78)/0.131688638449176_qp/
      data x(79)/0.373706088715419_qp/, a(79)/0.142096109318382_qp/
      data x(80)/0.227785851141645_qp/, a(80)/0.149172986472603_qp/
      data x(81)/0.076526521133497_qp/, a(81)/0.152753387130725_qp/
!**** n=24
      data x(82)/0.995187219997021_qp/, a(82)/0.012341229799987_qp/
      data x(83)/0.974728555971309_qp/, a(83)/0.028531388628933_qp/
      data x(84)/0.938274552002732_qp/, a(84)/0.044277438817419_qp/
      data x(85)/0.886415527004401_qp/, a(85)/0.059298584915436_qp/
      data x(86)/0.820001985973902_qp/, a(86)/0.073346481411080_qp/
      data x(87)/0.740124191578554_qp/, a(87)/0.086190161531953_qp/
      data x(88)/0.648093651936975_qp/, a(88)/0.097618652104113_qp/
      data x(89)/0.545421471388839_qp/, a(89)/0.107444270115965_qp/
      data x(90)/0.433793507626045_qp/, a(90)/0.115505668053725_qp/
      data x(91)/0.315042679696163_qp/, a(91)/0.121670472927803_qp/
      data x(92)/0.191118867473616_qp/, a(92)/0.125837456346828_qp/
      data x(93)/0.064056892862605_qp/, a(93)/0.127938195346752_qp/
!**** n=32
      data x(94)/0.997263861849481_qp/, a(94)/0.007018610009470_qp/
      data x(95)/0.985611511545268_qp/, a(95)/0.016274394730905_qp/
      data x(96)/0.964762255587506_qp/, a(96)/0.025392065309262_qp/
      data x(97)/0.934906075937739_qp/, a(97)/0.034273862913021_qp/
      data x(98)/0.896321155766052_qp/, a(98)/0.042835898022226_qp/
      data x(99)/0.849367613732569_qp/, a(99)/0.050998059262376_qp/
      data x(100)/0.794483795967942_qp/, a(100)/0.058684093478535_qp/
      data x(101)/0.732182118740289_qp/, a(101)/0.065822222776361_qp/
      data x(102)/0.663044266930215_qp/, a(102)/0.072345794108848_qp/
      data x(103)/0.587715757240762_qp/, a(103)/0.078193895787070_qp/
      data x(104)/0.506899908932229_qp/, a(104)/0.083311924226946_qp/
      data x(105)/0.421351276130635_qp/, a(105)/0.087652093004403_qp/
      data x(106)/0.331868602282127_qp/, a(106)/0.091173878695763_qp/
      data x(107)/0.239287362252137_qp/, a(107)/0.093844399080804_qp/
      data x(108)/0.144471961582796_qp/, a(108)/0.095638720079274_qp/
      data x(109)/0.048307665687738_qp/, a(109)/0.096540088514727_qp/
!**** n=40
      data x(110)/0.998237709710559_qp/, a(110)/0.004521277098533_qp/
      data x(111)/0.990726238699457_qp/, a(111)/0.010498284531152_qp/
      data x(112)/0.977259949983774_qp/, a(112)/0.016421058381907_qp/
      data x(113)/0.957916819213791_qp/, a(113)/0.022245849194166_qp/
      data x(114)/0.932812808278676_qp/, a(114)/0.027937006980023_qp/
      data x(115)/0.902098806968874_qp/, a(115)/0.033460195282547_qp/
      data x(116)/0.865959503212259_qp/, a(116)/0.038782167974472_qp/
      data x(117)/0.824612230833311_qp/, a(117)/0.043870908185673_qp/
      data x(118)/0.778305651426519_qp/, a(118)/0.048695807635072_qp/
      data x(119)/0.727318255189927_qp/, a(119)/0.053227846983936_qp/
      data x(120)/0.671956684614179_qp/, a(120)/0.057439769099391_qp/
      data x(121)/0.612553889667980_qp/, a(121)/0.061306242492928_qp/
      data x(122)/0.549467125095128_qp/, a(122)/0.064804013456601_qp/
      data x(123)/0.483075801686178_qp/, a(123)/0.067912045815233_qp/
      data x(124)/0.413779204371605_qp/, a(124)/0.070611647391286_qp/
      data x(125)/0.341994090825758_qp/, a(125)/0.072886582395804_qp/
      data x(126)/0.268152185007253_qp/, a(126)/0.074723169057968_qp/
      data x(127)/0.192697580701371_qp/, a(127)/0.076110361900626_qp/
      data x(128)/0.116084070675255_qp/, a(128)/0.077039818164247_qp/
      data x(129)/0.038772417506050_qp/, a(129)/0.077505947978424_qp/
!**** n=48
      data x(130)/0.998771007252426_qp/, a(130)/0.003153346052305_qp/
      data x(131)/0.993530172266350_qp/, a(131)/0.007327553901276_qp/
      data x(132)/0.984124583722826_qp/, a(132)/0.011477234579234_qp/
      data x(133)/0.970591592546247_qp/, a(133)/0.015579315722943_qp/
      data x(134)/0.952987703160430_qp/, a(134)/0.019616160457355_qp/
      data x(135)/0.931386690706554_qp/, a(135)/0.023570760839324_qp/
      data x(136)/0.905879136715569_qp/, a(136)/0.027426509708356_qp/
      data x(137)/0.876572020274247_qp/, a(137)/0.031167227832798_qp/
      data x(138)/0.843588261624393_qp/, a(138)/0.034777222564770_qp/
      data x(139)/0.807066204029442_qp/, a(139)/0.038241351065830_qp/
      data x(140)/0.767159032515740_qp/, a(140)/0.041545082943464_qp/
      data x(141)/0.724034130923814_qp/, a(141)/0.044674560856694_qp/
      data x(142)/0.677872379632663_qp/, a(142)/0.047616658492490_qp/
      data x(143)/0.628867396776513_qp/, a(143)/0.050359035553854_qp/
      data x(144)/0.577224726083972_qp/, a(144)/0.052890189485193_qp/
      data x(145)/0.523160974722233_qp/, a(145)/0.055199503699984_qp/
      data x(146)/0.466902904750958_qp/, a(146)/0.057277292100403_qp/
      data x(147)/0.408686481990716_qp/, a(147)/0.059114839698395_qp/
      data x(148)/0.348755886292160_qp/, a(148)/0.060704439165893_qp/
      data x(149)/0.287362487355455_qp/, a(149)/0.062039423159892_qp/
      data x(150)/0.224763790394689_qp/, a(150)/0.063114192286254_qp/
      data x(151)/0.161222356068891_qp/, a(151)/0.063924238584648_qp/
      data x(152)/0.097004699209462_qp/, a(152)/0.064466164435950_qp/
      data x(153)/0.032380170962869_qp/, a(153)/0.064737696812683_qp/
!**** n=64
      data x(154)/0.999305041735772_qp/, a(154)/0.001783280721696_qp/
      data x(155)/0.996340116771955_qp/, a(155)/0.004147033260562_qp/
      data x(156)/0.991013371476744_qp/, a(156)/0.006504457968978_qp/
      data x(157)/0.983336253884625_qp/, a(157)/0.008846759826363_qp/
      data x(158)/0.973326827789910_qp/, a(158)/0.011168139460131_qp/
      data x(159)/0.961008799652053_qp/, a(159)/0.013463047896718_qp/
      data x(160)/0.946411374858402_qp/, a(160)/0.015726030476024_qp/
      data x(161)/0.929569172131939_qp/, a(161)/0.017951715775697_qp/
      data x(162)/0.910522137078502_qp/, a(162)/0.020134823153530_qp/
      data x(163)/0.889315445995114_qp/, a(163)/0.022270173808383_qp/
      data x(164)/0.865999398154092_qp/, a(164)/0.024352702568710_qp/
      data x(165)/0.840629296252580_qp/, a(165)/0.026377469715054_qp/
      data x(166)/0.813265315122797_qp/, a(166)/0.028339672614259_qp/
      data x(167)/0.783972358943341_qp/, a(167)/0.030234657072402_qp/
      data x(168)/0.752819907260531_qp/, a(168)/0.032057928354851_qp/
      data x(169)/0.719881850171610_qp/, a(169)/0.033805161837141_qp/
      data x(170)/0.685236313054233_qp/, a(170)/0.035472213256882_qp/
      data x(171)/0.648965471254657_qp/, a(171)/0.037055128540240_qp/
      data x(172)/0.611155355172393_qp/, a(172)/0.038550153178615_qp/
      data x(173)/0.571895646202634_qp/, a(173)/0.039953741132720_qp/
      data x(174)/0.531279464019894_qp/, a(174)/0.041262563242623_qp/
      data x(175)/0.489403145707052_qp/, a(175)/0.042473515123653_qp/
      data x(176)/0.446366017253464_qp/, a(176)/0.043583724529323_qp/
      data x(177)/0.402270157963991_qp/, a(177)/0.044590558163756_qp/
      data x(178)/0.357220158337668_qp/, a(178)/0.045491627927418_qp/
      data x(179)/0.311322871990210_qp/, a(179)/0.046284796581314_qp/
      data x(180)/0.264687162208767_qp/, a(180)/0.046968182816210_qp/
      data x(181)/0.217423643740007_qp/, a(181)/0.047540165714830_qp/
      data x(182)/0.169644420423992_qp/, a(182)/0.047999388596458_qp/
      data x(183)/0.121462819296120_qp/, a(183)/0.048344762234802_qp/
      data x(184)/0.072993121787799_qp/, a(184)/0.048575467441503_qp/
      data x(185)/0.024350292663424_qp/, a(185)/0.048690957009139_qp/
!**** n=80
      data x(186)/0.999553822651630_qp/, a(186)/0.001144950003186_qp/
      data x(187)/0.997649864398237_qp/, a(187)/0.002663533589512_qp/
      data x(188)/0.994227540965688_qp/, a(188)/0.004180313124694_qp/
      data x(189)/0.989291302499755_qp/, a(189)/0.005690922451403_qp/
      data x(190)/0.982848572738629_qp/, a(190)/0.007192904768117_qp/
      data x(191)/0.974909140585727_qp/, a(191)/0.008683945269260_qp/
      data x(192)/0.965485089043799_qp/, a(192)/0.010161766041103_qp/
      data x(193)/0.954590766343634_qp/, a(193)/0.011624114120797_qp/
      data x(194)/0.942242761309872_qp/, a(194)/0.013068761592401_qp/
      data x(195)/0.928459877172445_qp/, a(195)/0.014493508040509_qp/
      data x(196)/0.913263102571757_qp/, a(196)/0.015896183583725_qp/
      data x(197)/0.896675579438770_qp/, a(197)/0.017274652056269_qp/
      data x(198)/0.878722567678213_qp/, a(198)/0.018626814208299_qp/
      data x(199)/0.859431406663111_qp/, a(199)/0.019950610878141_qp/
      data x(200)/0.838831473580255_qp/, a(200)/0.021244026115782_qp/
      data x(201)/0.816954138681463_qp/, a(201)/0.022505090246332_qp/
      data x(202)/0.793832717504605_qp/, a(202)/0.023731882865930_qp/
      data x(203)/0.769502420135041_qp/, a(203)/0.024922535764115_qp/
      data x(204)/0.744000297583597_qp/, a(204)/0.026075235767565_qp/
      data x(205)/0.717365185362099_qp/, a(205)/0.027188227500486_qp/
      data x(206)/0.689637644342027_qp/, a(206)/0.028259816057276_qp/
      data x(207)/0.660859898986119_qp/, a(207)/0.029288369583267_qp/
      data x(208)/0.631075773046871_qp/, a(208)/0.030272321759557_qp/
      data x(209)/0.600330622829751_qp/, a(209)/0.031210174188114_qp/
      data x(210)/0.568671268122709_qp/, a(210)/0.032100498673487_qp/
      data x(211)/0.536145920897131_qp/, a(211)/0.032941939397645_qp/
      data x(212)/0.502804111888784_qp/, a(212)/0.033733214984611_qp/
      data x(213)/0.468696615170544_qp/, a(213)/0.034473120451753_qp/
      data x(214)/0.433875370831756_qp/, a(214)/0.035160529044747_qp/
      data x(215)/0.398393405881969_qp/, a(215)/0.035794393953416_qp/
      data x(216)/0.362304753499487_qp/, a(216)/0.036373749905835_qp/
      data x(217)/0.325664370747701_qp/, a(217)/0.036897714638276_qp/
      data x(218)/0.288528054884511_qp/, a(218)/0.037365490238730_qp/
      data x(219)/0.250952358392272_qp/, a(219)/0.037776364362001_qp/
      data x(220)/0.212994502857666_qp/, a(220)/0.038129711314477_qp/
      data x(221)/0.174712291832646_qp/, a(221)/0.038424993006959_qp/
      data x(222)/0.136164022809143_qp/, a(222)/0.038661759774076_qp/
      data x(223)/0.097408398441584_qp/, a(223)/0.038839651059051_qp/
      data x(224)/0.058504437152420_qp/, a(224)/0.038958395962769_qp/
      data x(225)/0.019511383256793_qp/, a(225)/0.039017813656306_qp/
!**** n=96
      data x(226)/0.999689503883230_qp/, a(226)/0.000796792065552_qp/
      data x(227)/0.998364375863181_qp/, a(227)/0.001853960788946_qp/
      data x(228)/0.995981842987209_qp/, a(228)/0.002910731817934_qp/
      data x(229)/0.992543900323762_qp/, a(229)/0.003964554338444_qp/
      data x(230)/0.988054126329623_qp/, a(230)/0.005014202742927_qp/
      data x(231)/0.982517263563014_qp/, a(231)/0.006058545504235_qp/
      data x(232)/0.975939174585136_qp/, a(232)/0.007096470791153_qp/
      data x(233)/0.968326828463264_qp/, a(233)/0.008126876925698_qp/
      data x(234)/0.959688291448742_qp/, a(234)/0.009148671230783_qp/
      data x(235)/0.950032717784437_qp/, a(235)/0.010160770535008_qp/
      data x(236)/0.939370339752755_qp/, a(236)/0.011162102099838_qp/
      data x(237)/0.927712456722308_qp/, a(237)/0.012151604671088_qp/
      data x(238)/0.915071423120898_qp/, a(238)/0.013128229566961_qp/
      data x(239)/0.901460635315852_qp/, a(239)/0.014090941772314_qp/
      data x(240)/0.886894517402420_qp/, a(240)/0.015038721026994_qp/
      data x(241)/0.871388505909296_qp/, a(241)/0.015970562902562_qp/
      data x(242)/0.854959033434601_qp/, a(242)/0.016885479864245_qp/
      data x(243)/0.837623511228187_qp/, a(243)/0.017782502316045_qp/
      data x(244)/0.819400310737931_qp/, a(244)/0.018660679627411_qp/
      data x(245)/0.800308744139140_qp/, a(245)/0.019519081140145_qp/
      data x(246)/0.780369043867433_qp/, a(246)/0.020356797154333_qp/
      data x(247)/0.759602341176647_qp/, a(247)/0.021172939892191_qp/
      data x(248)/0.738030643744400_qp/, a(248)/0.021966644438744_qp/
      data x(249)/0.715676812348967_qp/, a(249)/0.022737069658329_qp/
      data x(250)/0.692564536642171_qp/, a(250)/0.023483399085926_qp/
      data x(251)/0.668718310043916_qp/, a(251)/0.024204841792364_qp/
      data x(252)/0.644163403784967_qp/, a(252)/0.024900633222483_qp/
      data x(253)/0.618925840125468_qp/, a(253)/0.025570036005349_qp/
      data x(254)/0.593032364777572_qp/, a(254)/0.026212340735672_qp/
      data x(255)/0.566510418561397_qp/, a(255)/0.026826866725591_qp/
      data x(256)/0.539388108324357_qp/, a(256)/0.027412962726029_qp/
      data x(257)/0.511694177154667_qp/, a(257)/0.027970007616848_qp/
      data x(258)/0.483457973920596_qp/, a(258)/0.028497411065085_qp/
      data x(259)/0.454709422167743_qp/, a(259)/0.028994614150555_qp/
      data x(260)/0.425478988407300_qp/, a(260)/0.029461089958167_qp/
      data x(261)/0.395797649828908_qp/, a(261)/0.029896344136328_qp/
      data x(262)/0.365696861472313_qp/, a(262)/0.030299915420827_qp/
      data x(263)/0.335208522892625_qp/, a(263)/0.030671376123669_qp/
      data x(264)/0.304364944354496_qp/, a(264)/0.031010332586313_qp/
      data x(265)/0.273198812591049_qp/, a(265)/0.031316425596861_qp/
      data x(266)/0.241743156163840_qp/, a(266)/0.031589330770727_qp/
      data x(267)/0.210031310460567_qp/, a(267)/0.031828758894411_qp/
      data x(268)/0.178096882367618_qp/, a(268)/0.032034456231992_qp/
      data x(269)/0.145973714654896_qp/, a(269)/0.032206204794030_qp/
      data x(270)/0.113695850110665_qp/, a(270)/0.032343822568575_qp/
      data x(271)/0.081297495464425_qp/, a(271)/0.032447163714064_qp/
      data x(272)/0.048812985136049_qp/, a(272)/0.032516118713868_qp/
      data x(273)/0.016276744849602_qp/, a(273)/0.032550614492363_qp/
!
!
!-----test n
      alpha=0.5_qp*(ax+bx)
      beta=0.5_qp*(bx-ax)
      if( n.lt.1 .or. n.gt.96 ) go to 100
      if(n.ne.1) go to 1
      z(1)=alpha
      w(1)=bx-ax
      return
!
    1 if (n.le.16) go to 3
      if (n.gt.24) go to 4
      n=4*(n/4)
      go to 3
    4 if (n.gt.48) go to 5
      n=8*(n/8)
      go to 3
    5 n=16*(n/16)
!
!----- set k equal to initial subscript and store results
    3 k=ktab(n)
      m=n/2
      do 2 j=1,m
      jtab=k-1+j
      wtemp=beta*a(jtab)
      delta=beta*x(jtab)
      z(j)=alpha-delta
      w(j)=wtemp
      jp=n+1-j
      z(jp)=alpha+delta
      w(jp)=wtemp
    2 continue
      if((n-m-m).eq.0) return
      z(m+1)=alpha
      jmid=k+m
      w(m+1)=beta*a(jmid)
      return
!
  100 zn=n
      write(6,200) zn
  200 format(/////' error in gset. n has the non-permissible value',
     &e11.3/' execution terminated.')
      stop
      end
! name:    dgelg
!        programmbibliothek rhrz bonn        02/02/81       dgelg
!                                            fortran iv     ibm 370/168
!
! purpose:
!
! to solve a general system of simultaneous linear equations.
!
! usage:   call dgelg(r,a,m,n,eps,ier)
!
! parameters:
!
! r:       double precision m by n right hand side matrix
!          (destroyed). on return r contains the solutions
!          of the equations.
!
! a:       double precision m by m coefficient matrix
!          (destroyed).
!
! m:       the number of equations in the system.
!
! n:       the number of right hand side vectors.
!
! eps:     single precision input constant which is used as
!          relative tolerance for test on loss of
!          significance.
!
! ier:     resulting error parameter coded as follows
!           ier=0  - no error,
!           ier=-1 - no result because of m less than 1 or
!                   pivot element at any elimination step
!                   equal to 0,
!           ier=k  - warning due to possible loss of signifi-
!                   cance indicated at elimination step k+1,
!                   where pivot element was less than or
!                   equal to the internal tolerance eps times
!                   absolutely greatest element of matrix a.
!
! remarks: (1) input matrices r and a are assumed to be stored
!              columnwise in m*n resp. m*m successive storage
!              locations. on return solution matrix r is stored
!              columnwise too.
!          (2) the procedure gives results if the number of equations m
!              is greater than 0 and pivot elements at all elimination
!              steps are different from 0. however warning ier=k - if
!              given indicates possible loss of significance. in case
!              of a well scaled matrix a and appropriate tolerance eps,
!              ier=k may be interpreted that matrix a has the rank k.
!              no warning is given in case m=1.
!
! method:
!
! solution is done by means of gauss-elimination with
! complete pivoting.
!
! programs required:
!          none
!
! access:
!
! load module:    sys3.fortlib(dgelg)
! source module:  sys3.symlib.fortran(dgelg)
! description:    sys3.infolib(dgelg)
!
! author:         ibm, ssp iii
! installation:   ibm 370/168, mvs-jes2, fortran iv (h ext. enh.)
!
!**********************************************************************
      subroutine dgelg_qp(r,a,m,n,eps,ier)
!
!
      use iso_fortran_env
      implicit complex*16 (a-h,o-z)
      integer,parameter:: qp=real64
      dimension a(*),r(*)
      real*16 eps
!
!
!
!
      if(m)23,23,1
!
!     search for greatest element in matrix a
    1 ier=0
      piv=0._qp
      mm=m*m
      nm=n*m
      do 3 l=1,mm
      tb=abs(a(l))
      if(real(tb-piv))3,3,2
    2 piv=tb
      i=l
    3 continue
      tol=eps*piv
!     a(i) is pivot element. piv contains the absolute value of a(i).
!
!
!     start elimination loop
      lst=1
      do 17 k=1,m
!
!     test on singularity
      if(real(piv))23,23,4
    4 if(real(ier))7,5,7
    5 if(real(piv-tol))6,6,7
    6 ier=k-1
    7 pivi=1._qp/a(i)
      j=(i-1)/m
      i=i-j*m-k
      j=j+1-k
!     i+k is row-index, j+k column-index of pivot element
!
!     pivot row reduction and row interchange in right hand side r
      do 8 l=k,nm,m
      ll=l+i
      tb=pivi*r(ll)
      r(ll)=r(l)
    8 r(l)=tb
!
!     is elimination terminated
      if(k-m)9,18,18
!
!     column interchange in matrix a
    9 lend=lst+m-k
      if(j)12,12,10
   10 ii=j*m
      do 11 l=lst,lend
      tb=a(l)
      ll=l+ii
      a(l)=a(ll)
   11 a(ll)=tb
!
!     row interchange and pivot row reduction in matrix a
   12 do 13 l=lst,mm,m
      ll=l+i
      tb=pivi*a(ll)
      a(ll)=a(l)
   13 a(l)=tb
!
!     save column interchange information
      a(lst)=j
!
!     element reduction and next pivot search
      piv=0._qp
      lst=lst+1
      j=0
      do 16 ii=lst,lend
      pivi=-a(ii)
      ist=ii+m
      j=j+1
      do 15 l=ist,mm,m
      ll=l-j
      a(l)=a(l)+pivi*a(ll)
      tb=abs(a(l))
      if(real(tb-piv))15,15,14
   14 piv=tb
      i=l
   15 continue
      do 16 l=k,nm,m
      ll=l+j
   16 r(ll)=r(ll)+pivi*r(l)
   17 lst=lst+m
!     end of elimination loop
!
!
!     back substitution and back interchange
   18 if(m-1)23,22,19
   19 ist=mm+m
      lst=m+1
      do 21 i=2,m
      ii=lst-i
      ist=ist-lst
      l=ist-m
      l=a(l)+.5_qp
      do 21 j=ii,nm,m
      tb=r(j)
      ll=j
      do 20 k=ist,mm,m
      ll=ll+1
   20 tb=tb-a(k)*r(ll)
      k=j+l
      r(j)=r(k)
   21 r(k)=tb
   22 return
!
!
!     error return
   23 ier=-1
      return
      end
!**********************************************************************
!********************************************************************
! name:    dminv
!        programmbibliothek rhrz bonn        28/11/78       dminv
!                                            fortran iv     ibm 370/168
!
! purpose:
!
! invert a matrix
!
! usage:   call dminv (a,n,d,l,m)
!
! parameters:
!
! a:       input matrix, destroyed in computation and replaced by
!          resultant inverse.
!          double precision required.
!
! n:       order of matrix a
!
! d:       resultant determinant
!          double precision required.
!
! l:       work vector of length n
!
! m:       work vector of length n
!
! remarks: matrix a must be a general matrix
!
! method:
!
! the standard gauss-jordan method is used. the determinant
! is also calculated. a determinant of zero indicates that
! the matrix is singular.
!
! programs required:
!          none
!
! author:  ibm, ssp iii
!
!**********************************************************************
      subroutine dminv_qp (a,n,d,l,m)
       use iso_fortran_env

      implicit complex*16 (a-h,o-z)
      integer,parameter:: qp=real64
      dimension a(*),l(*),m(*)

!
!
!        search for largest element
!
      d=1._qp
      nk=-n
      do 80 k=1,n
      nk=nk+n
      l(k)=k
      m(k)=k
      kk=nk+k
      biga=a(kk)
      do 20 j=k,n
      iz=n*(j-1)
      do 20 i=k,n
      ij=iz+i
   10 if (abs(biga)-abs(a(ij)))  15,20,20
   15 biga=a(ij)
      l(k)=i
      m(k)=j
   20 continue
!
!        interchange rows
!
      j=l(k)
      if(j-k) 35,35,25
   25 ki=k-n
      do 30 i=1,n
      ki=ki+n
      hold=-a(ki)
      ji=ki-k+j
      a(ki)=a(ji)
   30 a(ji) =hold
!
!        interchange columns
!
   35 i=m(k)
      if(i-k) 45,45,38
   38 jp=n*(i-1)
      do 40 j=1,n
      jk=nk+j
      ji=jp+j
      hold=-a(jk)
      a(jk)=a(ji)
   40 a(ji) =hold
!
!        divide column by minus pivot (value of pivot element is
!        contained in biga)
!
   45 if(real(biga)) 48,46,48
   46 d=0._qp
      return
   48 do 55 i=1,n
      if(i-k) 50,55,50
   50 ik=nk+i
      a(ik)=a(ik)/(-biga)
   55 continue
!
!        reduce matrix
!
      do 65 i=1,n
      ik=nk+i
      hold=a(ik)
      ij=i-n
      do 65 j=1,n
      ij=ij+n
      if(i-k) 60,65,60
   60 if(j-k) 62,65,62
   62 kj=ij-i+k
      a(ij)=hold*a(kj)+a(ij)
   65 continue
!
!        divide row by pivot
!
      kj=k-n
      do 75 j=1,n
      kj=kj+n
      if(j-k) 70,75,70
   70 a(kj)=a(kj)/biga
   75 continue
!
!        product of pivots
!
      d=d*biga
!
!        replace pivot by reciprocal
!
      a(kk)=1._qp/biga
   80 continue
!
!        final row and column interchange
!
      k=n
  100 k=(k-1)
      if(k) 150,150,105
  105 i=l(k)
      if(i-k) 120,120,108
  108 jq=n*(k-1)
      jr=n*(i-1)
      do 110 j=1,n
      jk=jq+j
      hold=a(jk)
      ji=jr+j
      a(jk)=-a(ji)
  110 a(ji) =hold
  120 j=m(k)
      if(j-k) 100,100,125
  125 ki=k-n
      do 130 i=1,n
      ki=ki+n
      hold=a(ki)
      ji=ki-k+j
      a(ki)=-a(ji)
  130 a(ji) =hold
      go to 100
  150 return
      end
!****************this is the end of the program n3lo **********************
      end module n3lo_qp_c_module
