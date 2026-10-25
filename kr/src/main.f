C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .    P R O G R A M   kr3e_opt.for                                   .
C .    DOUBLE IONIZATION (e,3e) OF Kr 4p -> Kr2+ (4p^4, 3P)           .
C .    MODEL: fast projectile in the first Born approximation (plane  .
C .    waves), the two slow ejected electrons in a BBK final state    .
C .    (two Coulomb waves with variable charge Z(r) + the e-e factor  .
C .    C12, Fourier transformed over p). Optimized version of         .
C .    DI-Kr-4p-BBK-zv-3P.for (same physics, same routines tcs/tcp).  .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .
C .  WHAT CHANGED WITH RESPECT TO DI-Kr-4p-BBK-zv-3P.for
C .  ----------------------------------------------------
C .  1. Wave functions: Clementi tables (iwf=1, tables 4.1-4.3) or the
C .     data of the original code (iwf=2), in one routine (KRWF). Every
C .     orbital is normalised at run time (norms printed). Table 4.1,
C .     3S column, 9th row: 0.02180 -> 0.00218 (misprint: with it the
C .     3s norm is 0.990 and <3s|4s> = 0.019 ; corrected : 1.000 , 0).
C .     The 2S column of table 4.1 is not orthonormal (norm 0.982,
C .     <1s|2s> = 0.019): it only enters Z(r), through its density, and
C .     is renormalised.
C .  2. Z(r) = 36 - sum occ * screening of each orbital density, for
C .     any Slater basis (ZKR). The original az2p/az3p/az4p dropped
C .     the 2p-4p cross terms (av2p4p summed into av2s4s).
C .  3. The ionized 4p radial function and the analytic part use the
C .     same coefficients (the original Fr4p/qkr4p had 0.003172
C .     instead of 0.03172 : 4p norm 1.0041).
C .  4. Analytic part: TCP3 gives m = -1, 0, 1 in one call (the
C .     q-functions are shared). Identical to tcp (checked).
C .  5. Z(r) correction (tdN of the original, which used ONE point,
C .     nlg=nlg1=1): exact partial-wave method (PWCOR). The difference
C .     D = conj(phi_Z(r) - phi_Z=2) depends on r and k.r only:
C .        D = sum_l d_l(r) P_l(k^.r^)
C .     and the Fourier integral with R4p(r) Y1m is
C .        I_m(q) = (2pi)^-1.5 sum_{l,L,M} 4pi i^L J_lL(q)
C .                 4pi/(2l+1) Y*_l,mu(k^) Y*_LM(q^) <Y_LM Y_1m Y_l,mu>
C .        J_lL(q) = int_0^dz d_l(r) R4p(r) j_L(qr) r^2 dr
C .     (mu = -m-M, L = l-1, l+1). J_lL is tabulated once on a |q| grid;
C .     the angular integration is exact.
C .  6. Electron 2 has a fixed momentum: its amplitudes are computed
C .     once on the whole p grid, not for every angle.
C .  7. OpenMP over the p grid ; no 960 MB static arrays.
C .
C .  VALIDATION (scripts in kr/tests)
C .   - TCP3 = tcp for m=-1,0,1 (all n, general vectors).
C .   - PWCOR = brute-force 3D integral of conj(phi_Z - phi_2) R4p Y1m
C .     exp(i q.r) (same 1F1 routine).
C .   - izv=0, iwf=2 : identical TDCS to the original program with
C .     the Z(r) term removed (and 0.03172), on a reduced p grid.
C .
C .  COMPILATION:
C .   gfortran -O3 -march=native -fopenmp -ffixed-line-length-132
C .            kr3e_opt.for -o kr3e_opt
C .   ifx -O3 -xHOST -qopenmp -extend-source 132 kr3e_opt.for
C .            -o kr3e_opt
C .  RUN:  export OMP_NUM_THREADS=8 ; echo "1 73" | ./kr3e_opt
C .  Input "ni n" = first and last index of theta_1 (theta_1 =
C .  dth1*(index-1)). Output: kr3e_3P.dat (theta_1 , TDCS).
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      program kr3e
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      parameter (lmx=30,nqx=4000,nrx=400)
      common /wfn/ nb(3),nn(3,12),zt(3,12),co(3,5,12),nocc(3,5)
     $ ,norb(3)
      dimension dxp(400),dwp(400),dxtp(200),dwtp(200),dxhp(200)
     $ ,dwhp(200)
      double precision rgr(nrx),wgr(nrx),r4g(nrx),xx(400),ww(400)
      double precision ZKR,R4P,GAMOWF,qmax
      complex*16, allocatable :: e2t(:,:),e2s(:,:)
      double precision, allocatable :: ppx(:),ppy(:),ppz(:),pwt(:)
      complex*16, allocatable :: dl1(:,:),dl2(:,:),jt1(:,:,:)
     $ ,jt2(:,:,:),cc1(:,:,:,:),cc2(:,:,:,:)
      complex*16 q1t(-1:1),q1s(-1:1),q2t(-1:1),q2s(-1:1)
      complex*16 a1(-1:1),a2(-1:1)
      complex*16 t3,t5,t8,s3,s5,s8,tk12

      api=dacos(-1.d0)
      xi=(0.d0,1.d0)
      acon=2.d0/27.21165d0
      arad=api/180.d0

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  SWITCHES
C .   iwf = 1 : Clementi tables 4.1-4.3 ; 2 : data of the original code
C .   izv = 1 : Z(r) correction (PWCOR) ; 0 : constant charge Z = 2
C .   lmax, dz, ndr x nlr : partial waves and radial grid of PWCOR
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      iwf=1
      izv=1
      lmax=16
      dz=8.d0
      ndr=8
      nlr=12

C .  KINEMATICS (as DI-Kr-4p-BBK-zv-3P.for)
      enioni=38.d0
      ee1=7.d0
      ee2=7.d0
      ees=5500.d0
      eei=ees+ee1+ee2+enioni
      ets=1.d0*arad
      his=180.d0*arad
      et2=104.d0*arad
      hi2=0.d0
      dth1=5.d0
C .  asymptotic charge seen by the ejected electrons ; C12 constants
      az=2.d0
      abeta=0.005d0
      az12=0.5d0

C .  p GRID (as the original : |p| in 4 intervals , 80 x 80 angles)
      dp0=0.d0 ; dp1=0.7d0 ; dp2=1.6d0 ; dp3=2.1d0 ; dpf=10.d0
      np1=20 ; np2=120 ; np3=40 ; np4=5
      nt=80
      nh=80

      aki=dsqrt(eei*acon)
      aks=dsqrt(ees*acon)
      ak1=dsqrt(ee1*acon)
      ak2=dsqrt(ee2*acon)
      akix=0.d0 ; akiy=0.d0 ; akiz=aki
      aksx=aks*dsin(ets)*dcos(his)
      aksy=aks*dsin(ets)*dsin(his)
      aksz=aks*dcos(ets)
      ak2x=ak2*dsin(et2)*dcos(hi2)
      ak2y=ak2*dsin(et2)*dsin(hi2)
      ak2z=ak2*dcos(et2)
      adx=akix-aksx ; ady=akiy-aksy ; adz=akiz-aksz
      akk=adx*adx+ady*ady+adz*adz

      allocate(dl1(nrx,0:lmx),dl2(nrx,0:lmx),jt1(0:nqx,0:lmx,2)
     $ ,jt2(0:nqx,0:lmx,2),cc1(-1:1,0:lmx,2,-lmx-1:lmx+1)
     $ ,cc2(-1:1,0:lmx,2,-lmx-1:lmx+1))

C .  WAVE FUNCTIONS , Z(r)
      call KRWF(iwf)
      write(*,'(A,F10.5,A,F10.5,A,F10.5)') ' Z(0)=',ZKR(1.d-8),
     $ '  Z(4)=',ZKR(4.d0),'  Z(dz)=',ZKR(dz)

C .  p grid as flat arrays (index ip)
      np=0
      call PINT(dp0,dp1,np1,xx,ww,np,dxp,dwp)
      call PINT(dp1,dp2,np2,xx,ww,np,dxp,dwp)
      call PINT(dp2,dp3,np3,xx,ww,np,dxp,dwp)
      call PINT(dp3,dpf,np4,xx,ww,np,dxp,dwp)
      call GAULEG(0.d0,api,dxtp,dwtp,nt,2d-16)
      call GAULEG(0.d0,2.d0*api,dxhp,dwhp,nh,2d-16)
      npt=np*nt*nh
      allocate(ppx(npt),ppy(npt),ppz(npt),pwt(npt))
      allocate(e2t(-1:1,npt),e2s(-1:1,npt))
      ip=0
      do j1=1,np
      do j2=1,nh
      do j3=1,nt
         ip=ip+1
         ppx(ip)=dxp(j1)*dsin(dxtp(j3))*dcos(dxhp(j2))
         ppy(ip)=dxp(j1)*dsin(dxtp(j3))*dsin(dxhp(j2))
         ppz(ip)=dxp(j1)*dcos(dxtp(j3))
         pwt(ip)=dwp(j1)*dxp(j1)*dxp(j1)*dwhp(j2)*dwtp(j3)
     $        *dsin(dxtp(j3))
      enddo
      enddo
      enddo
      write(*,'(A,I9)') ' p grid points: ',npt

C .  Z(r) CORRECTION : radial grid , 4p on it , d_l(r) , J_lL(q)
      qmax=dpf+dsqrt(akk)+0.5d0
      nq=min(nqx,nint(qmax/0.005d0))
      hq=qmax/nq
      if (izv.eq.1) then
         nr=0
         do id=1,ndr
            call GAULEG((id-1)*dz/ndr,id*dz/ndr,xx,ww,nlr,1d-15)
            do i=1,nlr
               nr=nr+1
               rgr(nr)=xx(i)
               wgr(nr)=ww(i)
               r4g(nr)=R4P(rgr(nr))
            enddo
         enddo
         call PWDL(ak1,az,nr,rgr,lmax,dl1)
         call JTAB(nr,rgr,wgr,r4g,dl1,lmax,nq,hq,jt1)
         if (dabs(ak2-ak1).lt.1d-12) then
            dl2=dl1
            jt2=jt1
         else
            call PWDL(ak2,az,nr,rgr,lmax,dl2)
            call JTAB(nr,rgr,wgr,r4g,dl2,lmax,nq,hq,jt2)
         endif
         call CCOEF(ak2x,ak2y,ak2z,lmax,cc2)
      endif

C .  ELECTRON 2 (fixed k2) ON THE WHOLE p GRID :
C .   e2t = <k2| at K-p ,  e2s = <k2| at -p   (analytic + Z(r) part)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,64) PRIVATE(ip,a1,a2,m)
      do ip=1,npt
         call AN4P(az,adx-ppx(ip),ady-ppy(ip),adz-ppz(ip)
     $    ,ak2x,ak2y,ak2z,a1)
         call AN4P(az,-ppx(ip),-ppy(ip),-ppz(ip),ak2x,ak2y,ak2z,a2)
         if (izv.eq.1) then
            call PWQ3(adx-ppx(ip),ady-ppy(ip),adz-ppz(ip)
     $       ,lmax,cc2,jt2,nq,hq,a1)
            call PWQ3(-ppx(ip),-ppy(ip),-ppz(ip),lmax,cc2,jt2,nq,hq,a2)
         endif
         do m=-1,1
            e2t(m,ip)=a1(m)
            e2s(m,ip)=a2(m)
         enddo
      enddo
!$OMP END PARALLEL DO

C .  GAMOW FACTORS OF THE EJECTED ELECTRONS
      anc1=GAMOWF(az/ak1)
      anc2=GAMOWF(az/ak2)

C .  ANGLES
      Read(*,*,end=999) ni,n
      goto 998
999   ni=1
      n=73
998   continue
      open(unit=4,file='kr3e_3P.dat',status='unknown')

      do 800 i1=ni,n
      at1=dth1*(i1-1)
      et1=at1*arad
      ak1x=ak1*dsin(et1)
      ak1y=0.d0
      ak1z=ak1*dcos(et1)
      ak12x=(ak1x-ak2x)/2.d0
      ak12y=(ak1y-ak2y)/2.d0
      ak12z=(ak1z-ak2z)/2.d0
      ak12=dsqrt(ak12x**2+ak12y**2+ak12z**2)
      an1=aks*ak1*ak2/(aki*akk*akk)
      anc3=GAMOWF(-az12/ak12)
      if (izv.eq.1) call CCOEF(ak1x,ak1y,ak1z,lmax,cc1)

      s3=(0.d0,0.d0) ; s5=(0.d0,0.d0) ; s8=(0.d0,0.d0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,64) PRIVATE(ip,q1t,q1s,q2t,q2s,m
!$OMP& ,tk12,t3,t5,t8) REDUCTION(+:s3,s5,s8)
      do ip=1,npt
         call AN4P(az,adx+ppx(ip),ady+ppy(ip),adz+ppz(ip)
     $    ,ak1x,ak1y,ak1z,q1t)
         call AN4P(az,ppx(ip),ppy(ip),ppz(ip),ak1x,ak1y,ak1z,q1s)
         if (izv.eq.1) then
            call PWQ3(adx+ppx(ip),ady+ppy(ip),adz+ppz(ip)
     $       ,lmax,cc1,jt1,nq,hq,q1t)
            call PWQ3(ppx(ip),ppy(ip),ppz(ip),lmax,cc1,jt1,nq,hq,q1s)
         endif
         do m=-1,1
            q2t(m)=e2t(m,ip)
            q2s(m)=e2s(m,ip)
         enddo
         tk12=tcs(abeta,-az12,ak12x-ppx(ip),ak12y-ppy(ip)
     $    ,ak12z-ppz(ip),ak12x,ak12y,ak12z,1)
C .   3P components : t(m1,m2) = q1t q2s + q1s q2t - 2 q1s q2s
         t3=TPAIR(-1,0,q1t,q1s,q2t,q2s)-TPAIR(0,-1,q1t,q1s,q2t,q2s)
         t5=TPAIR(0,1,q1t,q1s,q2t,q2s)-TPAIR(1,0,q1t,q1s,q2t,q2s)
         t8=TPAIR(-1,1,q1t,q1s,q2t,q2s)-TPAIR(1,-1,q1t,q1s,q2t,q2s)
         s3=s3+t3*tk12*pwt(ip)
         s5=s5+t5*tk12*pwt(ip)
         s8=s8+t8*tk12*pwt(ip)
      enddo
!$OMP END PARALLEL DO
      s3=s3/dsqrt(2.d0) ; s5=s5/dsqrt(2.d0) ; s8=s8/dsqrt(2.d0)
      atfi=3.d0*(cdabs(s3)**2+cdabs(s5)**2+cdabs(s8)**2)
      atfi=atfi*an1*16.d0*api*anc1*anc2*anc3
      write(4,'(F10.3,2X,E18.10)') at1,atfi
      call flush(4)
      write(*,'(A,F8.2,A,E16.8)') ' theta1 = ',at1,'  TDCS = ',atfi
 800  continue
      close(unit=4)
      end
