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
C .  Parameters: defaults below, or file kr3e.inp in the run directory
C .  (namelist, e.g.  &input ee1=15, ee2=15, degt2=104 /).
C .  Input "ni n" = first and last index of theta_1 (theta_1 =
C .  dth1*(index-1)) ; default 0..360 deg. Output: file fout
C .  (default kr3e_3P.dat) , columns theta_1 (deg) , TDCS.
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
      double precision ZKR,R4P,GAMOWF,qmax,plo
      double precision degts,degps,degt2,degp2
      character*80 fout
      logical lex
      namelist /input/ ee1,ee2,ees,enioni,degts,degps,degt2,degp2,dth1
     $ ,iwf,izv,lmax,dz,ndr,nlr,ilog,plo,dp1,dp2,dp3,dpf,np1,np2,np3
     $ ,np4,nt,nh,fout
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
C .   (interval edges dz*(i/ndr)**2 : dense near the nucleus, where
C .   Z(r) falls from 36)
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
C .  angles in degrees : scattered (theta_s , phi_s) , electron 2
C .  (theta_2 , phi_2) ; electron 1 is scanned in the plane phi_1 = 0
C .  with the step dth1
      degts=1.d0
      degps=180.d0
      degt2=104.d0
      degp2=0.d0
      dth1=5.d0
      fout='kr3e_3P.dat'
C .  asymptotic charge seen by the ejected electrons ; C12 constants
      az=2.d0
      abeta=0.005d0
      az12=0.5d0

C .  p GRID : |p| in 4 intervals , nt x nh angles. ilog=1 : the first
C .  interval [0,dp1] is Gauss-Legendre in ln(p) on [plo,dp1] (+ 4 points
C .  on [0,plo]) : the Fourier transform of C12 oscillates as
C .  exp(2i eta12 ln p) when p -> 0 (ilog=0 : uniform, as the original)
      ilog=1
      plo=1.d-4
      dp0=0.d0 ; dp1=0.7d0 ; dp2=1.6d0 ; dp3=2.1d0 ; dpf=10.d0
      np1=160 ; np2=120 ; np3=40 ; np4=5
      nt=80
      nh=80

C .  OPTIONAL INPUT FILE kr3e.inp (namelist &input ... /) : any of the
C .  variables above can be changed there without recompiling
      inquire(file='kr3e.inp',exist=lex)
      if (lex) then
         open(unit=7,file='kr3e.inp',status='old')
         read(7,nml=input)
         close(7)
         write(*,*) 'Parameters read from kr3e.inp'
      endif
      if (np1+np2+np3+np4+4.gt.400) stop 'np1+np2+np3+np4 > 396'
      if (nt.gt.200.or.nh.gt.200) stop 'nt or nh > 200'
      if (lmax.gt.lmx-1) stop 'lmax > 29'
      eei=ees+ee1+ee2+enioni
      ets=degts*arad
      his=degps*arad
      et2=degt2*arad
      hi2=degp2*arad
      write(*,'(A,4F10.3,A)') ' Ei, Es, E1, E2 = ',eei,ees,ee1,ee2,' eV'
      write(*,'(A,4F9.2,A)') ' theta_s, phi_s, theta_2, phi_2 = '
     $ ,degts,degps,degt2,degp2,' deg'

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
      if (ilog.eq.1) then
         call PINT(dp0,plo,4,xx,ww,np,dxp,dwp)
         call PINTL(plo,dp1,np1,xx,ww,np,dxp,dwp)
      else
         call PINT(dp0,dp1,np1,xx,ww,np,dxp,dwp)
      endif
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
            call GAULEG(dz*(dble(id-1)/ndr)**2,dz*(dble(id)/ndr)**2
     $       ,xx,ww,nlr,1d-15)
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
      n=nint(360.d0/dth1)+1
998   continue
      open(unit=4,file=trim(fout),status='unknown')

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
      write(4,'(F10.3,2X,ES18.10E3)') at1,atfi
      call flush(4)
      write(*,'(A,F8.2,A,ES16.8E3)') ' theta1 = ',at1,'  TDCS = ',atfi
 800  continue
      close(unit=4)
      end
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  KRWF : Slater data of Kr (ion Kr2+ 4p^4 occupations) in /wfn/  .
C .   l = 1 (s) , 2 (p) , 3 (d) ; for each l : nb basis functions    .
C .   (n , zeta) and the coefficients co(l,io,j) of norb orbitals.   .
C .   iwf = 1 : Clementi tables 4.1-4.3 ; 2 : data of the original   .
C .   code. Every orbital is normalised (norm printed).               .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      SUBROUTINE KRWF(iwf)
      implicit none
      integer iwf,nb,nn,nocc,norb,l,io,i,j
      double precision zt,co,s,SOVL
      common /wfn/ nb(3),nn(3,12),zt(3,12),co(3,5,12),nocc(3,5)
     $ ,norb(3)
      character*3 lab(3,5)
      integer ns(10),np(9),nd(5)
      double precision zs1(10),cs1(4,10),zp1(9),cp1(3,9),zd1(5),cd1(5)
      double precision zs2(10),cs2(4,10),zp2(9),cp2(3,9),zd2(5),cd2(5)
      data lab /'1s','2p','3d','2s','3p','  ','3s','4p','  ','4s'
     $ ,'  ','  ','  ','  ','  '/
      data ns /1,1,2,2,3,3,3,4,4,4/
      data np /2,2,3,3,3,4,4,4,4/
      data nd /3,3,3,3,3/
C .  Clementi tables 4.1 (s), 4.2 (p), 4.3 (d) ; 3s row 9 : 0.00218
      data zs1 /32.83510d0,40.94470d0,27.45800d0,16.06600d0,14.29620d0
     $ ,9.10937d0,6.37181d0,3.84546d0,2.57902d0,1.77192d0/
      data cs1 /0.71521d0,0.38139d0,-0.14543d0,-0.04349d0
     $ ,0.29911d0,-0.00123d0,0.00181d0,-0.00148d0
     $ ,-0.01854d0,0.17175d0,-0.09037d0,-0.03219d0
     $ ,0.00897d0,-1.07106d0,0.49528d0,0.16451d0
     $ ,-0.00464d0,-0.14131d0,0.25451d0,0.08852d0
     $ ,0.00190d0,-0.01520d0,-0.48504d0,-0.16671d0
     $ ,-0.00088d0,0.00401d0,-0.75593d0,-0.33291d0
     $ ,0.00026d0,-0.00122d0,-0.01203d0,0.46913d0
     $ ,-0.00018d0,0.00562d0,0.00218d0,0.55106d0
     $ ,0.00006d0,-0.00031d0,-0.00085d0,0.13572d0/
      data zp1 /17.03660d0,26.04380d0,15.51000d0,9.49403d0,6.57275d0
     $ ,5.38507d0,3.15603d0,2.02966d0,1.42733d0/
      data cp1 /0.72322d0,0.30185d0,0.08488d0
     $ ,0.06774d0,0.02508d0,0.00571d0
     $ ,0.22056d0,0.15903d0,0.04169d0
     $ ,0.04478d0,-0.28475d0,-0.07425d0
     $ ,-0.01672d0,-0.76440d0,-0.26866d0
     $ ,0.00609d0,-0.10700d0,0.01341d0
     $ ,-0.00195d0,-0.00562d0,0.51241d0
     $ ,0.00111d0,0.00137d0,0.42557d0
     $ ,-0.00040d0,-0.00053d0,0.18141d0/
      data zd1 /5.30650d0,3.36240d0,7.94963d0,10.35430d0,17.11420d0/
      data cd1 /0.50854d0,0.11070d0,0.24778d0,0.20584d0,0.02863d0/
C .  data of DI-Kr-4p-BBK-zv-3P.for (az1s..az4s, az2p..az4p, az3d)
      data zs2 /34.2726d0,40.40d0,26.8157d0,16.4472d0,13.9631d0
     $ ,8.89765d0,6.22432d0,3.83098d0,2.56514d0,1.73148d0/
      data cs2 /0.73837d0,-0.4169d0,-0.1555d0,-0.04719d0
     $ ,0.25206d0,0.0649d0,0.01678d0,0.0035d0
     $ ,0.02023d0,-0.21848d0,-0.11584d0,-0.03974d0
     $ ,-0.00881d0,1.08043d0,0.50977d0,0.16816d0
     $ ,0.00593d0,0.18819d0,0.27066d0,0.09647d0
     $ ,-0.00379d0,0.00226d0,-0.5713d0,-0.20412d0
     $ ,0.00191d0,0.00214d0,-0.68418d0,-0.30194d0
     $ ,-0.00064d0,-0.00052d0,-0.00538d0,0.47727d0
     $ ,0.00044d0,0.00028d0,-0.00082d0,0.55546d0
     $ ,-0.00015d0,-0.00010d0,0.00002d0,0.12492d0/
      data zp2 /16.3363d0,24.8543d0,15.4451d0,9.18326d0,6.4205d0
     $ ,5.25899d0,3.07835d0,1.9621d0,1.422d0/
      data cp2 /0.74427d0,-0.316d0,0.09018d0
     $ ,0.09648d0,-0.03474d0,0.00779d0
     $ ,0.16587d0,-0.12839d0,0.03172d0
     $ ,0.04848d0,0.3318d0,-0.08481d0
     $ ,-0.0215d0,0.73159d0,-0.26625d0
     $ ,0.00775d0,0.0878d0,0.0349d0
     $ ,-0.00253d0,0.00399d0,0.54322d0
     $ ,0.00152d0,-0.00065d0,0.40037d0
     $ ,-0.0006d0,0.00026d0,0.16488d0/
      data zd2 /5.19134d0,3.34983d0,7.68953d0,10.3162d0,16.7379d0/
      data cd2 /0.47749d0,0.10325d0,0.28079d0,0.20824d0,0.03159d0/

      nb(1)=10 ; nb(2)=9 ; nb(3)=5
      norb(1)=4 ; norb(2)=3 ; norb(3)=1
C .  occupations of the ion Kr2+ : 1s2 2s2 3s2 4s2 2p6 3p6 4p4 3d10
      nocc(1,1)=2 ; nocc(1,2)=2 ; nocc(1,3)=2 ; nocc(1,4)=2
      nocc(2,1)=6 ; nocc(2,2)=6 ; nocc(2,3)=4
      nocc(3,1)=10
      do j=1,10
         nn(1,j)=ns(j)
      enddo
      do j=1,9
         nn(2,j)=np(j)
      enddo
      do j=1,5
         nn(3,j)=nd(j)
      enddo
      if (iwf.eq.1) then
         write(*,*) 'Wave functions: Clementi tables 4.1-4.3'
         do j=1,10
            zt(1,j)=zs1(j)
            do io=1,4
               co(1,io,j)=cs1(io,j)
            enddo
         enddo
         do j=1,9
            zt(2,j)=zp1(j)
            do io=1,3
               co(2,io,j)=cp1(io,j)
            enddo
         enddo
         do j=1,5
            zt(3,j)=zd1(j)
            co(3,1,j)=cd1(j)
         enddo
      else
         write(*,*) 'Wave functions: data of the original code'
         do j=1,10
            zt(1,j)=zs2(j)
            do io=1,4
               co(1,io,j)=cs2(io,j)
            enddo
         enddo
         do j=1,9
            zt(2,j)=zp2(j)
            do io=1,3
               co(2,io,j)=cp2(io,j)
            enddo
         enddo
         do j=1,5
            zt(3,j)=zd2(j)
            co(3,1,j)=cd2(j)
         enddo
      endif
C .  normalisation of every orbital
      do l=1,3
      do io=1,norb(l)
         s=0.d0
         do i=1,nb(l)
         do j=1,nb(l)
            s=s+co(l,io,i)*co(l,io,j)
     $        *SOVL(nn(l,i),zt(l,i),nn(l,j),zt(l,j))
         enddo
         enddo
         write(*,'(A,A3,A,F9.5)') '   norm of ',lab(l,io),' = ',s
         do j=1,nb(l)
            co(l,io,j)=co(l,io,j)/dsqrt(s)
         enddo
      enddo
      enddo
      return
      end

C .  overlap of two normalised radial Slater functions
      FUNCTION SOVL(n1,z1,n2,z2)
      implicit none
      integer n1,n2
      double precision SOVL,z1,z2,fa,fac
      SOVL=fa(z1,n1-1)*fa(z2,n2-1)*fac(n1+n2)/(z1+z2)**(n1+n2+1)
      end

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  ZKR : Z(r) = 36 - sum_orbitals occ * r*V_Hartree of the orbital  .
C .  density (monopole) ; pair (i,j) of Slater functions of the same .
C .  orbital : p!/c^(p+1) - J_p(c,r) + r J_(p-1)(c,r) , p = ni+nj ,   .
C .  c = zi+zj (as az1s..az3d of the original code).                 .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      FUNCTION ZKR(ax)
      implicit none
      integer nb,nn,nocc,norb,l,io,i,j,ip
      double precision ZKR,ax,zt,co,bz,bv,c,fa,fac,avJn
      common /wfn/ nb(3),nn(3,12),zt(3,12),co(3,5,12),nocc(3,5)
     $ ,norb(3)
      bz=0.d0
      do l=1,3
      do io=1,norb(l)
         do i=1,nb(l)
         do j=1,nb(l)
            ip=nn(l,i)+nn(l,j)
            c=zt(l,i)+zt(l,j)
            bv=fac(ip)/(c**(ip+1))-avJn(c,ax,ip)+ax*avJn(c,ax,ip-1)
            bz=bz+nocc(l,io)*bv*co(l,io,i)*co(l,io,j)
     $        *fa(zt(l,i),nn(l,i)-1)*fa(zt(l,j),nn(l,j)-1)
         enddo
         enddo
      enddo
      enddo
      ZKR=36.d0-bz
      end

C .  R4P : radial function of the ionized 4p orbital
      FUNCTION R4P(ar)
      implicit none
      integer nb,nn,nocc,norb,j
      double precision R4P,ar,zt,co,fa
      common /wfn/ nb(3),nn(3,12),zt(3,12),co(3,5,12),nocc(3,5)
     $ ,norb(3)
      R4P=0.d0
      do j=1,nb(2)
         R4P=R4P+co(2,3,j)*fa(zt(2,j),nn(2,j)-1)*ar**(nn(2,j)-1)
     $    *dexp(-zt(2,j)*ar)
      enddo
      end

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  AN4P : analytic part for the 4p orbital, m = -1, 0, 1           .
C .  (= qkr4p of the original code, with the coefficients of /wfn/) .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      SUBROUTINE AN4P(az,adx,ady,adz,akx,aky,akz,res)
      implicit none
      integer nb,nn,nocc,norb,j,m
      double precision az,adx,ady,adz,akx,aky,akz,zt,co
      complex*16 res(-1:1),up(-1:1)
      common /wfn/ nb(3),nn(3,12),zt(3,12),co(3,5,12),nocc(3,5)
     $ ,norb(3)
      do m=-1,1
         res(m)=(0.d0,0.d0)
      enddo
      do j=1,nb(2)
         call TCP3(zt(2,j),az,adx,ady,adz,akx,aky,akz,nn(2,j),up)
         do m=-1,1
            res(m)=res(m)+co(2,3,j)*up(m)
         enddo
      enddo
      return
      end

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  TCP3 : tcp of the original code for c = -1, 0, 1 in one call.   .
C .  The products of q-functions do not depend on c : they are       .
C .  computed once ; only va2, vb2 and app change with c.            .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      SUBROUTINE TCP3(a,az,adx,ady,adz,akx,aky,akz,n,res)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      integer n,m
      complex*16 res(-1:1)
      double precision a,az,adx,ady,adz,akx,aky,akz
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      al1=akx*akx+aky*aky+akz*akz
      ak1=dsqrt(al1)
      aq=adx*adx+ady*ady+adz*adz
      al3=aq+al1-2.*(adx*akx+ady*aky+adz*akz)
      qz=az*xi/ak1
      qz5=2.*xi*ak1
      r=a
      x=r-xi*ak1
      qa=al3+r*r
      qab=aq+r*r-al1-2.*xi*r*ak1
      qb=qab-qa
      aqkz=adz-akz
      q01x=akx+xi*aky
      qk1x=-adx-xi*ady+akx+xi*aky
      q01y=-akx+xi*aky
      qk1y=adx-xi*ady-akx+xi*aky
C .  c-independent products
      if (n.eq.2) then
         up20=q1(qa,qz-2.,r)*q0(qab,-qz-1.)
         up21=q1(qab,-qz-1.,x)*q0(qa,qz-2.)
         w0=up20+up21
         w1=q0(qa,qz-2.)*q0(qab,-qz-1.)
      else if (n.eq.3) then
         up30=q2(qa,qz-2.,r)*q0(qab,-qz-1.)
         up31=q2(qab,-qz-1.,x)*q0(qa,qz-2.)
         up32=2.*q1(qa,qz-2.,r)*q1(qab,-qz-1.,x)
         w0=up30+up31+up32
         up34=q1(qa,qz-2.,r)*q0(qab,-qz-1.)
         up35=q1(qab,-qz-1.,x)*q0(qa,qz-2.)
         w1=(up34+up35)*2.
         w2=q0(qa,qz-2.)*q0(qab,-qz-1.)
      else
         up40=q3(qa,qz-2.,r)*q0(qab,-qz-1.)
         up41=q3(qab,-qz-1.,x)*q0(qa,qz-2.)
         up42=3.*q2(qa,qz-2.,r)*q1(qab,-qz-1.,x)
         up43=3.*q2(qab,-qz-1.,x)*q1(qa,qz-2.,r)
         w0=up40+up41+up42+up43
         up45=q2(qa,qz-2.,r)*q0(qab,-qz-1.)
         up46=q2(qab,-qz-1.,x)*q0(qa,qz-2.)
         up47=2.*q1(qa,qz-2.,r)*q1(qab,-qz-1.,x)
         w1=(up45+up46+up47)*3.
         up490=q1(qa,qz-2.,r)*q0(qab,-qz-1.)
         up491=q1(qab,-qz-1.,x)*q0(qa,qz-2.)
         w2=(up490+up491)*3.
      endif
      do m=-1,1
         if (m.eq.0) then
            app=6.**.5/api
            va2=aqkz
            vb2=akz
         else if (m.eq.1) then
            app=3.**.5/api
            va2=qk1x
            vb2=-q01x
         else
            app=3.**.5/api
            va2=qk1y
            vb2=-q01y
         endif
         if (n.eq.2) then
            res(m)=-(w0*tp0(qa,qb,qz,va2,vb2)
     $       +w1*tp1(qz,va2,vb2,r,qz5))*xi*app
         else if (n.eq.3) then
            res(m)=(w0*tp0(qa,qb,qz,va2,vb2)+w1*tp1(qz,va2,vb2,r,qz5)
     $       +w2*tp2(qz,va2,vb2))*app*xi
         else
            res(m)=-(w0*tp0(qa,qb,qz,va2,vb2)+w1*tp1(qz,va2,vb2,r,qz5)
     $       +w2*tp2(qz,va2,vb2))*xi*app
         endif
         res(m)=res(m)*fa(a,n-1)
      enddo
      return
      end

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  PWDL : Legendre components d_l(r) of D(r,u) = conj(phi_Z(r) -    .
C .  phi_Z=az) , u = cos(k^,r^) , phi as qf1 of the original code :   .
C .  phi_Z(r) = exp(ikru) exp(pi Z/2k) G(1+iZ/k) 1F1(-iZ/k,1,-ikr(1+u))
C .             / [exp(pi az/2k) G(1+i az/k)]                         .
C .  phi_az   = exp(ikru) 1F1(-i az/k,1,-ikr(1+u))                    .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      SUBROUTINE PWDL(ak,az,nr,rgr,lmax,dl)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      parameter (lmx=30,nux=100)
      integer nr,lmax,l,i,iu,nu
      double precision rgr(nr),du(nux),dwu(nux),pl(0:lmx),ZKR
      complex*16 dl(400,0:lmx),qgr,qdd
      api=dacos(-1.d0)
      xi=(0.d0,1.d0)
      nu=min(nux,2*lmax+24)
      call GAULEG(-1.d0,1.d0,du,dwu,nu,1d-15)
      alp=az/ak
      sq=dexp(api*alp/2.d0)*uc(1.d0+xi*alp)
      do i=1,nr
         ar=rgr(i)
         alr=ZKR(ar)/ak
         qgr=dexp(api*alr/2.d0)*uc(1.d0+xi*alr)/sq
         do l=0,lmax
            dl(i,l)=(0.d0,0.d0)
         enddo
         do iu=1,nu
            bu=du(iu)
            zc=-xi*ak*ar*(1.d0+bu)
            qe=cdexp(xi*ak*ar*bu)
            qdd=dconjg(qe*(qgr*uf(-xi*alr,(1.d0,0.d0),zc)
     $        -uf(-xi*alp,(1.d0,0.d0),zc)))
            pl(0)=1.d0
            pl(1)=bu
            do l=2,lmax
               pl(l)=((2*l-1)*bu*pl(l-1)-(l-1)*pl(l-2))/l
            enddo
            do l=0,lmax
               dl(i,l)=dl(i,l)+0.5d0*(2*l+1)*dwu(iu)*pl(l)*qdd
            enddo
         enddo
      enddo
      return
      end

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  JTAB : J_lL(q) = int d_l(r) R4p(r) j_L(qr) r^2 dr , L = l-1 (s=1) .
C .  and l+1 (s=2) , on q = iq*hq , iq = 0..nq                        .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      SUBROUTINE JTAB(nr,rgr,wgr,r4g,dl,lmax,nq,hq,jt)
      implicit none
      integer lmx,nr,lmax,nq,iq,i,l
      parameter (lmx=30)
      double precision rgr(nr),wgr(nr),r4g(nr),hq,bj(0:lmx+1),wr
      complex*16 dl(400,0:lmx),jt(0:4000,0:lmx,2)
      do iq=0,nq
         do l=0,lmax
            jt(iq,l,1)=(0.d0,0.d0)
            jt(iq,l,2)=(0.d0,0.d0)
         enddo
         do i=1,nr
            call SBJ(iq*hq*rgr(i),lmax+1,bj)
            wr=wgr(i)*r4g(i)*rgr(i)*rgr(i)
            do l=0,lmax
               if (l.ge.1) jt(iq,l,1)=jt(iq,l,1)+wr*dl(i,l)*bj(l-1)
               jt(iq,l,2)=jt(iq,l,2)+wr*dl(i,l)*bj(l+1)
            enddo
         enddo
      enddo
      return
      end

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  CCOEF : angular coefficients for one electron direction k^ :    .
C .  cc(m,l,s,M) = (2pi)^-1.5 (4pi)^2/(2l+1) i^L Y*_l,mu(k^)         .
C .               <Y_LM Y_1m Y_l,mu> ,  mu = -m-M , L = l-/+1        .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      SUBROUTINE CCOEF(akx,aky,akz,lmax,cc)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      parameter (lmx=30)
      integer lmax,l,is,ll,mm,m,mu
      double precision W3J
      complex*16 cc(-1:1,0:lmx,2,-lmx-1:lmx+1)
      api=dacos(-1.d0)
      xi=(0.d0,1.d0)
      ak=dsqrt(akx**2+aky**2+akz**2)
      ctk=akz/ak
      bphk=datan2(aky,akx)
      cc=(0.d0,0.d0)
      do l=0,lmax
      do is=1,2
         ll=l-1
         if (is.eq.2) ll=l+1
         if (ll.lt.0) cycle
         do m=-1,1
         do mm=-ll,ll
            mu=-m-mm
            if (abs(mu).gt.l) cycle
            g=dsqrt((2*ll+1)*3.d0*(2*l+1)/(4.d0*api))
     $       *W3J(ll,1,l,0,0,0)*W3J(ll,1,l,mm,m,mu)
            if (g.eq.0.d0) cycle
            cc(m,l,is,mm)=(4.d0*api)**2/(2*l+1)*xi**ll
     $       *dconjg(ylm(l,dble(mu),ctk,bphk))*g/dsqrt((2.d0*api)**3)
         enddo
         enddo
      enddo
      enddo
      return
      end

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  PWQ3 : Z(r) correction for m = -1,0,1 at the vector q :          .
C .  sum_{l,s,M} cc(m,l,s,M) J_lL(|q|) Y*_LM(q^)   (added to res)    .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      SUBROUTINE PWQ3(qx,qy,qz,lmax,cc,jt,nq,hq,res)
      implicit none
      integer lmx,lmax,nq,l,is,ll,mm,m,i0,k
      parameter (lmx=30)
      double precision qx,qy,qz,hq,aq,ct,ph,t,cl(4)
      complex*16 cc(-1:1,0:lmx,2,-lmx-1:lmx+1),jt(0:4000,0:lmx,2)
      complex*16 res(-1:1),y(0:lmx+1,-lmx-1:lmx+1),ajq,yc
      aq=dsqrt(qx*qx+qy*qy+qz*qz)
      if (aq.gt.1d-12) then
         ct=qz/aq
         ph=datan2(qy,qx)
      else
         ct=1.d0
         ph=0.d0
      endif
      call YALL(ct,ph,lmax+1,y)
C .  4-point Lagrange interpolation in |q|
      t=aq/hq
      i0=int(t)-1
      if (i0.lt.0) i0=0
      if (i0.gt.nq-3) i0=nq-3
      t=t-i0
      cl(1)=-(t-1.d0)*(t-2.d0)*(t-3.d0)/6.d0
      cl(2)=t*(t-2.d0)*(t-3.d0)/2.d0
      cl(3)=-t*(t-1.d0)*(t-3.d0)/2.d0
      cl(4)=t*(t-1.d0)*(t-2.d0)/6.d0
      do l=0,lmax
      do is=1,2
         ll=l-1
         if (is.eq.2) ll=l+1
         if (ll.lt.0) cycle
         ajq=(0.d0,0.d0)
         do k=1,4
            ajq=ajq+cl(k)*jt(i0+k-1,l,is)
         enddo
         do mm=-ll,ll
            yc=dconjg(y(ll,mm))*ajq
            do m=-1,1
               res(m)=res(m)+cc(m,l,is,mm)*yc
            enddo
         enddo
      enddo
      enddo
      return
      end

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  YALL : Y_LM(theta,phi) for all L <= lmax (Condon-Shortley phase, .
C .  same convention as ylm) , normalised Legendre recurrences        .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      SUBROUTINE YALL(ct,ph,lmax,y)
      implicit none
      integer lmx,lmax,l,m
      parameter (lmx=30)
      double precision ct,ph,st,p(0:lmx+1,0:lmx+1),api,a,b
      complex*16 y(0:lmx+1,-lmx-1:lmx+1),e
      api=dacos(-1.d0)
      st=dsqrt(dmax1(0.d0,1.d0-ct*ct))
      p(0,0)=1.d0/dsqrt(4.d0*api)
      do m=1,lmax
         p(m,m)=-dsqrt((2*m+1)/(2.d0*m))*st*p(m-1,m-1)
      enddo
      do m=0,lmax-1
         p(m+1,m)=dsqrt(2.d0*m+3.d0)*ct*p(m,m)
      enddo
      do m=0,lmax
         do l=m+2,lmax
            a=dsqrt((4.d0*l*l-1.d0)/(dble(l)*l-dble(m)*m))
            b=dsqrt((dble(l-1)**2-dble(m)*m)/(4.d0*(l-1)**2-1.d0))
            p(l,m)=a*(ct*p(l-1,m)-b*p(l-2,m))
         enddo
      enddo
      do l=0,lmax
         do m=0,l
            e=cdexp(dcmplx(0.d0,m*ph))
            y(l,m)=p(l,m)*e
            if (m.gt.0) y(l,-m)=(-1)**m*dconjg(y(l,m))
         enddo
      enddo
      return
      end

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  SBJ : spherical Bessel functions j_0..j_lmax(x)                  .
C .  series (x<0.05) , upward (x>lmax) , Miller downward otherwise    .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      SUBROUTINE SBJ(x,lmax,bj)
      implicit none
      integer lmax,l,ls
      double precision x,bj(0:lmax),df,t,s,b0,b1,b2,bt(0:200)
      if (x.lt.0.05d0) then
         df=1.d0
         t=1.d0
         do l=0,lmax
            if (l.gt.0) then
               df=df*(2*l+1)
               t=t*x
            endif
            bj(l)=t/df*(1.d0-x*x/(2.d0*(2*l+3)))
         enddo
         return
      endif
      if (x.gt.dble(lmax)) then
         bj(0)=dsin(x)/x
         if (lmax.ge.1) bj(1)=dsin(x)/(x*x)-dcos(x)/x
         do l=2,lmax
            bj(l)=(2*l-1)/x*bj(l-1)-bj(l-2)
         enddo
         return
      endif
      ls=lmax+20+int(x)
      if (ls.gt.200) ls=200
      b2=0.d0
      b1=1.d-30
      do l=ls,0,-1
         b0=(2*l+3)/x*b1-b2
         if (l.le.lmax) bt(l)=b0
         if (l.eq.0) exit
         b2=b1
         b1=b0
         if (dabs(b1).gt.1.d250) then
            b1=b1*1.d-250
            b2=b2*1.d-250
            bt(l:lmax)=bt(l:lmax)*1.d-250
         endif
      enddo
C .  last b0 is j_0 (unnormalised) ; bt(l) holds j_l
      s=dsin(x)/x
      t=bt(0)
      if (dabs(s).lt.1d-3) then
         s=dsin(x)/(x*x)-dcos(x)/x
         t=bt(1)
      endif
      do l=0,lmax
         bj(l)=bt(l)*s/t
      enddo
      return
      end

C .  Wigner 3j symbol (Racah formula), integer arguments
      FUNCTION W3J(j1,j2,j3,m1,m2,m3)
      implicit none
      integer j1,j2,j3,m1,m2,m3,k,kmin,kmax
      double precision W3J,fac,s,t
      W3J=0.d0
      if (m1+m2+m3.ne.0) return
      if (j3.lt.abs(j1-j2).or.j3.gt.j1+j2) return
      if (abs(m1).gt.j1.or.abs(m2).gt.j2.or.abs(m3).gt.j3) return
      kmin=max(0,j2-j3-m1,j1-j3+m2)
      kmax=min(j1+j2-j3,j1-m1,j2+m2)
      s=0.d0
      do k=kmin,kmax
         t=fac(k)*fac(j1+j2-j3-k)*fac(j1-m1-k)*fac(j2+m2-k)
     $    *fac(j3-j2+m1+k)*fac(j3-j1-m2+k)
         s=s+(-1)**k/t
      enddo
      W3J=(-1)**(j1-j2-m3)*dsqrt(fac(j1+j2-j3)*fac(j1-j2+j3)
     $ *fac(-j1+j2+j3)/fac(j1+j2+j3+1)*fac(j1+m1)*fac(j1-m1)
     $ *fac(j2+m2)*fac(j2-m2)*fac(j3+m3)*fac(j3-m3))*s
      end

C .  t(m1,m2) = q1t(m1) q2s(m2) + q1s(m1) q2t(m2) - 2 q1s(m1) q2s(m2)
      FUNCTION TPAIR(m1,m2,q1t,q1s,q2t,q2s)
      implicit none
      integer m1,m2
      complex*16 TPAIR,q1t(-1:1),q1s(-1:1),q2t(-1:1),q2s(-1:1)
      TPAIR=q1t(m1)*q2s(m2)+q1s(m1)*q2t(m2)-2.d0*q1s(m1)*q2s(m2)
      end

C .  Gamow factor 2 pi eta / (1 - exp(-2 pi eta)) , eta = Z/k
      FUNCTION GAMOWF(e)
      implicit none
      double precision GAMOWF,e,api
      api=dacos(-1.d0)
      if (dabs(e).lt.1d-7) then
         GAMOWF=1.d0
      else
         GAMOWF=2.d0*api*e/(1.d0-dexp(-2.d0*api*e))
      endif
      end

C .  append n Gauss-Legendre points of [a,b] to (dxp,dwp)
      SUBROUTINE PINT(a,b,n,xx,ww,np,dxp,dwp)
      implicit none
      integer n,np,i
      double precision a,b,xx(*),ww(*),dxp(*),dwp(*)
      call GAULEG(a,b,xx,ww,n,2d-16)
      do i=1,n
         dxp(np+i)=xx(i)
         dwp(np+i)=ww(i)
      enddo
      np=np+n
      return
      end

C .  append n Gauss-Legendre points in ln(p) on [a,b] (a > 0)
      SUBROUTINE PINTL(a,b,n,xx,ww,np,dxp,dwp)
      implicit none
      integer n,np,i
      double precision a,b,xx(*),ww(*),dxp(*),dwp(*)
      call GAULEG(dlog(a),dlog(b),xx,ww,n,2d-16)
      do i=1,n
         dxp(np+i)=dexp(xx(i))
         dwp(np+i)=ww(i)*dexp(xx(i))
      enddo
      np=np+n
      return
      end
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  ROUTINES OF DI-Kr-4p-BBK-zv-3P.for (unchanged)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      function avJn(al,ax,n)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
	alx=al*ax
	as=0.d0
	do m=0,n
	as=as+(alx**m)/fac(m)
	enddo
	avJn=fac(n)*dexp(-alx)*as/(al**(n+1))
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function tcs(a,az,adx,ady,adz,akx,aky,akz,n)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      dimension us(2)
c NON slater
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      aps=1./(2.**.5*api)

      al1=akx*akx+aky*aky+akz*akz
      ak1=dsqrt(al1)
      aq=adx*adx+ady*ady+adz*adz
      al3=aq+al1-2.*(adx*akx+ady*aky+adz*akz)
      qz=az*xi/ak1
      sz=1.-qz

      r=a
      x=r-xi*ak1
      qa=al3+r*r
      qab=aq+r*r-al1-2.*xi*r*ak1
      qb=qab-qa
c      ucd=uc(sz)*cdexp(api*qz/(2.*xi))
       ucd=(1.d0,0.d0)


c 1s
      us10=q1(qa,qz-1.,r)*q0(qab,-qz)
      us11=q0(qa,qz-1.)*q1(qab,-qz,x)
      us(1)=-(us10+us11)*aps*ucd

c 2s
      us20=q2(qa,qz-1.,r)*q0(qab,-qz)
      us21=q2(qab,-qz,x)*q0(qa,qz-1.)
      us22=2.*q1(qa,qz-1.,r)*q1(qab,-qz,x)
      us(2)=(us20+us21+us22)*aps*ucd
      tcs=us(n)*fa(a,n-1)
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function tcp(a,az,adx,ady,adz,akx,aky,akz,n,c)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      dimension up(4)

      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      

      al1=akx*akx+aky*aky+akz*akz
      ak1=dsqrt(al1)
      aq=adx*adx+ady*ady+adz*adz
      al3=aq+al1-2.*(adx*akx+ady*aky+adz*akz)
      qz=az*xi/ak1
      sz=1.-qz

      qz5=2.*xi*ak1
      r=a
      x=r-xi*ak1
      qa=al3+r*r
      qab=aq+r*r-al1-2.*xi*r*ak1
      qb=qab-qa
c      ucd=uc(sz)*cdexp(api*qz/(2.*xi))
       ucd=(1.d0,0.d0)
      aqkz=adz-akz
      q01x=akx+xi*aky
      qk1x=-adx-xi*ady+akx+xi*aky
      q01y=-akx+xi*aky
      qk1y=adx-xi*ady-akx+xi*aky
      akp2=aq-al3-al1

     
      if(c.eq.0.) then
      app=6.**.5/api
      va2=aqkz
      vb2=akz
      goto 272   
      endif  
      if(c.eq.1.) then 
      app=3.**.5/api
      va2=qk1x
      vb2=-q01x
      goto 272   
      endif  
      if(c.eq.(-1.)) then 
      app=3.**.5/api
      va2=qk1y
      vb2=-q01y
      goto 272   
      endif  
 
 272  if (n.eq.2) goto 991
      if (n.eq.3) goto 992
      if (n.eq.4) goto 993

c 2p    
 991  up20=q1(qa,qz-2.,r)*q0(qab,-qz-1.)
      up21=q1(qab,-qz-1.,x)*q0(qa,qz-2.)
      up22=(up20+up21)*tp0(qa,qb,qz,va2,vb2)
      up23=q0(qa,qz-2.)*q0(qab,-qz-1.)
      up24=up23*tp1(qz,va2,vb2,r,qz5)
      up(2)=-(up22+up24)*xi*app*ucd
	goto 9999
c 3p
 992     up30=q2(qa,qz-2.,r)*q0(qab,-qz-1.)
      up31=q2(qab,-qz-1.,x)*q0(qa,qz-2.)
      up32=2.*q1(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up33=(up30+up31+up32)*tp0(qa,qb,qz,va2,vb2)
      up34=q1(qa,qz-2.,r)*q0(qab,-qz-1.)
      up35=q1(qab,-qz-1.,x)*q0(qa,qz-2.)
      up36=(up34+up35)*2.*tp1(qz,va2,vb2,r,qz5)
      up37=q0(qa,qz-2.)*q0(qab,-qz-1.)*tp2(qz,va2,vb2)
      up(3)=(up33+up36+up37)*app*xi*ucd
	goto 9999

c 4p
993      up40=q3(qa,qz-2.,r)*q0(qab,-qz-1.)
      up41=q3(qab,-qz-1.,x)*q0(qa,qz-2.)
      up42=3.*q2(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up43=3.*q2(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up44=(up40+up41+up42+up43)*tp0(qa,qb,qz,va2,vb2)
      up45=q2(qa,qz-2.,r)*q0(qab,-qz-1.)
      up46=q2(qab,-qz-1.,x)*q0(qa,qz-2.)
      up47=2.*q1(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up48=(up45+up46+up47)*3.*tp1(qz,va2,vb2,r,qz5)
      up490=q1(qa,qz-2.,r)*q0(qab,-qz-1.)
      up491=q1(qab,-qz-1.,x)*q0(qa,qz-2.)
      up49=(up490+up491)*3.*tp2(qz,va2,vb2)
      up(4)=-(up44+up48+up49)*xi*app*ucd
	goto 9999

9999	 continue
      tcp=up(n)*fa(a,n-1)
	end
c fonctions annexes pour calculer les transformees fouriercoulomb
      function q0(x,y)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      q0=x**y
      end
      function q1(x,y,z)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      q1=2.*z*y*(x**(y-1.))
      end
      function q2(x,y,z)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      q2=2.*y*(x**(y-2.))*(x+2.*z*z*(y-1.))
      end
      function q3(x,y,z)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      q3=z*y*(y-1.)*x**(y-3.)*(8.*z*z*(y-2.)+12.*x)
      end
      function q4(x,y,z)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      q4=3.*x*x+12.*z*z*(y-2.)*x+4.*z**4.*(y-2.)*(y-3.)
      q4=q4*4.*y*(y-1.)*x**(y-4.)
      end
      function tp0(x,y,z1,z2,z3)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tp0=(x+y*(1.-z1))*z2+x*z1*z3
      end
      function tp1(z1,z2,z3,z4,z5)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tp1=(2.*z4-z5*(1.-z1))*z2+2.*z4*z1*z3
      end
      function tp2(z1,z2,z3)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tp2=2.*(z2+z1*z3)
      end


c fonction de normation des slater
      function fa(a,l1)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      fa=(a+a)**(1.5+l1)/dsqrt(fac(2+l1+l1))
      end

      function fac(n)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      fac=1.
      do 700 i=2,n
      fac=fac*i
 700  continue
      end

      function uf(za,zc,z)
c fonction 1f1  pour z inf a 1
c handbook p.504
c et pour z grand
       implicit double precision (a-h)
       implicit complex*16 (o-z)
       dimension u(291)
       z1=1.
       i00=0
       ic=int(zc)
       ia=int(za)
       ica=ic-ia
       aua=dreal(za)
       bua=dimag(za)
       bva=dabs(bua)
       ina=nint(aua)
       au4=dabs(aua-ina)
       aee=.0001 
cas divergence zc entier et za entier et c-a sup 0
       if (zc.eq.ic.and.ic.le.i00.and.za.eq.ia.and.ica.gt.i00)then
       write(*,*) "impossible,c=-m et a ou b neq -n"
       stop
       endif



cas divergence zc entier negatif ou nul 
       if (zc.eq.ic.and.ic.le.i00) then
       write(*,*) "impossible,c=-m et a ou b neq -n"
       stop
       endif


       if (zc.eq.za) then
       uf=cdexp(z)
       goto 9133
       endif
cas z grand (sup.20)  
       anorm=cdabs(z)
       acut=22.
       if (anorm.gt.acut) then
       uf=usr(za,zc,z)
       goto 9133
       endif 

       zs=0. 
       u(1)=za*z/zc
c       do 9131 i=1,4
       do 9131 i=1,290
       u(i+1)=u(i)*(za+i)*z/((zc+i)*(i+1.))
       zs=zs+u(i)
       am=dreal(zs/1000000.)
       au=dreal(u(i+1))
       bm=dimag(zs/1000000.)
       bu=dimag(u(i+1))
       if (dabs(au).le.dabs(am).and.dabs(bu).le.dabs(bm)) then
       goto 9132
       endif
 9131  continue
 9132  uf=zs+1.
        if (i.ge.90) then
       write(*,*)"attention 90 points utilises pour serie uf",i 
       endif
 9133  end


      function ulog(z)
c fonction log avec complexes
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.00d00)
      api=3.141592654d00
      bx=dreal(z)
      by=dimag(z)
      cn=cdabs(z)
      au0=0.
      if(bx.eq.au0) then
       if(by.lt.au0) then
        ulog=-0.5*xi*api+dlog(cn)
        goto 9301
        else 
        ulog=0.5*xi*api+dlog(cn)
        goto 9301
        endif
        else 
      atet=datan(by/bx)
      if(bx.lt.au0.and.by.ge.au0) then
      atet=api+atet
      goto 9300
      endif
      if(bx.lt.au0.and.by.lt.au0) then
      atet=-(api-atet)
      endif
 9300 ulog=dlog(cn)+xi*atet
 9301 endif
      end

     
	function uga(sz)
cfonction gamma incomplete pour reel de z positif
c handbook p.256 6.1.34  
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      dimension a(19)
      api=3.141592654d00 
      au=2.0
      a(1)=1.0
      a(2)=.5772156649
      a(3)=-.6558780715
      a(4)=-.042002635
      a(5)=.1665386113
      a(6)=-.0421977345
      a(7)=-.0096219715
      a(8)=.0072189432
      a(9)=-.0011651675
      a(10)=-.0002152416
      a(11)=.0001280502
      a(12)=-.0000201348
      a(13)=-.0000012504
      a(14)=.0000011330
      a(15)=-.0000002056
      a(16)=.0000000061
      a(17)=.000000005
      a(18)=-.0000000011
      a(19)=.0000000001
      ann=cdabs(sz)
      if (ann.le.au) then
      us=0.
      do 9010 i=1,19
      ui=i+.0
      sn=a(i)*(sz**ui)
      us=us+sn
 9010 continue
      uga=1./us
      else
      uga=stir(sz)
      endif
      end

      function stir(sz)
c fonction gamma stirling handbook p.257 6.1.37
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00 
      s2=1./(12.*sz)
      s3=1./(288.*sz*sz)
      s4=-139./(51840.*sz*sz*sz)
      s5=-571./(2488320.*sz*sz*sz*sz)
      ss=(1.+s2+s3+s4+s5)*dsqrt(2.*api)
      stir=ss*cdexp((sz-.5)*ulog(sz))*cdexp(-sz)
      end

      function uc(sz)
c fonction gamma complete handbook p.256 6.1.17
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      au=dreal(sz)
      bu=dimag(sz)
      inz=nint(au)
      au8=dabs(au-inz)
      bv=dabs(bu)
      ae8=0.0000001
      i00=0
      au0=0.
      if (bv.lt.ae8.and.au8.lt.ae8.and.inz.le.i00) then
      write(*,*)"impossible pour gamma entier negatif ou nul"
      stop
      endif
      if(au.eq.au0.and.bu.lt.au0) then
      bum=-bu
      ug1=2.*(api/bum)/(dexp(api*bum)-dexp(-api*bum))
      sz2=-sz      
      uc=ug1/uga(sz2)
      go to 9900
      endif

      if(au.lt.au0) then
      ugg=api/(cdsin(api*sz))
      sz1=1.-sz
      uc=ugg/uga(sz1)
      else
      uc=uga(sz)
      endif
 9900 end 

      function qcj(x)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      qcj=x*dconjg(x)
      end


	SUBROUTINE GAULEG(X1,X2,X,W,N,EPS)
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .   P R O G R A M                                                   .
C .                  TO CALCULATE THE NODES AND WEIGHTS OF THE        .
C .                  GAUSS-LEGENDRE QUADRATURE. PUBLISHED IN:         .
C .                  W.H. PRESS, B.F. FLANERY, S.A. TEUKOLSKY,        .
C .                  AND W.T.VETTERLEY, NUMERICAL RECIPES: THE        .
C .                  ART OF SCIENTIFIC COMPUTING (CAMBRIDGE           .
C .                  UNIVERSITY PRESS, CAMBRIDGE, 1986), CHAP. 4      .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
       IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N),W(N)
      DATA DZERO/0.D0/, QUAT/0.25D0/, HALF/0.5D0/, ONE/1.D0/, TWO/2.D0/
      PI = DACOS(-ONE)
      M=(N+1)/2
      XM=HALF*(X2+X1)
      XL=HALF*(X2-X1)
      DO 12 I = 1 , M
        Z=DCOS(PI*(I-QUAT)/(N+HALF))
    1   CONTINUE
          P1=ONE
          P2=DZERO
          DO 11 J = 1 , N
            P3=P2
            P2=P1
            P1=((TWO*J-ONE)*Z*P2-(J-ONE)*P3)/J
   11     CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-ONE)
          Z1=Z
          Z=Z1-P1/PP
        IF (ABS(Z-Z1).GT.EPS) GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=TWO*XL/((ONE-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
   12 CONTINUE
      RETURN				
      END
C
      function usr(za,zc,z)
c fonction 1f1 incomplete pour z grand
c handbook p.508 13.5.1
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      dimension u(101),v(101)


      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      api2=api*.5
      api2m=-api*.5
      api32=3.*api2
      api32m=-api32
      au0=0.
      arz=dreal(z)
      aiz=dimag(z)
      if (aiz.gt.au0.and.arz.eq.au0) then
      atat=api/2.
      goto 9520
      endif
      if (aiz.lt.au0.and.arz.eq.au0) then
      atat=-api/2.
      goto 9520
      endif

      atat=datan(aiz/arz)
      
 9520 zrr=0.


      zss=0.
      u(1)=za*(1.+za-zc)/(-z)
      v(1)=(zc-za)*(1.-za)/(z)

      do 9521 i=1,100
      u(i+1)=u(i)*(za+i)*(1.+za-zc+i)/((i+1.)*(-z))
      zrr=zrr+u(i)
      arm=dreal(zrr/10000.)
      aru=dreal(u(i+1))
      brm=dimag(zrr/10000.)
      bru=dimag(u(i+1))
      if (dabs(aru).le.dabs(arm).and.dabs(bru).le.dabs(brm)) then
      goto 9522
      endif
 9521 continue
 9522 usr1=1.+zrr

      if (i.ge.90) then
      write (*,*) "attention 90 pts pour serie usr1",i 
      endif

      do 9531 j=1,100
      v(j+1)=v(j)*(zc-za+j)*(1.-za+j)/((j+1.)*z)
      zss=zss+v(j)
      asm=dreal(zss/10000.)
      asu=dreal(v(j+1))
      bsm=dimag(zss/10000.)
      bsu=dimag(v(j+1))
      if (dabs(asu).le.dabs(asm).and.dabs(bsu).le.dabs(bsm)) then
      goto 9532
      endif
 9531 continue
 9532 usr2=1.+zss



      if (j.ge.110) then
      write (*,*) "attention 90 pts pour serie usr2",j 
      endif

      if (atat.gt.api2m.and.atat.le.api32) then
      usr3=usr1*(cdexp(xi*api*za))*z**(-za)*uc(zc)/(uc(zc-za))
      goto 9533
      endif
      if (atat.gt.api32m.and.atat.le.api2m) then
      usr3=usr1*(cdexp(-xi*api*za))*z**(-za)*uc(zc)/(uc(zc-za))
      goto 9533
      endif
      if (atat.eq.api32m.or.atat.eq.api32) then
      write (*,*) "attention arctang non defini pour usr"
      endif
 9533 usr6=usr2*uc(zc)*(cdexp(z))*z**(za-zc)/(uc(za))
      usr=usr3+usr6
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       function algndr(l,m,ax)
c fonction polynomes de legendre associ�s aux harmoniques sph�riques  
c ax=cos (theta)     
      implicit double precision (a-h)
      bx=dabs(ax)
      a1=1.d00
c      if(m.lt.0.or.m.gt.l.or.bx.gt.a1) then
c      write(6,*) 'PB sur le polynome de legendre l,m,ax'
c      endif
      amm=1.d00
      if(m.gt.0) then 
      asomx2=dsqrt((1.d00-ax)*(1.d00+ax))
      fact=1.d00
      do 11 i=1,m
      amm=amm*fact*asomx2
      fact=fact+2.d00
  11  continue
      endif
      if(l.eq.m) then
      algndr=amm
      else
      ammp1=ax*(2*m+1)*amm
      if(l.eq.m+1) then
      algndr=ammp1
      else
      do 12 ll=m+2,l
      all=(ax*(2*ll-1)*ammp1-(ll+m-1)*amm)/(ll-m)
      amm=ammp1
      ammp1=all
  12  continue
      algndr=all
      endif
      endif
      return
      end

      function ylm(l,am,ax,ay)
c fonction y(l,m) harmoniques sph�riques 
c ax=cos (theta) 
c ay=phi      
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00 
      xi=(0.0d00,1.00d00) 
      m=nint(am)
c     write (*,*) m,am
c     stop
      if(l.eq.0) then
      if(m.ne.0) then
      ylm=0.d00
      else
c     write (*,*) m,am
      ylm=1./dsqrt(4.*api)
c     write(*,*)  1./dsqrt(4.*api)
      endif
      endif
      if(am.ge.0) then
      ap=dsqrt(fac(l-m)/fac(l+m))
      an=((-1.)**m)*dsqrt((2.*l+1)/(4.*api))
      ylm=ap*an*cdexp(xi*m*ay)*algndr(l,m,ax)
      else
      m=-m
      ap=dsqrt(fac(l-m)/fac(l+m))
      an=dsqrt(((2.*l+1)/(4.*api)))
      ylm=ap*an*cdexp(-xi*m*ay)*algndr(l,m,ax)
      endif 
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
