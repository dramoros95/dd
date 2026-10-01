C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .    P R O G R A M   nh3_opt_claude.for                             .
C .                  MODEL BBK3CW-Z(VARIABLE)                         .
C .                  WITHOUT EXCHANGE EFFECTS                         .
C .                  NH3 (3a1 , 1e OR 2a1 orbital)                    .
C .                  one-center wave functions of Moccia (1964)       .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .
C .  WHAT THIS PROGRAM COMPUTES
C .  ---------------------------
C .  The triple differential cross section (TDCS) of NH3 (e,2e) in the
C .  BBK3CW-Z model with variable charges. It is ch4_opt_claude.for
C .  (same physics routines, same integration grids, same kinematics,
C .  same geometries igeom=1/2, same fast organisation) with the CH4
C .  target replaced by NH3:
C .
C .  TARGET (SUBROUTINE TARGET , R. Moccia, J. Chem. Phys. 40, 2176
C .  (1964), TABLE I , pyramidal NH3 , calculated equilibrium):
C .    N nucleus (Z=7) at the origin, the 3 H nuclei smeared on the
C .    sphere r = R(N-H) = 1.928 a.u. (charge 3), 10 electrons:
C .       1a1 (2e)  2a1 (2e)  3a1 (2e)  1e = 1ex+1ey (4e)
C .    every MO = sum of Slater functions r**(n-1)exp(-zeta r) S(l,m)
C .    (all coefficients of Table I, norms checked at run time).
C .  IONIZED ORBITAL (switch iorb below):
C .    iorb=1 : 3a1 (HOMO, lone pair) 96.4% p , 3.3% s , 0.04% d
C .    iorb=2 : 1e  (one of 1ex/1ey)  93.6% p , 5.8% d
C .    iorb=3 : 2a1                   95.8% s , 2.8% p , 0.5% d
C .    (iorb=0 : CH4 1t2 p part = exactly the target of
C .     ch4_opt_claude.for , kept as a regression test)
C .  The f parts (l=3) are NOT ionized (no analytic term for l=3):
C .  0.23% (3a1), 0.53% (1e), 0.89% (2a1) of the orbital. They are
C .  included in the screening Z(r).
C .
C .  TDCS (randomly oriented molecules, one-center expansion):
C .    the ionized orbital is  phi = sum_g R_g(r) S(l_g,m_g)
C .    (g = one (l,m) part of the orbital, R_g = its Slater sum).
C .    The orientation average of |T|^2 makes the parts g incoherent:
C .      <|T|^2> = sum_g  1/(2l_g+1) * sum_{mu=-l..l} |T[R_g Y(l,mu)]|^2
C .    and, as in ch4_opt_claude.for (factor 2 = 2 electrons of the
C .    orbital):
C .      TDCS = 2*aNs*anc1*anc0*anc01 * <|T|^2>
C .    For CH4 1t2 (p only) this is exactly 2*aNs*(|qs0|^2+|qs1|^2+
C .    |qs2|^2)/3*... of ch4_opt_claude.for.
C .    T[R Y(l,mu)] = analytic part (Z=1 Coulomb wave, TFGA: s, p, d,
C .    any n) + numerical part (TECWN: Z(r)-wave minus Z=1 wave, r<d).
C .    NOTE 1e : the TDCS is for ONE orbital (1ex , 2 electrons) ;
C .    the whole 1e level (1ex+1ey, 4 electrons) is 2 x TDCS.
C .
C .  WHAT CHANGED WITH RESPECT TO ch4_opt_claude.for
C .  ------------------------------------------------
C .  1. Z(r) (azs, aze : ion ; azi : neutral) computed from the Slater
C .     data of SUBROUTINE TARGET (ZSCR), any molecule XHn.
C .  2. Analytic part : TFP3 (CH4 p only) replaced by TFGB/TFGA (s, p,
C .     d , any n , any zeta), from BBK3CWZ_CH4_opt.f.
C .     TFGB (part depending only on |p| and |ke|) once per point p.
C .  3. TECWN : the 3 p columns (m=0,+1,-1) become ncm columns , one
C .     per ionized (l,mu). The mirror symmetry (igeom=1) is
C .     generalised to any l (header C of ch4_opt_claude.for):
C .       c_mu(g') = (-1)**mu c_-mu(g)  (g' = mirror y->-y of g) ,
C .       A_mu = sum_half (c_mu + (-1)**mu c_-mu) E+  (column  +mu)
C .       B_mu = sum_half (c_mu - (-1)**mu c_-mu) E-  (column  -mu)
C .       S_mu = (A_mu+B_mu)/2 , S_-mu = (-1)**mu (A_mu-B_mu)/2 , S_0=A_0
C .     (for l=1 these are exactly the formulas of ch4_opt_claude.for)
C .  4. Ion charge and Eio set by TARGET for the chosen orbital.
C .
C .  VALIDATION (done before delivery)
C .  ---------------------------------
C .   - iorb=0 (CH4 data) reproduces ch4_opt_claude.for : identical
C .     TDCS (10 digits) at 48, 54, 60, 66 deg.
C .   - TFGA (s n=1,2 ; p n=2,3 ; d n=3,4 ; all mu) vs direct
C .     numerical integration of the same Fourier-Coulomb integral
C .     (general directions of p and ke) : agreement 1e-7 (accuracy
C .     of the numerical 1F1).
C .   - Mirror-symmetry path vs full-grid path for 1e (p and d
C .     channels, mu up to 2) : identical TDCS (10 digits).
C .   - Norms of all the MOs of Table I = 1 (+-1e-5), Z(r) -> 1 (ion)
C .     and 0 (neutral) far away.
C .
C .  DIFFERENCES WITH NH3_test.for (reference program)
C .  -------------------------------------------------
C .   - NH3_test: the analytic part TFP still used the CH4 1t2 Slater
C .     data (zeta 1.373/2.95, c 1.25998...) while TECWN used NH3 ->
C .     here TFGA uses the NH3 data of the ionized MO.
C .   - NH3_test: azs (scattered) used the NEUTRAL Z(r) and azi
C .     (incident) the ION Z(r) (swapped) -> here as in Ch4.for:
C .     scattered/ejected = ion , incident = neutral.
C .   - NH3_test: only m=0 computed, then divided by 3 -> here all mu
C .     (orientation average , header above).
C .   - NH3_test: coefficients partly different from Table I (2a1,
C .     1a1, radial function of 3a1), the 2s zeta=1.75 function taken
C .     as 1s, p part of 2a1 and f(-3) parts omitted, R(N-H)=1.9124
C .     -> here every coefficient of Table I, R(N-H)=1.928 (Moccia).
C .   - NH3_test: s and d parts of 3a1 not ionized -> included here.
C .
C .  COMPILATION:
C .   Intel oneAPI (as ch4_opt_claude.for):
C .   ifx -O3 -xHOST -qopenmp -extend-source 132 -no-prec-div
C .       -fp-model fast=2 nh3_opt_claude.for -o nh3_opt_claude.exe
C .   gfortran:
C .   gfortran -O3 -march=native -fopenmp -ffixed-line-length-132
C .       nh3_opt_claude.for -o nh3_opt_claude
C .   export OMP_NUM_THREADS=16
C .   ./nh3_opt_claude.exe <<< "1 61"
C .  Input "ni n" = first and last angle index (angle = 6 deg *
C .  (index-1)); without input: ni=1, n=61. With ni > 1 the results
C .  are APPENDED to the existing output files.
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      IMPLICIT DOUBLEPRECISION (a-h)
      IMPLICIT COMPLEX*16 (o-z)
      character*100 fname,rname
      PARAMETER (nleg=10,nd=4,md=1,np=5,nt=3,nf=3)
      PARAMETER (nlp=60,nlt=40,nlf=40)
C . ntec = points of the TECWN grid (ejected electron, md interval)
C . nts  = points of the TSCWN grid (incident/scattered, nd intervals)
      PARAMETER (ntec=md*nleg**3,nts=nd*nleg**3)
C . max number of (l,m) parts of the target, of Slater functions,
C . of charged shells, of ionized (l,mu)
      PARAMETER (norbx=40,nstox=400,nshx=4,ncmx=16)
      PARAMETER (api=3.141592654d00,acon=0.07349793195d00
     $,arad=0.0174532925d00)
      PARAMETER (xi=(0.0d00,1.0d00))
      DIMENSION azs(nd,nleg),aze(md,nleg),azi(nd,nleg),dx(nleg),dw(nleg)
     $  ,cp(10),ct(10),cf(4),qffs(nd,nleg,nleg,nleg)
     $  ,qffi(nd,nleg,nleg,nleg),q0ffi(nd,nleg,nleg,nleg)
     $  ,q1ffs(nd,nleg,nleg,nleg),qffe(md,nleg,nleg,nleg)
     $  ,dxp(500),dxt(500),dxf(500),dwp(500),dwt(500),dwf(500)
      DIMENSION cd(nd+1)

C . Geometry / orbital switches and symmetry flag
      INTEGER igeom,iorb,nteff
      LOGICAL lsym

C . THE TARGET (SUBROUTINE TARGET): nuclear charge, charged shells,
C . its (l,m) parts (name, l, m, occupations neutral/ion, ionized?,
C . Slater functions n , zeta , c)
      CHARACTER*8 alab(norbx)
      CHARACTER*8 aname
      DOUBLE PRECISION aqsh(nshx),arsh(nshx),aocn(norbx),aoci(norbx)
     $  ,esz(nstox),csz(nstox)
      INTEGER lorb(norbx),morb(norbx),iono(norbx),ist(norbx),nst(norbx)
     $  ,nsn(nstox),nsh,norb
      DOUBLE PRECISION RADF,ORBNRM,ZSCR

C . THE IONIZED PARTS: nio parts (jorb), their ncm (l,mu) channels
C . (iocm = part, lcm = l, mcm = mu, wcm = 1/(2l+1), ipcm = channel
C . -mu), radial function at the ejected nodes (fre), their nsl Slater
C . terms for the analytic part (l, n, zeta, cfv, part)
      INTEGER jorb(norbx),iocm(ncmx),lcm(ncmx),mcm(ncmx),ipcm(ncmx)
     $  ,icpm(ncmx),ltv(nstox),nsv(nstox),iov(nstox),nio,ncm,nsl,kmx
      DOUBLE PRECISION wcm(ncmx),fre(md,nleg,norbx),esv(nstox)
     $  ,cfv(nstox)
      INTEGER icm,je,mm,lq,mq,io,jo,ln,icol,ip,im

C . Precomputed grids (REAL, explicitly declared) and TSCWN coefficients
C . split in real/imaginary parts (for vectorisation)
      DOUBLE PRECISION r_tec_x(ntec),r_tec_y(ntec),r_tec_z(ntec)
      DOUBLE PRECISION r_ts_x(nts),r_ts_y(nts),r_ts_z(nts)
      DOUBLE PRECISION cts1r(nts),cts1i(nts),cts2r(nts),cts2i(nts)
      COMPLEX*16 qf,wbase,base_cf,qbr,qy,qym,qcc
      DOUBLE PRECISION wgt,sqrt2pi3,sgm,ake2
      INTEGER g,nang,ia,ntask,ndone,nprint,nmine,nkh

C . Per-angle data (index ia = 1..nang)
C .   ctr/cti(g,icm,ia) : TECWN coefficients (real/imag), ncm columns
C .   akexa..akeza      : ejected electron momentum ke
C .   ak01xa..ak01za    : K01=(Ks-Ke)/2 ,  anc01a : its Gamow factor
C .   qsacc(icm,ia)     : the amplitude T of the (l,mu) icm
      DOUBLE PRECISION, ALLOCATABLE :: ctr(:,:,:),cti(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: akexa(:),akeya(:),akeza(:)
      DOUBLE PRECISION, ALLOCATABLE :: at11a(:)
      DOUBLE PRECISION, ALLOCATABLE :: ak01xa(:),ak01ya(:),ak01za(:)
      DOUBLE PRECISION, ALLOCATABLE :: anc01a(:)
      COMPLEX*16, ALLOCATABLE :: qsacc(:,:),qsl(:,:)

C . Thread-private work arrays for one block of nlf points (all phi_p
C . values of one (p,theta_p,kf) line): TECWN phases E+ (epp) and E-
C . (emp) ; TFGB data of the nlf points (qbv, pbv, qbk) ; the sums of
C . the columns (qcol) and the TECWN amplitudes (qsm) ; TFGA (tfv)
      DOUBLE PRECISION, ALLOCATABLE :: eppr(:,:),eppi(:,:),empr(:,:)
     $  ,empi(:,:)
      COMPLEX*16, ALLOCATABLE :: qbv(:,:),pbv(:,:),qbk(:,:,:)
     $  ,qcol(:,:),qsm(:,:),tfv(:)
      DOUBLE PRECISION apxb(nlf),apyb(nlf),apzb(nlf),chb(nlf)
      COMPLEX*16 qa1(nlf),qa2(nlf)

C . Scalars used inside the parallel region (all declared explicitly)
      DOUBLE PRECISION ap,atp,afp,apx,apy,apz,apex,apey,apez
      DOUBLE PRECISION akpx,akpy,akpz,akp01x,akp01y,akp01z,aqp
      DOUBLE PRECISION cp_del,ct_del,cf_del,adr,dcs,dsn
      DOUBLE PRECISION ew1r,ew1i,ew2r,ew2i
      DOUBLE PRECISION e1r,e1i,e2r,e2i,bwr,bwi,bcs,bsn
      DOUBLE PRECISION acc1r,acc1i,acc2r,acc2i
      COMPLEX*16 QT1,QC2,QC1,QT2,qdiff
      COMPLEX*16 q0s,tcs1
      INTEGER kp,i1,kt,j,kf,k,ig,kk

      if (mod(nlf,2).ne.0) stop 'nlf must be even'
      if (mod(nleg,2).ne.0) stop 'nleg must be even (mirror symmetry)'

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  GEOMETRY:  1 = coplanar (theta_e scanned in the x-z plane)
C .             2 = Ch4.for  (Ete fixed, azimuth phi_e scanned)
C .  ORBITAL :  1 = 3a1 , 2 = 1e , 3 = 2a1  (0 = CH4 1t2 test)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      igeom=1
      iorb=1
      Ete=api/2.d0
C . The mirror-symmetry shortcut is exact only when ke_y=0 for
C . every angle, i.e. in the coplanar geometry.
      lsym=(igeom.eq.1)

C .  THE TARGET (Eio = ionization energy of the chosen orbital)
      CALL TARGET(iorb,aname,acen,Eio,nsh,aqsh,arsh,norb,alab,lorb
     $,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox,nshx)

C .  FOR INTEGRAL OVER R INSIDE SUBROUTINES/FUNCTIONS
      CALL GAULEG(-1.d0,1.d0,dx,dw,nleg,2d-16)

C .  FOR INTEGRAL OVER P THETA PHI
      CALL GAULEG(-1.d0,1.d0,dxp,dwp,nlp,2d-16)
      CALL GAULEG(-1.d0,1.d0,dxt,dwt,nlt,2d-16)
      CALL GAULEG(-1.d0,1.d0,dxf,dwf,nlf,2d-16)

C .  CONSTANTS / ENERGIES / ANGLES (same meaning as in Ch4.for)
      d=4.d00 ; az1=1.d00 ; az= 1.d00 ; alpha=0.005d00 ; almda=0.005d00
      az01=-0.5d00 ; eps=2.d00*arad

C . ANGLES & ENERGIES (same kinematics as ch4_opt_claude.for :
C . Ei=500 eV , ejected Ee=74 eV , scattering angle 6 deg).
C . Eio comes from TARGET (vertical IP of the chosen NH3 orbital).
C . CHANGE THEM HERE FOR YOUR EXPERIMENT.
      ets=6.d00*arad ; afs=180.d00*arad
      Ee=74.d00 ; Ei=500.d00
      Es=Ei-Ee-Eio
      if (Es.le.0.d0) stop 'Es = Ei-Ee-Eio <= 0 : check the energies'
      aki=dsqrt(Ei*acon)
      aks=dsqrt(Es*acon)
      ake=dsqrt(Ee*acon)

C .  Ki along z ; Ks in the (x,z) plane (afs=180 deg)
      akix=0.d0
      akiy=0.d0
      akiz=aki

      aksx=aks*dsin(ets)*dcos(afs)
      aksy=aks*dsin(ets)*dsin(afs)
      aksz=aks*dcos(ets)

C .  Transfer K = Ki-Ks, its magnitude and direction (for the grid)
      adx=akix-aksx
      ady=akiy-aksy
      adz=akiz-aksz
      ad=dsqrt(adx**2.+ady**2.+adz**2.)
      atd=dacos(adz/ad)
      atk=api-atd

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  CHECKS OF THE TARGET : norm of every MO (sum of its (l,m) parts
C .  with the same name), charge far away of the ion (1) and of the
C .  neutral molecule (0)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      write(*,'(A,A,A,F8.3,A)') ' Target ',trim(aname),
     $' , Eio = ',Eio,' eV'
      do io=1,norb
         ln=1
         do jo=1,io-1
            if (alab(jo).eq.alab(io)) ln=0
         enddo
         if (ln.eq.1) then
            bn=0.d0
            do jo=io,norb
               if (alab(jo).eq.alab(io))
     $            bn=bn+ORBNRM(jo,ist,nst,nsn,esz,csz)
            enddo
            write(*,'(A,A8,A,F9.5)') '   norm of ',alab(io),' = ',bn
            if (dabs(bn-1.d0).gt.1.d-3) write(*,*) ' WARNING : NORM'
         endif
      enddo
      bzi=ZSCR(1.d3,2,acen,nsh,aqsh,arsh,norb,aocn,aoci,ist,nst
     $,nsn,esz,csz)
      bzn=ZSCR(1.d3,1,acen,nsh,aqsh,arsh,norb,aocn,aoci,ist,nst
     $,nsn,esz,csz)
      write(*,'(A,F9.5,A,F9.5)') '   Z(r=1000) ion =',bzi,
     $' , neutral =',bzn
      if (dabs(bzi-az).gt.1.d-2.or.dabs(bzn).gt.1.d-2)
     $  write(*,*) ' WARNING : CHARGES OF THE TARGET'

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  THE IONIZED PARTS, THEIR (l,mu) CHANNELS (mu = 0,1,-1,2,-2) AND
C .  THEIR SLATER TERMS FOR TFGA
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      nio=0
      ncm=0
      nsl=0
      kmx=0
      do io=1,norb
      if (iono(io).eq.1) then
         lq=lorb(io)
         if (lq.gt.2) stop 'ionized part with l > 2 : not possible'
         nio=nio+1
         jorb(nio)=io
         do mm=0,2*lq
            ncm=ncm+1
            if (ncm.gt.ncmx) stop 'too many ionized (l,mu): ncmx'
            iocm(ncm)=nio
            lcm(ncm)=lq
            mcm(ncm)=((mm+1)/2)*(1-2*mod(mm+1,2))
            wcm(ncm)=1.d0/dble(2*lq+1)
         enddo
         do i=ist(io),ist(io)+nst(io)-1
            nsl=nsl+1
            ltv(nsl)=lq
            nsv(nsl)=nsn(i)
            esv(nsl)=esz(i)
C .  cfv = (-1)**(n-l) * (n-l)! * c * N (see TFGA)
            cfv(nsl)=((-1)**(nsn(i)-lq))*fac(nsn(i)-lq)*csz(i)
     $           *fa(esz(i),nsn(i)-1)
            iov(nsl)=nio
            kmx=max(kmx,nsn(i)-lq)
         enddo
      endif
      enddo
      if (kmx.gt.12) stop 'ionized Slater function with n-l > 12'
C .  channel -mu of each channel (ipcm) ; column type icpm : 1 = E+ ,
C .  2 = E- (mirror symmetry: column +mu -> A_mu on E+ , column -mu ->
C .  B_mu on E- ; full grid: every column on E+ = exp(-i p.r))
      do icm=1,ncm
         ipcm(icm)=icm
         do jo=1,ncm
            if (iocm(jo).eq.iocm(icm).and.mcm(jo).eq.-mcm(icm))
     $         ipcm(icm)=jo
         enddo
         icpm(icm)=1
         if (lsym.and.mcm(icm).lt.0) icpm(icm)=2
      enddo
      write(*,'(A,I0,A,I0,A)') '   ionized : ',nio,' (l,m) parts , ',
     $ncm,' (l,mu) channels'
      do je=1,nio
         io=jorb(je)
         write(*,'(A,A8,A,I0,A,I0,A,F9.5)') '     ',alab(io),' l=',
     $   lorb(io),' m=',morb(io),' , weight in the MO =',
     $   ORBNRM(io,ist,nst,nsn,esz,csz)
      enddo
      call flush(6)

C .   VARIABLE CHARGE Z (scattered, ejected : ion ; incident : neutral)
      CALL ZTAB(nd,d,dx,nleg,azs,2,acen,nsh,aqsh,arsh,norb,aocn,aoci
     $,ist,nst,nsn,esz,csz)
      CALL ZTAB(md,d,dx,nleg,aze,2,acen,nsh,aqsh,arsh,norb,aocn,aoci
     $,ist,nst,nsn,esz,csz)
      CALL ZTAB(nd,d,dx,nleg,azi,1,acen,nsh,aqsh,arsh,norb,aocn,aoci
     $,ist,nst,nsn,esz,csz)

C .   Radial functions of the ionized parts at the ejected nodes
      do kd=1,md+1
         cd(kd)=(kd-1.d0)*d/md
      enddo
      do je=1,nio
      do id=1,md
      do i=1,nleg
         ar=(dx(i)*(cd(id+1)-cd(id))+cd(id+1)+cd(id))*0.5d00
         fre(id,i,je)=RADF(ar,jorb(je),ist,nst,nsn,esz,csz)
      enddo
      enddo
      enddo

C .   Gamow factors anc1 (ejected) and anc0 (scattered)
      azke=az/ake
      if (azke.ge.-1d-7.and.azke.le.1d-7) then
         anc1=1.d0
      else
         anc1=2.d0*api*azke/(1.d0-dexp(-2.d0*api*azke))
      endif

      azks=az1/aks
      if (azks.ge.-1d-7.and.azks.le.1d-7) then
         anc0=1.d0
      else
         anc0=2.d0*api*azks/(1.d0-dexp(-2.d0*api*azks))
      endif

C .   Coulomb normalisations
      ucds=uc(1.d00+xi*azks)*dexp(api*azks/2.d00)
      ucde=uc(1.d00+xi*azke)*dexp(api*azke/2.d00)
      ucdi=1.d0+0.d0*xi

C .   Division of the integration intervals over p / theta_p / phi_p
C .   (finer around |p|=K and around the direction of -K)
      cp(1)=0.d0 ; cp(2)=0.05d0 ; cp(3)=0.2d0 ; cp(4)=0.98*ad
      cp(5)=1.02d0*ad ; cp(6)=4.d0
      ct(1)=0.d0 ; ct(2)=atk-eps ; ct(3)=atk+eps ; ct(4)=api
      cf(1)=0.d0 ; cf(2)=api-eps ; cf(3)=api+eps ; cf(4)=2.d0*api

C .   COULOMB WAVE FUNCTIONS ON THE r GRID (angle independent)
C .      qffs  : scattered electron, variable charge (ion)
C .      qffi  : incident  electron, variable charge (neutral)
C .      q1ffs : scattered electron, Z=1
C .      q0ffi : incident  electron, Z=0 (plane wave)
      CALL QSCWF(qffs,azs,aksx,aksy,aksz,dx,nleg,nleg,nleg,d,nd,ucds)
      CALL QICWF(qffi,azi,akix,akiy,akiz,dx,nleg,nleg,nleg,d,nd,ucdi)
      CALL Q1SCWF(q1ffs,aksx,aksy,aksz,dx,nleg,nleg,nleg,d,nd)
      CALL Q0SCWF(q0ffi,akix,akiy,akiz,dx,nleg,nleg,nleg,d,nd)

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  PRECOMPUTE THE TSCWN GRID AND COEFFICIENTS (as ch4_opt_claude)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      do kd=1,nd+1
         cd(kd)=(kd-1.)*d/nd
      enddo
      ay00=1.d0/dsqrt(4.d0*api)
      sqrt2pi3=(2.d0*api)*dsqrt(2.d0*api)

      g=0
      do id=1,nd
      do i=1,nleg
      do j=1,nleg
      do k=1,nleg
         g=g+1
         ar=(dx(i)*(cd(id+1)-cd(id))+cd(id+1)+cd(id))/2.d00
         et=api*(dx(j)+1.d00)/2.d00
         ef=api*(dx(k)+1.d00)
         r_ts_x(g)=ar*dsin(et)*dcos(ef)
         r_ts_y(g)=ar*dsin(et)*dsin(ef)
         r_ts_z(g)=ar*dcos(et)
         qf=dconjg(qffs(id,i,j,k))*qffi(id,i,j,k)
     $     -dconjg(q1ffs(id,i,j,k))*q0ffi(id,i,j,k)
         wbase=qf*dexp(-alpha*ar)*dw(i)*dw(j)*dw(k)*dsin(et)
     $        *(api**2)*(cd(id+1)-cd(id))/4.d0/sqrt2pi3
         cts1r(g)=dreal(wbase*ar)
         cts1i(g)=dimag(wbase*ar)
         cts2r(g)=dreal(wbase*ar*ar*ay00)
         cts2i(g)=dimag(wbase*ar*ar*ay00)
      enddo
      enddo
      enddo
      enddo

C .  Angle indices to compute (angle = 6 deg * (index-1))
      Read(*,*,end=999) ni,n
      goto 998
999   ni=1
      n=61
998   continue
      nang=n-ni+1
      if (nang.lt.1) stop 'nothing to do (n < ni)'

C .  Output files: nh3_<orbital>_<geometry>_<theta_s>_<Ee>ev.dat and a
C .  results file
      if (igeom.eq.1) then
         write(fname,'(A,A,A,I0,A,I0,A)') 'nh3_',trim(aname),
     $     '_coplanar_',nint(ets/arad),'_',nint(Ee),'ev.dat'
         rname='results_nh3_coplanar.dat'
      else
         write(fname,'(A,A,A,I0,A,I0,A)') 'nh3_',trim(aname),
     $     '_phie_',nint(ets/arad),'_',nint(Ee),'ev.dat'
         rname='results_nh3_phie.dat'
      endif
      if (ni.eq.1) then
         open(unit=2,file=trim(fname),status='replace')
         open(unit=3,file=trim(rname),status='replace')
      else
        open(unit=2,file=trim(fname),status='unknown',position='append')
        open(unit=3,file=trim(rname),status='unknown',position='append')
      endif
      nh=30
      h=api/nh

C .  Size of the TECWN grid actually used: half grid with the mirror
C .  symmetry, full grid otherwise
      if (lsym) then
         nkh=nleg/2
      else
         nkh=nleg
      endif
      nteff=md*nleg*nleg*nkh

      allocate(ctr(nteff,ncm,nang),cti(nteff,ncm,nang))
      allocate(akexa(nang),akeya(nang),akeza(nang),at11a(nang))
      allocate(ak01xa(nang),ak01ya(nang),ak01za(nang),anc01a(nang))
      allocate(qsacc(ncm,nang))

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  PRECOMPUTE ALL ANGLE-DEPENDENT DATA (once per angle, cheap)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      do ia=1,nang
      i_ang=ni+ia-1
      ate=h*(i_ang-1)
      at11a(ia)=ate/arad

C .   Ke: igeom=1 theta_e=ate in the x-z plane ; igeom=2 phi_e=ate
      if (igeom.eq.1) then
         akex=ake*dsin(ate)
         akey=0.d00
         akez=ake*dcos(ate)
      else
         akex=ake*dsin(Ete)*dcos(ate)
         akey=ake*dsin(Ete)*dsin(ate)
         akez=ake*dcos(Ete)
      endif
      akexa(ia)=akex
      akeya(ia)=akey
      akeza(ia)=akez

C .   Ejected electron Coulomb function (minus its Z=1 part) on r grid
      CALL QECWF(qffe,aze,akex,akey,akez,dx,nleg,nleg,nleg,d,md,ucde)

C .   TECWN coefficients. Original function (for l, mu):
C .     TECWN = 1/(2pi)^1.5 * sum_g conj(qffe(g))*R(r)*Y(l,mu)(g)
C .             * exp(i pe.r_g) * r^2 * weights ,   pe = -p
C .   c_mu(g) = everything except the phase. Columns stored:
C .     lsym  : column +mu : c_mu + (-1)**mu c_-mu  (A_mu , on E+)
C .             column -mu : c_mu - (-1)**mu c_-mu  (B_mu , on E-)
C .             on the half grid (k=1..nleg/2 ; k'=nleg+1-k is the
C .             mirror y -> -y , qffe(k')=qffe(k) because ke_y=0)
C .     else  : column mu : c_mu on the full grid
      do kd=1,md+1
         cd(kd)=(kd-1.)*d/md
      enddo
      if (lsym.and.akey.ne.0.d0) stop 'mirror symmetry needs akey=0'
      g=0
      do id=1,md
      do i=1,nleg
      do j=1,nleg
      do k=1,nkh
         g=g+1
         ar=(dx(i)*(cd(id+1)-cd(id))+cd(id+1)+cd(id))*0.5d0
         et=api*(dx(j)+1.d00)*0.5d0
         ef=api*(dx(k)+1.d00)
         if (ia.eq.1) then
            r_tec_x(g)=ar*dsin(et)*dcos(ef)
            r_tec_y(g)=ar*dsin(et)*dsin(ef)
            r_tec_z(g)=ar*dcos(et)
         endif
         cth=dcos(et)
         wgt=dw(i)*dw(j)*dw(k)*dsin(et)*ar*ar*(api**2)*(cd(id+1)-cd(id))
     $     /4.d0
         base_cf=dconjg(qffe(id,i,j,k))*wgt/sqrt2pi3
         do icm=1,ncm
            lq=lcm(icm)
            mq=mcm(icm)
            qbr=base_cf*fre(id,i,iocm(icm))
            qy=ylm(lq,dble(mq),cth,ef)
            if (lsym.and.mq.ne.0) then
               qym=ylm(lq,dble(-mq),cth,ef)
               sgm=dble(1-2*mod(abs(mq),2))
               if (mq.gt.0) then
                  qcc=qbr*(qy+sgm*qym)
               else
                  qcc=qbr*(qym-sgm*qy)
               endif
            else
               qcc=qbr*qy
            endif
            ctr(g,icm,ia)=dreal(qcc)
            cti(g,icm,ia)=dimag(qcc)
         enddo
      enddo
      enddo
      enddo
      enddo

C .   K01 = (Ks-Ke)/2 and its Gamow factor anc01
      ak01xa(ia)=(aksx-akex)/2.d00
      ak01ya(ia)=(aksy-akey)/2.d00
      ak01za(ia)=(aksz-akez)/2.d00
      ak01=dsqrt(ak01xa(ia)**2.+ak01ya(ia)**2.+ak01za(ia)**2.)

      azk01=az01/ak01
      if (azk01.ge.-1d-7.and.azk01.le.1d-7) then
         anc01a(ia)=1.d0
      else
         anc01a(ia)=2.d0*api*azk01/(1.d0-dexp(-2.d0*api*azk01))
      endif
      enddo

C .  |ke|**2 for TFGB (the same for all the angles)
      ake2=ake*ake

      qsacc=(0.d0,0.d0)
      ntask=np*nlp*nt
      ndone=0
      nprint=max(1,ntask/20)
      write(*,'(A,I0,A,I0,A,I0,A)') ' Computing angles ',ni,' to ',n,
     $' (',nang,' angles) in a single pass over the (p,theta,phi) grid'
      if (lsym) then
         write(*,'(A)') ' coplanar geometry: mirror-symmetric half grid'
      else
         write(*,'(A)') ' general geometry: full TECWN grid'
      endif
      call flush(6)

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  MAIN INTEGRATION over p=(ap,atp,afp): every point contributes to
C .  every ejection angle. Each thread accumulates its own partial sums
C .  qsl(ncm,nang); they are added together at the end (CRITICAL).
C .  Parallel tasks = (kp,i1,kt) = 5*60*3 = 900, dynamic scheduling.
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!$OMP PARALLEL DEFAULT(NONE)
!$OMP& SHARED(cp,ct,cf,dxp,dxt,dxf,dwp,dwt,dwf,
!$OMP& r_tec_x,r_tec_y,r_tec_z,ctr,cti,lsym,nteff,
!$OMP& r_ts_x,r_ts_y,r_ts_z,cts1r,cts1i,cts2r,cts2i,
!$OMP& akexa,akeya,akeza,ak01xa,ak01ya,ak01za,ake2,
!$OMP& nang,aki,az,az01,alpha,almda,aksx,aksy,aksz,
!$OMP& ncm,nio,nsl,kmx,iocm,lcm,mcm,ipcm,icpm,ltv,nsv,esv,cfv,iov,
!$OMP& qsacc,ntask,ndone,nprint)
!$OMP& PRIVATE(kp,i1,kt,j,kf,k,ig,ia,ap,atp,afp,apx,apy,apz,
!$OMP& apex,apey,apez,akpx,akpy,akpz,akp01x,akp01y,akp01z,aqp,
!$OMP& cp_del,ct_del,cf_del,adr,dcs,dsn,ew1r,ew1i,ew2r,ew2i,
!$OMP& e1r,e1i,e2r,e2i,bwr,bwi,bcs,bsn,acc1r,acc1i,acc2r,acc2i,
!$OMP& QT1,QC2,QC1,QT2,qdiff,icm,icol,ip,im,kk,sgm,
!$OMP& eppr,eppi,empr,empi,qbv,pbv,qbk,qcol,qsm,tfv,
!$OMP& apxb,apyb,apzb,chb,qa1,qa2,qsl,nmine)
      allocate(eppr(nteff,nlf),eppi(nteff,nlf))
      allocate(empr(nteff,nlf),empi(nteff,nlf))
      allocate(qbv(nsl,nlf),pbv(nsl,nlf),qbk(0:kmx,nsl,nlf))
      allocate(qcol(ncm,nlf),qsm(ncm,nlf),tfv(ncm))
      allocate(qsl(ncm,nang))
      qsl=(0.d0,0.d0)

!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC,1)
      do kp=1,np
      do i1=1,nlp
      do kt=1,nt
         cp_del=cp(kp+1)-cp(kp)
         ap=(dxp(i1)*cp_del+cp(kp+1)+cp(kp))/2.d0
         ct_del=ct(kt+1)-ct(kt)

      do j=1,nlt
         atp=(dxt(j)*ct_del+ct(kt+1)+ct(kt))/2.d0

      do kf=1,nf
         cf_del=cf(kf+1)-cf(kf)

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C . (a) ANGLE-INDEPENDENT PART, once per point, for a block of nlf
C .     points (all phi_p of this line)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      do k=1,nlf
         afp=(dxf(k)*cf_del+cf(kf+1)+cf(kf))/2.d0
         apx=ap*dsin(atp)*dcos(afp)
         apy=ap*dsin(atp)*dsin(afp)
         apz=ap*dcos(atp)
         apxb(k)=apx
         apyb(k)=apy
         apzb(k)=apz
C .    integration weight ("chandler" in Ch4.for)
         chb(k)=dwp(i1)*dwt(j)*dwf(k)*dsin(atp)*ap*ap
     $    *cp_del*ct_del*cf_del/8.d00

C .    TSCWN sums: sw1 = TSCWN(nn=2,mn=1), sw2 = TSCWN(nn=1,mn=0)
         ew1r=0.d0 ; ew1i=0.d0 ; ew2r=0.d0 ; ew2i=0.d0
!$OMP SIMD REDUCTION(+:ew1r,ew1i,ew2r,ew2i) PRIVATE(adr,dcs,dsn)
         do ig=1,nts
            adr=apx*r_ts_x(ig)+apy*r_ts_y(ig)+apz*r_ts_z(ig)
            dcs=dcos(adr)
            dsn=dsin(adr)
            ew1r=ew1r+cts2r(ig)*dcs-cts2i(ig)*dsn
            ew1i=ew1i+cts2r(ig)*dsn+cts2i(ig)*dcs
            ew2r=ew2r+cts1r(ig)*dcs-cts1i(ig)*dsn
            ew2i=ew2i+cts1r(ig)*dsn+cts1i(ig)*dcs
         enddo

C .    Ki+p ; QT1 and QC2 do not depend on the ejection angle
         akpx=apx
         akpy=apy
         akpz=aki+apz
         QT1=tcs1(alpha,az,akpx,akpy,akpz,aksx,aksy,aksz)
         QC2=q0s(alpha,az,akpx,akpy,akpz,aksx,aksy,aksz)
C .    qa1 = QT1+QW1 (2nd term), qa2 = QC2+QW2 (3rd term) of Ch4.for
         qa1(k)=QT1+dcmplx(ew1r,ew1i)
         qa2(k)=QC2+dcmplx(ew2r,ew2i)

C .    Analytic part of the direct term: the part which depends only
C .    on |p| and |ke| (TFGB), once per point
         aqp=apx*apx+apy*apy+apz*apz
         CALL TFGB(az,aqp,ake2,nsl,ltv,nsv,esv,kmx,qbv(1,k),pbv(1,k)
     $   ,qbk(0,1,k))

C .    TECWN phases, e = exp(i pe.r) with pe=-p
         if (lsym) then
C .    half grid point g and its mirror g' (y -> -y):
C .       E+ = e+e' = 2cos(py y) exp(-i(px x+pz z))
C .       E- = e-e' = -2i sin(py y) exp(-i(px x+pz z))
!$OMP SIMD PRIVATE(adr,dcs,dsn,bcs,bsn)
         do ig=1,nteff
            adr=-(apx*r_tec_x(ig)+apz*r_tec_z(ig))
            dcs=dcos(adr)
            dsn=dsin(adr)
            bcs=2.d0*dcos(apy*r_tec_y(ig))
            bsn=2.d0*dsin(apy*r_tec_y(ig))
            eppr(ig,k)=bcs*dcs
            eppi(ig,k)=bcs*dsn
            empr(ig,k)=bsn*dsn
            empi(ig,k)=-bsn*dcs
         enddo
         else
C .    full grid: the same phase e for all the columns (E+ only)
!$OMP SIMD PRIVATE(adr,dcs,dsn)
         do ig=1,nteff
            adr=-(apx*r_tec_x(ig)+apy*r_tec_y(ig)+apz*r_tec_z(ig))
            eppr(ig,k)=dcos(adr)
            eppi(ig,k)=dsin(adr)
         enddo
         endif
      enddo

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C . (b) ANGLE-DEPENDENT PART: loop over all ejection angles
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      do ia=1,nang

C .  TECWN : the ncm column sums, two points (k,k+1) at a time (the
C .  coefficients of a column are loaded once for both points)
      do k=1,nlf,2
         do icol=1,ncm
            acc1r=0.d0 ; acc1i=0.d0 ; acc2r=0.d0 ; acc2i=0.d0
            if (icpm(icol).eq.1) then
!$OMP SIMD REDUCTION(+:acc1r,acc1i,acc2r,acc2i)
!$OMP& PRIVATE(e1r,e1i,e2r,e2i,bwr,bwi)
            do ig=1,nteff
               bwr=ctr(ig,icol,ia)
               bwi=cti(ig,icol,ia)
               e1r=eppr(ig,k)
               e1i=eppi(ig,k)
               e2r=eppr(ig,k+1)
               e2i=eppi(ig,k+1)
               acc1r=acc1r+bwr*e1r-bwi*e1i
               acc1i=acc1i+bwr*e1i+bwi*e1r
               acc2r=acc2r+bwr*e2r-bwi*e2i
               acc2i=acc2i+bwr*e2i+bwi*e2r
            enddo
            else
!$OMP SIMD REDUCTION(+:acc1r,acc1i,acc2r,acc2i)
!$OMP& PRIVATE(e1r,e1i,e2r,e2i,bwr,bwi)
            do ig=1,nteff
               bwr=ctr(ig,icol,ia)
               bwi=cti(ig,icol,ia)
               e1r=empr(ig,k)
               e1i=empi(ig,k)
               e2r=empr(ig,k+1)
               e2i=empi(ig,k+1)
               acc1r=acc1r+bwr*e1r-bwi*e1i
               acc1i=acc1i+bwr*e1i+bwi*e1r
               acc2r=acc2r+bwr*e2r-bwi*e2i
               acc2i=acc2i+bwr*e2i+bwi*e2r
            enddo
            endif
            qcol(icol,k)=dcmplx(acc1r,acc1i)
            qcol(icol,k+1)=dcmplx(acc2r,acc2i)
         enddo
C .  TECWN amplitudes S_mu of the channels
         do kk=k,k+1
         do icm=1,ncm
            if (.not.lsym.or.mcm(icm).eq.0) then
               qsm(icm,kk)=qcol(icm,kk)
            else if (mcm(icm).gt.0) then
               ip=icm
               im=ipcm(icm)
               qsm(icm,kk)=0.5d0*(qcol(ip,kk)+qcol(im,kk))
            else
               ip=ipcm(icm)
               im=icm
               sgm=dble(1-2*mod(abs(mcm(icm)),2))
               qsm(icm,kk)=(0.5d0*sgm)*(qcol(ip,kk)-qcol(im,kk))
            endif
         enddo
         enddo
      enddo

      do k=1,nlf
         apex=-apxb(k)
         apey=-apyb(k)
         apez=-apzb(k)

C .  DIRECT TERM: analytic part of all the channels (TFGA) + TECWN
         CALL TFGA(az,apex,apey,apez,akexa(ia),akeya(ia),akeza(ia)
     $   ,nsl,ltv,nsv,esv,cfv,iov,kmx,qbv(1,k),pbv(1,k),qbk(0,1,k)
     $   ,nio,ncm,iocm,lcm,mcm,tfv)

C .  2ND AND 3RD TERMS: TQ1-TQ2 = QC1*(QT1+QW1) - QT2*(QC2+QW2)
C .  (QC1 and QT2 depend on the angle through K01)
         akp01x=ak01xa(ia)-apxb(k)
         akp01y=ak01ya(ia)-apyb(k)
         akp01z=ak01za(ia)-apzb(k)
         QC1=q0s(almda,az01,akp01x,akp01y,akp01z,
     $           ak01xa(ia),ak01ya(ia),ak01za(ia))
         QT2=tcs1(almda,az01,akp01x,akp01y,akp01z,
     $           ak01xa(ia),ak01ya(ia),ak01za(ia))

C .  INTEGRATION SUMMATION (one amplitude per (l,mu) channel)
         qdiff=(QC1*qa1(k)-QT2*qa2(k))*chb(k)
         do icm=1,ncm
            qsl(icm,ia)=qsl(icm,ia)+qdiff*(tfv(icm)+qsm(icm,k))
         enddo
      enddo
      enddo

      enddo
      enddo

C .  progress counter (ATOMIC CAPTURE: no race on ndone)
!$OMP ATOMIC CAPTURE
      ndone=ndone+1
      nmine=ndone
!$OMP END ATOMIC
      if (mod(nmine,nprint).eq.0) then
         write(*,'(A,I3,A)') ' progress: ',nint(100.*nmine/ntask),' %'
         call flush(6)
      endif
      enddo
      enddo
      enddo
!$OMP END DO

C .  add the partial sums of this thread to the total
!$OMP CRITICAL
      qsacc=qsacc+qsl
!$OMP END CRITICAL
      deallocate(eppr,eppi,empr,empi,qbv,pbv,qbk,qcol,qsm,tfv,qsl)
!$OMP END PARALLEL

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  CROSS SECTION : 2*aNs*Gamow * sum_channels 1/(2l+1)*|T|^2
C .  (= 2*aNs*(|qs0|^2+|qs1|^2+|qs2|^2)/3*... of ch4_opt_claude.for
C .  for a p orbital)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      aNs=ake*aks/(aki*api)
      do ia=1,nang
         i_ang=ni+ia-1
         asigp=0.d0
         do icm=1,ncm
            asigp=asigp+wcm(icm)*dreal(qsacc(icm,ia)
     $           *dconjg(qsacc(icm,ia)))
         enddo
         asigmaD=asigp*anc1*anc0*anc01a(ia)
         asigma=2.0*aNs*(asigmaD)

         write(*,'(A,I3,A,I3,A,F7.1,A,E14.6)') ' Angle ', i_ang, ' / ',
     $    n, ' : angle = ', at11a(ia), ' deg -> TDCS = ', asigma
         write(2,'(F12.4,4X,E18.10)') at11a(ia), asigma
         write(3,'(F12.4,4X,E18.10)') at11a(ia), asigma
      enddo

      close(unit=2)
      close(unit=3)
      end


C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  TARGET : NH3 , ONE-CENTER SCF MO'S OF R. MOCCIA , J. CHEM. PHYS.  .
C .  40 , 2176 (1964) , TABLE I (pyramidal NH3 , calculated           .
C .  equilibrium : R(N-H) = 1.9280 a.u. , theta = 108.90 deg).        .
C .  Every MO is given as a list of (l,m) parts ; each part is a sum  .
C .  of normalised Slater functions c * N * r**(n-1) * exp(-zeta*r)   .
C .  times the real spherical harmonic S(l,m) (same basis for all    .
C .  the MOs of the same symmetry , coefficients of Table I).         .
C .  For each part : name of its MO , l , m , occupation in the      .
C .  neutral molecule and in the ion , 1 if it is ionized.           .
C .  The ion has the hole in the ionized MO. The 3 protons are a    .
C .  charge 3 on the sphere r = 1.928 a.u.                           .
C .  iorb = 1 : 3a1 , 2 : 1e (1ex) , 3 : 2a1 ;                       .
C .  iorb = 0 : CH4 1t2 p part (target of ch4_opt_claude.for , test) .
C .  Eio : vertical ionization energies (photoelectron spectra) ,     .
C .  CHANGE THEM IF YOUR REFERENCE USES OTHER VALUES.                .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      SUBROUTINE TARGET(iorb,aname,acen,Eio,nsh,aqsh,arsh,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox,nshx)
      IMPLICIT NONE
      INTEGER iorb,nsh,norb,norbx,nstox,nshx
      CHARACTER*(*) aname
      CHARACTER*8 alab(norbx)
      DOUBLE PRECISION acen,Eio,aqsh(nshx),arsh(nshx),aocn(norbx)
     $,aoci(norbx),esz(nstox),csz(nstox)
      INTEGER lorb(norbx),morb(norbx),iono(norbx),ist(norbx),nst(norbx)
     $,nsn(nstox)
C .  NH3 : basis of the a1 MOs (s , d0 , p0 , f0 , f-3) and of the e
C .  MOs (d , p , f) : n and zeta
      INTEGER ns(5),nd(2),np(3),nf(1)
      DOUBLE PRECISION es(5),ed(2),ep(3),ef(1)
C .  coefficients (Table I) : a1 MOs (s , d0 , p0 , f0 , f-3)
      DOUBLE PRECISION c1s(5),c1d(2),c1p(3),c1f0(1),c1f3(1)
      DOUBLE PRECISION c2s(5),c2d(2),c2p(3),c2f0(1),c2f3(1)
      DOUBLE PRECISION c3s(5),c3d(2),c3p(3),c3f0(1),c3f3(1)
C .  1ex : d(m=1) , d(m=-2) , p(m=1) , f(m=1) , f(m=-2)
C .  1ey : d(m=2) , d(m=-1) , p(m=-1) , f(m=2) , f(m=-1)
      DOUBLE PRECISION cxd1(2),cxd2(2),cxp(3),cxf1(1),cxf2(1)
      DOUBLE PRECISION cyd2(2),cyd1(2),cyp(3),cyf2(1),cyf1(1)
C .  occupations (neutral , ion) and ionized flags of the 5 MOs
      DOUBLE PRECISION on(5),oi(5)
      INTEGER iz(5)
C .  CH4 (ch4_opt_claude.for , Moccia 1964 paper I)
      INTEGER mss(5),msp(3),msd(2),msf(1)
      DOUBLE PRECISION mes(5),mep(3),med(2),mef(1)
      DOUBLE PRECISION m1s(5),m1f(1),m2s(5),m2f(1),mtp(3),mtd(2),mtf(1)
      INTEGER i

      DATA ns /1,1,2,2,2/
      DATA es /11.000d0,6.400d0,1.750d0,1.280d0,2.560d0/
      DATA nd /3,3/
      DATA ed /1.600d0,2.350d0/
      DATA np /2,2,2/
      DATA ep /1.340d0,1.990d0,2.900d0/
      DATA nf /4/
      DATA ef /2.000d0/
C .  1a1
      DATA c1s /0.06572d0,0.93704d0,-0.01261d0,0.00524d0,0.01545d0/
      DATA c1d /0.00002d0,-0.00006d0/
      DATA c1p /-0.00164d0,0.00393d0,-0.00355d0/
      DATA c1f0 /0.00011d0/
      DATA c1f3 /-0.00020d0/
C .  2a1
      DATA c2s /0.01157d0,-0.23268d0,0.75114d0,0.12576d0,0.14793d0/
      DATA c2d /-0.07830d0,0.00659d0/
      DATA c2p /-0.14357d0,-0.01826d0,-0.00938d0/
      DATA c2f0 /0.04992d0/
      DATA c2f3 /-0.08013d0/
C .  3a1
      DATA c3s /0.00605d0,-0.06461d0,0.24313d0,-0.14177d0,0.07510d0/
      DATA c3d /-0.01440d0,-0.00699d0/
      DATA c3p /0.95405d0,-0.29504d0,0.40188d0/
      DATA c3f0 /-0.04098d0/
      DATA c3f3 /0.02420d0/
C .  1ex
      DATA cxd1 /-0.18794d0,0.03710d0/
      DATA cxd2 /-0.22929d0,0.05282d0/
      DATA cxp /1.00304d0,-0.28579d0,0.31169d0/
      DATA cxf1 /-0.04008d0/
      DATA cxf2 /0.06080d0/
C .  1ey
      DATA cyd2 /-0.22929d0,0.05282d0/
      DATA cyd1 /-0.18794d0,0.03710d0/
      DATA cyp /1.00304d0,-0.28579d0,0.31170d0/
      DATA cyf2 /0.06080d0/
      DATA cyf1 /-0.04008d0/
C .  CH4 (same data as ch4_opt_claude.for)
      DATA mss /1,1,2,4,4/
      DATA mes /9.5d0,5.5d0,1.5d0,2.d0,3.d0/
      DATA msp /2,3,4/
      DATA mep /1.373d0,2.950d0,2.95d0/
      DATA msd /4,4/
      DATA med /2.4d0,1.9d0/
      DATA msf /7/
      DATA mef /2.9d0/
      DATA m1s /0.05838d0,0.93837d0,0.0715d0,-0.0331d0,-0.03118d0/
      DATA m1f /0.00039d0/
      DATA m2s /.00877d00,-.21248d00,.98204d0,.05076d00,-.01799d00/
      DATA m2f /0.14254d0/
      DATA mtp /1.25998d00,-.05762d00,-.26738d00/
      DATA mtd /-.06691d00,.32775d00/
      DATA mtf /-.08695d00/

      norb=0

      if (iorb.eq.0) then
C . . . CH4 1t2 (regression test : ch4_opt_claude.for) . . . . . . . .
         aname='ch4test'
         acen=6.d0
         Eio=12.6d0
         nsh=1
         aqsh(1)=4.d0
         arsh(1)=2.08d0
         CALL ADDORB('1a1',0,0,2.d0,2.d0,0,5,mss,mes,m1s,norb,alab
     $   ,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
         CALL ADDORB('1a1',3,0,2.d0,2.d0,0,1,msf,mef,m1f,norb,alab
     $   ,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
         CALL ADDORB('2a1',0,0,2.d0,2.d0,0,5,mss,mes,m2s,norb,alab
     $   ,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
         CALL ADDORB('2a1',3,0,2.d0,2.d0,0,1,msf,mef,m2f,norb,alab
     $   ,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
C .  1t2 : 6 electrons (5 in the ion), ionized : its p part
         CALL ADDORB('1t2',1,0,6.d0,5.d0,1,3,msp,mep,mtp,norb,alab
     $   ,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
         CALL ADDORB('1t2',2,0,6.d0,5.d0,0,2,msd,med,mtd,norb,alab
     $   ,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
         CALL ADDORB('1t2',3,0,6.d0,5.d0,0,1,msf,mef,mtf,norb,alab
     $   ,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
         return
      endif

C . . . NH3 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      acen=7.d0
      nsh=1
      aqsh(1)=3.d0
      arsh(1)=1.928d0
C .  MOs : 1 = 1a1 , 2 = 2a1 , 3 = 3a1 , 4 = 1ex , 5 = 1ey
      do i=1,5
         on(i)=2.d0
         oi(i)=2.d0
         iz(i)=0
      enddo
      if (iorb.eq.1) then
C .  3a1 (HOMO) : vertical IP 10.85 eV
         aname='3a1'
         Eio=10.85d0
         oi(3)=1.d0
         iz(3)=1
      else if (iorb.eq.2) then
C .  1e (1ex) : vertical IP 16.4 eV
         aname='1e'
         Eio=16.4d0
         oi(4)=1.d0
         iz(4)=1
      else if (iorb.eq.3) then
C .  2a1 : vertical IP 27.0 eV
         aname='2a1'
         Eio=27.0d0
         oi(2)=1.d0
         iz(2)=1
      else
         stop 'TARGET : iorb must be 0, 1, 2 or 3'
      endif

C .  1a1
      CALL ADDORB('1a1',0,0,on(1),oi(1),iz(1),5,ns,es,c1s,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('1a1',2,0,on(1),oi(1),iz(1),2,nd,ed,c1d,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('1a1',1,0,on(1),oi(1),iz(1),3,np,ep,c1p,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('1a1',3,0,on(1),oi(1),0,1,nf,ef,c1f0,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('1a1',3,-3,on(1),oi(1),0,1,nf,ef,c1f3,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
C .  2a1
      CALL ADDORB('2a1',0,0,on(2),oi(2),iz(2),5,ns,es,c2s,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('2a1',2,0,on(2),oi(2),iz(2),2,nd,ed,c2d,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('2a1',1,0,on(2),oi(2),iz(2),3,np,ep,c2p,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('2a1',3,0,on(2),oi(2),0,1,nf,ef,c2f0,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('2a1',3,-3,on(2),oi(2),0,1,nf,ef,c2f3,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
C .  3a1
      CALL ADDORB('3a1',0,0,on(3),oi(3),iz(3),5,ns,es,c3s,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('3a1',2,0,on(3),oi(3),iz(3),2,nd,ed,c3d,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('3a1',1,0,on(3),oi(3),iz(3),3,np,ep,c3p,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('3a1',3,0,on(3),oi(3),0,1,nf,ef,c3f0,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('3a1',3,-3,on(3),oi(3),0,1,nf,ef,c3f3,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
C .  1ex
      CALL ADDORB('1ex',2,1,on(4),oi(4),iz(4),2,nd,ed,cxd1,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('1ex',2,-2,on(4),oi(4),iz(4),2,nd,ed,cxd2,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('1ex',1,1,on(4),oi(4),iz(4),3,np,ep,cxp,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('1ex',3,1,on(4),oi(4),0,1,nf,ef,cxf1,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('1ex',3,-2,on(4),oi(4),0,1,nf,ef,cxf2,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
C .  1ey
      CALL ADDORB('1ey',2,2,on(5),oi(5),iz(5),2,nd,ed,cyd2,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('1ey',2,-1,on(5),oi(5),iz(5),2,nd,ed,cyd1,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('1ey',1,-1,on(5),oi(5),iz(5),3,np,ep,cyp,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('1ey',3,2,on(5),oi(5),0,1,nf,ef,cyf2,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      CALL ADDORB('1ey',3,-1,on(5),oi(5),0,1,nf,ef,cyf1,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox)
      return
      end


C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .   ADDS ONE (l,m) PART (name of its MO , l , m , occupations       .
C .   neutral/ion , 1 if ionized , ns Slater functions nv , ev , cv)  .
C .   TO THE TABLES OF THE TARGET                                     .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      SUBROUTINE ADDORB(aname,l,m,aon,aoi,ionz,ns,nv,ev,cv
     $,norb,alab,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz
     $,norbx,nstox)
      IMPLICIT NONE
      CHARACTER*(*) aname
      INTEGER l,m,ionz,ns,nv(ns),norb,norbx,nstox
      DOUBLE PRECISION aon,aoi,ev(ns),cv(ns)
      CHARACTER*8 alab(norbx)
      INTEGER lorb(norbx),morb(norbx),iono(norbx),ist(norbx),nst(norbx)
     $,nsn(nstox)
      DOUBLE PRECISION aocn(norbx),aoci(norbx),esz(nstox),csz(nstox)
      INTEGER i
      norb=norb+1
      if (norb.gt.norbx) stop 'TARGET : too many (l,m) parts : norbx'
      if (norb.eq.1) then
         ist(norb)=1
      else
         ist(norb)=ist(norb-1)+nst(norb-1)
      endif
      if (ist(norb)+ns-1.gt.nstox) stop 'TARGET : too many Slater'
      alab(norb)=aname
      lorb(norb)=l
      morb(norb)=m
      aocn(norb)=aon
      aoci(norb)=aoi
      iono(norb)=ionz
      nst(norb)=ns
      do i=1,ns
         nsn(ist(norb)+i-1)=nv(i)
         esz(ist(norb)+i-1)=ev(i)
         csz(ist(norb)+i-1)=cv(i)
      enddo
      return
      end


C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .   Z(r) AT THE nd*nleg RADIAL NODES OF [0,d] (Azv/Azvi OF THE      .
C .   ORIGINAL) : iocc=1 NEUTRAL MOLECULE , iocc=2 ION                 .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      SUBROUTINE ZTAB(nd,d,dx,nleg,azvr,iocc,acen,nsh,aqsh,arsh,norb
     $,aocn,aoci,ist,nst,nsn,esz,csz)
      IMPLICIT DOUBLEPRECISION (a-h)
      IMPLICIT COMPLEX*16 (o-z)
      DOUBLEPRECISION ZSCR
      DIMENSION dx(nleg),azvr(nd,nleg),cd(nd+1),aqsh(*),arsh(*)
     $,aocn(*),aoci(*),ist(*),nst(*),nsn(*),esz(*),csz(*)
      do kd=1,nd+1
         cd(kd)=(kd-1.d0)*d/nd
      enddo
      do id=1,nd
      do i=1,nleg
         ar=(dx(i)*(cd(id+1)-cd(id))+cd(id+1)+cd(id))*0.5d00
         azvr(id,i)=ZSCR(ar,iocc,acen,nsh,aqsh,arsh,norb,aocn,aoci
     $   ,ist,nst,nsn,esz,csz)
      enddo
      enddo
      return
      end


C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .   Z(r) = CHARGE OF THE NUCLEI - SCREENING BY THE ELECTRONS AT ax  .
C .   (azw/azwi , az1a1.. OF Ch4.for FOR ANY TARGET) : FOR THE        .
C .   DENSITY r**N*exp(-c*r) (N=n1+n2 , c=zeta1+zeta2) OF A PRODUCT   .
C .   OF 2 SLATER FUNCTIONS OF THE SAME (l,m) PART , THE SCREENING IS .
C .      N!/c**(N+1) - J(N,c,r) + r*J(N-1,c,r)                         .
C .   J(n,c,r) = INTEGRAL FROM r TO INFINITY OF t**n*exp(-c*t)         .
C .   (THE PARTS WITH DIFFERENT (l,m) DO NOT MIX IN THE SPHERICAL     .
C .   AVERAGE). A SHELL OF CHARGE q AT THE RADIUS R GIVES q*min(r/R,1).
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      FUNCTION ZSCR(ax,iocc,acen,nsh,aqsh,arsh,norb,aocn,aoci
     $,ist,nst,nsn,esz,csz)
      IMPLICIT DOUBLEPRECISION (a-h)
      IMPLICIT COMPLEX*16 (o-z)
      DOUBLEPRECISION ZSCR
      DIMENSION aqsh(*),arsh(*),aocn(*),aoci(*),ist(*),nst(*),nsn(*)
     $,esz(*),csz(*)
      bzn=acen
      do is=1,nsh
         bzn=bzn+aqsh(is)*dmin1(ax/arsh(is),1.d0)
      enddo
      bze=0.d0
      do io=1,norb
         if (iocc.eq.1) then
            aoc=aocn(io)
         else
            aoc=aoci(io)
         endif
         if (aoc.ne.0.d0) then
            bor=0.d0
            do i=ist(io),ist(io)+nst(io)-1
            do j=ist(io),ist(io)+nst(io)-1
               nn=nsn(i)+nsn(j)
               c=esz(i)+esz(j)
               bv=fac(nn)/(c**(nn+1))-avJn(c,ax,nn)
               bv=bv+ax*avJn(c,ax,nn-1)
               bor=bor+bv*(csz(i)*fa(esz(i),nsn(i)-1))
     $            *(csz(j)*fa(esz(j),nsn(j)-1))
            enddo
            enddo
            bze=bze+aoc*bor
         endif
      enddo
      ZSCR=bzn-bze
      end


C .   R(ar) OF THE (l,m) PART io (Frwf OF THE ORIGINAL)
      FUNCTION RADF(ar,io,ist,nst,nsn,esz,csz)
      IMPLICIT DOUBLEPRECISION (a-h)
      IMPLICIT COMPLEX*16 (o-z)
      DOUBLEPRECISION RADF
      DIMENSION ist(*),nst(*),nsn(*),esz(*),csz(*)
      bf=0.d0
      do i=ist(io),ist(io)+nst(io)-1
         bf=bf+csz(i)*(ar**(nsn(i)-1))*dexp(-esz(i)*ar)
     $     *fa(esz(i),nsn(i)-1)
      enddo
      RADF=bf
      end


C .   NORM OF THE (l,m) PART io : SUM c1*c2*N1*N2*(n1+n2)!/
C .   (zeta1+zeta2)**(n1+n2+1)
      FUNCTION ORBNRM(io,ist,nst,nsn,esz,csz)
      IMPLICIT DOUBLEPRECISION (a-h)
      IMPLICIT COMPLEX*16 (o-z)
      DOUBLEPRECISION ORBNRM
      DIMENSION ist(*),nst(*),nsn(*),esz(*),csz(*)
      bn=0.d0
      do i=ist(io),ist(io)+nst(io)-1
      do j=ist(io),ist(io)+nst(io)-1
         nn=nsn(i)+nsn(j)
         bn=bn+csz(i)*fa(esz(i),nsn(i)-1)*csz(j)*fa(esz(j),nsn(j)-1)
     $     *fac(nn)/((esz(i)+esz(j))**(nn+1))
      enddo
      enddo
      ORBNRM=bn
      end


C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .   ANALYTIC PART OF THE EJECTED ELECTRON TERM FOR ANY s , p OR d    .
C .   SLATER FUNCTION r**(n-1)*exp(-zeta*r) : tcs(zeta,az,ad,k,n)      .
C .   (s) , tcp(zeta,az,ad,k,n,m) (p) AND tcd(zeta,az,ad,k,n,m) (d)    .
C .   OF THE ORIGINAL , ANY n (from BBK3CWZ_CH4_opt.f , with the      .
C .   intrinsic exp/log/sin/cos/atan2).                                .
C .   THE TERMS ARE LEIBNIZ SUMS OF THE DERIVATIVES (WITH RESPECT TO  .
C .   zeta) OF qa**ya AND qab**yb :                                    .
C .      L(k) = SUM(j=0..k) C(k,j) d^j(qa**ya) d^(k-j)(qab**yb)       .
C .   s : tcs = (-1)**n * aps * L(n)       ya=qz-1 , yb=-qz            .
C .   p : tcp = (-1)**(n-1)*xi*app*SUM(s=0..2) C(n-1,s)*L(n-1-s)*tp_s .
C .                                        ya=qz-2 , yb=-qz-1          .
C .   d : tcd = (-1)**(n-1)*apd*SUM(s=0..4) C(n-2,s)*L(n-2-s)*td_s    .
C .                                        ya=qz-3 , yb=-qz-2          .
C .   X**y WITH X = u**2+CONST (dX/dzeta = 2u) : THE SCALED RATIOS    .
C .   a(k) = d^k(X**y)/(k!*X**y)  FOLLOW                              .
C .   (k+1)*X*a(k+1) = 2u*(y-k)*a(k) + (2*y-k+1)*a(k-1)               .
C .   (u=zeta FOR qa , u=zeta-i*|k| FOR qab) , SO L(k) = k!*SUM a*b , .
C .   THE k! AND THE SIGN ARE IN cfv = (-1)**(n-l)*(n-l)!*c*N(n,zeta). .
C .   TFGB : THE PART WHICH DEPENDS ONLY ON |ad| AND |k| ;             .
C .   TFGA : THE REST , FOR ad AND k.                                  .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

C .   TFGB : FOR THE nts SLATER TERMS (l , n , zeta) , aq=|ad|**2 ,
C .   al1=|k|**2 : qbv=qab , pbv=qab**yb , qbk(k)=d^k(qab**yb)/
C .   (k!*qab**yb)
      SUBROUTINE TFGB(az,aq,al1,nts,ltv,nsv,esv,kmx,qbv,pbv,qbk)
      IMPLICIT DOUBLEPRECISION (a-h)
      IMPLICIT COMPLEX*16 (o-z)
      DIMENSION ltv(nts),nsv(nts),esv(nts),qbv(nts),pbv(nts)
     $,qbk(0:kmx,nts)
      xi=(0.0d00,1.0d00)
      ak1=dsqrt(al1)
      qz=az*xi/ak1
      eta=dimag(qz)
      do it=1,nts
         r=esv(it)
         x=r-xi*ak1
         qab=aq+r*r-al1-2.*xi*r*ak1
         qbv(it)=qab
         bre=dreal(qab)
         bim=dimag(qab)
         ab2=bre*bre+bim*bim
C .  qab**(-qz) = exp(eta*arg(qab))*exp(-i*eta*log|qab|)
         aev=dexp(eta*datan2(bim,bre))
         alv=-eta*0.5d0*dlog(ab2)
         qib1=dconjg(qab)/ab2
C .  qab**(-qz-l) (s : l=0 , p : l=1 , d : l=2)
         l=ltv(it)
         yb=-qz-l
         qpb=aev*dcmplx(dcos(alv),dsin(alv))
         do k=1,l
            qpb=qpb*qib1
         enddo
         pbv(it)=qpb
         ncap=nsv(it)-l
         qbk(0,it)=(1.d0,0.d0)
         if (ncap.ge.1) qbk(1,it)=2.*x*yb*qib1
         do k=1,ncap-1
            qbk(k+1,it)=(2.d0*x*(yb-k)*qbk(k,it)
     $      +(2.d0*yb-(k-1))*qbk(k-1,it))*(qib1/(k+1))
         enddo
      enddo
      return
      end


C .   TFGA : THE ANALYTIC PART tfv(icm) OF THE (l,m) icm OF THE ncm
C .   (l,m) OF THE nio IONIZED PARTS (s : m=0 , p : m=0,1,-1 ,
C .   d : m=0,1,-1,2,-2) FOR ad AND k , qbv/pbv/qbk FROM TFGB (SAME
C .   |ad| AND |k|) , cfv = (-1)**(n-l)*(n-l)!*c*N(n,zeta) , iov =
C .   IONIZED PART OF THE TERM.
C .   p : tcp = xi*app*(va2*qau+vb2*qz*qav) (va2 , vb2 , app : m)
C .   d : tcd = -apd*(zd1*qd1+zd2*qd2+zd3*qd3) (zd1 , zd2 , zd3 : m)
      SUBROUTINE TFGA(az,adx,ady,adz,akx,aky,akz,nts,ltv,nsv,esv,cfv
     $,iov,kmx,qbv,pbv,qbk,nio,ncm,iocm,lcm,mcm,tfv)
      IMPLICIT DOUBLEPRECISION (a-h)
      IMPLICIT COMPLEX*16 (o-z)
      DIMENSION ltv(nts),nsv(nts),esv(nts),cfv(nts),iov(nts),qbv(nts)
     $,pbv(nts),qbk(0:kmx,nts),iocm(ncm),lcm(ncm),mcm(ncm),tfv(ncm)
C .  LOCAL ARRAYS OF FIXED SIZE (niox = norbx OF THE MAIN PROGRAM ,
C .  n-l <= 12)
      PARAMETER (niox=40)
      DIMENSION qak(0:12)
     $,qau(niox),qav(niox),qss(niox),qc1(12,0:2),qc2(12,0:2),qya(0:2)
     $,ark(12),qd1(niox),qd2(niox),qd3(niox),qls(0:4)
C .  ark(k) = 1/k
      DATA ark /1.d0,0.5d0,0.333333333333333333333d0,0.25d0,0.2d0
     $,0.166666666666666666667d0,0.142857142857142857143d0,0.125d0
     $,0.111111111111111111111d0,0.1d0,0.0909090909090909090909d0
     $,0.0833333333333333333333d0/
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
C .  6.**.5 , 3.**.5 , 2.**.5 , 10.**.5 ... : SINGLE PRECISION AS IN
C .  THE ORIGINAL tcs/tcp/tcd
      app0=6.**.5/api
      app1=3.**.5/api
      aps=1./(2.**.5*api)
      apd0=10.**.5/api
      apd1=60.**.5/api
      apd2=15.**.5/api
      al1=akx*akx+aky*aky+akz*akz
      ak1=dsqrt(al1)
      aq=adx*adx+ady*ady+adz*adz
      al3=aq+al1-2.*(adx*akx+ady*aky+adz*akz)
      qz=az*xi/ak1
      eta=dimag(qz)
      qz5=2.*xi*ak1
      qz1=1.d0-qz
      qz51=qz5*qz1
C .  COEFFICIENTS OF THE RECURRENCE (qa SIDE) : ya=qz-1-l
C .  a(k+1) = (zeta*qc1(k)*a(k) + qc2(k)*a(k-1))/qa
      do l=0,2
         yal=qz-1.d0-l
         qya(l)=yal
         do k=1,kmx-1
            qc1(k,l)=(2.d0*ark(k+1))*(yal-k)
            qc2(k,l)=ark(k+1)*(2.d0*yal-(k-1))
         enddo
      enddo
      do je=1,nio
         qau(je)=(0.d0,0.d0)
         qav(je)=(0.d0,0.d0)
         qss(je)=(0.d0,0.d0)
         qd1(je)=(0.d0,0.d0)
         qd2(je)=(0.d0,0.d0)
         qd3(je)=(0.d0,0.d0)
      enddo
      do it=1,nts
         ae=esv(it)
         l=ltv(it)
         ncap=nsv(it)-l
         je=iov(it)
C .  qa (REAL) AND exp(i*eta*log(qa))
         aqa=al3+ae*ae
         alv=eta*dlog(aqa)
         ca1=1.d0/aqa
C .  cfv*P = cfv * qa**ya * qab**yb (qa**(-1-l) , ya=qz-1-l)
         if (l.eq.0) then
            qpw=(cfv(it)*ca1)*dcmplx(dcos(alv),dsin(alv))*pbv(it)
         else if (l.eq.1) then
            qpw=(cfv(it)*ca1*ca1)*dcmplx(dcos(alv),dsin(alv))*pbv(it)
         else
            qpw=(cfv(it)*ca1*ca1*ca1)*dcmplx(dcos(alv),dsin(alv))
     $         *pbv(it)
         endif
C .  a(k) = d^k(qa**ya)/(k!*qa**ya) , a(0)=1
         qak(1)=(2.d0*ae*ca1)*qya(l)
         if (ncap.ge.2) qak(2)=((ae*qc1(1,l))*qak(1)+qc2(1,l))*ca1
         do k=2,ncap-1
            qak(k+1)=((ae*qc1(k,l))*qak(k)+qc2(k,l)*qak(k-1))*ca1
         enddo
C .  L(m)/m! = SUM(j=0..m) a(j)*b(m-j) = a(m)+b(m)+SUM(j=1..m-1)
C .  (m >= 1 , a(0)=b(0)=1)
         ql0=qak(ncap)+qbk(ncap,it)
         do j=1,ncap-1
            ql0=ql0+qak(j)*qbk(ncap-j,it)
         enddo
         if (l.eq.0) then
C .  s : (-1)**n * L(n)
            qss(je)=qss(je)+qpw*ql0
         else if (l.eq.2) then
C .  d : SUM(s=0..4) L(n-2-s)/(n-2-s)! * td_s/s! , qls(s)=L(n-2-s)/..
            qls(0)=ql0
            do is=1,4
               j=ncap-is
               if (j.lt.0) then
                  qls(is)=(0.d0,0.d0)
               else if (j.eq.0) then
                  qls(is)=(1.d0,0.d0)
               else
                  ql=qak(j)+qbk(j,it)
                  do i=1,j-1
                     ql=ql+qak(i)*qbk(j-i,it)
                  enddo
                  qls(is)=ql
               endif
            enddo
C .  td_s/s! = A_s*zd1 + B_s*zd2 + C_s*zd3 (x=qa , y=qb=qab-qa ,
C .  z4=zeta , z5=qz5 IN td0..td4 OF THE ORIGINAL)
            qb=qbv(it)-aqa
            qd1(je)=qd1(je)+qpw*(qls(0)*(aqa*aqa)+qls(1)*(4.d0*aqa*ae)
     $      +qls(2)*(2.d0*aqa+4.d0*ae*ae)+qls(3)*(4.d0*ae)+qls(4))
            qd2(je)=qd2(je)+qpw*(-qls(0)*(aqa*qb)
     $      +qls(1)*(aqa*qz5-2.d0*ae*qb)+qls(2)*(2.d0*ae*qz5-qb)
     $      +qls(3)*qz5)
            qd3(je)=qd3(je)+qpw*(qls(0)*(qb*qb)-qls(1)*(2.d0*qb*qz5)
     $      +qls(2)*(qz5*qz5))
         else
C .  p : (-1)**(n-1)*( L(n-1)*tp0 + (n-1)*L(n-2)*tp1
C .      + (n-1)(n-2)/2*L(n-3)*tp2 ) , qt0b=qa*qz , qt1b=2*zeta*qz
C .      (qz IS FACTORISED IN qav)
            m1=ncap-1
            m2=ncap-2
            if (m1.eq.0) then
               ql1=(1.d0,0.d0)
               ql2=(0.d0,0.d0)
            else
               ql1=qak(m1)+qbk(m1,it)
               do j=1,m1-1
                  ql1=ql1+qak(j)*qbk(m1-j,it)
               enddo
               if (m2.eq.0) then
                  ql2=(1.d0,0.d0)
               else
                  ql2=qak(m2)+qbk(m2,it)
                  do j=1,m2-1
                     ql2=ql2+qak(j)*qbk(m2-j,it)
                  enddo
               endif
            endif
            qt0a=aqa+(qbv(it)-aqa)*qz1
            qt1a=2.d0*ae-qz51
            qu=ql0*qt0a+ql1*qt1a+ql2
            qv=ql0*aqa+ql1*(2.d0*ae)+ql2
            qau(je)=qau(je)+qpw*qu
            qav(je)=qav(je)+qpw*qv
         endif
      enddo
C .  THE (l,m) : s -> aps*SUM , p -> xi*app*(va2*qau+vb2*qz*qav)
C .  m=0 : va2=adz-akz , vb2=akz   m=1 : va2=qk1x , vb2=-q01x
C .  m=-1 : va2=qk1y , vb2=-q01y   (AS IN tcp)
      aqkz=adz-akz
      q01x=akx+xi*aky
      qk1x=-adx-xi*ady+akx+xi*aky
      q01y=-akx+xi*aky
      qk1y=adx-xi*ady-akx+xi*aky
      qz03=qz*qz-3.*qz+2.
      do icm=1,ncm
         je=iocm(icm)
         m=mcm(icm)
         if (lcm(icm).eq.0) then
            tfv(icm)=aps*qss(je)
         else if (lcm(icm).eq.2) then
C .  d : va3 , vb3 , vc3 , apd OF tcd FOR m
            if (m.eq.0) then
               akp2=aq-al3-al1
               va3=3.*(aqkz*aqkz)-al3
               vb3=6.*aqkz*akz-akp2
               vc3=3.*(akz*akz)-al1
               apd=apd0
            else if (m.eq.1) then
               va3=aqkz*qk1x
               vb3=-aqkz*q01x+akz*qk1x
               vc3=-akz*q01x
               apd=apd1
            else if (m.eq.-1) then
               va3=aqkz*qk1y
               vb3=-aqkz*q01y+akz*qk1y
               vc3=-akz*q01y
               apd=apd1
            else if (m.eq.2) then
               va3=qk1x*qk1x
               vb3=-2.*q01x*qk1x
               vc3=q01x*q01x
               apd=apd2
            else
               va3=qk1y*qk1y
               vb3=-2.*q01y*qk1y
               vc3=q01y*q01y
               apd=apd2
            endif
            zd1=2.*va3+2.*qz*vb3+vc3*qz*(qz+1.)
            zd2=(qz-1.)*(4.*va3+qz*vb3)
            zd3=qz03*va3
            tfv(icm)=-apd*(zd1*qd1(je)+zd2*qd2(je)+zd3*qd3(je))
         else
            if (m.eq.0) then
               qva=adz-akz
               qvb=akz
               app=app0
            else if (m.eq.1) then
               qva=qk1x
               qvb=-q01x
               app=app1
            else
               qva=qk1y
               qvb=-q01y
               app=app1
            endif
            tfv(icm)=xi*app*(qva*qau(je)+qvb*(qav(je)*qz))
         endif
      enddo
      return
      end


C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .  tcs1 = tcs(a,az,...,n=1). tcs computed the 12 orders us(1..12)    .
C .  (~150 complex powers) at every call and returned us(n); it is     .
C .  only called with n=1. Same trick as TFP3 for the powers:          .
C .    us(1) = -(q1(qa,qz-1,r)*q0(qab,-qz) + q0(qa,qz-1)*q1(qab,-qz,x)) .
C .            * aps                                                   .
C .  Verified against tcs(...,1): max relative difference 3e-15.       .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	FUNCTION tcs1(a,az,adx,ady,adz,akx,aky,akz)
	IMPLICIT NONE
	COMPLEX*16 tcs1
	DOUBLE PRECISION a,az,adx,ady,adz,akx,aky,akz
	DOUBLE PRECISION api,aps,al1,ak1,aq,al3,r,cz,qa,lqa,lr,li
	COMPLEX*16 xi,qz,x,qab,pa,pb,us10,us11

	api=3.141592654d00
	xi=(0.0d00,1.0d00)
	aps=1./(2.**.5*api)

	al1=akx*akx+aky*aky+akz*akz
	ak1=dsqrt(al1)
	aq=adx*adx+ady*ady+adz*adz
	al3=aq+al1-2.*(adx*akx+ady*aky+adz*akz)
	cz=az/ak1
	qz=cz*xi
	r=a
	x=r-xi*ak1
	qa=al3+r*r
	qab=aq+r*r-al1-2.*xi*r*ak1

C . us10=q1(qa,qz-1,r)*q0(qab,-qz) ; us11=q0(qa,qz-1)*q1(qab,-qz,x)
C . pa=qa**(qz-1) (qa real>0) ; pb=qab**(-qz) (principal log)
	lqa=dlog(qa)
	pa=dcmplx(dcos(cz*lqa),dsin(cz*lqa))/qa
	lr=0.5d0*dlog(dreal(qab)**2+dimag(qab)**2)
	li=datan2(dimag(qab),dreal(qab))
	pb=dexp(cz*li)*dcmplx(dcos(cz*lr),-dsin(cz*lr))
	us10=2.*r*(qz-1.)*(pa/qa)*pb
	us11=pa*2.*x*(-qz)*(pb/qab)
	tcs1=-(us10+us11)*aps
	end

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   LIBRARY OF ch4_opt_claude.for / Ch4.for (unchanged). The CH4    .
C .   specific routines (Azv, Azvi, Frwfp/d/f, azw, azwi, az1a1,      .
C .   az2a1, az1t2z, TFP, TFD, TFF, TFP3) are replaced by TARGET,     .
C .   ZTAB/ZSCR, RADF and TFGB/TFGA above.                            .
C .   tcp, tcs, TECWN, TSCWN are kept for reference/validation but    .
C .   are not called by the program above.                            .
C .                                                                   .
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

C .																	.
C .	Calculating Qffs scattered electron Colomb wave function		.
C .																	.
	SUBROUTINE QSCWF(qf,azc,akx,aky,akz,dx,n1,n2,n3,d,nd,ucd1)
	IMPLICIT DOUBLEPRECISION (a-h)
	IMPLICIT COMPLEX*16 (o-z)
	DIMENSION dx(n1),qf(nd,n1,n2,n3),azc(nd,n1),cd(nd+1)



	api=dacos(-1.d00)
      xi=(0.0d00,1.0d00)
      ak=dsqrt(akx**2+aky**2+akz**2)

      do kd=1,nd+1
	cd(kd)=(kd-1.d0)*d/nd
	enddo

	do id=1,nd
	   do i=1,n1
	   do j=1,n2
	   do k=1,n3
	   ar=(dx(i)*(cd(id+1)-cd(id))+cd(id+1)+cd(id))*0.5d00
	   et=(api/2.d0)*(dx(j)+1.d0)
	   ef=api*(dx(k)+1.d0)
C .																	.
C . Les Arguments de r (rx,ry,rz)										.
C .																	.
	arx=ar*dsin(et)*dcos(ef)
	ary=ar*dsin(et)*dsin(ef)
	arz=ar*dcos(et)
C .																	.
C .   akr ===> ke.r (Vector)								        	.
C .  			    adr ===> k.r  (Vector)								.
C .                  amkr ===> ke.r (Scalar)							.
C .																	. 
      qkx=-xi*(akx*arx+aky*ary+akz*arz+ak*ar)
      qrx=xi*(akx*arx+aky*ary+akz*arz)



C . ARGUMENTS OF UC AND UF										    .
	alp=azc(id,i)/ak											
      qzk=-xi*alp
C .																	.
C .																	.
      q1=1.d00+xi*0.d00
      szz=1.d00-qzk
C .
C .  COLOMB FUNCTION
C .

	qf(id,i,j,k)=cdexp(qrx+api*alp/2.d00)*uf(qzk,q1,qkx)*uc(szz)/ucd1
	    enddo
	    enddo
	    enddo
	enddo
	return
	end
C .																	.
C .	Calculating QffI Incident electron Colomb wave function		.
C .																	.
	SUBROUTINE QICWF(qf,azc,akx,aky,akz,dx,n1,n2,n3,d,nd,ucd1)
	IMPLICIT DOUBLEPRECISION (a-h)
	IMPLICIT COMPLEX*16 (o-z)
	DIMENSION dx(n1),qf(nd,n1,n2,n3),azc(nd,n1),cd(nd+1)



	api=dacos(-1.d00)
      xi=(0.0d00,1.0d00)
      ak=dsqrt(akx**2+aky**2+akz**2)

      do kd=1,nd+1
	cd(kd)=(kd-1.d0)*d/nd
	enddo

	do id=1,nd
	   do i=1,n1
	   do j=1,n2
	   do k=1,n3
	   ar=(dx(i)*(cd(id+1)-cd(id))+cd(id+1)+cd(id))*0.5d00
	   et=(api/2.d0)*(dx(j)+1.d0)
	   ef=api*(dx(k)+1.d0)
C .																	.
C . Les Arguments de r (rx,ry,rz)										.
C .																	.
	arx=ar*dsin(et)*dcos(ef)
	ary=ar*dsin(et)*dsin(ef)
	arz=ar*dcos(et)
C .																	.
C .   akr ===> ke.r (Vector)								        	.
C .  			    adr ===> k.r  (Vector)								.
C .                  amkr ===> ke.r (Scalar)							.
C .																	. 
      qkx=-xi*(akx*arx+aky*ary+akz*arz-ak*ar)
      qrx=xi*(akx*arx+aky*ary+akz*arz)



C . ARGUMENTS OF UC AND UF										    .
	alp=-azc(id,i)/ak											
      qzk=-xi*alp
C .																	.
C .																	.
      q1=1.d00+xi*0.d00
	szz=1.d00+xi*alp

C .
C .  COLOMB FUNCTION
C .

	qf(id,i,j,k)=cdexp(qrx-api*alp/2.d00)*uf(qzk,q1,qkx)*uc(szz)/ucd1
	    enddo
	    enddo
	    enddo
	enddo
	return
	end	
		
C .																	.
C .	Calculating Qffs scattered electron Colomb wave function Z=1	.
C .																	.
	SUBROUTINE Q1SCWF(qf,akx,aky,akz,dx,n1,n2,n3,d,nd)
	IMPLICIT DOUBLEPRECISION (a-h)
	IMPLICIT COMPLEX*16 (o-z)
	DIMENSION dx(n1),qf(nd,n1,n2,n3),cd(nd+1)

	api=dacos(-1.d00)
      xi=(0.0d00,1.0d00)
      ak=dsqrt(akx**2+aky**2+akz**2)

      do kd=1,nd+1
	cd(kd)=(kd-1.d0)*d/nd
	enddo

	do id=1,nd
	   do i=1,n1
	   do j=1,n2
	   do k=1,n3
	   ar=(dx(i)*(cd(id+1)-cd(id))+cd(id+1)+cd(id))*0.5d00
	   et=(api/2.d0)*(dx(j)+1.d0)
	   ef=api*(dx(k)+1.d0)
C .																	.
C . Les Arguments de r (rx,ry,rz)										.
C .																	.
	arx=ar*dsin(et)*dcos(ef)
	ary=ar*dsin(et)*dsin(ef)
	arz=ar*dcos(et)
C .																	.
C .   akr ===> ke.r (Vector)								        	.
C .  			    adr ===> k.r  (Vector)								.
C .                  amkr ===> ke.r (Scalar)							.
C .																	. 
      qkx=-xi*(akx*arx+aky*ary+akz*arz+ak*ar)
      qrx=xi*(akx*arx+aky*ary+akz*arz)



C . ARGUMENTS OF  UF										    	   .
	alp=1.d0/ak											
      qzk=-xi*alp																
      q1=1.d00+xi*0.d00
C .
C .  COLOMB FUNCTION
C .

	qf(id,i,j,k)=cdexp(qrx)*uf(qzk,q1,qkx)
	    enddo
	    enddo
	    enddo
	enddo
	return
	end

C .																	.
C .	Calculating Q0ffi Incident electron Colomb wave function		.
C .																	.
	SUBROUTINE Q0SCWF(qf,akx,aky,akz,dx,n1,n2,n3,d,nd)
	IMPLICIT DOUBLEPRECISION (a-h)
	IMPLICIT COMPLEX*16 (o-z)
	DIMENSION dx(n1),qf(nd,n1,n2,n3),cd(nd+1)



	api=dacos(-1.d00)
      xi=(0.0d00,1.0d00)
      ak=dsqrt(akx**2+aky**2+akz**2)

      do kd=1,nd+1
	cd(kd)=(kd-1.d0)*d/nd
	enddo

	do id=1,nd
	   do i=1,n1
	   do j=1,n2
	   do k=1,n3
	   ar=(dx(i)*(cd(id+1)-cd(id))+cd(id+1)+cd(id))*0.5d00
	   et=(api/2.d0)*(dx(j)+1.d0)
	   ef=api*(dx(k)+1.d0)
C .																	.
C . Les Arguments de r (rx,ry,rz)										.
C .																	.
	arx=ar*dsin(et)*dcos(ef)
	ary=ar*dsin(et)*dsin(ef)
	arz=ar*dcos(et)

      qrx=xi*(akx*arx+aky*ary+akz*arz)
C .
C .  COLOMB FUNCTION
C .
	qf(id,i,j,k)=cdexp(qrx)
	    enddo
	    enddo
	    enddo
	enddo
	return
	end
C .																	.
C .	Calculating QEff scattered electron Colomb wave function		.
C .																	.
      SUBROUTINE QECWF(qff,azc,akx,aky,akz,dx,n1,n2,n3,d,nd,ucd1)
      IMPLICIT DOUBLEPRECISION (a-h)            
      IMPLICIT COMPLEX*16 (o-z)
      DIMENSION dx(n1),qff(nd,n1,n2,n3),azc(nd,n1),cd(nd+1)



      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      ak=dsqrt(akx**2+aky**2+akz**2)


      do kd=1,nd+1
	cd(kd)=(kd-1.d0)*d/nd
	enddo

      do id=1,nd
         do i=1,n1
         do j=1,n2
         do k=1,n3

         ar=(dx(i)*(cd(id+1)-cd(id))+cd(id+1)+cd(id))/2.d00
	   et=(api/2.d0)*(dx(j)+1.d0)
	   ef=api*(dx(k)+1.d0)
     
C .																	.
C . Les Arguments de r (rx,ry,rz)										.
C .																	.
	arx=ar*dsin(et)*dcos(ef)
	ary=ar*dsin(et)*dsin(ef)
	arz=ar*dcos(et)


      qkx=-xi*(akx*arx+aky*ary+akz*arz+ak*ar)
      qrx=xi*(akx*arx+aky*ary+akz*arz)

      q1=1.d00+xi*0.d00
c
      alp=azc(id,i)/ak
      qzk=-xi*alp
      sz2=1.d00-qzk
c

      alw=1.d0/ak
      qzw=-xi*alw

      qf2=cdexp(qrx+api*alp/2.d00)*uf(qzk,q1,qkx)*uc(sz2)/ucd1
      qf1=cdexp(qrx)*uf(qzw,q1,qkx)


      qff(id,i,j,k)=qf2-qf1

      enddo
      enddo
      enddo
      enddo

      return
      end

C .																	.
C .                						                            .
C .  																	.

      FUNCTION TECWN(qff,fg,adx,ady,adz,n1,n2,n3,dx,dw,d,nd,l,am)
      IMPLICIT DOUBLEPRECISION (a-h)            
      IMPLICIT COMPLEX*16 (o-z)

      DIMENSION dx(n1),dw(n1),cd(nd+1)
     $         ,qff(nd,n1,n2,n3)
     $		 ,fg(nd,n1)
      xi=(0.0d00,1.d00)


      api=3.141592654d00

      do kd=1,nd+1
 	cd(kd)=(kd-1.)*d/nd
      enddo 

      qs=(0.00d0,0.00d0)
      do id=1,nd
      do i=1,n1
      do j=1,n2
      do k=1,n3

      ar=(dx(i)*(cd(id+1)-cd(id))+cd(id+1)+cd(id))*0.5d0
      et=api*(dx(j)+1.d00)*0.5d0
      ef=api*(dx(k)+1.d00)

	arx=ar*dsin(et)*dcos(ef)
	ary=ar*dsin(et)*dsin(ef)
	arz=ar*dcos(et)

	adr=adx*arx+ady*ary+adz*arz
	his=fg(id,i)
	atk=dcos(et)
	qylm=ylm(l,am,atk,ef)

C ..
C ..            ax,ay ==> cos(theta), phi respectivement
C ..			  qylm ==> Les harmoniques Spheriques
C ..

	tf=dconjg(qff(id,i,j,k))*his*cdexp(xi*adr)*qylm
C  .. L'INTEGRALE
	qs=qs+dw(i)*dw(j)*dw(k)*dsin(et)*ar*ar
     $*(api**2)*tf*(cd(id+1)-cd(id))/4.d0	 
      enddo
      enddo
      enddo
      enddo

      TECWN=qs/dsqrt((2.d00*api)**3)
      return
      end
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  . .
      FUNCTION TSCWN(alp,qffi,qffs,q0ffi,q1ffs,dx,dw
     $,adx,ady,adz,n1,n2,n3,d,nd,nn,mn)
      IMPLICIT DOUBLEPRECISION (a-h)            
      IMPLICIT COMPLEX*16 (o-z)

      DIMENSION dx(n1),dw(n1),cd(nd+1),qffi(nd,n1,n2,n3)
     $         ,qffs(nd,n1,n2,n3),q0ffi(nd,n1,n2,n3),q1ffs(nd,n1,n2,n3)



      api=3.141592654d00
      xi=(0.0d00,1.d00)

      do kd=1,nd+1
 	cd(kd)=(kd-1.)*d/nd
      enddo 

	qs=(0.d0,0.d0)
      do id=1,nd
	  do i=1,n1
	  do j=1,n2
	  do k=1,n3

		ar=(dx(i)*(cd(id+1)-cd(id))+cd(id+1)+cd(id))/2.d00
		et=api*(dx(j)+1.d00)/2.d00
		ef=api*(dx(k)+1.d00)

      arx=ar*dsin(et)*dcos(ef)
      ary=ar*dsin(et)*dsin(ef)
      arz=ar*dcos(et)

      adr=adx*arx+ady*ary+adz*arz

      qdr=cdexp(xi*adr)

	qf=dconjg(qffs(id,i,j,k))*qffi(id,i,j,k)
     $  -dconjg(q1ffs(id,i,j,k))*q0ffi(id,i,j,k)

      qs=qs+qf*qdr*dexp(-alp*ar)
     $*dw(i)*dw(j)*dw(k)*dsin(et)*(ar**nn)
     $*(api**2)*(cd(id+1)-cd(id))/4.d0

      enddo
      enddo
      enddo
      enddo
 	ay00=1.d0/dsqrt(4.d0*api)
      TSCWN=(ay00**mn)*qs/dsqrt((2.d00*api)**3)
      return
      end

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








      SUBROUTINE GAULEG(X1,X2,X,W,N,EPS)
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C .                                                                   .
C .   P R O G R A M                                                   .
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
      SUBROUTINE GAULAG(X,W,N,ALF)
      INTEGER N,MAXIT
      DOUBLE PRECISION ALF,W(N),X(N)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.D-14,MAXIT=10)
CU    USES GAMMLN
      INTEGER I,ITS,J
      DOUBLE PRECISION AI,GAMMLN
      DOUBLE PRECISION P1,P2,P3,PP,Z,Z1
      DO 13 I=1,N
        IF(I.EQ.1)THEN
          Z=(1.D0+ALF)*(3.D0+.92D0*ALF)/(1.D0+2.4D0*N+1.8D0*ALF)
        ELSE IF(I.EQ.2)THEN
          Z=Z+(15.D0+6.25D0*ALF)/(1.D0+.9D0*ALF+2.5D0*N)
        ELSE
          AI=I-2
          Z=Z+((1.D0+2.55D0*AI)/(1.9D0*AI)+1.26D0*AI*ALF/(1.D0+3.5D0*AI)
     *)*
     *(Z-X(I-2))/(1.D0+.3D0*ALF)
        ENDIF
        DO 12 ITS=1,MAXIT
          P1=1.D0
          P2=0.D0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2*J-1+ALF-Z)*P2-(J-1+ALF)*P3)/J
11        CONTINUE
          PP=(N*P1-(N+ALF)*P2)/Z
          Z1=Z
          Z=Z1-P1/PP
          IF(ABS(Z-Z1).LE.EPS)GOTO 1
12      CONTINUE
        PAUSE 'TOO MANY ITERATIONS IN GAULAG'
1       X(I)=Z
        W(I)=-EXP(GAMMLN(ALF+N)-GAMMLN( DBLE(N)))/(PP*N*P2)
13    CONTINUE
      RETURN
      END
C
      FUNCTION GAMMLN(XX)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER,XX,GAMMLN
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c fonction factorielle
      function fac(n)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      dimension factab(0:24)
      data factab /1.d0,1.d0,2.d0,6.d0,24.d0,120.d0,720.d0,5040.d0,
     $ 40320.d0,362880.d0,3628800.d0,39916800.d0,479001600.d0,
     $ 6227020800.d0,87178291200.d0,1307674368000.d0,
     $ 20922789888000.d0,355687428096000.d0,6402373705728000.d0,
     $ 121645100408832000.d0,2432902008176640000.d0,
     $ 51090942171709440000.d0,1124000727777607680000.d0,
     $ 25852016738884976640000.d0,620448401733239439360000.d0/
      if(n.le.24) then
        fac=factab(n)
      else
        fac=1.d0
        do i=2,n
          fac=fac*i
        enddo
      endif
      end

c fonction de normation des slater
      function fa(a,l1)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      fa=(a+a)**(1.5+l1)/dsqrt(fac(2+l1+l1))
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function alm(aj1,aj2,aj3,am1,am2,am3)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      parameter (api=3.141592654d00)
      parameter (ak0=0.d0)
c	api=3.141592654d00
c	xi=(0.0d00,1.0d00)
c	ak0=0.d0
      alm1=(2.d0*aj1+1.d0)*(2.d0*aj2+1.d0)*(2.d0*aj3+1.d0)/4.d0/api
      alm1=dsqrt(alm1)
      alm=alm1*a3j(aj1,aj2,aj3,ak0,ak0,ak0)*a3j(aj1,aj2,aj3,am1,am2,am3)
c	print*,'alm=',alm,88888
      return
      end



ccccccccc       calcul des symboles 3j          cccccccccccccccccccccccccccccccccc
      function a3j(aj1,aj2,aj3,am1,am2,am3) 
      implicit double precision(a-h,o-z)
      dimension am(9) 
      amm=am1+am2+am3
      if(amm.ne.0.) then
      a3j=0.
      goto 189
      else 
      goto 188
      endif

 188  a=(2*aj3+1) 
      am(1)=(aj1+aj2-aj3) 
      am(2)=(aj1-aj2+aj3) 
      am(3)=(aj2-aj1+aj3) 
      am(4)=aj1+am1 
      am(5)=aj1-am1 
      am(6)=aj2+am2 
      am(7)=aj2-am2 
      am(8)=aj3+am3 
      am(9)=aj3-am3 
      dnm=aj1+aj2+aj3+1 
      afact=1 
      do 1005 i=1,9 
      n=am(i) 
      call at(n,fact,if) 
      afact=afact*fact 
 1005  continue 
      n=dnm 
      call at(n,fact,IF) 
      ay=fact 
      ax=dsqrt(a*afact/ay) 
      asum=0.0 
      do 1006 kk=1,20 
      k=kk-1 
      am(1)=k 
      am(2)=aj1+aj2-aj3-k 
      am(3)=aj1-am1-k 
      am(4)=aj2+am2-k 
      am(5)=aj3-aj2+am1+k 
      am(6)=aj3-aj1-am2+k 
      afact=1 
      do 1007 I=1,6 
      n=am(I) 
      If(n) 1006,1008,1008 
 1008 call at(n,fact,if) 
      afact=afact*fact 
 1007 continue 
      ay=afact 
      asum=asum+(-1)**k/ay 
 1006 continue 
      clebsh=ax*asum 
      a3j=clebsh*((-1)**(aj1-aj2+am3))/dsqrt(2.d0*aj3+1.d0)
 189  return 
      end 
      
      subroutine at(n,fact,if) 
      implicit double precision(a-h,o-z)
      if=0 
      if(n) 1001,1002,1003 
 1001 if=1 
      return 
 1002 fact=1 
      return 
 1003 fact=1 
      do 1004 i=1,n 
      fact=fact*i 
 1004 continue 
      return 
      end 


      function uyl(atk,ate,fik,fie,l)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      ll=2*l+1
      uy=(0.00d0,0.00d0)
      do j=1,ll
      am=-l+j-1
      uy=uy+dconjg(ylm(l,am,ate,fie))*ylm(l,am,atk,fik)
      enddo
      uyl=uy
      return
      end


      function uy1(atk,ate,fik,fie)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      am0=0.d0
      am1=1.d0
      am11=-1.d0
    
      uy1=ylm(1,am0,ate,fie)*ylm(1,am0,atk,fik)
     $-ylm(1,am1,ate,fie)*ylm(1,am11,atk,fik)
     $-ylm(1,am11,ate,fie)*ylm(1,am1,atk,fik)
      return
      end

      function uy2(atk,ate,fik,fie)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      am0=0.d0
      am1=1.d0
      am11=-1.d0
      am2=2.d0
      am22=-2.d0
    
      uy2=ylm(2,am0,ate,fie)*ylm(2,am0,atk,fik)
     $-ylm(2,am1,ate,fie)*ylm(2,am11,atk,fik)
     $-ylm(2,am11,ate,fie)*ylm(2,am1,atk,fik)
     $+ylm(2,am2,ate,fie)*ylm(2,am22,atk,fik)
     $+ylm(2,am22,ate,fie)*ylm(2,am2,atk,fik)
      return
      end

      function uy3(atk,ate,fik,fie)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      am0=0.d0
      am1=1.d0
      am11=-1.d0
      am2=2.d0
      am22=-2.d0
      am3=3.d0
      am33=-3.d0
    
      uy3=ylm(3,am0,ate,fie)*ylm(3,am0,atk,fik)
     $-ylm(3,am1,ate,fie)*ylm(3,am11,atk,fik)
     $-ylm(3,am11,ate,fie)*ylm(3,am1,atk,fik)
     $+ylm(3,am2,ate,fie)*ylm(3,am22,atk,fik)
     $+ylm(3,am22,ate,fie)*ylm(3,am2,atk,fik)
     $-ylm(3,am3,ate,fie)*ylm(3,am33,atk,fik)
     $-ylm(3,am33,ate,fie)*ylm(3,am3,atk,fik)
      return
      end


      function uy4(atk,ate,fik,fie)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      am0=0.d0
      am1=1.d0
      am11=-1.d0
      am2=2.d0
      am22=-2.d0
      am3=3.d0
      am33=-3.d0
      am4=4.d0
      am44=-4.d0
    
      uy4=ylm(4,am0,ate,fie)*ylm(4,am0,atk,fik)
     $-ylm(4,am1,ate,fie)*ylm(4,am11,atk,fik)
     $-ylm(4,am11,ate,fie)*ylm(4,am1,atk,fik)
     $+ylm(4,am2,ate,fie)*ylm(4,am22,atk,fik)
     $+ylm(4,am22,ate,fie)*ylm(4,am2,atk,fik)
     $-ylm(4,am3,ate,fie)*ylm(4,am33,atk,fik)
     $-ylm(4,am33,ate,fie)*ylm(4,am3,atk,fik)
     $+ylm(4,am4,ate,fie)*ylm(4,am44,atk,fik)
     $+ylm(4,am44,ate,fie)*ylm(4,am4,atk,fik)
      return
      end

      function uy5(atk,ate,fik,fie)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      am0=0.d0
      am1=1.d0
      am11=-1.d0
      am2=2.d0
      am22=-2.d0
      am3=3.d0
      am33=-3.d0
      am4=4.d0
      am44=-4.d0
      am5=5.d0
      am55=-5.d0
    
      uy5=ylm(5,am0,ate,fie)*ylm(5,am0,atk,fik)
     $-ylm(5,am1,ate,fie)*ylm(5,am11,atk,fik)
     $-ylm(5,am11,ate,fie)*ylm(5,am1,atk,fik)
     $+ylm(5,am2,ate,fie)*ylm(5,am22,atk,fik)
     $+ylm(5,am22,ate,fie)*ylm(5,am2,atk,fik)
     $-ylm(5,am3,ate,fie)*ylm(5,am33,atk,fik)
     $-ylm(5,am33,ate,fie)*ylm(5,am3,atk,fik)
     $+ylm(5,am4,ate,fie)*ylm(5,am44,atk,fik)
     $+ylm(5,am44,ate,fie)*ylm(5,am4,atk,fik)
     $-ylm(5,am5,ate,fie)*ylm(5,am55,atk,fik)
     $-ylm(5,am55,ate,fie)*ylm(5,am5,atk,fik)
      return
      end

      function uy6(atk,ate,fik,fie)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      am0=0.d0
      am1=1.d0
      am11=-1.d0
      am2=2.d0
      am22=-2.d0
      am3=3.d0
      am33=-3.d0
      am4=4.d0
      am44=-4.d0
      am5=5.d0
      am55=-5.d0
      am6=6.d0
      am66=-6.d0
    
      uy6=ylm(6,am0,ate,fie)*ylm(6,am0,atk,fik)
     $-ylm(6,am1,ate,fie)*ylm(6,am11,atk,fik)
     $-ylm(6,am11,ate,fie)*ylm(6,am1,atk,fik)
     $+ylm(6,am2,ate,fie)*ylm(6,am22,atk,fik)
     $+ylm(6,am22,ate,fie)*ylm(6,am2,atk,fik)
     $-ylm(6,am3,ate,fie)*ylm(6,am33,atk,fik)
     $-ylm(6,am33,ate,fie)*ylm(6,am3,atk,fik)
     $+ylm(6,am4,ate,fie)*ylm(6,am44,atk,fik)
     $+ylm(6,am44,ate,fie)*ylm(6,am4,atk,fik)
     $-ylm(6,am5,ate,fie)*ylm(6,am55,atk,fik)
     $-ylm(6,am55,ate,fie)*ylm(6,am5,atk,fik)
     $+ylm(6,am6,ate,fie)*ylm(6,am66,atk,fik)
     $+ylm(6,am66,ate,fie)*ylm(6,am6,atk,fik)
      return
      end

      function uy7(atk,ate,fik,fie)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      am0=0.d0
      am1=1.d0
      am11=-1.d0
      am2=2.d0
      am22=-2.d0
      am3=3.d0
      am33=-3.d0
      am4=4.d0
      am44=-4.d0
      am5=5.d0
      am55=-5.d0
      am6=6.d0
      am66=-6.d0
      am7=7.d0
      am77=-7.d0
    
      uy7=ylm(7,am0,ate,fie)*ylm(7,am0,atk,fik)
     $-ylm(7,am1,ate,fie)*ylm(7,am11,atk,fik)
     $-ylm(7,am11,ate,fie)*ylm(7,am1,atk,fik)
     $+ylm(7,am2,ate,fie)*ylm(7,am22,atk,fik)
     $+ylm(7,am22,ate,fie)*ylm(7,am2,atk,fik)
     $-ylm(7,am3,ate,fie)*ylm(7,am33,atk,fik)
     $-ylm(7,am33,ate,fie)*ylm(7,am3,atk,fik)
     $+ylm(7,am4,ate,fie)*ylm(7,am44,atk,fik)
     $+ylm(7,am44,ate,fie)*ylm(7,am4,atk,fik)
     $-ylm(7,am5,ate,fie)*ylm(7,am55,atk,fik)
     $-ylm(7,am55,ate,fie)*ylm(7,am5,atk,fik)
     $+ylm(7,am6,ate,fie)*ylm(7,am66,atk,fik)
     $+ylm(7,am66,ate,fie)*ylm(7,am6,atk,fik)
     $-ylm(7,am7,ate,fie)*ylm(7,am77,atk,fik)
     $-ylm(7,am77,ate,fie)*ylm(7,am7,atk,fik)
      return
      end


      function uy8(atk,ate,fik,fie)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      am0=0.d0
      am1=1.d0
      am11=-1.d0
      am2=2.d0
      am22=-2.d0
      am3=3.d0
      am33=-3.d0
      am4=4.d0
      am44=-4.d0
      am5=5.d0
      am55=-5.d0
      am6=6.d0
      am66=-6.d0
      am7=7.d0
      am77=-7.d0
      am8=8.d0
      am88=-8.d0
    
      uy8=ylm(8,am0,ate,fie)*ylm(8,am0,atk,fik)
     $-ylm(8,am1,ate,fie)*ylm(8,am11,atk,fik)
     $-ylm(8,am11,ate,fie)*ylm(8,am1,atk,fik)
     $+ylm(8,am2,ate,fie)*ylm(8,am22,atk,fik)
     $+ylm(8,am22,ate,fie)*ylm(8,am2,atk,fik)
     $-ylm(8,am3,ate,fie)*ylm(8,am33,atk,fik)
     $-ylm(8,am33,ate,fie)*ylm(8,am3,atk,fik)
     $+ylm(8,am4,ate,fie)*ylm(8,am44,atk,fik)
     $+ylm(8,am44,ate,fie)*ylm(8,am4,atk,fik)
     $-ylm(8,am5,ate,fie)*ylm(8,am55,atk,fik)
     $-ylm(8,am55,ate,fie)*ylm(8,am5,atk,fik)
     $+ylm(8,am6,ate,fie)*ylm(8,am66,atk,fik)
     $+ylm(8,am66,ate,fie)*ylm(8,am6,atk,fik)
     $-ylm(8,am7,ate,fie)*ylm(8,am77,atk,fik)
     $-ylm(8,am77,ate,fie)*ylm(8,am7,atk,fik)
     $+ylm(8,am8,ate,fie)*ylm(8,am88,atk,fik)
     $+ylm(8,am88,ate,fie)*ylm(8,am8,atk,fik)
      return
      end








ccccc  argument l+xi*az/ake cccccc
      function xe(m,qiy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)

      xe=m+1.+qiy
      return
      end
ccc   argument l-xi*az/ak 
      function xee(m,qiy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xee=m+1.-qiy
      return
      end



ccc   argument 2l+2 de fonction hypergeometrique
      function ql(l)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      ql=l+xi*0.d0
      return
      end


ccccc  dephasage coulombien appelé cs cccccc



      function cs(xe)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      cs=imag(uc(xe))/real(uc(xe))
      cs=datan(cs)
      return
      end


ccccc  norme de la fonction gamma  ccc
      function agl(xe)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      bgl=uc(xe)*dconjg(uc(xe))
      agl=dsqrt(bgl)
      return
      end
cccccccccccccccccccccccccccccccccccccc


c      function delt0(a1s,ak,d)
      function delt0(a1s,ak)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
c      dimension ax(10),dx(10),dw(10)
      dimension ax(10),dxl(10),dwl(10)
cccccccccccccccccc points et poids de gauss-Legendre ccccccccccccccccccccc
c      data dx /-0.973906528517,-0.865063366689,-0.679409568299,
c     $-0.433395394129,-0.148874338982,0.973906528517,0.865063366689,
c     $0.679409568299,0.433395394129,0.148874338982/
     
c      data dw /0.0666713443087,0.149451349151,0.219086362516,
c     $0.26926671931,0.295524224715,0.0666713443087,0.149451349151,
c     $0.219086362516,0.26926671931,0.295524224715/
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccc points et poids de gauss-Laguerre cccccccc
      data dxl /0.13779347054,0.729454549503,1.80834290174,
     $3.40143369785,5.55249614006,8.33015274676,11.8437858379,
     $16.2792578314,21.996585812,29.9206970123/
     
      data dwl /0.354009738607,0.831902301044,1.33028856175,
     $1.86306390311,2.45025555808,3.12276415514,3.9341526956,
     $4.99241487226,6.57220248513,9.78469584034/
      asom=0.d0

      do i=1,10
      a=0.001
c	ax(i)=(dx(i)*(d-a)+d+a)/2.d0
      ax(i)=(dxl(i)-a)
      akx=ax(i)*ak
      aj0=dsin(akx)/akx 
      avstat=-1.-(1.d0+a1s*ax(i))*dexp(-2.*a1s*ax(i))

      ft=aj0**2*avstat*ax(i)
c	aint=ft*dw(i)*(d-a)/2.
      aint=ft*dwl(i)
      asom=asom+ft
      enddo
      delt0=-2.*asom*ak
      return
      end


c      function delt1(a1s,ak,d)
      function delt1(a1s,ak)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
c      dimension ax(10),dx(10),dw(10)
      dimension ax(10),dxl(10),dwl(10)
cccccccccccccccccc points et poids de gauss-Legendre ccccccccccccccccccccc
c      data dx /-0.973906528517,-0.865063366689,-0.679409568299,
c     $-0.433395394129,-0.148874338982,0.973906528517,0.865063366689,
c     $0.679409568299,0.433395394129,0.148874338982/
     
c      data dw /0.0666713443087,0.149451349151,0.219086362516,
c     $0.26926671931,0.295524224715,0.0666713443087,0.149451349151,
c     $0.219086362516,0.26926671931,0.295524224715/
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccc points et poids de gauss-Laguerre cccccccc
      data dxl /0.13779347054,0.729454549503,1.80834290174,
     $3.40143369785,5.55249614006,8.33015274676,11.8437858379,
     $16.2792578314,21.996585812,29.9206970123/
     
      data dwl /0.354009738607,0.831902301044,1.33028856175,
     $1.86306390311,2.45025555808,3.12276415514,3.9341526956,
     $4.99241487226,6.57220248513,9.78469584034/
cccccccccccccccccccccccccccccccccccccccccccccccccccc
	
     	asom=0.d0
      do i=1,10
      a=0.001
c	ax(i)=(dx(i)*(d-a)+d+a)/2.d0
      ax(i)=(dxl(i)-a)
      akx=ax(i)*ak
      aj1=dsin(akx)/akx/akx-dcos(akx)/akx 
      avstat=-1.-(1.d0+a1s*ax(i))*dexp(-2.*a1s*ax(i))

      ft=aj1**2*avstat*ax(i)
c	aint=ft*dw(i)*(d-a)/2.
      aint=ft*dwl(i)
      asom=asom+ft
      enddo
      delt1=-2.*asom*ak
      return
      end


c      function delt2(a1s,ak,d)
      function delt2(a1s,ak)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
c      dimension ax(10),dx(10),dw(10)
      dimension ax(10),dxl(10),dwl(10)
cccccccccccccccccc points et poids de gauss-Legendre ccccccccccccccccccccc
c      data dx /-0.973906528517,-0.865063366689,-0.679409568299,
c     $-0.433395394129,-0.148874338982,0.973906528517,0.865063366689,
c     $0.679409568299,0.433395394129,0.148874338982/
     
c      data dw /0.0666713443087,0.149451349151,0.219086362516,
c     $0.26926671931,0.295524224715,0.0666713443087,0.149451349151,
c     $0.219086362516,0.26926671931,0.295524224715/
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccc points et poids de gauss-Laguerre cccccccc
      data dxl /0.13779347054,0.729454549503,1.80834290174,
     $3.40143369785,5.55249614006,8.33015274676,11.8437858379,
     $16.2792578314,21.996585812,29.9206970123/
     
      data dwl /0.354009738607,0.831902301044,1.33028856175,
     $1.86306390311,2.45025555808,3.12276415514,3.9341526956,
     $4.99241487226,6.57220248513,9.78469584034/
cccccccccccccccccccccccccccccccccccccccccccccccccccc
	
     	asom=0.d0
      do i=1,10
      a=0.001
c	ax(i)=(dx(i)*(d-a)+d+a)/2.d0
      ax(i)=(dxl(i)-a)
      akx=ax(i)*ak
      aj2=3.*dsin(akx)/(akx**3)-dsin(akx)/akx-3.*dcos(akx)/(akx**2) 
      avstat=-1.-(1.d0+a1s*ax(i))*dexp(-2.*a1s*ax(i))
      ft=aj2**2*avstat*ax(i)
c	aint=ft*dw(i)*(d-a)/2.
      aint=ft*dwl(i)
      asom=asom+ft
      enddo
      delt2=-2.*asom*ak

      return
      end



cccc
      function ay1(cz)
c      function ay1(cz,ak)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      ay1=dsqrt(3.d0/(4.d0*api))*cz
      return
      end

      function ry11(cx,cy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry11=-dsqrt(3.d0/(8.d0*api))*rdxy
      return
      end

      function ry111(cx,cy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry11=-dsqrt(3.d0/(8.d0*api))*rdxy
      ry111=-dconjg(ry11)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccc
      function ay2(cz,ak2)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      ay2=dsqrt(5./(16.d0*api))*(3.*cz**2-ak2**2)
      return
      end

      function ry21(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry21=-3.*dsqrt(5.d0/api/6.d0)*rdxy*cz/2.
      return
      end

      function ry211(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry21=-3.*dsqrt(5.d0/api/6.d0)*rdxy*cz/2.
      ry211=-dconjg(ry21)
      return
      end

      function ry22(cx,cy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
c      ak=dsqrt(cz**2+cx**2+cy**2)
      ry22=dsqrt(15.d0/api/32.d0)*rdxy**2
      return
      end

      function ry222(cx,cy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry22=dsqrt(15.d0/api/32.d0)*rdxy**2
      ry222=dconjg(ry22)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function ay3(cz,ak)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      ay3=-dsqrt(7.d0/(16.d0*api))*(3.*cz*ak**2-5.*cz**3)
      return
      end			 
	
     
      function ry31(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry31=3.*dsqrt(7.d0/api/3.d0)*rdxy*(ak**2-5.*cz**2)/8.
      return
      end

      function ry311(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry31=3.*dsqrt(7.d0/api/3.d0)*rdxy*(ak**2-5.*cz**2)/8.
      ry311=-dconjg(ry31)
      return
      end

      function ry32(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry32=dsqrt(105.d0/api/2.d0)*rdxy**2*cz/4.
      return
      end

      function ry322(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry32=dsqrt(105.d0/api/2.d0)*rdxy**2*cz/4.
      ry322=dconjg(ry32)
      return
      end

      function ry33(cx,cy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry33=-dsqrt(35.d0/api/64.d0)*rdxy**3
      return
      end

      function ry333(cx,cy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)

      api=3.141592654d00
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry33=-dsqrt(35.d0/api/64.d0)*rdxy**3
      ry333=-dconjg(ry33)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function ay4(ab,ak)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      ay4=(3./16.)*dsqrt(1.d0/api)*(3.*ak**4-30.*(ak*ab)**2
     $+35.*ab**4)
      return
      end			 

      function ry41(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry41=-(15./8.)*dsqrt(1.d0/5./api)*(7.*cz**3-3.*cz*ak**2)*rdxy
      return
      end

      function ry411(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry41=-(15./8.)*dsqrt(1.d0/5./api)*(7.*cz**3-3.*cz*ak**2)*rdxy
      ry411=-dconjg(ry41)
      return
      end

      function ry42(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry42=(15./8.)*dsqrt(1.d0/10./api)*(7.*cz**2-ak**2)*rdxy**2
      return
      end

      function ry422(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry42=(15./8.)*dsqrt(1.d0/10./api)*(7.*cz**2-ak**2)*rdxy**2
      ry422=dconjg(ry42)
      return
      end

      function ry43(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry43=-(105./8.)*dsqrt(1.d0/35./api)*cz*rdxy**3
      return
      end

      function ry433(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry43=-(105./8.)*dsqrt(1.d0/35./api)*cz*rdxy**3
      ry433=-dconjg(ry43)
      return
      end


      function ry44(cx,cy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)
      api=3.141592654d00
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry44=(105./16.)*dsqrt(1.d0/70./api)*rdxy**4
      return
      end

      function ry444(cx,cy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      xi=(0.0d00,1.0d00)

      api=3.141592654d00
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry44=(105./16.)*dsqrt(1.d0/70./api)*rdxy**4
      ry444=dconjg(ry44)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function ay5(cz,ak)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      ay5=-(1./16.)*dsqrt(11.d0/api)*(-63.*cz**5+70.*cz**3*ak**2
     $-15.*cz*ak**4)
      return
      end			 


      function ry51(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry51=dsqrt(165.d0/8.d0/api)*rdxy*(14.*cz**2*ak**2-21.*cz**4
     $-ak**4)/8.
      return
      end

      function ry511(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry51=dsqrt(165.d0/8.d0/api)*rdxy*(14.*cz**2*ak**2-21.*cz**4
     $-ak**4)/8.
      ry511=-dconjg(ry51)
      return
      end

      function ry52(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry52=dsqrt(1155.d0/api/2.d0)*rdxy**2*cz*(-ak**2+3.*cz**2)/8.
      return
      end

      function ry522(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry52=dsqrt(1155.d0/api/2.d0)*rdxy**2*cz*(-ak**2+3.*cz**2)/8.
      ry522=dconjg(ry52)
      return
      end

      function ry53(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry53=dsqrt(385.d0/api)*(ak**2-9.*cz**2)*rdxy**3/32.d0
      return
      end


      function ry533(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry53=dsqrt(385.d0/api)*(ak**2-9.*cz**2)*rdxy**3/32.d0
      ry533=-dconjg(ry53)
      return
      end

      function ry54(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry54=3.*dsqrt(385.d0/api/2.d0)*cz*rdxy**4/16.
      return
      end

      function ry544(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry54=3.*dsqrt(385.d0/api/2.d0)*cz*rdxy**4/16.
      ry544=dconjg(ry54)
      return
      end


      function ry55(cx,cy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry55=-3.*dsqrt(77.d0/api)*rdxy**5/32.
      return
      end

      function ry555(cx,cy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry55=-3.*dsqrt(77.d0/api)*rdxy**5/32.
      ry555=-dconjg(ry55)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function ay6(cz,ak)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      ay6=(1./32.)*dsqrt(13.d0/api)*(231.*cz**6-315.*cz**4*ak**2
     $+105.*cz**2*ak**4-ak**6)
      return
      end			 

      function ry61(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry61=-dsqrt(273.d0/2.d0/api)*rdxy*(33.*cz**5-30.*cz**3*ak**2
     $+5.*cz*ak**4)/16.
      return
      end

      function ry611(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry61=-dsqrt(273.d0/2.d0/api)*rdxy*(33.*cz**5-30.*cz**3*ak**2
     $+5.*cz*ak**4)/16.
      ry611=-dconjg(ry61)
      return
      end

      function ry62(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry62=dsqrt(1365.d0/api)*rdxy**2*(33.*cz**4-18.*cz**2*ak**2
     $+ak**4)/64.
      return
      end

      function ry622(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry62=dsqrt(1365.d0/api)*rdxy**2*(33.*cz**4-18.*cz**2*ak**2
     $+ak**4)/64.
      ry622=dconjg(ry62)
      return
      end


      function ry63(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry63=-dsqrt(1365.d0/api)*rdxy**3*(11.*cz**3-3.*cz*ak**2)/32.
      return
      end

      function ry633(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry63=-dsqrt(1365.d0/api)*rdxy**3*(11.*cz**3-3.*cz*ak**2)/32.
      ry633=-dconjg(ry63)
      return
      end

      function ry64(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry64=3.*dsqrt(91.d0/2./api)*rdxy**4*(11.*cz**2-ak**2)/32.
      return
      end

      function ry644(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry64=3.*dsqrt(91.d0/2./api)*rdxy**4*(11.*cz**2-ak**2)/32.
      ry644=dconjg(ry64)
      return
      end


      function ry65(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry65=-3.*dsqrt(1001.d0/api)*rdxy**5*ak/32.
      return
      end

      function ry655(cx,cy,cz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
      ak=dsqrt(cz**2+cx**2+cy**2)
      ry65=-3.*dsqrt(1001.d0/api)*rdxy**5*ak/32.
      ry655=-dconjg(ry65)
      return
      end

      function ry66(cx,cy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry66=dsqrt(3003.d0/api)*rdxy**6/63.
      return
      end

      function ry666(cx,cy)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      rdxy=cx+xi*cy
c	ak=dsqrt(cz**2+cx**2+cy**2)
      ry66=dsqrt(3003.d0/api)*rdxy**6/63.
      ry666=dconjg(ry66)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function v1f1(x1,x2,x3)
c fonction log avec complexes
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      v1=x1*x3/x2+x1*(x1+1.)*x3**2/x2/(x2+1.)/2.
     $+x1*(x1+1.)*(x1+2.)*x3**3/x2/(x2+1.)/(x2+2.)/6.
     $+x1*(x1+1.)*(x1+2.)*(x1+3.)*x3**4/x2/(x2+1.)/(x2+2.)/(x2+3.)/24.

      v1f1=1.+v1
      return
      end


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

      function uf(za,zc,z)
c fonction 1f1  pour z inf a 1
c handbook p.504
c et pour z grand
       implicit double precision (a-h)
       implicit complex*16 (o-z)
       dimension u(1000)
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


      function usr(za,zc,z)
c fonction 1f1 incomplete pour z grand
c handbook p.508 13.5.1
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      dimension u(700),v(700)


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

      if (j.ge.90) then
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




       function algndr(l,m,ax)
c fonction polynomes de legendre associés aux harmoniques sphériques  
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
c fonction y(l,m) harmoniques sphériques 
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
      function q5(x,y,z)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      q5=15.*x*x+20.*z*z*(y-3.)*x+4.*z**4.*(y-3.)*(y-4.)
      q5=q5*8.*z*y*(y-1.)*(y-2.)*x**(y-5.)
      end
      function q6(x,y,z)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      q6=15.*x*x*x+90.*z*z*x*x*(y-3.)+60.*x*z**4.*(y-3.)*(y-4.)
      q6=q6+8.*z**6.*(y-3.)*(y-4.)*(y-5.)
      q6=q6*8.*y*(y-1.)*(y-2.)*x**(y-6.)
      end
      function q7(x,y,z)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      q7=105.*x*x*x+210.*z*z*(y-4.)*x*x+84.*z**4.*(y-4.)*(y-5.)*x
      q7=q7+8.*z**6.*(y-4.)*(y-5.)*(y-6.)
      q7=q7*16.*y*(y-1.)*(y-2.)*(y-3.)*z*x**(y-7.)
      end
      function q8(x,y,z)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      q8=105.*(x**4.+8.*z*z*(y-4.)*x*x*x+8.*z**4.*(y-4.)*(y-5.)*x*x)
      q8=q8+224.*z**6.*(y-4)*(y-5.)*(y-6.)*x
      q8=q8+16.*z**8.*(y-4.)*(y-5.)*(y-6.)*(y-7.)
      q8=q8*16.*y*(y-1.)*(y-2.)*(y-3.)*x**(y-8.)
      end
      function q9(x,y,z)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      q9=x*x*(945.*x*x+2520.*z*z*(y-5.)*x+1512.*z**4.*(y-5.)*(y-6.))
      q9=q9+288.*z**6.*(y-5.)*(y-6.)*(y-7.)*x
      q9=q9+16.*z**8.*(y-5.)*(y-6.)*(y-7.)*(y-8.)
      q9=q9*32.*y*(y-1.)*(y-2.)*(y-3.)*(y-4.)*z*x**(y-9.)
      end
      function q10(x,y,z)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      q10=945.*x**5.+9450.*z*z*(y-5.)*x**4.
      q10=q10+12600.*z**4.*(y-5.)*(y-6.)*x**3.
      q10=q10+5040.*z**6.*(y-5.)*(y-6.)*(y-7.)*x*x
      q10=q10+720.*z**8.*(y-5.)*(y-6.)*(y-7.)*(y-8.)*x
      q10=q10+32.*z**10.*(y-5.)*(y-6.)*(y-7.)*(y-8.)*(y-9.)
      q10=q10*32.*y*(y-1.)*(y-2.)*(y-3.)*(y-4.)*x**(y-10.)
      end
      function q11(x,y,z) 
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      u11=10395.*x**5.+34650.*z*z*(y-6.)*x**4.
      u12=27720.*z**4.*(y-6.)*(y-7.)*x**3.
      u13=7920.*z**6.*(y-6.)*(y-7.)*(y-8.)*x*x
      u14=880.*z**8.*(y-6.)*(y-7.)*(y-8.)*(y-9.)*x
      u15=32.*z**10.*(y-6.)*(y-7.)*(y-8.)*(y-9.)*(y-10.)
      q11=(u11+u12+u13+u14+u15)*64.*z*y*(y-1.)*(y-2.)
     1     *(y-3.)*(y-4.)*(y-5.)*x**(y-11.)
      end
      function q12(x,y,z)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      u11=10395.*x**6.+124740.*z*z*(y-6.)*x**5.
      u12=207900.*z**4.*(y-6.)*(y-7.)*x**4.
      u13=110880.*z**6.*(y-6.)*(y-7.)*(y-8.)*x**3.
      u14=23760.*z**8.*(y-6.)*(y-7.)*(y-8.)*(y-9.)*x*x
      u15=2112.*z**10.*(y-6.)*(y-7.)*(y-8.)*(y-9.)*(y-10.)*x
      u16=64.*z**12.*(y-6.)*(y-7.)*(y-8.)*(y-9.)*(y-10.)*(y-11.)
      q12=(u11+u12+u13+u14+u15+u16)*64.*y*(y-1.)*(y-2.)
     1     *(y-3.)*(y-4.)*(y-5.)*x**(y-12.)
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

      function td0(x,y,z1,z2,z3)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      td0=x*x*z1-x*y*z2+y*y*z3
      end
      function td1(x,y,z1,z2,z3,z4,z5)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      td1=4.*x*z1*z4-2.*y*z2*z4+x*z2*z5-2.*y*z3*z5
      end
      function td2(x,y,z1,z2,z3,z4,z5)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      td2=2.*(z1*(2.*x+4.*z4*z4)+z2*(2.*z4*z5-y)+z5*z5*z3)
      end
      function td3(z1,z2,z4,z5)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      td3=6.*(4.*z4*z1+z5*z2)
      end
      function td4(z1)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      td4=24.*z1
      end

      function tf0(x,y,z1,z2,z3,z4)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tf0=x*x*x*z1-x*x*y*z2+x*y*y*z3-y*y*y*z4
      end
      function tf1(x,y,z1,z2,z3,z4,z5,z6)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tf1=6.*x*x*z1*z5-z2*(4.*x*y*z5-z6*x*x)+z3*2.*(y*y*z5-y*x*z6)
      tf1=tf1+z4*3.*y*y*z6
      end
      function tf2(x,y,z1,z2,z3,z4,z5,z6)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tf2=6.*z1*(x*x+4.*x*z5*z5)-z2*4.*(2.*z5*z5*y+x*y-2.*z5*z6*x)
      tf2=tf2+2.*z3*(y*y-4.*z5*z6*y+z6*z6*x)-6.*z4*y*z6*z6
      end
      function tf3(x,y,z1,z2,z3,z4,z5,z6)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tf3=24.*z1*z5*(3.*x+2.*z5*z5)-12.*z2*(2.*y*z5-2.*z5*z5*z6-z6*x)
      tf3=tf3+12.*z3*z6*(-y+z6*z5)+6.*z4*z6*z6*z6
      end
      function tf4(x,y,z1,z2,z3,z5,z6)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tf4=72.*z1*(x+4.*z5*z5)-24.*z2*(y-4.*z5*z6)
      tf4=tf4+24.*z3*z6*z6
      end
      function tf5(z1,z2,z5,z6)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tf5=720.*z5*z1+120.*z6*z2
      end
      function tf6(z1)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tf6=720.*z1
      end


      function tg0(x,y,z1,z2,z3,z4,z5)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tg0=x*x*x*x*z1-x*x*x*y*z2+x*x*y*y*z3-x*y*y*y*z4+y*y*y*y*z5
      end
      function tg1(x,y,z1,z2,z3,z4,z5,z6,z7)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tg1=8.*x*x*x*z1*z6-z2*x*x*(6.*z6*y-z7*x)
      tg1=tg1+z3*2.*x*y*(2.*z6*y-z7*x)
      tg1=tg1-z4*y*y*(2.*z6*y-3.*z7*x)-4.*z7*y*y*y*z5
      end
      function tg2(x,y,z1,z2,z3,z4,z5,z6,z7)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tg2=8.*z1*x*x*(x+6.*z6*z6)-6.*x*z2*(x*y+4.*z6*z6*y-2.*z7*z6*x)
      tg2=tg2+2.*z3*(2.*x*y*y+4.*z6*z6*y*y-8.*z7*z6*x*y+z7*z7*x*x)
      tg2=tg2-2.*z4*y*(y*y-6.*z7*z6*y+3.*z7*z7*x)+12.*z7*z7*y*y*z5
      end
      function tg3(x,y,z1,z2,z3,z4,z5,z6,z7)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tg30=48.*z1*x*z6*(3.*x+4.*z6*z6)
      tg31=-6.*z2*(12.*x*y*z6-3.*x*x*z7+8.*z6**3.*y-12.*z7*z6*z6*x)
      tg32=24.*z3*(z6*y*y-z7*x*y-2.*z7*z6*z6*y+z7*z7*z6*x)
      tg33=-6.*z4*z7*(-3.*y*y+6.*z7*z6*y-z7*z7*x)
      tg3=tg30+tg31+tg32+tg33-24.*z7**3.*y*z5
      end
      function tg4(x,y,z1,z2,z3,z4,z5,z6,z7)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tg40=z1*(144.*x*x+1152.*x*z6*z6+384.*z6**4.)
      tg41=-z2*(72.*x*y+288.*z6*z6*y-288.*z7*z6*x-192.*z7*z6*z6*z6)
      tg42=z3*(24.*y*y-192.*z7*z6*y+48.*z7*z7*x+96.*z7*z7*z6*z6)
      tg43=-z4*(72.*z7*z7*y-48.*z6*z7*z7*z7)
      tg4=tg40+tg41+tg42+tg43+24.*z7**4.*z5
      end
      function tg5(x,y,z1,z2,z3,z4,z6,z7)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tg50=z1*(2880.*x*z6+3840.*z6**3.)
      tg51=-z2*(720.*y*z6-360.*z7*x-1440.*z6*z6*z7)
      tg52=z3*(-240.*y*z7+480.*z6*z7*z7)
      tg53=z4*120.*z7**3.
      tg5=tg50+tg51+tg52+tg53
      end
      function tg6(x,y,z1,z2,z3,z6,z7)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tg60=z1*(2880.*x+17280.*z6*z6)
      tg61=-z2*(720.*y-4320.*z7*z6)
      tg62=z3*720.*z7*z7
      tg6=tg60+tg61+tg62
      end
      function tg7(z1,z2,z6,z7)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tg70=z1*40320.*z6
      tg71=z2*5040.*z7
      tg7=tg70+tg71
      end
      function tg8(z1)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tg8=z1*40320.
      end
  
      function th0(x,y,z1,z2,z3,z4,z5,z6)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      th0=x*x*x*x*x*z1-x*x*x*x*y*z2+x*x*x*y*y*z3-x*x*y*y*y*z4
     1 +x*y*y*y*y*z5-y*y*y*y*y*z6
      end
      function th1(x,y,z1,z2,z3,z4,z5,z6,z7,z8)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      th1=10.*x*x*x*x*z1*z7-z2*x*x*x*(8.*z7*y-z8*x)
      th1=th1+z3*x*x*y*(6.*z7*y-2.*z8*x)-z4*x*y*y*(4.*z7*y-3.*z8*x)
     1 +z5*y*y*y*(2.*z7*y-4.*z8*x)+5.*z6*z8*y*y*y*y
      end
      function th2(x,y,z1,z2,z3,z4,z5,z6,z7,z8)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      th2=x*x*x*z1*(10.*x+80.*z7*z7)
     1  -z2*x*x*(8.*x*y+48.*z7*z7*y-16.*z7*z8*x)
      th2=th2+z3*x*(6.*x*y*y+24.*z7*z7*y*y-24.*z7*z8*x*y+2.*z8*z8*x*x)
     1 -z4*y*(4.*x*y*y+8.*z7*z7*y*y-24.*z7*z8*x*y+6.*z8*z8*x*x)
     1 +z5*y*y*(2.*y*y-16.*z7*z8*y+12.*z8*z8*x)-z6*20.*z8*z8*y*y*y
      end
      function th3(x,y,z1,z2,z3,z4,z5,z6,z7,z8)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      th3=240.*x*x*z1*z7*(x+2.*z7*z7)
     1  -z2*x*(144.*z7*x*y-24.*z8*x*x+192.*z7*z7*z7*y-144.*z8*z7*z7*x)
     1 +z3*(72.*z7*x*y*y-36.*z8*x*x*y+48.*z7*z7*z7*y*y
     1              -144.*z8*z7*z7*x*y+36.*z7*z8*z8*x*x)
     1 -z4*(24.*z7*y*y*y-36.*z8*x*y*y-72.*z7*z7*z8*y*y
     1               +72.*z8*z8*z7*x*y-6.*z8*z8*z8*x*x)
     1 +z5*y*(-24.*z8*y*y+72.*z7*z8*z8*y-24.*z8*z8*z8*x)
     1  +z6*60.*z8*z8*z8*y*y
      end
      function th4(x,y,z1,z2,z3,z4,z5,z6,z7,z8)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      th4=240.*x*z1*(x*x+12.*z7*z7*x+8.*z7*z7*z7*z7)
     1  -z2*(144.*x*x*y+1152.*x*y*z7*z7-576.*z7*z8*x*x
     1          -768.*z8*z7*z7*z7*x+384.*z7*z7*z7*z7*y)
     1 +z3*(72.*x*y*y+288.*z7*z7*y*y-576.*z7*z8*x*y+72.*z8*z8*x*x
     1            -384.*z8*z7*z7*z7*y+288.*z8*z8*z7*z7*x)
     1 -z4*(24.*y*y*y-288.*z7*z8*y*y+144.*z8*z8*x*y+288.*z7*z7*z8*z8*y
     1               -96.*z8*z8*z8*z7*x)
     1 +z5*(144.*z8*z8*y*y-192.*z7*z8*z8*z8*y+24.*z8*z8*z8*z8*x)
     1  -z6*120.*z8*z8*z8*z8*y
      end
      function th5(x,y,z1,z2,z3,z4,z5,z6,z7,z8)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      th5=z1*(7200.*z7*x*x+19200.*z7*z7*z7*x+3840.*z7*z7*z7*z7*z7)
     1  -z2*(2880.*z7*x*y-720.*z8*x*x+3840.*z7*z7*z7*y
     1       -5760.*z8*z7*z7*x-1920.*z7*z7*z7*z7*z8)
     1 +z3*(720.*z7*y*y-720.*z8*x*y-2880.*z8*z7*z7*y
     1        +1440.*z8*z8*z7*x+960.*z8*z8*z7*z7*z7)
     1 -z4*(-360.*z8*y*y+1440.*z7*z8*z8*y-240.*z8*z8*z8*x
     1       -480.*z7*z7*z8*z8*z8)
     1 +z5*240.*z8*z8*z8*(-2.*y+z7*z8)+z6*120.*z8*z8*z8*z8*z8
      end
      function th6(x,y,z1,z2,z3,z4,z5,z7,z8)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      th6=z1*(7200.*x*x+86400.*z7*z7*x+57600.*z7*z7*z7*z7)
     1  -z2*(2880.*x*y+17280.*z7*z7*y-17280.*z8*z7*x
     1       -23040.*z8*z7*z7*z7)
     1 +z3*(720.*y*y-8640.*z8*z7*y+2160.*z8*z8*x+8640.*z8*z8*z7*z7)
     1 -z4*(2160.*z8*z8*y-2880.*z7*z8*z8*z8)
     1 +z5*720.*z8*z8*z8*z8
      end
      function th7(x,y,z1,z2,z3,z4,z7,z8)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      th7=z1*(201600.*z7*x+403200.*z7*z7*z7)
     1  -z2*(40320.*z7*y-20160.*z8*x-120960.*z8*z7*z7)
     1 +z3*(-10080.*z8*y+30240.*z8*z8*z7)
     1 +z4*5040.*z8*z8*z8
      end


      function ti0(x,y,z1,z2,z3,z4,z5,z6,z7)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      ti0=x**6.*z1-x**5.*y*z2+x**4.*y*y*z3-(x*y)**3.*z4+x*x*y**4.*z5
      ti0=ti0-x*y**5.*z6+y**6.*z7
      end
      function ti1(x,y,z1,z2,z3,z4,z5,z6,z7,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      ti1=12.*x**5.*z1*ze-z2*x**4.*(10.*ze*y-zb*x)
      ti1=ti1+z3*2.*x**3.*y*(4.*ze*y-zb*x)
      ti1=ti1-3.*z4*x*x*y*y*(2.*ze*y-zb*x)
      ti1=ti1+4.*x*y*y*y*z5*(y*ze-zb*x)
      ti1=ti1-z6*y**4.*(2.*ze*y-5.*zb*x)-6.*z7*zb*y**5.
      end
      function ti2(x,y,z1,z2,z3,z4,z5,z6,z7,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      ti2=12.*z1*x**4.*(x+10.*ze*ze)
     1 -10.*z2*x**3.*(x*y+8.*ze*ze*y-2.*ze*zb*x)
     1 +2.*z3*x*x*(4.*x*y*y+24.*(ze*y)**2.-16.*ze*zb*x*y+x*x*zb*zb)
     1 -6.*z4*x*y*(x*y*y+4.*ze*ze*y*y-6.*ze*zb*x*y+zb*zb*x*x)
     1 +4.*z5*y*y*(x*y*y+2.*ze*ze*y*y-8.*ze*zb*x*y+3.*(zb*x)**2.)
     1 -2.*z6*y**3.*(y*y-10.*ze*zb*y+10.*zb*zb*x)+30.*z7*zb*zb*y**4.
       end
      function ti3(x,y,z1,z2,z3,z4,z5,z6,z7,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      ti3=120.*z1*x**3.*ze*(3.*x+8.*ze*ze)
     1 -10.*z2*x*x*(24.*ze*x*y-3.*zb*x*x+48.*ze**3.*y-24.*ze*ze*zb*x)
     1 +48.*z3*x*(3.*ze*x*y*y-zb*x*x*y+4.*ze**3.*y*y-6.*ze*ze*zb*x*y
     1            +ze*(zb*x)**2.)
     1 -6.*z4*(12.*ze*x*y**3.+8.*(ze*y)**3.-9.*zb*(x*y)**2.
     1         -36.*ze*ze*zb*x*y*y+18.*ze*zb*zb*x*x*y-(zb*x)**3.)
     1 +24.*z5*y*(ze*y**3.-2.*zb*x*y*y-4.*ze*ze*zb*y*y
     1 +6.*ze*zb*zb*x*y-zb**3.*x*x)
     1 -10.*z6*y**2.*(-3.*zb*y*y+12.*ze*zb*zb*y-6.*zb**3.*x)
     1 -120.*z7*(zb*y)**3.
       end
      function ti4(x,y,z1,z2,z3,z4,z5,z6,z7,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      ti4=360.*z1*x*x*(x*x+16.*ze*ze*x+16.*ze**4.)
     1 -240.*z2*x*(x*x*y+12.*ze*ze*x*y-4.*ze*zb*x*x+8.*ze**4.*y
     1         -8.*ze**3.*zb*x)
     1 +48.*z3*(3.*(x*y)**2.+24.*(ze*y)**2.*x-24.*ze*zb*x*x*y
     1 +2.*zb*zb*x**3.+8.*ze**4.*y*y-32.*ze**3.*zb*x*y
     1 +12.*(ze*zb*x)**2.)
     1 -72.*z4*(y**3.*(x+4.*ze*ze)-12.*ze*zb*x*y*y-8.*ze**3.*zb*y*y
     1     +3.*(zb*x)**2.*y+12.*(ze*zb)**2.*x*y-2.*ze*zb**3.*x*x)
     1 +24.*z5*(y**3.*(y-16.*ze*zb)+12.*(zb*y)**2.*x+24.*(ze*zb*y)**2.
     1 -16.*ze*zb**3.*x*y+zb**4.*x*x)
     1 -240.*z6*y*zb*zb*(y*y-2.*ze*zb*y+x*.5*zb*zb)+360.*z7*zb**4.*y*y
       end
      function ti5(x,y,z1,z2,z3,z4,z5,z6,z7,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      ti5=2880.*z1*x*ze*(5.*x*x+20.*ze*ze*x+8.*ze**4.)
     1 -240.*z2*(30.*ze*x*x*y-5.*zb*x**3.+80.*ze**3.*x*y
     1 -60.*ze*ze*zb*x*x+16.*ze**5.*y-40.*ze**4.*zb*x)
     1 +480.*z3*(6.*ze*x*y*y-3.*zb*x*x*y+8.*ze**3.*y*y
     1     -24.*ze*ze*zb*x*y+6.*ze*zb*zb*x*x-8.*ze**4.*zb*y
     1          +8.*ze**3.*zb*zb*x)
     1 -360.*z4*(2.*ze*y**3.-3.*zb*x*y*y-12.*ze*ze*zb*y*y
     1      +12.*ze*zb*zb*x*y-zb**3.*x*x-4.*ze*ze*zb**3.*x
     1       +8.*ze**3.*zb*zb*y)
     1 +240.*z5*(-2.*zb*y**3.+12.*ze*zb*zb*y*y-8.*ze*ze*zb**3.*y
     1        -4.*zb**3.*x*y+2.*ze*zb**4.*x)
     1 -120.*z6*zb*(-10.*zb*zb*y*y+10.*ze*zb**3.*y-zb**4.*x)
     1 -720.*z7*zb**5.*y
       end
      function ti6(x,y,z1,z2,z3,z4,z5,z6,z7,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      ti6=2880.*z1*(5.*x**3.+90.*(ze*x)**2.+120.*x*ze**4.+16.*ze**6.)
     1 -1440.*z2*(5.*x*x*y+60.*ze*ze*x*y-30.*ze*zb*x**2.+40.*ze**4.*y
     1            -80.*ze**3.*zb*x-16.*ze**5.*zb)
     1 +2880.*z3*(x*y*y+6.*(ze*y)**2.-12.*ze*zb*x*y+1.5*(zb*x)**2.
     1-16.*ze**3.*zb*y+12.*x*(ze*zb)**2.+4.*ze**4.*zb*zb)
     1 -720.*z4*(y**3.-18.*ze*zb*y*y+9.*zb*zb*x*y+36.*y*(ze*zb)**2.
     1  -12.*ze*x*zb**3.-8.*(ze*zb)**3.)
     1 +5760.*z5*(.75*(zb*y)**2.-2.*ze*zb**3.*y+.5*ze*ze*zb**4.
     1     +.25*zb**4.*x)
     1 -720.*z6*zb**4.*(5.*y-2.*ze*zb)+720.*z7*zb**6.
       end


      function tj0(x,y,z1,z2,z3,z4,z5,z6,z7,z8)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tj0=x**7.*z1-x**6.*y*z2+x**5.*y*y*z3-x**4.*y**3.*z4
     1 +x**3.*y**4.*z5-x*x*y**5.*z6+x*y**6.*z7-y**7.*z8
      end
      function tj1(x,y,z1,z2,z3,z4,z5,z6,z7,z8,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tj1=14.*x**6.*z1*ze-z2*x**5.*(12.*ze*y-zb*x)
     1 +z3*2.*x**4.*y*(5.*ze*y-zb*x)-z4*x**3.*y*y*(8.*ze*y-3.*zb*x)
     1 +x*x*y*y*y*z5*(6.*y*ze-4.*zb*x)-z6*x*y**4.*(4.*ze*y-5.*zb*x)
     1 +z7*y**5.*(2.*ze*y-6.*zb*x)+z8*7.*zb*y**6.
      end
      function tj2(x,y,z1,z2,z3,z4,z5,z6,z7,z8,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tj2=14.*x**5.*z1*(x+12.*ze*ze)
     1  -12.*z2*x**4.*(x*y+10.*ze*ze*y-2.*zb*ze*x)
     1 +z3*x**3.*(10.*x*y*y+80.*ze*ze*y*y-40.*ze*zb*x*y+2.*zb*zb*x*x)
     1 -z4*x**2.*y*(8.*x*y*y+48.*ze*ze*y*y-48.*ze*zb*x*y+6.*zb*zb*x*x)
     1 +x*y*y*z5*(6.*x*y*y+24.*ze*ze*y*y-48.*ze*zb*x*y+12.*zb*zb*x*x)
     1 -z6*(4.*x*y**5.+8.*ze*ze*y**5.-40.*ze*zb*x*y**4.
     1 +20.*zb*zb*x*x*y**3.)
     1 +z7*y**4.*(2.*y*y-24.*ze*zb*y+30.*zb*zb*x)-z8*42.*zb*zb*y**5.
      end
      function tj3(x,y,z1,z2,z3,z4,z5,z6,z7,z8,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tj3=168.*z1*ze*x**4.*(3.*x+10.*ze*ze)-24.*z2*(15.*ze*x**4.*y
     1 -1.5*zb*x**5.-15.*ze*ze*zb*x**4.+40.*ze**3.*x**3.*y)
     1 +120.*z3*(2.*ze*x**3.*y*y-.5*zb*x**4.*y+4.*ze**3.*(x*y)**2.
     1 -4.*ze*ze*zb*x**3.*y+.5*ze*zb*zb*x**4.)
     1 -48.*z4*(3.*ze*x*x*y**3.-1.5*zb*x**3.*y*y+4.*ze**3.*x*y**3.
     1 -9.*ze*ze*zb*(x*y)**2.+3.*ze*zb*zb*x**3.*y-.125*zb**3.*x**4.)
     1 +24.*z5*(3.*ze*x*y**4.-3.*zb*x**2.*y**3.+2.*ze**3.*y**4.
     1 -12.*ze*ze*zb*x*y**3.+9.*ze*zb*zb*(x*y)**2.
     1 -zb**3.*x**3.*y)
     1 -24.*z6*(ze*y**5.-2.5*zb*x*y**4.-5.*ze*ze*zb*y**4.
     1 +10.*ze*zb*zb*x*y**3.-2.5*zb**3.*(x*y)**2.)
     1 +24.*z7*(-1.5*zb*y**5.+7.5*ze*zb*zb*y**4.-5.*zb**3.*x*y**3.)
     1 +210.*z8*zb**3.*y**4.
      end
      function tj4(x,y,z1,z2,z3,z4,z5,z6,z7,z8,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tj4=168.*z1*(3.*x**5.+60.*ze*ze*x**4.+80.*ze**4.*x**3.)
     1-120.*z2*(3.*x**4.*y+48.*ze*ze*x**3.*y-12.*ze*zb*x**4.
     1 -32.*ze**3.*zb*x**3.+48.*ze**4.*x*x*y)
     1+240.*z3*(x**3.*y*y+12.*(ze*x*y)**2.-8.*ze*zb*x**3.*y
     1 +.5*zb*zb*x**4.+8.*ze**4.*x*y*y-16.*ze**3.*zb*x*x*y
     1 +4.*(ze*zb)**2.*x**3.)     
     1 -48.*z4*(3.*x*x*y**3.+24.*ze*ze*x*y**3.-36.*ze*zb*(x*y)**2.
     1 +6.*zb*zb*x**3.*y+8.*ze**4.*y**3.-48.*ze**3.*zb*x*y*y
     1 +36.*(ze*zb*x)**2.*y-4.*ze*(zb*x)**3.)
     1 +24.*z5*(3.*x*y**4.+12.*ze*ze*y**4.-48.*ze*zb*x*y**3.
     1  +18.*(zb*x*y)**2.-32.*ze**3.*zb*y**3.+72.*(ze*zb)**2.*x*y*y
     1  -24.*ze*zb**3.*x*x*y+zb**4.*x**3.)
     1 -24.*z6*(y**5.-20.*ze*zb*y**4.+40.*(ze*zb)**2.*y**3.
     1 +20.*zb*zb*x*y**3.-40.*ze*zb**3.*x*y*y+5.*zb**4.*x*x*y)
     1 +240.*z7*(1.5*zb*zb*y**4.-4.*ze*(zb*y)**3.+1.5*zb**4.*x*y*y)
     1 -840.*z8*zb**4.*y**3.
      end
      function tj5(x,y,z1,z2,z3,z4,z5,z6,z7,z8,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tj5=1680.*z1*(15.*ze*x**4.+80.*(ze*x)**3.+48.*ze**5.*x*x)
     1 -360.*z2*(40.*ze*x**3.*y-5.*zb*x**4.+160.*ze**3.*x*x*y
     1 -80.*ze*ze*zb*x**3.-80.*ze**4.*zb*x*x+64.*ze**5.*x*y)
     1 +240.*z3*(30.*ze*(x*y)**2.-10.*zb*x**3.*y+80.*ze**3.*x*y*y
     1 -120.*ze*ze*zb*x*x*y+20.*ze*zb*zb*x**3.+16.*ze**5.*y*y
     1 -80.*ze**4.*zb*x*y+40.*ze**3.*(zb*x)**2.)
     1 -48.*z4*(60.*ze*x*y**3.-45.*zb*(x*y)**2.+48.*(ze*y)**3.
     1 -360.*ze*ze*zb*x*y*y+180.*ze*zb*zb*x*x*y-10.*(zb*x)**3.
     1 +32.*(ze*y)**3.-120.*ze**4.*zb*y*y+240.*ze**3.*zb*zb*x*y
     1 -60.*ze*ze*zb**3.*x*x)
     1 +24.*z5*(30.*ze*y**4.-60.*zb*x*y**3.-240.*ze*ze*zb*y**3.
     1 +360.*ze*zb*zb*x*y*y-60.*zb**3.*x*x*y+240.*ze**3.*zb*zb*y*y
     1 -240.*ze*ze*zb**3.*x*y+30.*ze*zb**4.*x*x)
     1-24.*z6*(-25.*zb*y**4.+200.*ze*zb*zb*y**3.-200.*ze*ze*zb**3.*y*y
     1 -100.*zb**3.*x*y*y+100.*ze*zb**4.*x*y-5.*zb**5.*x*x)
     1 +240.*z7*(-10.*(zb*y)**3.+15.*ze*zb**4.*y*y-3.*zb**5.*x*y)
     1 +840.*z8*3.*zb**5.*y*y
      end
      function tj6(x,y,z1,z2,z3,z4,z5,z6,z7,z8,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tj6=5040.*z1*(5.*x**4.+120.*ze*ze*x**3.+240.*ze**4.*x*x
     1 +64.*ze**6.*x)
     1 -2880.*z2*(5.*x**3.*y+90.*ze*ze*x*x*y-30.*ze*zb*x**3.
     1+120.*ze**4.*x*y-120.*ze**3.*zb*x*x-48.*ze**5.*zb*x+16.*ze**6.*y)
     1 +1440.*z3*(5.*(x*y)**2.+60.*ze*ze*x*y*y-60.*ze*zb*x*x*y
     1 +5.*zb*zb*x**3.+40.*ze**4.*y*y-160.*ze**3.*zb*x*y
     1 +60.*(ze*zb*x)**2.-32.*ze**5.*zb*y+40.*ze**4.*zb*zb*x)
     1 -192.*z4*(15.*x*y**3.+66.*ze*ze*y**3.-270.*ze*zb*x*y*y
     1 +.5*135.*zb*zb*x*x*y-360.*ze**3.*zb*y*y+540.*ze*ze*zb*zb*x*y
     1 -90.*ze*zb**3.*x*x+24.*ze*ze*y**3.+180.*ze**4.*zb*zb*y
     1 -120.*(ze*zb)**3.*x)
     1 +144.*z5*(5.*y**4.-120.*ze*zb*y**3.+360.*(ze*zb*y)**2.
     1 +90.*zb*zb*x*y*y-240.*ze*zb**3.*x*y+15.*zb**4.*x*x
     1 -160.*(ze*zb)**3.*y+60.*ze*ze*zb**4.*x)
     1 -192.*z6*(37.5*zb*zb*y**3.-150.*ze*zb**3.*y*y
     1 +75.*ze*ze*zb**4.*y+37.5*zb**4.*x*y-15.*ze*zb**5.*x)
     1 +720.*z7*(10.*zb**4.*y*y+5.*zb**4.*y*y-12.*ze*zb**5.*y
     1 +zb**6.*x)
     1 -5040.*z8*zb**6.*y
      end

      function tk0(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tk0=x**8.*z1-x**7.*y*z2+x**6.*y*y*z3-x**5.*y**3.*z4
     1 +(x*y)**4.*z5-x**3.*y**5.*z6+x*x*y**6.*z7-x*y**7.*z8+y**8.*z9
      end
      function tk1(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tk1=16.*x**7.*z1*ze-z2*x**6.*(14.*ze*y-zb*x)
     1 +z3*2.*x**5.*y*(6.*ze*y-zb*x)-z4*x**4.*y*y*(10.*ze*y-3.*zb*x)
     1 +(x*y)**3.*z5*(8.*y*ze-4.*zb*x)-z6*x*x*y**4.*(6.*ze*y-5.*zb*x)
     1 +z7*x*y**5.*(4.*ze*y-6.*zb*x)-z8*y**6.*(2.*ze*y-7.*zb*x)
     1 -8.*z9*zb*y**7.
      end
      function tk2(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tk2=16.*z1*(x**7.+14.*ze*ze*x**6.)
     1 -14.*z2*(x**6.*y+12.*ze*ze*x**5.*y-2.*ze*zb*x**6.)
     1 +2.*z3*(6.*x**5.*y*y+60.*ze*ze*x**4.*y*y-24.*ze*zb*x**5.*y
     1 +zb*zb*x**6.)
     1 -2.*z4*(5.*x**4.*y**3.+40.*ze*ze*(x*y)**3.-30.*ze*zb*x**4.*y*y
     1 +3.*zb*zb*x**5.*y)
     1 +8.*z5*(x**3.*y**4.+6.*ze*ze*x*x*y**4.-8.*ze*zb*(x*y)**3.
     1 +1.5*zb*zb*x**4.*y*y)
     1 -2.*z6*(3.*x*x*y**5.+12.*ze*ze*x*y**5.-30.*ze*zb*x*x*y**4.
     1 +10.*zb*zb*(x*y)**3.)
     1 +2.*z7*(2.*x*y**6.+4.*ze*ze*y**6.-24.*ze*zb*x*y**5.
     1 +15.*zb*zb*x*x*y**4.)
     1 -2.*z8*(y**7.-14.*ze*zb*y**6.+21.*zb*zb*x*y**5.)
     1 +56.*z9*zb*zb*y**6.
      end
      function tk3(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tk3=672.*z1*(ze*x**6.+4.*ze**3.*x**5.)
     1 -42.*z2*(12.*ze*x**5.*y-zb*x**6.+40.*ze**3.*x**4.*y
     1 -12.*ze*ze*zb*x**5.)
     1 +24.*z3*(15.*ze*x**4.*y*y-3.*zb*x**5.*y+40.*ze**3.*x**3.*y*y
     1        -30.*ze*ze*zb*x**4.*y+3.*ze*zb*zb*x**5.)
     1 -6.*z4*(40.*ze*(x*y)**3.-15.*zb*x**4.*y*y+80.*ze**3.*x*x*y**3.
     1    -120.*ze*ze*zb*x**3.*y*y+30.*ze*zb*zb*x**4.*y-zb**3.*x**5.)
     1 +24.*z5*(6.*ze*x*x*y**4.-4.*zb*(x*y)**3.+8.*ze**3.*x*y**4.
     1 -24.*ze*ze*zb*x*x*y**3.+12.*ze*zb*zb*x**3.*y*y-zb**3.*x**4.*y)
     1 -12.*z6*(6.*ze*x*y**5.-7.5*zb*x*x*y**4.+4.*ze**3.*y**5.
     1 -30.*ze*ze*zb*x*y**4.+30.*ze*zb*zb*x*x*y**3.-5.*(zb*x)**3.*y*y)
     1 +24.*z7*(ze*y**6.-3.*zb*x*y**5.-6.*ze*ze*zb*y**5.
     1        +15.*ze*zb*zb*x*y**4.-5.*zb**3.*x*x*y**3.)
     1 -84.*z8*(-.5*zb*y**6.+3.*ze*zb*zb*y**5.-2.5*zb**3.*x*y**4.)
     1 -336.*z9*zb**3.*y**5.
      end
      function tk4(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tk4=672.*z1*(x**6.+24.*ze**2.*x**5.+40.*(ze*x)**4.)
     1 -42.*z2*(12.*x**5.*y+240.*ze*ze*x**4.*y-24.*ze*zb*x**5.
     1 +320.*ze**4.*x**3.*y-160.*ze**3.*zb*x**4.-24.*ze*zb*x**5.)
     1 +24.*z3*(15.*x**4.*y*y+240.*ze*ze*x**3.*y*y-120.*ze*zb*x**4.*y
     1 +6.*zb*zb*x**5.+240.*ze**4.*(x*y)**2.-320.*ze**3.*zb*x**3.*y
     1 +60.*ze*ze*zb*zb*x**4.)
     1 -6.*z4*(40.*(x*y)**3.+480.*ze*ze*x*x*y**3.-480.*ze*zb*x**3.*y*y
     1  +60.*zb*zb*x**4.*y+320.*ze**4.*x*y**3.-960.*ze**3.*zb*(x*y)**2.
     1  +480.*(ze*zb)**2.*x**3.*y-40.*ze*zb**3.*x**4.)
     1 +24.*z5*(6.*x*x*y**4.+48.*ze*ze*x*y**4.-96.*ze*zb*x*x*y**3.
     1 +24.*zb*zb*x**3.*y*y+16.*(ze*y)**4.-128.*ze**3.*zb*x*y**3.
     1 +144.*(ze*zb*x*y)**2.-32.*ze*(zb*x)**3.*y+(zb*x)**4.)
     1 -12.*z6*(6.*x*y**5.+24.*ze*ze*y**5.-120.*ze*zb*x*y**4.
     1+60.*zb*zb*x*x*y**3.-80.*ze**3.*zb*y**4.+240.*ze*ze*zb*zb*x*y**3.
     1 -120.*ze*zb**3.*(x*y)**2.+10.*zb**4.*x**3.*y)
     1 +24.*z7*(y**6.-24.*ze*zb*y**5.+30.*zb*zb*x*y**4.
     1 +60.*(ze*zb)**2.*y**4.-80.*ze*zb**3.*x*y**3.
     1 +15.*zb**4.*(x*y)**2.)
     1 -84.*z8*(6.*zb*zb*y**5.-20.*ze*zb**3.*y**4.+10.*zb**4.*x*y**3.)
     1 +1680.*z9*(zb*y)**4.
      end
      function tk5(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tk5=13440.*z1*(3.*ze*x**5.+20.*ze**3.*x**4.+16.*ze**5.*x**3.)
     1 -840.*z2*(30.*ze*x**4.*y-3.*zb*x**5.+160.*(ze*x)**3.*y
     1 -60.*ze*ze*zb*x**4.+96.*ze**5.*x*x*y-80.*ze**4.*zb*x**3.)
     1 +240.*z3*(60.*ze*x**3.*y*y-15.*zb*x**4.*y+240.*ze**3.*(x*y)**2.
     1 -240.*ze*ze*zb*x**3.*y+30.*ze*zb*zb*x**4.+96.*ze**5.*x*y*y
     1 -240.*ze**4.*zb*x*x*y+80.*ze**3.*zb*zb*x**3.)
     1 -240.*z4*(30.*ze*x*x*y**3.-15.*zb*x**3.*y*y+80.*ze**3.*x*y**3.
     1 -180.*ze*ze*zb*(x*y)**2.+60.*ze*zb*zb*x**3.*y-2.5*zb**3.*x**4.
     1 +16.*ze**5.*y**3.-120.*ze**4.*zb*x*y*y+120.*ze**3.*zb*zb*x*x*y
     1 -20.*ze*ze*(zb*x)**3.)
     1 +960.*z5*(3.*ze*x*y**4.-3.*zb*x*x*y**3.+4.*ze**3.*y**4.
     1 -24.*ze*ze*zb*x*y**3.+18.*ze*zb*zb*(x*y)**2.-2.*zb**3.*x**3.*y
     1-8.*ze**4.*zb*y**3.+24.*ze**3.*zb*zb*x*y*y-12.*ze*ze*zb**3.*x*x*y
     1 +ze*zb**4.*x**3.)
     1 -120.*z6*(6.*ze*y**5.-15.*zb*x*y**4.-60.*ze*ze*zb*y**4.
     1 +120.*ze*zb*zb*x*y**3.-30.*zb**3.*(x*y)**2.+80.*(ze*y)**3.*zb*zb
     1 -120.*ze*ze*zb**3.*x*y*y+30.*ze*zb**4.*x*x*y-zb**5.*x**3.)
     1 +48.*z7*(-15.*zb*y**5.+150.*ze*zb*zb*y**4.-100.*(zb*y)**3.*x
     1 -200.*ze*ze*(zb*y)**3.+150.*ze*zb**4.*x*y*y-15.*zb**5.*x*x*y)
     1 -840.*z8*(-5.*zb**3.*y**4.+10.*ze*zb**4.*y**3.-3.*zb**5.*x*y*y)
     1 -6720.*z9*zb**5.*y**3.
      end
      function tk6(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tk6=40320.*z1*(x**5.+30.*ze*ze*x**4.+80.*ze**4.*x**3.
     1  +32.*ze**6.*x*x)
     1 -5040.*z2*(5.*x**4.*y+120.*ze*ze*x**3.*y-30.*ze*zb*x**4.
     1 +240.*ze**4.*x*x*y-160.*(ze*x)**3.*zb+64.*ze**6.*x*y
     1 -96.*ze**5.*zb*x*x)
     1 +720.*z3*(20.*x**3.*y*y+360.*(ze*x*y)**2.-240.*ze*zb*x**3.*y
     1 +15.*zb*zb*x**4.+480.*ze**4.*x*y*y-960.*ze**3.*zb*x*x*y
     1 +240.*ze*ze*zb*zb*x**3.+64.*ze**6.*y*y-384.*ze**5.*zb*x*y
     1 +240.*ze**4.*zb*zb*x*x)
     1 -1440.*z4*(5.*x*x*y**3.+60.*ze*ze*x*y**3.-90.*ze*zb*x*x*y*y
     1 +15.*zb*zb*x**3.*y+40.*ze**4.*y**3.-240.*ze**3.*zb*x*y*y
     1 +180.*ze*ze*zb*zb*x*x*y-20.*ze*(zb*x)**3.-48.*ze**5.*zb*y*y
     1 +120.*ze**4.*zb*zb*x*y-40.*(ze*zb)**3.*x*x)
     1 +2880.*z5*(x*y**4.+6.*ze*ze*y**4.-24.*ze*zb*x*y**3.
     1 +9.*zb*zb*x*x*y*y+72.*ze*ze*zb*zb*x*y*y-32.*(ze*y)**3.*zb
     1 -24.*ze*zb**3.*x*x*y+zb**4.*x**3.+24.*ze**4.*zb*zb*y*y
     1 -32.*(ze*zb)**3.*x*y+6.*ze*ze*zb**4.*x*x)
     1 -720.*z6*(y**5.-30.*ze*zb*y**4.+30.*zb*zb*x*y**3.
     1 +120.*ze*ze*zb*zb*y**3.-120.*ze*zb**3.*x*y*y+15.*zb**4.*x*x*y
     1 -80.*(ze*zb)**3.*y*y+60.*ze*ze*zb**4.*x*y-6.*ze*zb**5.*x*x)
     1 +720.*z7*(15.*zb*zb*y**4.-80.*ze*(zb*y)**3.+30.*zb**4.*x*y*y
     1 +60.*ze*ze*zb**4.*y*y-24.*ze*zb**5.*x*y+zb**6.*x*x)
     1 -5040.*z8*(5.*zb**4.*y**3.-6.*ze*zb**5.*y*y+zb**6.*x*y)
     1 +20160.*z9*zb**6.*y*y
      end

      function tl0(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tl0=x**9.*z1-x**8.*y*z2+x**7.*y*y*z3-x**6.*y**3.*z4
     1 +x**5.*y**4.*z5-x**4.*y**5.*z6+x**3.*y**6.*z7-x*x*y**7.*z8
     1 +x*y**8.*z9-y**9.*z10
      end
      function tl1(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tl1=18.*ze*x**8.*z1-z2*x**7.*(16.*ze*y-zb*x)
     1 +z3*2.*x**6.*y*(7.*ze*y-zb*x)
     1 -z4*3.*x**5.*y*y*(4.*ze*y-zb*x)
     1 +z5*2.*x**4.*y**3.*(5.*ze*y-2.*zb*x)
     1 -z6*x**3.*y**4.*(8.*ze*y-5.*zb*x)
     1 +z7*6.*x*x*y**5.*(ze*y-zb*x)
     1 -z8*x*y**6.*(4.*ze*y-7.*zb*x)
     1 +z9*2.*y**7.*(ze*y-4.*zb*x)+9.*z10*zb*y**8.
      end
      function tl2(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tl2=18.*x**7.*z1*(x+16.*ze*ze)
     1 -16.*z2*x**6.*(x*y+14.*ze*ze*y-2.*ze*zb*x)
     1 +z3*2.*x**5.*(7.*x*y*y+84.*ze*ze*y*y-28.*ze*zb*x*y+zb*zb*x*x)
     1 -z4*6.*x**4.*y*(2.*x*y*y+20.*ze*ze*y*y-12.*ze*zb*x*y+zb*zb*x*x)
     1 +z5*2.*x**3.*y*y*(5.*x*y*y+40.*ze*ze*y*y-40.*ze*zb*x*y
     1 +6.*zb*zb*x*x)
     1 -4.*z6*x*x*y**3.*(2.*x*y*y+12.*ze*ze*y*y-20.*ze*zb*x*y
     1 +5.*zb*zb*x*x)
     1 +z7*6.*x*y**4.*(x*y*y+4.*ze*ze*y*y-12.*ze*zb*x*y+5.*zb*zb*x*x)
     1 -2.*z8*y**5.*(2.*x*y*y+4.*ze*ze*y*y-28.*ze*zb*x*y+21.*zb*zb*x*x)
     1 +z9*2.*y**6.*(y*y-16.*ze*zb*y+28.*zb*zb*x)-72.*z10*zb*zb*y**7.
      end
      function tl3(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tl3=288.*x**6.*z1*ze*(3.*x+14.*ze*ze)
     1 -48.*z2*x**5.*(14.*ze*x*y-zb*x*x+56.*ze**3.*y-14.*ze*ze*zb*x)
     1 +84.*z3*x**4.*(6.*ze*x*y*y-zb*x*x*y+20.*ze**3.*y*y
     1 -12.*ze*ze*zb*x*y+ze*zb*zb*x*x)
     1 -z4*6.*x**3.*(60.*ze*x*y**3.-18.*zb*x*x*y*y+160.*(ze*y)**3.
     1 -180.*ze*ze*zb*x*y*y+36.*ze*zb*zb*x*x*y-(zb*x)**3.)
     1 +z5*24.*x*x*y*(10.*ze*x*y**3.-5.*zb*x*x*y*y+20.*ze**3.*y*y*y
     1 -40.*ze*ze*zb*x*y*y+15.*ze*zb*zb*x*x*y-(zb*x)**3.)
     1 -12.*z6*x*y*y*(12.*ze*x*y*y*y-10.*zb*x*x*y*y+16.*(ze*y)**3.
     1 -60.*ze*ze*zb*x*y*y+40.*ze*zb*zb*x*x*y-5.*zb**3.*x**3.)
     1 +12.*z7*y**3.*(6.*ze*x*y**3.-9.*zb*x*x*y*y+4.*(ze*y)**3.
     1 -36.*ze*ze*zb*x*y*y+45.*ze*zb*zb*x*x*y-10.*(zb*x)**3.)
     1 -6.*z8*y**4.*(4.*ze*y**3.-14.*zb*x*y*y-28.*ze*ze*zb*y*y
     1 +84.*ze*zb*zb*x*y-35.*zb**3.*x*x)
     1 +48.*z9*y**5.*(-zb*y*y+7.*ze*zb*zb*y-7.*zb**3.*x)
     1 +504.*z10*zb*zb*zb*y**6.
      end
      function tl4(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tl4=864.*x**5.*z1*(x*x+28.*ze*ze*x+56.*ze**4.)
     1 -672.*z2*x**4.*(x*x*y+24.*ze*ze*x*y-4.*ze*zb*x*x+40.*ze**4.*y
     1 -16.*ze**3.*zb*x)
     1 +84.*z3*x**3.*(6.*(x*y)**2.+120.*ze*ze*x*y*y-48.*ze*zb*x*x*y
     1 +zb*zb*x**3.+160.*ze**4.*y*y-160.*ze**3.*zb*x*y+24.*(ze*zb*x)**2.
     1 +zb*zb*x**3.)
     1 -72.*z4*x**2.*(5.*x*x*y**3.+80.*ze*ze*x*y**3.-60.*ze*zb*x*x*y*y
     1 +6.*zb*zb*x**3.*y+80.*ze**4.*y**3.-160.*ze**3.*zb*x*y*y
     1 +60.*ze*ze*zb*zb*x*x*y-4.*ze*(zb*x)**3.)
     1 +24.*z5*(10.*x**3.*y**4.+120.*ze*ze*x*x*y**4.
     1 -160.*ze*zb*(x*y)**3.+30.*zb*zb*x**4.*y*y+80.*ze**4.*x*y**4.
     1 -320.*ze**3.*zb*x*x*y**3.+240.*ze*ze*zb*zb*x**3.*y*y
     1 -40.*ze*zb**3.*x**4.*y+zb**4.*x**5.)
     1 -24.*z6*(6.*x*x*y**5.+48.*ze*ze*x*y**5.-120.*ze*zb*x*x*y**4.
     1 +40.*zb*zb*(x*y)**3.+16.*ze**4.*y**5.-160.*ze**3.*zb*x*y**4.
     1 +240.*ze*ze*zb*zb*x*x*y**3.-80.*ze*(zb*x)**3.*y*y
     1 +5.*(zb*x)**4.*y)
     1 +24.*z7*(3.*x*y**6.+12.*ze*ze*y**6.-72.*ze*zb*x*y**5.
     1+45.*zb*zb*x*x*y**4.-48.*ze**3.*zb*y**5.+180.*ze*ze*zb*zb*x*y**4.
     1 -120.*ze*(zb*y)**3.*x*x+15.*zb**4.*x**3.*y*y)
     1 -24.*z8*(y**7.-28.*ze*zb*y**6.+42.*zb*zb*x*y**5.
     1 +84.*ze*ze*zb*zb*y**5.-140.*ze*zb**3.*x*y**4.
     1 +35.*zb**4.*x*x*y**3.)
     1 +336.*z9*(2.*zb*zb*y**6.-8.*ze*zb**3.*y**5.+5.*(zb*y)**4.*x)
     1 -3024.*z10*zb**4.*y**5.
      end
      function tl5(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tl5=60480.*z1*(ze*x**6.+8.*ze**3.*x**5.+8.*ze**5.*x**4.)
     1 -3360.*z2*(12.*ze*x**5.*y+80.*ze**3.*x**4.*y-zb*x**6.
     1 -24.*ze*ze*zb*x**5.+64.*ze**5.*x**3.*y-40.*(ze*x)**4.*zb)
     1 +1680.*z3*(15.*ze*x**4.*y*y-3.*zb*x**5.*y+80.*ze**3.*x**3.*y*y
     1 -60.*ze*ze*zb*x**4.*y+6.*ze*zb*zb*x**5.+48.*ze**5.*x*x*y*y
     1 -80.*ze**4.*zb*x**3.*y+20.*ze**3.*zb*zb*x**4.)
     1 -360.*z4*(40.*ze*(x*y)**3.-15.*zb*x**4.*y*y+160.*(ze*y)**3.*x*x
     1 -240.*ze*ze*zb*x**3.*y*y+60.*ze*zb*zb*x**4.*y-2.*zb**3.*x**5.
     1 +64.*ze**5.*x*y**3.-240.*ze**4.*zb*x*x*y*y
     1 +160.*(ze*x)**3.*zb*zb*y-20.*ze*ze*zb**3.*x**4.)
     1 +240.*z5*(30.*ze*x*x*y**4.-20.*zb*(x*y)**3.+80.*ze**3.*x*y**4.
     1 -240.*ze*ze*zb*x*x*y**3.+120.*ze*zb*zb*x**3.*y*y
     1 -10.*zb**3.*x**4.*y+16.*ze**5.*y**4.-160.*ze**4.*zb*x*y**3.
     1 +240.*ze**3.*zb*zb*x*x*y*y-80.*ze*ze*zb**3.*x**3.*y
     1 +5.*ze*zb**4.*x**4.)
     1 -120.*z6*(24.*ze*x*y**5.-30.*zb*x*x*y**4.-240.*ze*ze*zb*x*y**4.
     1 +32.*ze**3.*y**5.+240.*ze*zb*zb*x*x*y**3.-40.*(zb*x)**3.*y*y
     1 -80.*(ze*y)**4.*zb+320.*(ze*y)**3.*zb*zb*x
     1 -240.*ze*ze*zb**3.*x*x*y*y+40.*ze*zb**4.*x**3.*y-zb**5.*x**4.)
     1 +144.*z7*(5.*ze*y**6.-15.*zb*x*y**5.-60.*ze*ze*zb*y**5.
     1 +150.*ze*zb*zb*x*y**4.-50.*zb**3.*x*x*y**3.
     1 +100.*ze**3.*zb*zb*y**4.-200.*ze*ze*zb**3.*x*y**3.
     1 +75.*ze*zb**4.*x*x*y*y-5.*zb**5.*x**3.*y)
     1 -168.*z8*(-5.*zb*y**6.+60.*ze*zb*zb*y**5.-50.*zb**3.*x*y**4.
     1 -100.*ze*ze*zb**3.*y**4.+100.*ze*zb**4.*x*y**3.
     1 -15.*zb**5.*x*x*y*y)
     1 +3360.*z9*(-2.*zb**3.*y**5.+5.*ze*(zb*y)**4.-2.*zb**5.*x*y**3.)
     1 +15120.*z10*zb**5.*y**4.
      end
      function tl6(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tl6=60480.*z1*x**3.*(x**3.+36.*ze*ze*x*x+120.*ze**4*x+64.*ze**6.)
     1 -40320.*z2*(x**5.*y+30.*ze*ze*x**4.*y-6.*ze*zb*x**5.
     1 +80.*ze**4.*x**3.*y-40.*ze**3.*zb*x**4.+32.*ze**6.*x*x*y
     1 -32.*ze**5.*zb*x**3.)
     1 +5040.*z3*(5.*x**4.*y*y+120.*ze*ze*x**3.*y*y-60.*ze*zb*x**4.*y
     1 +60.*ze*ze*zb*zb*x**4.+3.*zb*zb*x**5.+240.*ze**4.*x*x*y*y
     1 -320.*(ze*x)**3.*zb*y+64.*ze**6.*x*y*y-192.*ze**5.*zb*x*x*y
     1 +80.*ze**4.*zb*zb*x**3.)
     1-720.*z4*(20.*(x*y)**3.+360.*ze*ze*x*x*y**3.-360.*ze*zb*x**3.*y*y
     1 +45.*zb*zb*x**4.*y+480.*ze**4.*x*y**3.-1440.*ze**3.*zb*x*x*y*y
     1 +720.*ze*ze*zb*zb*x**3.*y-60.*ze*zb**3.*x**4.+64.*ze**6.*y**3.
     1-576.*ze**5.*zb*x*y*y+720.*ze**4.*zb*zb*x*x*y-160.*(ze*zb*x)**3.)
     1 +720.*z5*(10.*x*x*y**4.+120.*ze*ze*x*y**4.-240.*ze*zb*x*x*y**3.
     1 +60.*zb*zb*x**3.*y*y+80.*ze**4.*y**4.-640.*(ze*y)**3.*zb*x
     1 +720.*ze*ze*zb*zb*x*x*y*y-160.*ze*(zb*x)**3.*y+5.*(zb*x)**4.
     1 -128.*ze**5.*zb*y**3.+480.*ze**4.*zb*zb*x*y*y
     1 -320.*(ze*zb)**3.*x*x*y+40.*ze*ze*zb**4.*x**3.)
     1 -2880.*z6*(x*y**5.+6.*ze*ze*y**5.-30.*ze*zb*x*y**4.
     1+15.*zb*zb*x*x*y**3.-40.*ze**3.*zb*y**4.+120.*ze*ze*zb*zb*x*y**3.
     1 -60.*ze*zb**3.*x*x*y*y+5.*zb**4.*x**3.*y+40.*ze**4.*zb*zb*y**3.
     1-80.*(ze*zb)**3.*x*y*y+30.*ze*ze*zb**4.*x*x*y-2.*ze*zb**5.*x**3.)
     1 +720.*z7*(y**6.-36.*ze*zb*y**5.+45.*zb*zb*x*y**4.
     1 +180.*ze*ze*zb*zb*y**4.-240.*ze*(zb*y)**3.*x+45.*zb**4.*x*x*y*y
     1-160.*(ze*zb*y)**3.+180.*ze*ze*zb**4.*x*y*y-36.*ze*zb**5.*x*x*y
     1 +zb**6.*x**3.)
     1 -5040.*z8*(3.*zb*zb*y**5.-20.*ze*zb**3.*y**4.+10.*zb**4.*x*y**3.
     1 +20.*ze*ze*zb**4.*y**3.-12.*ze*zb**5.*x*y*y+zb**6.*x*x*y)
     1 +10080.*z9*(5.*(zb*y)**4.-8.*ze*zb**5.*y**3.+2.*zb**6.*x*y*y)
     1 -60480.*z10*zb**6.*y**3.
      end

    
      function tm0(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tm0=x**10.*z1-x**9.*y*z2+x**8.*y*y*z3-x**7.*y**3.*z4
     1 +x**6.*y**4.*z5-x**5.*y**5.*z6+x**4.*y**6.*z7-x**3.*y**7.*z8
     1 +x*x*y**8.*z9-x*y**9.*z10+y**10.*z11
      end
      function tm1(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tm1=20.*ze*x**9.*z1-z2*x**8.*(18.*ze*y-zb*x)
     1 +z3*2.*x**7.*y*(8.*ze*y-zb*x)
     1 -z4*x**6.*y*y*(14.*ze*y-3.*zb*x)
     1 +z5*2.*x**5.*y**3.*(6.*ze*y-2.*zb*x)
     1 -z6*x**4.*y**4.*(10.*ze*y-5.*zb*x)
     1 +z7*x**3.*y**5.*(8.*ze*y-6.*zb*x)
     1 -z8*x*x*y**6.*(6.*ze*y-7.*zb*x)
     1 +z9*4.*x*y**7.*(ze*y-2.*zb*x)-z10*y**8.*(2.*ze*y-9.*zb*x)
     1 -10.*z11*zb*y**9.
      end
      function tm2(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tm2=20.*z1*x**8.*(x+18.*ze*ze)
     1 -z2*18.*x**7.*(x*y+16.*ze*ze*y-2.*zb*ze*x)
     1 +2.*z3*x**6.*(8.*x*y*y+112.*ze*ze*y*y-32.*ze*zb*x*y+zb*zb*x*x)
     1 -2.*z4*x**5.*y*(7.*x*y*y+84.*ze*ze*y*y-42.*ze*zb*x*y
     1 +3.*zb*zb*x*x)
     1 +12.*z5*x**4.*y**2.*(x*y*y+10.*ze*ze*y*y-8.*ze*zb*x*y+zb*zb*x*x)
     1-10.*z6*(x*y)**3.*(x*y*y+8.*ze*ze*y*y-10.*ze*zb*x*y+2.*zb*zb*x*x)
     1 +2.*z7*x*x*y**4.*(4.*x*y*y+24.*ze*ze*y*y-48.*ze*zb*x*y
     1 +15.*zb*zb*x*x)
     1 -6.*z8*x*y**5.*(x*y*y+4.*ze*ze*y*y-14.*ze*zb*x*y+7.*zb*zb*x*x)
     1 +4.*z9*y**6.*(x*y*y+2.*ze*ze*y*y-16.*ze*zb*x*y+14.*zb*zb*x*x)
     1 -2.*z10*y**7.*(y*y-18.*ze*zb*y+36.*zb*zb*x)
     1 +90.*z11*zb*zb*y**8.
      end
      function tm3(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tm3=360.*z1*x**7.*ze*(3.*x+16.*ze*ze)
     1-z2*18.*x**6.*(48.*ze*x*y-3.*zb*x*x+224.*ze**3.*y-48.*ze*ze*zb*x)
     1 +96.*z3*x**5.*(7.*ze*x*y*y-zb*x*x*y+28.*ze**3.*y*y
     1 -14.*ze*ze*zb*x*y+ze*zb*zb*x*x)
     1 -6.*z4*x**4.*(84.*ze*x*y**3.-21.*zb*x*x*y*y+280.*(ze*y)**3.
     1 -252.*ze*ze*zb*x*y*y+42.*ze*zb*zb*x*x*y-(zb*x)**3.)
     1 +24.*z5*x**3.*y*(15.*ze*x*y**3.-6.*zb*(x*y)**2.+40.*(ze*y)**3.
     1 -60.*ze*ze*zb*x*y*y+18.*ze*zb*zb*x*x*y-(zb*x)**3.)
     1-10.*z6*(x*y)**2.*(24.*ze*x*y**3.-15.*zb*x*x*y*y+48.*(ze*y)**3.
     1 -120.*ze*ze*zb*x*y*y+60.*ze*zb*zb*x*x*y-6.*(zb*x)**3.)
     1 +24.*z7*x*y**3.*(6.*ze*x*y**3.-6.*zb*x*x*y*y+8.*(ze*y)**3.
     1 -36.*ze*ze*zb*x*y*y+30.*ze*zb*zb*x*x*y-5.*(zb*x)**3.)
     1 -6.*z8*y**4.*(12.*ze*x*y**3.-21.*zb*x*x*y*y+8.*(ze*y)**3.
     1 -84.*ze*ze*zb*x*y*y+126.*ze*zb*zb*x*x*y-35.*(zb*x)**3.)
     1 +24.*z9*y**5.*(ze*y**3.-4.*zb*x*y*y-8.*ze*ze*zb*y*y
     1 +28.*ze*zb*zb*x*y-14.*zb**3.*x*x)
     1 -18.*z10*y**6.*(-3.*zb*y*y+24.*ze*zb*zb*y-28.*zb**3.*x)
     1 -720.*z11*zb*zb*zb*y**7.
      end
      function tm4(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tm4=360.*z1*x**6.*(3.*x*x+96.*ze*ze*x+224.*ze**4.)
     1-z2*288.*x**5.*(3.*x*x*y+84.*ze*ze*x*y-12.*zb*ze*x*x
     1 -56.*ze**3.*zb*x+168.*ze**4.*y)
     1 +96.*z3*x**4.*(7.*x*x*y*y+168.*ze*ze*x*y*y-56.*ze*zb*x*x*y
     1 +2.*zb*zb*x**3.+280.*ze**4.*y*y-224.*ze**3.*zb*x*y
     1 +28.*ze*ze*zb*zb*x*x)
     1 -24.*z4*x**3.*(21.*x*x*y**3.+420.*ze*ze*x*y**3.
     1 -252.*ze*zb*x*x*y*y+21.*zb*zb*x**3.*y+560.*ze**4.*y**3.
     1 -840.*ze**3.*zb*x*y*y+252.*ze*ze*zb*zb*x*x*y-14.*ze*(zb*x)**3.)
     1 +24.*z5*x*x*(15.*x*x*y**4.+240.*ze*ze*x*y**4.
     1 -240.*ze*zb*x*x*y**3.+36.*zb*zb*x**3.*y*y+240.*(ze*y)**4.
     1 -640.*(ze*y)**3.*zb*x+360.*(ze*zb*x*y)**2.-48.*ze*(zb*x)**3.*y
     1 +(zb*x)**4.) 
     1-40.*z6*x*y*(6.*x**2.*y**4.+72.*ze*ze*x*y**4.
     1 -120.*ze*zb*x*x*y**3.+30.*zb*zb*x**3.*y*y-240.*(ze*y)**3.*zb*x
     1 +48.*(ze*y)**4.+240.*(ze*zb*x*y)**2.-60.*ze*(zb*x)**3.*y
     1 +3.*(zb*x)**4.)
     1 +24.*z7*y*y*(6.*x*x*y**4.+48.*ze*ze*x*y**4.-144.*ze*zb*x*x*y**3.
     1 +60.*zb*zb*x**3.*y*y+16.*(ze*y)**4.-192.*(ze*y)**3.*zb*x
     1 +360.*(ze*zb*x*y)**2.-160.*ze*(zb*x)**3.*y+15.*(zb*x)**4.)
     1 -12.*z8*y**3.*(6.*x*y**4.+24.*ze*ze*y**4.-168.*ze*zb*x*y**3.
     1 +126.*(zb*x*y)**2.-112.*(ze*y)**3.*zb+504.*ze*ze*zb*zb*x*y*y
     1 -420.*ze*zb**3.*x*x*y+70.*zb**4.*x**3.)
     1 +24.*z9*y**4.*(y**4.-32.*ze*zb*y**3.+56.*zb*zb*x*y*y
     1 +112.*(ze*zb*y)**2.-224.*ze*zb**3.*x*y+70.*zb**4.*x*x)
     1 -144.*z10*y**5.*(6.*zb*zb*y*y-28.*ze*zb**3.*y+21.*zb**4.*x)
     1 +5040.*z11*zb**4.*y**6.
      end
      function tm5(x,y,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,ze,zb)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      tm5=5760.*z1*x**5.*ze*(15.*x*x+140.*ze*ze*x+168.*ze**4.)
     1-z2*4320.*x**4.*(14.*ze*x*x*y-zb*x**3.+112.*ze**3.*x*y
     1 -28.*ze*ze*zb*x*x-56.*ze**4.*zb*x+112.*ze**5.*y)
     1 +6720.*z3*x**3.*(6.*ze*x*x*y*y-zb*x**3.*y+40.*ze**3.*x*y*y
     1 -24.*ze*ze*zb*x*x*y+2.*ze*zb*zb*x**3.+32.*ze**5.*y*y
     1 -40.*ze**4.*zb*x*y+8.*ze**3.*zb*zb*x*x)
     1 -168.*z4*x*x*(150.*ze*x*x*y**3.-45.*zb*x**3.*y*y
     1 +800.*(ze*y)**3.*x-900.*ze*ze*zb*x*x*y*y+180.*ze*zb*zb*x**3.*y
     1 -5.*zb**3.*x**4.+480.*ze**5.*y**3.-1200.*ze**4.*zb*x*y*y
     1 +600.*ze**3.*zb*zb*x*x*y-60.*ze*ze*(zb*x)**3.)
     1 +1440.*z5*x*(10.*ze*x*x*y**4.-5.*zb*x**3.*y**3.
     1 +40.*ze**3.*x*y**4.-80.*ze*ze*zb*x*x*y**3.
     1 +30.*ze*zb*zb*x**3.*y*y-2.*zb**3.*x**4.*y+16.*ze**5.*y**4.
     1 -80.*ze**4.*zb*x*y**3.+80.*ze**3.*zb*zb*x*x*y*y
     1 -20.*ze*ze*(zb*x)**3.*y+ze*(zb*x)**4.) 
     1-120.*z6*(60.*ze*x*x*y**5.-50.*zb*x**3.*y**4.+160.*ze**3.*x*y**5.
     1 -600.*ze*ze*zb*x*x*y**4.+400.*ze*zb*zb*(x*y)**3.
     1 -50.*zb**3.*x**4.*y*y+800.*(ze*y)**3.*zb*zb*x*x
     1 -400.*(ze*y)**4.*zb*x+32.*(ze*y)**5.-400.*ze*ze*(zb*x)**3.*y*y
     1 +50.*ze*(zb*x)**4.*y-(zb*x)**5.)
     1+240.*z7*y*(12.*ze*x*y**5.-144.*ze*ze*zb*x*y**4.-18.*zb*x*x*y**4.
     1 +16.*ze**3.*y**5.+180.*ze*zb*zb*x*x*y**3.-40.*(zb*x)**3.*y*y
     1 -48.*(ze*y)**4.*zb+240.*(ze*y)**3.*zb*zb*x
     1 -240.*zb**3.*(ze*x*y)**2.+60.*ze*zb**4.*x**3.*y-3.*zb**5.*x**4.)
     1 -360.*z8*y**2.*(2.*ze*y**5.-7.*zb*x*y**4.-28.*ze*ze*zb*y**4.
     1 +84.*ze*zb*zb*x*y**3.-35.*zb**3.*x*x*y*y+56.*ze**3.*zb*zb*y**3.
     1 -140.*ze*ze*zb**3.*x*y*y+70.*ze*zb**4.*x*x*y-7.*zb**5.*x**3.)
     1 +960.*z9*y**3.*(-zb*y**4.+14.*ze*zb*zb*y**3.-14.*zb**3.*x*y*y
     1 -28.*ze*ze*zb**3.*y*y+35.*ze*zb**4.*x*y-7.*zb**5.*x*x)
     1 -5040.*z10*y**4.*(-2.*zb**3.*y*y+6.*ze*zb**4.*y-3.*zb**5.*x)
     1 -30240.*z11*(zb*y)**5.
      end




c integrale N1N2exp(-ar)exp(-br) r**n1 r**n2 r**2 dr
c non slater
      function fb(a,b,n1,n2)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
c      fb=fa(a,n1)*fa(b,n2)*fac(n1+n2+2)/((a+b)**(n1+n2+3))
      fb=fac(n1+n2+2)/((a+b)**(n1+n2+3))
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





      function tcs(a,az,adx,ady,adz,akx,aky,akz,n)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      dimension us(12)
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
      ucd=1.d0+xi*0.d0
      
      
c 1s
      us10=q1(qa,qz-1.,r)*q0(qab,-qz)
      us11=q0(qa,qz-1.)*q1(qab,-qz,x)
      us(1)=-(us10+us11)*aps*ucd

c 2s
      us20=q2(qa,qz-1.,r)*q0(qab,-qz)
      us21=q2(qab,-qz,x)*q0(qa,qz-1.)
      us22=2.*q1(qa,qz-1.,r)*q1(qab,-qz,x)
      us(2)=(us20+us21+us22)*aps*ucd

c 3s
      us30=q3(qa,qz-1.,r)*q0(qab,-qz)
      us31=q0(qa,qz-1.)*q3(qab,-qz,x)
      us32=3.*q2(qa,qz-1.,r)*q1(qab,-qz,x)
      us33=3.*q2(qab,-qz,x)*q1(qa,qz-1.,r)
      us(3)=-(us30+us31+us32+us33)*aps*ucd

c 4s
      us40=q4(qa,qz-1.,r)*q0(qab,-qz)
      us41=q4(qab,-qz,x)*q0(qa,qz-1.)
      us42=4.*q3(qa,qz-1.,r)*q1(qab,-qz,x)
      us43=4.*q3(qab,-qz,x)*q1(qa,qz-1.,r)
      us44=6.*q2(qa,qz-1.,r)*q2(qab,-qz,x)
      us(4)=(us40+us41+us42+us43+us44)*aps*ucd
     
c 5s
      us50=q5(qa,qz-1.,r)*q0(qab,-qz)
      us51=q5(qab,-qz,x)*q0(qa,qz-1.)
      us52=5.*q4(qa,qz-1.,r)*q1(qab,-qz,x)
      us53=5.*q4(qab,-qz,x)*q1(qa,qz-1.,r)
      us54=10.*q3(qa,qz-1.,r)*q2(qab,-qz,x)
      us55=10.*q3(qab,-qz,x)*q2(qa,qz-1.,r)
      us(5)=-(us50+us51+us52+us53+us54+us55)*aps*ucd
     
c 6s
      us60=q6(qa,qz-1.,r)*q0(qab,-qz)
      us61=q6(qab,-qz,x)*q0(qa,qz-1.)
      us62=6.*q5(qa,qz-1.,r)*q1(qab,-qz,x)
      us63=6.*q5(qab,-qz,x)*q1(qa,qz-1.,r)
      us64=15.*q4(qa,qz-1.,r)*q2(qab,-qz,x)
      us65=15.*q4(qab,-qz,x)*q2(qa,qz-1.,r)
      us66=20.*q3(qa,qz-1.,r)*q3(qab,-qz,x)
      us(6)=(us60+us61+us62+us63+us64+us65+us66)*aps*ucd
      
c 7s
      us70=q7(qa,qz-1.,r)*q0(qab,-qz)
      us71=q7(qab,-qz,x)*q0(qa,qz-1.)
      us72=7.*q6(qa,qz-1.,r)*q1(qab,-qz,x)
      us73=7.*q6(qab,-qz,x)*q1(qa,qz-1.,r)
      us74=21.*q5(qa,qz-1.,r)*q2(qab,-qz,x)
      us75=21.*q5(qab,-qz,x)*q2(qa,qz-1.,r)
      us76=35.*q4(qa,qz-1.,r)*q3(qab,-qz,x)
      us77=35.*q4(qab,-qz,x)*q3(qa,qz-1.,r)
      us(7)=-(us70+us71+us72+us73+us74+us75+us76+us77)*aps*ucd
     
c 8s
      us80=q8(qa,qz-1.,r)*q0(qab,-qz)
      us81=q8(qab,-qz,x)*q0(qa,qz-1.)
      us82=8.*q7(qa,qz-1.,r)*q1(qab,-qz,x)
      us83=8.*q7(qab,-qz,x)*q1(qa,qz-1.,r)
      us84=28.*q6(qa,qz-1.,r)*q2(qab,-qz,x)
      us85=28.*q6(qab,-qz,x)*q2(qa,qz-1.,r)
      us86=56.*q5(qa,qz-1.,r)*q3(qab,-qz,x)
      us87=56.*q5(qab,-qz,x)*q3(qa,qz-1.,r)
      us88=70.*q4(qa,qz-1.,r)*q4(qab,-qz,x)
      us(8)=(us80+us81+us82+us83+us84+us85+us86+us87+us88)
     1  *aps*ucd

c 9s
      us90=q9(qa,qz-1.,r)*q0(qab,-qz)
      us91=q9(qab,-qz,x)*q0(qa,qz-1.)
      us92=9.*q8(qa,qz-1.,r)*q1(qab,-qz,x)
      us93=9.*q8(qab,-qz,x)*q1(qa,qz-1.,r)
      us94=36.*q7(qa,qz-1.,r)*q2(qab,-qz,x)
      us95=36.*q7(qab,-qz,x)*q2(qa,qz-1.,r)
      us96=84.*q6(qa,qz-1.,r)*q3(qab,-qz,x)
      us97=84.*q6(qab,-qz,x)*q3(qa,qz-1.,r)
      us98=126.*q5(qa,qz-1.,r)*q4(qab,-qz,x)
      us99=126.*q5(qab,-qz,x)*q4(qa,qz-1.,r)
      us(9)=-(us90+us91+us92+us93+us94+us95+us96+us97+us98+us99)
     1  *aps*ucd

c 10s
      us100=q10(qa,qz-1.,r)*q0(qab,-qz)
      us101=q10(qab,-qz,x)*q0(qa,qz-1.)
      us102=10.*q9(qa,qz-1.,r)*q1(qab,-qz,x)
      us103=10.*q9(qab,-qz,x)*q1(qa,qz-1.,r)
      us104=45.*q8(qa,qz-1.,r)*q2(qab,-qz,x)
      us105=45.*q8(qab,-qz,x)*q2(qa,qz-1.,r)
      us106=120.*q7(qa,qz-1.,r)*q3(qab,-qz,x)
      us107=120.*q7(qab,-qz,x)*q3(qa,qz-1.,r)
      us108=210.*q6(qa,qz-1.,r)*q4(qab,-qz,x)
      us109=210.*q6(qab,-qz,x)*q4(qa,qz-1.,r)
      us1010=252.*q5(qa,qz-1.,r)*q5(qab,-qz,x)
      us(10)=(us100+us101+us102+us103+us104+us105+us106+us107+us108
     1         +us109+us1010)*aps*ucd
c 11s 
      us110=q11(qa,qz-1.,r)*q0(qab,-qz)
      us111=q11(qab,-qz,x)*q0(qa,qz-1.)
      us112=11.*q10(qa,qz-1.,r)*q1(qab,-qz,x)
      us113=11.*q10(qab,-qz,x)*q1(qa,qz-1.,r)
      us114=55.*q9(qa,qz-1.,r)*q2(qab,-qz,x)
      us115=55.*q9(qab,-qz,x)*q2(qa,qz-1.,r)
      us116=165.*q8(qa,qz-1.,r)*q3(qab,-qz,x)
      us117=165.*q8(qab,-qz,x)*q3(qa,qz-1.,r)
      us118=330.*q7(qa,qz-1.,r)*q4(qab,-qz,x)
      us119=330.*q7(qab,-qz,x)*q4(qa,qz-1.,r)
      us1110=462.*q6(qa,qz-1.,r)*q5(qab,-qz,x)
      us1111=462.*q5(qa,qz-1.,r)*q6(qab,-qz,x)
      us(11)=-(us110+us111+us112+us113+us114+us115+us116+us117+us118
     1         +us119+us1110+us1111)*aps*ucd
c 12s 
      us120=q12(qa,qz-1.,r)*q0(qab,-qz)
      us121=q12(qab,-qz,x)*q0(qa,qz-1.)
      us122=12.*q11(qa,qz-1.,r)*q1(qab,-qz,x)
      us123=12.*q11(qab,-qz,x)*q1(qa,qz-1.,r)
      us124=66.*q10(qa,qz-1.,r)*q2(qab,-qz,x)
      us125=66.*q10(qab,-qz,x)*q2(qa,qz-1.,r)
      us126=220.*q9(qa,qz-1.,r)*q3(qab,-qz,x)
      us127=220.*q9(qab,-qz,x)*q3(qa,qz-1.,r)
      us128=495.*q8(qa,qz-1.,r)*q4(qab,-qz,x)
      us129=495.*q8(qab,-qz,x)*q4(qa,qz-1.,r)
      us1210=792.*q7(qa,qz-1.,r)*q5(qab,-qz,x)
      us1211=792.*q5(qa,qz-1.,r)*q7(qab,-qz,x)
      us1212=924.*q6(qa,qz-1.,r)*q6(qab,-qz,x)
      us(12)=(us120+us121+us122+us123+us124+us125+us126+us127+us128
     1         +us129+us1210+us1211+us1212)*aps*ucd

      tcs=us(n)
      end

      function tcp(a,az,adx,ady,adz,akx,aky,akz,n,c)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      dimension up(12)
c NON slater
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
      ucd=1.d0+xi*0.d0
      
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

         
 272  if(n.eq.2) goto 2
      if(n.eq.3) goto 3 
      if(n.eq.4) goto 4 
      if(n.eq.5) goto 5 
      if(n.eq.6) goto 6 
      if(n.eq.7) goto 7 
      if(n.eq.8) goto 8 
      if(n.eq.9) goto 9 
      if(n.eq.10) goto 10     
      if(n.eq.11) goto 11     
      if(n.eq.12) goto 12     


c 2p
  2   up20=q1(qa,qz-2.,r)*q0(qab,-qz-1.)
      up21=q1(qab,-qz-1.,x)*q0(qa,qz-2.)
      up22=(up20+up21)*tp0(qa,qb,qz,va2,vb2)
      up23=q0(qa,qz-2.)*q0(qab,-qz-1.)
      up24=up23*tp1(qz,va2,vb2,r,qz5)
      up(2)=-(up22+up24)*xi*app*ucd
      goto 500

c 3p
  3   up30=q2(qa,qz-2.,r)*q0(qab,-qz-1.)
      up31=q2(qab,-qz-1.,x)*q0(qa,qz-2.)
      up32=2.*q1(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up33=(up30+up31+up32)*tp0(qa,qb,qz,va2,vb2)
      up34=q1(qa,qz-2.,r)*q0(qab,-qz-1.)
      up35=q1(qab,-qz-1.,x)*q0(qa,qz-2.)
      up36=(up34+up35)*2.*tp1(qz,va2,vb2,r,qz5)
      up37=q0(qa,qz-2.)*q0(qab,-qz-1.)*tp2(qz,va2,vb2)
      up(3)=(up33+up36+up37)*app*xi*ucd
      goto 500

c 4p
 4    up40=q3(qa,qz-2.,r)*q0(qab,-qz-1.)
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
      goto 500
c 5p      
  5   up500=q4(qa,qz-2.,r)*q0(qab,-qz-1.)
      up501=q4(qab,-qz-1.,x)*q0(qa,qz-2.)
      up502=4.*q3(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up503=4.*q3(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up504=6.*q2(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up50=(up500+up501+up502+up503+up504)*tp0(qa,qb,qz,va2,vb2)
      up510=q3(qa,qz-2.,r)*q0(qab,-qz-1.)
      up511=q3(qab,-qz-1.,x)*q0(qa,qz-2.)
      up512=3.*q2(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up513=3.*q2(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up51=(up510+up511+up512+up513)*4.*tp1(qz,va2,vb2,r,qz5)
      up520=q2(qa,qz-2.,r)*q0(qab,-qz-1.)
      up521=q2(qab,-qz-1.,x)*q0(qa,qz-2.)
      up522=2.*q1(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up52=(up520+up521+up522)*6.*tp2(qz,va2,vb2)
      up(5)=(up50+up51+up52)*xi*app*ucd
      goto 500
c 6p
  6   up600=q5(qa,qz-2.,r)*q0(qab,-qz-1.)
      up601=q5(qab,-qz-1.,x)*q0(qa,qz-2.)
      up602=5.*q4(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up603=5.*q4(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up604=10.*q3(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up605=10.*q3(qab,-qz-1.,x)*q2(qa,qz-2.,r)
      up60=(up600+up601+up602+up603+up604+up605)*tp0(qa,qb,qz,va2,vb2)
      up610=q4(qa,qz-2.,r)*q0(qab,-qz-1.)
      up611=+q4(qab,-qz-1.,x)*q0(qa,qz-2.)
      up612=4.*q3(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up613=4.*q3(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up614=6.*q2(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up61=(up610+up611+up612+up613+up614)*5.*tp1(qz,va2,vb2,r,qz5)
      up620=q3(qa,qz-2.,r)*q0(qab,-qz-1.)
      up621=q3(qab,-qz-1.,x)*q0(qa,qz-2.)
      up622=3.*q2(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up623=3.*q2(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up62=(up620+up621+up622+up623)*10.*tp2(qz,va2,vb2)
      up(6)=-(up60+up61+up62)*xi*app*ucd
      goto 500
c 7p 
  7   up700=q6(qa,qz-2.,r)*q0(qab,-qz-1.)
      up701=q6(qab,-qz-1.,x)*q0(qa,qz-2.)
      up702=6.*q5(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up703=6.*q5(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up704=15.*q4(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up705=15.*q4(qab,-qz-1.,x)*q2(qa,qz-2.,r)
      up706=20.*q3(qab,-qz-1.,x)*q3(qa,qz-2.,r)
      up70=(up700+up701+up702+up703+up704+up705+up706)
     1     *tp0(qa,qb,qz,va2,vb2)
      up710=6.*q5(qa,qz-2.,r)*q0(qab,-qz-1.)
      up711=6.*q5(qab,-qz-1.,x)*q0(qa,qz-2.)
      up712=30.*q4(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up713=30.*q4(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up714=60.*q3(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up715=60.*q2(qa,qz-2.,r)*q3(qab,-qz-1.,x)
      up71=(up710+up711+up712+up713+up714+up715)
     1    *tp1(qz,va2,vb2,r,qz5)
      up720=15.*q4(qa,qz-2.,r)*q0(qab,-qz-1.)
      up721=15.*q4(qab,-qz-1.,x)*q0(qa,qz-2.)
      up722=60.*q3(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up723=60.*q3(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up724=90.*q2(qab,-qz-1.,x)*q2(qa,qz-2.,r)
      up72=(up720+up721+up722+up723+up724)*tp2(qz,va2,vb2)
      up(7)=(up70+up71+up72)*xi*app*ucd
      goto 500
c 8p 
  8   up800=q7(qa,qz-2.,r)*q0(qab,-qz-1.)
      up801=q7(qab,-qz-1.,x)*q0(qa,qz-2.)
      up802=7.*q6(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up803=7.*q6(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up804=21.*q5(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up805=21.*q5(qab,-qz-1.,x)*q2(qa,qz-2.,r)
      up806=35.*q4(qab,-qz-1.,x)*q3(qa,qz-2.,r)
      up807=35.*q3(qab,-qz-1.,x)*q4(qa,qz-2.,r)
      up80=(up800+up801+up802+up803+up804+up805+up806+up807)
     1     *tp0(qa,qb,qz,va2,vb2)
      up810=7.*q6(qa,qz-2.,r)*q0(qab,-qz-1.)
      up811=7.*q6(qab,-qz-1.,x)*q0(qa,qz-2.)
      up812=42.*q5(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up813=42.*q5(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up814=105.*q4(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up815=105.*q2(qa,qz-2.,r)*q4(qab,-qz-1.,x)
      up816=140.*q3(qa,qz-2.,r)*q3(qab,-qz-1.,x)
      up81=(up810+up811+up812+up813+up814+up815+up816)
     1    *tp1(qz,va2,vb2,r,qz5)
      up820=21.*q5(qa,qz-2.,r)*q0(qab,-qz-1.)
      up821=21.*q5(qab,-qz-1.,x)*q0(qa,qz-2.)
      up822=105.*q4(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up823=105.*q4(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up824=210.*q3(qab,-qz-1.,x)*q2(qa,qz-2.,r)
      up825=210.*q2(qab,-qz-1.,x)*q3(qa,qz-2.,r)
      up82=(up820+up821+up822+up823+up824+up825)*tp2(qz,va2,vb2)
      up(8)=-(up80+up81+up82)*xi*app*ucd
      goto 500
c 9p 
  9   up900=q8(qa,qz-2.,r)*q0(qab,-qz-1.)
      up901=q8(qab,-qz-1.,x)*q0(qa,qz-2.)
      up902=8.*q7(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up903=8.*q7(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up904=28.*q6(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up905=28.*q6(qab,-qz-1.,x)*q2(qa,qz-2.,r)
      up906=56.*q5(qab,-qz-1.,x)*q3(qa,qz-2.,r)
      up907=56.*q3(qab,-qz-1.,x)*q5(qa,qz-2.,r)
      up908=70.*q4(qab,-qz-1.,x)*q4(qa,qz-2.,r)
      up90=(up900+up901+up902+up903+up904+up905+up906+up907+up908)
     1     *tp0(qa,qb,qz,va2,vb2)
      up910=8.*q7(qa,qz-2.,r)*q0(qab,-qz-1.)
      up911=8.*q7(qab,-qz-1.,x)*q0(qa,qz-2.)
      up912=56.*q6(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up913=56.*q6(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up914=168.*q5(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up915=168.*q2(qa,qz-2.,r)*q5(qab,-qz-1.,x)
      up916=280.*q4(qa,qz-2.,r)*q3(qab,-qz-1.,x)
      up917=280.*q3(qa,qz-2.,r)*q4(qab,-qz-1.,x)
      up91=(up910+up911+up912+up913+up914+up915+up916+up917)
     1    *tp1(qz,va2,vb2,r,qz5)
      up920=28.*q6(qa,qz-2.,r)*q0(qab,-qz-1.)
      up921=28.*q6(qab,-qz-1.,x)*q0(qa,qz-2.)
      up922=168.*q5(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up923=168.*q5(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up924=420.*q4(qab,-qz-1.,x)*q2(qa,qz-2.,r)
      up925=420.*q2(qab,-qz-1.,x)*q4(qa,qz-2.,r)
      up926=560.*q3(qab,-qz-1.,x)*q3(qa,qz-2.,r)
      up92=(up920+up921+up922+up923+up924+up925+up926)
     1   *tp2(qz,va2,vb2)
      up(9)=(up90+up91+up92)*xi*app*ucd
      goto 500

c 10p 
  10  up1000=q9(qa,qz-2.,r)*q0(qab,-qz-1.)
      up1001=q9(qab,-qz-1.,x)*q0(qa,qz-2.)
      up1002=9.*q8(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up1003=9.*q8(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up1004=36.*q7(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up1005=36.*q7(qab,-qz-1.,x)*q2(qa,qz-2.,r)
      up1006=84.*q6(qab,-qz-1.,x)*q3(qa,qz-2.,r)
      up1007=84.*q3(qab,-qz-1.,x)*q6(qa,qz-2.,r)
      up1008=126.*q5(qab,-qz-1.,x)*q4(qa,qz-2.,r)
      up1009=126.*q4(qab,-qz-1.,x)*q5(qa,qz-2.,r)
      up100=(up1000+up1001+up1002+up1003+up1004+up1005+up1006+up1007
     1 +up1008+up1009)*tp0(qa,qb,qz,va2,vb2)
      up1010=9.*q8(qa,qz-2.,r)*q0(qab,-qz-1.)
      up1011=9.*q8(qab,-qz-1.,x)*q0(qa,qz-2.)
      up1012=72.*q7(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up1013=72.*q7(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up1014=252.*q6(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up1015=252.*q2(qa,qz-2.,r)*q6(qab,-qz-1.,x)
      up1016=504.*q5(qa,qz-2.,r)*q3(qab,-qz-1.,x)
      up1017=504.*q3(qa,qz-2.,r)*q5(qab,-qz-1.,x)
      up1018=630.*q4(qa,qz-2.,r)*q4(qab,-qz-1.,x)
      up101=(up1010+up1011+up1012+up1013+up1014+up1015+up1016+up1017
     1  +up1018)*tp1(qz,va2,vb2,r,qz5)
      up1020=36.*q7(qa,qz-2.,r)*q0(qab,-qz-1.)
      up1021=36.*q7(qab,-qz-1.,x)*q0(qa,qz-2.)
      up1022=252.*q6(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up1023=252.*q6(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up1024=756.*q5(qab,-qz-1.,x)*q2(qa,qz-2.,r)
      up1025=756.*q2(qab,-qz-1.,x)*q5(qa,qz-2.,r)
      up1026=1260.*q4(qab,-qz-1.,x)*q3(qa,qz-2.,r)
      up1027=1260.*q3(qab,-qz-1.,x)*q4(qa,qz-2.,r)
      up102=(up1020+up1021+up1022+up1023+up1024+up1025+up1026+up1027)
     1   *tp2(qz,va2,vb2)
      up(10)=-(up100+up101+up102)*xi*app*ucd
      goto 500
c 11p 
  11  up1100=q10(qa,qz-2.,r)*q0(qab,-qz-1.)
      up1101=q10(qab,-qz-1.,x)*q0(qa,qz-2.)
      up1102=10.*q9(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up1103=10.*q9(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up1104=45.*q8(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up1105=45.*q8(qab,-qz-1.,x)*q2(qa,qz-2.,r)
      up1106=120.*q7(qab,-qz-1.,x)*q3(qa,qz-2.,r)
      up1107=120.*q3(qab,-qz-1.,x)*q7(qa,qz-2.,r)
      up1108=210.*q6(qab,-qz-1.,x)*q4(qa,qz-2.,r)
      up1109=210.*q4(qab,-qz-1.,x)*q6(qa,qz-2.,r)
      up1110=252.*q5(qab,-qz-1.,x)*q5(qa,qz-2.,r)
      up110=(up1100+up1101+up1102+up1103+up1104+up1105+up1106+up1107
     1 +up1108+up1109+up1110)*tp0(qa,qb,qz,va2,vb2)
      up1120=10.*q9(qa,qz-2.,r)*q0(qab,-qz-1.)
      up1121=10.*q9(qab,-qz-1.,x)*q0(qa,qz-2.)
      up1122=90.*q8(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up1123=90.*q8(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up1124=360.*q7(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up1125=360.*q2(qa,qz-2.,r)*q7(qab,-qz-1.,x)
      up1126=840.*q6(qa,qz-2.,r)*q3(qab,-qz-1.,x)
      up1127=840.*q3(qa,qz-2.,r)*q6(qab,-qz-1.,x)
      up1128=1260.*q5(qa,qz-2.,r)*q4(qab,-qz-1.,x)
      up1129=1260.*q4(qa,qz-2.,r)*q5(qab,-qz-1.,x)
      up111=(up1120+up1121+up1122+up1123+up1124+up1125+up1126+up1127
     1  +up1128+up1129)*tp1(qz,va2,vb2,r,qz5)
      up1140=45.*q8(qa,qz-2.,r)*q0(qab,-qz-1.)
      up1141=45.*q8(qab,-qz-1.,x)*q0(qa,qz-2.)
      up1142=360.*q7(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up1143=360.*q7(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up1144=1260.*q6(qab,-qz-1.,x)*q2(qa,qz-2.,r)
      up1145=1260.*q2(qab,-qz-1.,x)*q6(qa,qz-2.,r)
      up1146=2520.*q5(qab,-qz-1.,x)*q3(qa,qz-2.,r)
      up1147=2520.*q3(qab,-qz-1.,x)*q5(qa,qz-2.,r)
      up1148=3150.*q4(qab,-qz-1.,x)*q4(qa,qz-2.,r)
      up112=(up1140+up1141+up1142+up1143+up1144+up1145+up1146+up1147
     1 +up1148)*tp2(qz,va2,vb2)
      up(11)=(up110+up111+up112)*xi*app*ucd
      goto 500
c 12p  
  12  up1200=q11(qa,qz-2.,r)*q0(qab,-qz-1.)
      up1201=q11(qab,-qz-1.,x)*q0(qa,qz-2.)
      up1202=11.*q10(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up1203=11.*q10(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up1204=55.*q9(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up1205=55.*q9(qab,-qz-1.,x)*q2(qa,qz-2.,r)
      up1206=165.*q8(qab,-qz-1.,x)*q3(qa,qz-2.,r)
      up1207=165.*q3(qab,-qz-1.,x)*q8(qa,qz-2.,r)
      up1208=330.*q7(qab,-qz-1.,x)*q4(qa,qz-2.,r)
      up1209=330.*q4(qab,-qz-1.,x)*q7(qa,qz-2.,r)
      up1210=462.*q6(qab,-qz-1.,x)*q5(qa,qz-2.,r)
      up1211=462.*q5(qab,-qz-1.,x)*q6(qa,qz-2.,r)
      up120=(up1200+up1201+up1202+up1203+up1204+up1205+up1206+up1207
     1 +up1208+up1209+up1210+up1211)*tp0(qa,qb,qz,va2,vb2)
      up1220=11.*q10(qa,qz-2.,r)*q0(qab,-qz-1.)
      up1221=11.*q10(qab,-qz-1.,x)*q0(qa,qz-2.)
      up1222=110.*q9(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up1223=110.*q9(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up1224=495.*q8(qa,qz-2.,r)*q2(qab,-qz-1.,x)
      up1225=495.*q2(qa,qz-2.,r)*q8(qab,-qz-1.,x)
      up1226=1320.*q7(qa,qz-2.,r)*q3(qab,-qz-1.,x)
      up1227=1320.*q3(qa,qz-2.,r)*q7(qab,-qz-1.,x)
      up1228=2310.*q6(qa,qz-2.,r)*q4(qab,-qz-1.,x)
      up1229=2310.*q4(qa,qz-2.,r)*q6(qab,-qz-1.,x)
      up1230=2772.*q5(qa,qz-2.,r)*q5(qab,-qz-1.,x)
      up121=(up1220+up1221+up1222+up1223+up1224+up1225+up1226+up1227
     1  +up1228+up1229+up1230)*tp1(qz,va2,vb2,r,qz5)
      up1240=55.*q9(qa,qz-2.,r)*q0(qab,-qz-1.)
      up1241=55.*q9(qab,-qz-1.,x)*q0(qa,qz-2.)
      up1242=495.*q8(qa,qz-2.,r)*q1(qab,-qz-1.,x)
      up1243=495.*q8(qab,-qz-1.,x)*q1(qa,qz-2.,r)
      up1244=1980.*q7(qab,-qz-1.,x)*q2(qa,qz-2.,r)
      up1245=1980.*q2(qab,-qz-1.,x)*q7(qa,qz-2.,r)
      up1246=4620.*q6(qab,-qz-1.,x)*q3(qa,qz-2.,r)
      up1247=4620.*q3(qab,-qz-1.,x)*q6(qa,qz-2.,r)
      up1248=6930.*q5(qab,-qz-1.,x)*q4(qa,qz-2.,r)
      up1249=6930.*q4(qab,-qz-1.,x)*q5(qa,qz-2.,r)
      up122=(up1240+up1241+up1242+up1243+up1244+up1245+up1246+up1247
     1 +up1248+up1249)*tp2(qz,va2,vb2)
      up(12)=-(up120+up121+up122)*xi*app*ucd

 500  tcp=up(n)
      end

      function q0s(ae,az1,adx,ady,adz,akx,aky,akz)
c tfc de exp(-ae r)/r (pas de coef )
c dr exp(-ae r)/r exp(iu.r) exp(-ikd.r) exp(pi az1/2kd) 1F1(iaz1/kd,1,i(kdr+kd.r)/(2pi)**1.5
c gamma(1-iaz1/kd)
c adx,ady,adz composantes de u;  akx,aky,akz composantes de kd 
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00 
      xi=(0.0d00,1.00d00)
      al1=akx*akx+aky*aky+akz*akz
      ak=dsqrt(al1)
      akk=adx*adx+ady*ady+adz*adz
      aq=akk+al1-2.d00*(adx*akx+ady*aky+adz*akz)
      za=ae*ae+aq
      zab=ae*ae+akk-al1-2.d00*xi*ae*ak
      ya=az1/ak
      sz=1.d00-xi*ya 
      zf0=cdexp((-sz)*ulog(za))
      amod=zab*dconjg(zab)
      aer=.00001d00
      if  (amod.lt.aer) then 
      zg=1.d00
      go to 5120
      endif 
      zg=cdexp((-xi*ya)*ulog(zab))
 5120 z0s=2.d00**.5d00*zf0*zg/(api**.5d00)
c      q0s=uc(sz)*z0s*cdexp(api*ya/2.d00)
      q0s=z0s
      end



      function zfac(l,az,ak)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
	xi=(0.0d00,1.0d00)
	ay=az/ak
	zfac=1.d0+xi*0.d0
	do j=0,l
	zfac=zfac*(j*1.d00-xi*ay)
	enddo
	end


      function cl2(l,az,ak)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      ay=az/ak
      m=2*l+1

	qzk=xi*ay
	zz=l+1.d0-qzk
	ab=uc(zz)*dconjg(uc(zz))
c      ab=1.d0
c	cl2=((2.d0)**l)*dexp(api*ay/2.d0)/fac(m)
	cl2=dsqrt(ab)*dexp(api*ay/2.d0)*(2.d0**l)/(fac(m))
   	end



      function cl(l,az,ak)
cc   avec approximations
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      api=3.141592654d00
      xi=(0.0d00,1.0d00)
      ay=az/ak
      m=2*l+1

	qzk=xi*ay
	zz=l+1.d0-qzk
	if(ay.ge.400.d0) then
	dzfac=zfac(l,az,ak)*dconjg(zfac(l,az,ak))
	amd=dsqrt(2.d0*api/ay)*(1.d0+dexp(-2.d0*api*ay)/2.d0)
      cl=amd*dsqrt(dzfac)*(2.d0**l)/(fac(m))
	goto 410

	else
	ab=uc(zz)*dconjg(uc(zz))
	cl=dsqrt(ab)*dexp(api*ay/2.d0)*(2.d0**l)/(fac(m))
410   endif
	end

******************************************
c calcul des fonctions de Bessel spheriques (salim)
c******************************************
      function asphb(n,bz)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      asphb0 = dsin(bz)/bz
      asphb1 = dsin(bz)/(bz**2.d0)-dcos(bz)/bz
      asphb2 = dsin(bz)*(3.d0/(bz**3.d0)-1.d0/bz)
     $-3.d0*dcos(bz)/(bz**2.d0)

      nif=n-1
      if(nif)10,11,12
 10   asphb = asphb0
      return
 11   asphb = asphb1
      return
 12   if(n.eq.2)then
      asphb = asphb2
      return
      else
      do 20 i=3,n
      asphbr =(2.d0*i-1.d0)*asphb2/bz - asphb1
      asphb1 = asphb2
      asphb2 = asphbr
 20   continue
      asphb = asphbr
      return
      endif
      return
      end
c FIN Bessel

      function asph(n,bz)
c fonction spherique bessel pour bz PETIT
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      dimension av(101)
      api=3.141592654d00
      xi=(0.0d00,1.00d00)
      azf=0.d00
      av0=1.d00/(uc(n+1.5+0.*xi))
      av(1)=(-.25d00*bz*bz)/(uc(n+2.5+0.*xi))

       do 9131 i=2,100
       av(i)=(-.25d00*bz*bz)**i/(fac(i)*uc(n+1.5+i+0.*xi))
       azf=azf+av(i-1)
       am=azf/1000000.d00
       au=av(i)

       if (dabs(au).le.dabs(am)) then
       goto 9132
       endif
 9131  continue
 9132  avf=(azf+av0)*(.5d00*api)**.5d00*((.5d00)**(n+.5))
       asph=(bz**n)*avf
      end

      function asphv(n,bz)
c fonction spherique bessel pour bz, COMPLETE fait appel à asphb et asph
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      dimension a0(46)

      api=3.141592654d00
      xi=(0.0d00,1.00d00)

      if(n.eq.0) then
      if(bz.lt.0.01) then
      asphv=asph(0,bz)
      else
      asphv=asphb(0,bz)
      endif
      goto 81
      endif

      a0(1)=0.02
      a0(2)=0.04
      a0(3)=0.06
      a0(4)=0.12
      a0(5)=0.3
      a0(6)=0.5
      a0(7)=0.85
      a0(8)=1.2
      a0(9)=1.6
      a0(10)=1.9
      a0(11)=2.6
      a0(12)=3.4
      a0(13)=3.7
      a0(14)=4.4
      a0(15)=4.7
      a0(16)=5.2
      a0(17)=6.1
      a0(18)=6.9
      a0(19)=7.8
      a0(20)=8.5
      a0(21)=9.2
      a0(22)=9.6
      a0(23)=10.1
      a0(24)=10.8
      a0(25)=12.1
      a0(26)=12.9
      a0(27)=13.8
      a0(28)=14.3
      a0(29)=14.9
      a0(30)=15.9
      a0(31)=16.8
      a0(32)=17.5
      a0(33)=18.3
      a0(34)=18.3
      a0(35)=20.1
      a0(36)=20.1
      a0(37)=20.4
      a0(38)=21.7
      a0(39)=22.9
      a0(40)=23.8
      a0(41)=23.8
      a0(42)=25.9
      if(bz.lt.a0(n)) then
      asphv=asph(n,bz)
      else
      asphv=asphb(n,bz)
      endif
  81  end








