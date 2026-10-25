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
