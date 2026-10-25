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
