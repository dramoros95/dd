cccccccccccccccccccc DOUBLE ionisation de 3P de Ar avec configuration d'interaction ccccccccccccccc

      implicit double precision (a-h)
      implicit complex*16 (o-z)
	parameter (nlg=1,nlg1=1)
      dimension dxp(500),dwp(500),dxtp(200),dwtp(200),dxhp(200)
	1,dwhp(200)
      dimension dxp1(300),dwp1(300),dxp2(300),dwp2(300)
     $,dxp3(300),dwp3(300)
     $,dxp4(300),dwp4(300)
     $,dxp5(300),dwp5(300)	      
      dimension et1(1000),att1(1000),azkr(300),dw1(nlg),dw2(nlg1)
	1,dw3(nlg1)
     1,qxx1(nlg,nlg1,nlg1)
	1,dxr(300),dx2(200),dx3(200)
     1,fr4p(nlg),qxx2(nlg,nlg1,nlg1)
	dimension apa(200),etp(200),hip(200)
	dimension apx(200,200,200),apy(200,200,200),apz(200,200,200)
	dimension ampx(200,200,200),ampy(200,200,200),ampz(200,200,200)
	dimension ampdx(200,200,200),ampdy(200,200,200)
     $,ampdz(200,200,200)
	dimension apk12x(200,200,200),apk12y(200,200,200)
     $,apk12z(200,200,200)
	dimension apdx(200,200,200),apdy(200,200,200)
     $,apdz(200,200,200)
        

ccccccccc calcule de l'integrale	    
      api=dacos(-1.d0)
      di2=dsqrt(2.d0)
c     facteurs de conversion 
      acon=2./27.21165
      arad=api/180.
c      charge vue par l'electron ejecte
      az=2.d0 
	az1=az
	az2=az
c     constantes      
      xi=(0.0d00,1.0d00)
      ak0=0.d0
c
	am0=0.d0
	amm1=-1.d0
	am1=1.d0

ccccccccccccccccccccccccccccccccccccccccccc
	enioni=38.d0
      ak1=7.d00
      ak2=7.d00   
c	aki=5658.4d0   
c	aks=aki-ak1-ak2-enioni
	aks=5500.d0
	aki=aks+ak1+ak2+enioni
ccccccccccccccccccccccccccazzccccccccc
      aki=dsqrt(aki*acon)
      aks=dsqrt(aks*acon)
      ak1=dsqrt(ak1*acon)
      ak2=dsqrt(ak2*acon)
ccccccccccccccccccccccccccccccccccccccccccccccccc
      et2=104.d0*arad
	hi2=0.d0*arad
	
      ets=1.d0*arad
	his=180.d0*arad	

ccccccccccccccccccccccccccccccccccccccccccccccccc
	akix=0.d0
	akiy=0.d0
	akiz=aki
ccccccccccccccccccccccccccccccccccccccccccccccccc
	aksx=aks*dsin(ets)*dcos(his)
	aksy=aks*dsin(ets)*dsin(his)
	aksz=aks*dcos(ets)
ccccccccccccccccccccccccccccccccccccccccccccccccc
	ak2x=ak2*dsin(et2)*dcos(hi2)
 	ak2y=ak2*dsin(et2)*dsin(hi2)
	ak2z=ak2*dcos(et2)

c     moment de transfert      
      akk=aki*aki+aks*aks-2.*aki*aks*dcos(ets)
      ak=dsqrt(akk)
c  composantes du moment de transfert c
      adx=akix-aksx
      ady=akiy-aksy
      adz=akiz-aksz
c       direction du moment de transert   c
      ctantetak=dsqrt(adx*adx+ady*ady)/adz
	atetak=360.d0- datan(ctantetak)/arad
c	print*,ak,atetak		  
c	stop
c  interval d'intergration Pr 

	dpi=0.d0
	dp1=0.7d0
	dp2=1.6d0
	dp3=2.1d0
c	dp4=3.5d0
	dpf=10.d0
      np1=20
      np2=120
	np3=40
	np4=5
c	np5=10
      np=np1+np2+np3+np4
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	a=0.d0
	b=10.d0
	c=api
 	d=2.*api
	xi=(0.00d0,1.00d0)
		
	call gauleg(a,b,dxr,dw1,nlg,0.001d0)
	call gauleg(a,c,dx2,dw2,nlg1,0.001d0)
	call gauleg(a,d,dx3,dw3,nlg1,0.001d0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call GAULEG(dpi,dp1,dxp1,dwp1,np1,2d-16)
	do ip1=1,np1
	dxp(ip1)=dxp1(ip1)
	dwp(ip1)=dwp1(ip1)
	enddo

      call GAULEG(dp1,dp2,dxp2,dwp2,np2,2d-16)
	do ip2=1,np2
	dxp(ip2+np1)=dxp2(ip2)
	dwp(ip2+np1)=dwp2(ip2)
	enddo

      call GAULEG(dp2,dp3,dxp3,dwp3,np3,2d-16)
	do ip3=1,np3
	dxp(ip3+np1+np2)=dxp3(ip3)
	dwp(ip3+np1+np2)=dwp3(ip3)
	enddo 

      call GAULEG(dp3,dpf,dxp4,dwp4,np4,2d-16)
 	do ip4=1,np4
	dxp(ip4+np1+np2+np3)=dxp4(ip4)
	dwp(ip4+np1+np2+np3)=dwp4(ip4)
	enddo

c      call GAULEG(dp4,dpf,dxp5,dwp5,np5,2d-16)
c 	do ip5=1,np5
c	dxp(ip5+np1+np2+np3+np4)=dxp5(ip5)
c	dwp(ip5+np1+np2+np3+np4)=dwp5(ip5)
c	enddo

c  interval d'intergration Pt,Ph 

      dtf=api
	dhf=2.*api
	 
      nt=80
      nh=80
      call GAULEG(0.d0,dtf,dxtp,dwtp,nt,2d-16)
      call GAULEG(0.d0,dhf,dxhp,dwhp,nh,2d-16)
      n1=72
      h1=2.*api/n1
	call bzkr(dxr,nlg,azkr)
	call df4p(nlg,dxr,fr4p)
c	azk1=az1/ak1
c      azk2=az2/ak2	

c      ucd1=uc(1.d00+xi*azk1)*dexp(api*azk1/2.d00)
c	ucd2=uc(1.d00+xi*azk2)*dexp(api*azk2/2.d00)

	alp11=az1/ak1

	sz11=1.d0+xi*alp11
	sq1=dexp(api*alp11/2.d0)*uc(sz11)

	alp22=az2/ak2

	sz22=1.d0+xi*alp22
	sq2=dexp(api*alp22/2.d0)*uc(sz22)

      do 800 i1=1,n1+1
c	read(*,*) nii,nff
c      do 800 i1=nii,nff
      et1(i1)=h1*(i1-1)
	att1(i1)=et1(i1)/arad
	hi1=0.d0*arad

ccccccccccccccccccccccccccccccccccccccccccccccccc
	ak1x=ak1*dsin(et1(i1))*dcos(hi1)
 	ak1y=ak1*dsin(et1(i1))*dsin(hi1)
	ak1z=ak1*dcos(et1(i1))
c    copmosants k12
      ak12x=(ak1x-ak2x)/2.
      ak12y=(ak1y-ak2y)/2.
      ak12z=(ak1z-ak2z)/2.
cccccccccccccccccccccccccccccccccccccccccccccccc
      ak12=dsqrt(ak12x*ak12x+ak12y*ak12y+ak12z*ak12z)
c   constante intervenant dans la section efficace
      an1=aks*ak1*ak2/(aki*akk*akk) 
ccccccccccccccccc  ELEMENT DE MATRICES ccccccccccccccccccccccccccccccccccc
 	alp1=az1/ak1
      if (alp1.ge.-1d-7.and.alp1.le.1d-7) then
      anc1=1.d0
      goto 620
      else
      anc1=2.*api*alp1/(1.-dexp(-2.*api*alp1))
 620  endif
 	alp2=az2/ak2
      if (alp2.ge.-1d-7.and.alp2.le.1d-7) then

      anc2=1.d0
      goto 621
      else
      anc2=2.*api*alp2/(1.-dexp(-2.*api*alp2))
 621  endif
      call qf1(azkr,ak1x,ak1y,ak1z,nlg,nlg1,dxr,dx2
	1,dx3,sq1,qxx1) 

      call qf1(azkr,ak2x,ak2y,ak2z,nlg,nlg1,dxr,dx2
	1,dx3,sq2,qxx2) 
c	print*,qxx1,qxx2
c	stop

cccccccccccccccccccc   integration sur p   ccccccccccccccccccccccccccccc

c 	som=0.d0

      resj1a3=0.d0
      resj1a5=0.d0
      resj1a8=0.d0

      do 108 j1=1,np
      resj2a3=0.d0
      resj2a5=0.d0
      resj2a8=0.d0

      do 208 j2=1,nh
      resj3a3=0.d0
      resj3a5=0.d0
	resj3a8=0.d0
      do 308 j3=1,nt

      etp(j3)=dxtp(j3)
	hip(j2)=dxhp(j2)
      apa(j1)=dxp(j1)
c  vecteur  p
	apx(j1,j2,j3)=apa(j1)*dsin(etp(j3))*dcos(hip(j2))
	apy(j1,j2,j3)=apa(j1)*dsin(etp(j3))*dsin(hip(j2))
	apz(j1,j2,j3)=apa(j1)*dcos(etp(j3))
c    vecteur -p
	ampx(j1,j2,j3)=-apx(j1,j2,j3)
	ampy(j1,j2,j3)=-apy(j1,j2,j3)
	ampz(j1,j2,j3)=-apz(j1,j2,j3)
c    vecteur  k+p
	apdx(j1,j2,j3)=adx+apx(j1,j2,j3)
	apdy(j1,j2,j3)=ady+apy(j1,j2,j3)
	apdz(j1,j2,j3)=adz+apz(j1,j2,j3)
c  vecteur k-p 
	ampdx(j1,j2,j3)=adx-apx(j1,j2,j3)
 	ampdy(j1,j2,j3)=ady-apy(j1,j2,j3)
	ampdz(j1,j2,j3)=adz-apz(j1,j2,j3)
c   vecteur k12-p
	apk12x(j1,j2,j3)=ak12x-apx(j1,j2,j3)
	apk12y(j1,j2,j3)=ak12y-apy(j1,j2,j3)
	apk12z(j1,j2,j3)=ak12z-apz(j1,j2,j3)

ccccccccccc �lecton 1 pour les =/ m  ccccccccccccccccccccccccc
C m=-1
      qk1tm1=qkr4p(az,apdx(j1,j2,j3),apdy(j1,j2,j3),apdz(j1,j2,j3)
     $,ak1x,ak1y,ak1z,amm1)
	1+tdN(fr4p,az,apdx(j1,j2,j3),apdy(j1,j2,j3),apdz(j1,j2,j3),ak1
	1,dw1,dw2,dw3,nlg,nlg1,1,amm1,dxr,dx2,dx3,qxx1)

	qk1stm1=qkr4p(az,apx(j1,j2,j3),apy(j1,j2,j3),apz(j1,j2,j3)
     $,ak1x,ak1y,ak1z,amm1)
 	1+tdN(fr4p,az,apx(j1,j2,j3),apy(j1,j2,j3),apz(j1,j2,j3),ak1
	1,dw1,dw2,dw3,nlg,nlg1,1,amm1,dxr,dx2,dx3,qxx1)

C m=0
      qk1t0=qkr4p(az,apdx(j1,j2,j3),apdy(j1,j2,j3),apdz(j1,j2,j3)
     $,ak1x,ak1y,ak1z,am0)
 	1+tdN(fr4p,az,apdx(j1,j2,j3),apdy(j1,j2,j3),apdz(j1,j2,j3),ak1
	1,dw1,dw2,dw3,nlg,nlg1,1,am0,dxr,dx2,dx3,qxx1)

	qk1st0= qkr4p(az,apx(j1,j2,j3),apy(j1,j2,j3),apz(j1,j2,j3)
     $,ak1x,ak1y,ak1z,am0)
	1+tdN(fr4p,az,apx(j1,j2,j3),apy(j1,j2,j3),apz(j1,j2,j3),ak1
	1,dw1,dw2,dw3,nlg,nlg1,1,am0,dxr,dx2,dx3,qxx1)
C m=1
      qk1t1=qkr4p(az,apdx(j1,j2,j3),apdy(j1,j2,j3),apdz(j1,j2,j3)
     $,ak1x,ak1y,ak1z,am1)
	1+tdN(fr4p,az,apdx(j1,j2,j3),apdy(j1,j2,j3),apdz(j1,j2,j3),ak1
	1,dw1,dw2,dw3,nlg,nlg1,1,am1,dxr,dx2,dx3,qxx1)

	qk1st1=qkr4p(az,apx(j1,j2,j3),apy(j1,j2,j3),apz(j1,j2,j3)
     $ ,ak1x,ak1y,ak1z,am1)
	1+tdN(fr4p,az,apx(j1,j2,j3),apy(j1,j2,j3),apz(j1,j2,j3),ak1
	1,dw1,dw2,dw3,nlg,nlg1,1,am1,dxr,dx2,dx3,qxx1)


ccccccccccc �lecton 2 pour les =/ m  ccccccccccccccccccccccccc
C m=-1
      qk2tm1=qkr4p(az,ampdx(j1,j2,j3),ampdy(j1,j2,j3),ampdz(j1,j2,j3)
     $,ak2x,ak2y,ak2z,amm1)
	1+tdN(fr4p,az,ampdx(j1,j2,j3),ampdy(j1,j2,j3),ampdz(j1,j2,j3),ak2
	1,dw1,dw2,dw3,nlg,nlg1,1,amm1,dxr,dx2,dx3,qxx2)

	qk2stm1=qkr4p(az,ampx(j1,j2,j3),ampy(j1,j2,j3),ampz(j1,j2,j3)
     $,ak2x,ak2y,ak2z,amm1)
 	1+tdN(fr4p,az,ampx(j1,j2,j3),ampy(j1,j2,j3),ampz(j1,j2,j3),ak2
	1,dw1,dw2,dw3,nlg,nlg1,1,amm1,dxr,dx2,dx3,qxx2)

C m=0
      qk2t0=qkr4p(az,ampdx(j1,j2,j3),ampdy(j1,j2,j3),ampdz(j1,j2,j3)
     $,ak2x,ak2y,ak2z,am0)
 	1+tdN(fr4p,az,ampdx(j1,j2,j3),ampdy(j1,j2,j3),ampdz(j1,j2,j3),ak2
	1,dw1,dw2,dw3,nlg,nlg1,1,am0,dxr,dx2,dx3,qxx2)

	qk2st0= qkr4p(az,ampx(j1,j2,j3),ampy(j1,j2,j3),ampz(j1,j2,j3)
     $,ak2x,ak2y,ak2z,am0)
	1+tdN(fr4p,az,ampx(j1,j2,j3),ampy(j1,j2,j3),ampz(j1,j2,j3),ak2
	1,dw1,dw2,dw3,nlg,nlg1,1,am0,dxr,dx2,dx3,qxx2)
C m=1
      qk2t1=qkr4p(az,ampdx(j1,j2,j3),ampdy(j1,j2,j3),ampdz(j1,j2,j3)
     $,ak2x,ak2y,ak2z,am1)
	1+tdN(fr4p,az,ampdx(j1,j2,j3),ampdy(j1,j2,j3),ampdz(j1,j2,j3),ak2
	1,dw1,dw2,dw3,nlg,nlg1,1,am1,dxr,dx2,dx3,qxx2)

	qk2st1=qkr4p(az,ampx(j1,j2,j3),ampy(j1,j2,j3),ampz(j1,j2,j3)
     $,ak2x,ak2y,ak2z,am1)
	1+tdN(fr4p,az,ampx(j1,j2,j3),ampy(j1,j2,j3),ampz(j1,j2,j3),ak2
	1,dw1,dw2,dw3,nlg,nlg1,1,am1,dxr,dx2,dx3,qxx2)

ccccccccccc pour E3=(P-1,P0) ccccccccccccccccccccccccccccccccccc 
	
	tdm10=qk1tm1*qk2st0
      tem10=qk1stm1*qk2t0
	tnm10=qk1stm1*qk2st0

	tm10=tdm10+tem10-2.*tnm10
ccccccccccc pour E4=(P0,P-1) ccccccccccccccccccccccccccccccccccc 
	
	td0m1=qk1t0*qk2stm1
      te0m1=qk1st0*qk2tm1
	tn0m1=qk1st0*qk2stm1

	t0m1=td0m1+te0m1-2.*tn0m1 
ccccccccccc pour E5=(P0,P1) ccccccccccccccccccccccccccccccccccc 
	td01=qk1t0*qk2st1
      te01=qk1st0*qk2t1
	tn01=qk1st0*qk2st1

	t01=td01+te01-2.*tn01 

ccccccccccc pour E6=(P1,P0) ccccccccccccccccccccccccccccccccccc 
	td10=qk1t1*qk2st0
      te10=qk1st1*qk2t0
	tn10=qk1st1*qk2st0

	t10=td10+te10-2.*tn10

ccccccccccc pour E8=(P-1,P1) ccccccccccccccccccccccccccccccccccc 

	tdm11=qk1tm1*qk2st1
 	tem11=qk2t1*qk1stm1
	tnm11=qk1stm1*qk2st1

	tm11=tdm11+tem11-2.*tnm11
cccccccccccccccccc E9=(P1,P-1) cccccccccccccccccccccccccccc
      
	td1m1=qk1t1*qk2stm1
	te1m1=qk2tm1*qk1st1
	tn1m1= qk1st1*qk2stm1

	t1m1=td1m1+te1m1-2.*tn1m1
ccccccccccccccccccccc  la 3eme coulombienne   cccccccccccccccccccccccc
	abeta=0.005d0
      az12=0.5d00
      amz12=-az12
	tk12=tcs(abeta,amz12,apk12x(j1,j2,j3),apk12y(j1,j2,j3)
     $,apk12z(j1,j2,j3),ak12x,ak12y,ak12z,1)
	azz12=amz12/ak12
      if (azz12.ge.-1d-7.and.azz12.le.1d-7) then
      anc3=1.
      goto 616
      else
      anc3=2.*api*azz12/(1.-dexp(-2.*api*azz12))
 616  endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	tta3=(tm10-t0m1)*(1./dsqrt(2.d0))*tk12
	tta5=(t01-t10)*(1./dsqrt(2.d0))*tk12
	tta8=(tm11-t1m1)*(1./dsqrt(2.d0))*tk12

c	tt=2.*tta3+tta8
c     som=som+tt*dsin(etp)*dwtp(j3)*dwhp(j2)
      resj3a3=resj3a3+dwtp(j3)*tta3*dsin(etp(j3))
      resj3a5=resj3a5+dwtp(j3)*tta5*dsin(etp(j3))
      resj3a8=resj3a8+dwtp(j3)*tta8*dsin(etp(j3))
cc	print*,117,resj3
 308  continue

      resj2a3=resj2a3+dwhp(j2)*resj3a3
      resj2a5=resj2a5+dwhp(j2)*resj3a5
      resj2a8=resj2a8+dwhp(j2)*resj3a8


 208  continue

	resj1a3=resj1a3+dwp(j1)*resj2a3*apa(j1)*apa(j1)
	resj1a5=resj1a5+dwp(j1)*resj2a5*apa(j1)*apa(j1)
	resj1a8=resj1a8+dwp(j1)*resj2a8*apa(j1)*apa(j1)

c	som1=som*apa*apa
 108  continue
c	call atrapeze(som1,1000,0.d0,10.d0,asintr) 
	tma3int=resj1a3
	tma5int=resj1a5
	tma8int=resj1a8

	atma3=tma3int*dconjg(tma3int)
	atma5=tma5int*dconjg(tma5int)
	atma8=tma8int*dconjg(tma8int)


      atfi=(atma3+atma5+atma8)*3.
c	atfi=asintr*asintr
	const=16.*api

c      ares1d=ares1d*an1*gamov*gamov*anc1*anc2
      atfi=atfi*an1*const*anc1*anc2*anc3

      write(4,*) att1(i1),atfi
      print*,att1(i1),atfi
 800	 continue

      close (unit=4)
      stop
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine atrapeze(f,n,a,b,asintr)
	implicit double precision(a-h)
 	implicit complex*16(o-z)
	dimension f(5000)
	h=(b-a)/n 
	asintr=0
	do i=1,n
 	atr=(f(i)+f(i+1))*h/2.
	asintr=asintr+atr
	enddo
	return 
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function tcpzr(az,akx,aky,akz,akex,akey,akez,am)
 	implicit double precision (a-h)
      implicit complex*16 (o-z)

	tcpzr=qintz(akx,aky,akz,akex,akey,akez,am)
     $+qkr4p(az,akx,aky,akz,akex,akey,akez,am)
	return
	end

calcule de la fonction radiale Rnlm
      FUNCTION Fr4p(ar)
	implicit double precision (a-h)
	implicit complex*16 (o-z)
	dimension a4p(9),e4p(9)

	data a4p/0.09018d00,0.00779d0,0.003172,-0.08481d00,-0.26625d00,
     $0.0349d00,0.54322d00,0.40037d00,0.16488d00/
	data e4p/16.3363d00,24.8543d00,15.4451d00,9.18326d00,6.4205d00,
     $5.25899d00,3.07835d00,1.9621d00,1.422d00/


  	asomn2=0.d0
 	asomn3=0.d0
 	asomn4=0.d0

      do i=1,2
	f2n=fa(e4p(i),1)*a4p(i)*(ar)*dexp(-ar*e4p(i))
	asomn2=asomn2+f2n
		enddo

      do i=3,5
	f3n=fa(e4p(i),2)*a4p(i)*(ar**2)*dexp(-ar*e4p(i))
	asomn3=asomn3+f3n
     	enddo

      do i=6,9
	f4n=fa(e4p(i),3)*a4p(i)*(ar**3)*dexp(-ar*e4p(i))
 	asomn4=asomn4+f4n



     	enddo

      fr4p=asomn2+asomn3+asomn4
	end


cccccccccc calcul de tfc avec z(r) cccccccccccccccccccccccccc
	
	function qcw(az,ar,axx,ayy,azz,akex,akey,akez)
	implicit double precision (a-h)
	implicit complex*16(o-z) 

	api=dacos(-1.d0)
	xi=(0.00d0,1.00d0)
	ac=1./(dsqrt(2.*api)**3)
	ake=dsqrt(akex**2+akey**2+akez**2)
	
	
	al=-az/ake
	ac2=-api*al/2.d0
	apkr=akex*axx+akey*ayy+akez*azz
	za=xi*al
	zb=1.d0+xi*0.d00
	z=-xi*(ake*ar+apkr) 
	sz=1.-xi*al
	qcw=ac*cdexp(xi*apkr)*dexp(ac2)*uc(sz)*uf(za,zb,z)
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function qintz(akx,aky,akz,akex,akey,akez,am)
 	implicit double precision (a-h)
	implicit complex*16(o-z)
	dimension dxr(100),dxt(100),dxf(100),dwr(100),dwt(100),dwf(100)
     $,ar(100),at(100),afi(100),azw(100)
	api=dacos(-1.d0)
	a=0.d0
	b=10.d0
	c=api
 	d=2.*api
	n=20
	xi=(0.00d0,1.00d0)

	
	call gauleg(a,b,dxr,dwr,n,0.001d0)
	call gauleg(a,c,dxt,dwt,n,0.001d0)
	call gauleg(a,d,dxf,dwf,n,0.001d0)
	
	som=(0.00d0,0.00d0)

	do i=1,n
	do j=1,n
	do k=1,n
	ar(i)=dxr(i)
	at(j)=dxt(j)
	afi(k)=dxf(k)

	axx=ar(i)*dsin(at(j))*dcos(afi(k))
	ayy=ar(i)*dsin(at(j))*dsin(afi(k))
	azz=ar(i)*dcos(at(j))

	apkr=akx*axx+aky*ayy+akz*azz  
 	ak=dsqrt(akx**2+aky**2+akz**2)

	azw(i)=36.d00-(2.d00*az1s(ar(i))+2.d00*az2s(ar(i))
     $+2.d00*az3s(ar(i))+2.d00*az4s(ar(i))+6.d00*az2p(ar(i))
     $+6.d00*az3p(ar(i))+4.d00*az4p(ar(i))+10.d00*az3d(ar(i)))

	q1=qcw(azw(i),ar(i),axx,ayy,azz,akex,akey,akez)
	q2=qcw(2.d0,ar(i),axx,ayy,azz,akex,akey,akez)

	som=som+((q1-q2)*cdexp(xi*apkr)
     $*fr4p(ar(i))*ylm(1,am,dcos(at(j)),afi(k))
     $*ar(i)*ar(i)*dwr(i)*dsin(at(j))*dwt(j)*dwf(k))


	enddo
	enddo
	enddo

	qintz=som
	return 
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bzkr(dxr,nlg,azkr)
      implicit double precision (a-h)		 
      implicit complex*16 (o-z)
	dimension azkr(nlg),dxr(nlg)
      do k=1,nlg
	ax=dxr(k)
      azkr(k)=36.d00-(2.d00*az1s(ax)+2.d00*az2s(ax)+2.d00*az3s(ax)
     $+2.d00*az4s(ax)+6.d00*az2p(ax)+6.d00*az3p(ax)+4.d00*az4p(ax)
     $+10.d00*az3d(ax))

      enddo
      return 
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	SUBROUTINE df4p(nlg,dxr,fr4p)
      implicit double precision (a-h)
      implicit complex*16 (o-z)
	dimension fr4p(nlg),dxr(nlg)
	dimension a4p(9),e4p(9)

	data a4p/0.09018d00,0.00779d0,0.003172,-0.08481d00,-0.26625d00,
     $0.0349d00,0.54322d00,0.40037d00,0.16488d00/
	data e4p/16.3363d00,24.8543d00,15.4451d00,9.18326d00,6.4205d00,
     $5.25899d00,3.07835d00,1.9621d00,1.422d00/

	do ii=1,nlg
	ar=dxr(ii)

  	asomn2=0.d0
 	asomn3=0.d0
 	asomn4=0.d0

      do i=1,2
	f2n=fa(e4p(i),1)*a4p(i)*(ar)*dexp(-ar*e4p(i))
	asomn2=asomn2+f2n
		enddo

      do i=3,5
	f3n=fa(e4p(i),2)*a4p(i)*(ar**2)*dexp(-ar*e4p(i))
	asomn3=asomn3+f3n
     	enddo

      do i=6,9
	f4n=fa(e4p(i),3)*a4p(i)*(ar**3)*dexp(-ar*e4p(i))
 	asomn4=asomn4+f4n
     	enddo

      fr4p(ii)=asomn2+asomn3+asomn4
	enddo
	return 
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function azw(ax)
      implicit double precision (a-h)		 
      implicit complex*16 (o-z)  
      azw=36.d00-(2.d00*az1s(ax)+2.d00*az2s(ax)+2.d00*az3s(ax)
     $+2.d00*az4s(ax)+6.d00*az2p(ax)+6.d00*az3p(ax)+4.d00*az4p(ax)
     $+10.d00*az3d(ax))
      end


	function az1s(ax)
      implicit double precision (a-h)		 
      implicit complex*16 (o-z)  
	       
	dimension e1s(2),a1s(2),e2s(2),a2s(2)
	dimension e3s(3),a3s(3),e4s(3),a4s(3)

      data e1s /34.2726d00,40.40d00/
      data a1s /0.73837d00,0.25206d00/

      data e2s /26.8157d00,16.4472d00/
      data a2s /0.02023d00,-0.00881d00/

      data e3s /13.9631d00,8.89765d00,6.22432d00/
      data a3s /0.00593d00,-0.00379d00,0.00191d00/

      data e4s /3.83098d00,2.56514d00,1.73148d00/
      data a4s /-0.00064d00,0.00044d00,-0.00015d00/


      av1s1s=0.d0
      do k1=1,2
      do k2=1,2
	c11=e1s(k1)+e1s(k2)
	bv=2.d0/(c11**3)-avJn(c11,ax,2)
	bv=bv+ax*avJn(c11,ax,1)
	bv=bv*a1s(k1)*a1s(k2)*fa(e1s(k1),0)*fa(e1s(k2),0)
	av1s1s=av1s1s+bv
	enddo
	enddo

	av2s2s=0.d0
      do j1=1,2
      do j2=1,2
	c22=e2s(j1)+e2s(j2)
	bv=fac(4)/(c22**5)-avJn(c22,ax,4)
	bv=bv+ax*avJn(c22,ax,3)
	bv=bv*a2s(j1)*a2s(j2)*fa(e2s(j1),1)*fa(e2s(j2),1)
	av2s2s=av2s2s+bv
	enddo
	enddo

	av3s3s=0.d0
      do j1=1,3
      do j2=1,3
	c33=e3s(j1)+e3s(j2)
	bv=fac(6)/(c33**7)-avJn(c33,ax,6)
	bv=bv+ax*avJn(c33,ax,5)
	bv=bv*a3s(j1)*a3s(j2)*fa(e3s(j1),2)*fa(e3s(j2),2)
	av3s3s=av3s3s+bv
	enddo
	enddo

	av4s4s=0.d0
      do j1=1,3
      do j2=1,3
	c44=e4s(j1)+e4s(j2)
	bv=fac(8)/(c44**9)-avJn(c44,ax,8)
	bv=bv+ax*avJn(c44,ax,7)
	bv=bv*a4s(j1)*a4s(j2)*fa(e4s(j1),3)*fa(e4s(j2),3)
	av4s4s=av4s4s+bv
	enddo
	enddo

	av1s2s=0.d0
      do i1=1,2
      do i2=1,2
	c12=e1s(i1)+e2s(i2)
	bv=fac(3)/(c12**4)-avJn(c12,ax,3)
	bv=bv+ax*avJn(c12,ax,2)
	bv=bv*a1s(i1)*a2s(i2)*fa(e1s(i1),0)*fa(e2s(i2),1)
	av1s2s=av1s2s+bv
	enddo
	enddo

	av1s3s=0.d0
      do i1=1,2
      do i2=1,3
	c13=e1s(i1)+e3s(i2)
	bv=fac(4)/(c13**5)-avJn(c13,ax,4)
	bv=bv+ax*avJn(c13,ax,3)
	bv=bv*a1s(i1)*a3s(i2)*fa(e1s(i1),0)*fa(e3s(i2),2)
	av1s3s=av1s3s+bv
	enddo
	enddo

	av1s4s=0.d0
      do i1=1,2
      do i2=1,3
	c14=e1s(i1)+e4s(i2)
	bv=fac(5)/(c14**6)-avJn(c14,ax,5)
	bv=bv+ax*avJn(c14,ax,4)
	bv=bv*a1s(i1)*a4s(i2)*fa(e1s(i1),0)*fa(e4s(i2),3)
	av1s4s=av1s4s+bv
	enddo
	enddo

	av2s3s=0.d0
      do i1=1,2
      do i2=1,3
	c23=e2s(i1)+e3s(i2)
	bv=fac(5)/(c23**6)-avJn(c23,ax,5)
	bv=bv+ax*avJn(c23,ax,4)
	bv=bv*a2s(i1)*a3s(i2)*fa(e2s(i1),1)*fa(e3s(i2),2)
	av2s3s=av2s3s+bv
	enddo
	enddo

	av2s4s=0.d0
      do i1=1,2
      do i2=1,3
	c24=e2s(i1)+e4s(i2)
	bv=fac(6)/(c24**7)-avJn(c24,ax,6)
	bv=bv+ax*avJn(c24,ax,5)
	bv=bv*a2s(i1)*a4s(i2)*fa(e2s(i1),1)*fa(e4s(i2),3)
	av2s4s=av2s4s+bv
	enddo
	enddo

	av3s4s=0.d0
      do i1=1,3
      do i2=1,3
	c34=e3s(i1)+e4s(i2)
	bv=fac(7)/(c34**8)-avJn(c34,ax,7)
	bv=bv+ax*avJn(c34,ax,6)
	bv=bv*a3s(i1)*a4s(i2)*fa(e3s(i1),2)*fa(e4s(i2),3)
	av3s4s=av3s4s+bv
	enddo
	enddo


	az1s=av1s1s+av2s2s+av3s3s+av4s4s+2.d0*av1s2s+2.d0*av1s3s
     $+2.d0*av1s4s+2.d0*av2s3s+2.d0*av2s4s+2.d0*av3s4s
	end

	function az2s(ax)
      implicit double precision (a-h)		 
      implicit complex*16 (o-z)       
	dimension e1s(2),a1s(2),e2s(2),a2s(2)
	dimension e3s(3),a3s(3),e4s(3),a4s(3)
      data e1s /34.2726d00,40.40d00/
      data a1s /-0.4169d00,0.0649d00/
      data e2s /26.8157d00,16.4472d00/
      data a2s /-0.21848d00,1.08043d00/
      data e3s /13.9631d00,8.89765d00,6.22432d00/  
	  data a3s /0.18819d00,0.00226d00,0.00214d00/
	        data e4s /3.83098d00,2.56514d00,1.73148d00/
      data a4s /-0.00052d00,0.00028d00,-0.00010d00/


      av1s1s=0.d0
      do k1=1,2
      do k2=1,2
	c11=e1s(k1)+e1s(k2)
	bv=2.d0/(c11**3)-avJn(c11,ax,2)
	bv=bv+ax*avJn(c11,ax,1)
	bv=bv*a1s(k1)*a1s(k2)*fa(e1s(k1),0)*fa(e1s(k2),0)
	av1s1s=av1s1s+bv
	enddo
	enddo

	av2s2s=0.d0
      do j1=1,2
      do j2=1,2
	c22=e2s(j1)+e2s(j2)
	bv=fac(4)/(c22**5)-avJn(c22,ax,4)
	bv=bv+ax*avJn(c22,ax,3)
 	bv=bv*a2s(j1)*a2s(j2)*fa(e2s(j1),1)*fa(e2s(j2),1)
	av2s2s=av2s2s+bv
	enddo
	enddo

	av3s3s=0.d0
      do j1=1,3
      do j2=1,3
	c33=e3s(j1)+e3s(j2)
	bv=fac(6)/(c33**7)-avJn(c33,ax,6)
	bv=bv+ax*avJn(c33,ax,5)
	bv=bv*a3s(j1)*a3s(j2)*fa(e3s(j1),2)*fa(e3s(j2),2)
	av3s3s=av3s3s+bv
	enddo
	enddo

	av4s4s=0.d0
      do j1=1,3
      do j2=1,3
	c44=e4s(j1)+e4s(j2)
	bv=fac(8)/(c44**9)-avJn(c44,ax,8)
	bv=bv+ax*avJn(c44,ax,7)
	bv=bv*a4s(j1)*a4s(j2)*fa(e4s(j1),3)*fa(e4s(j2),3)
	av4s4s=av4s4s+bv
	enddo
	enddo

	av1s2s=0.d0
      do i1=1,2
      do i2=1,2
	c12=e1s(i1)+e2s(i2)
	bv=fac(3)/(c12**4)-avJn(c12,ax,3)
	bv=bv+ax*avJn(c12,ax,2)
	bv=bv*a1s(i1)*a2s(i2)*fa(e1s(i1),0)*fa(e2s(i2),1)
	av1s2s=av1s2s+bv
	enddo
	enddo

	av1s3s=0.d0
      do i1=1,2
      do i2=1,3
	c13=e1s(i1)+e3s(i2)
	bv=fac(4)/(c13**5)-avJn(c13,ax,4)
	bv=bv+ax*avJn(c13,ax,3)
	bv=bv*a1s(i1)*a3s(i2)*fa(e1s(i1),0)*fa(e3s(i2),2)
	av1s3s=av1s3s+bv
	enddo
	enddo

	av1s4s=0.d0
      do i1=1,2
      do i2=1,3
	c14=e1s(i1)+e4s(i2)
	bv=fac(5)/(c14**6)-avJn(c14,ax,5)
	bv=bv+ax*avJn(c14,ax,4)
	bv=bv*a1s(i1)*a4s(i2)*fa(e1s(i1),0)*fa(e4s(i2),3)
	av1s4s=av1s4s+bv
	enddo
	enddo

	av2s3s=0.d0
      do i1=1,2
      do i2=1,3
	c23=e2s(i1)+e3s(i2)
	bv=fac(5)/(c23**6)-avJn(c23,ax,5)
	bv=bv+ax*avJn(c23,ax,4)
	bv=bv*a2s(i1)*a3s(i2)*fa(e2s(i1),1)*fa(e3s(i2),2)
	av2s3s=av2s3s+bv
	enddo
	enddo

	av2s4s=0.d0
      do i1=1,2
      do i2=1,3
	c24=e2s(i1)+e4s(i2)
	bv=fac(6)/(c24**7)-avJn(c24,ax,6)
	bv=bv+ax*avJn(c24,ax,5)
	bv=bv*a2s(i1)*a4s(i2)*fa(e2s(i1),1)*fa(e4s(i2),3)
	av2s4s=av2s4s+bv
	enddo
	enddo

	av3s4s=0.d0
      do i1=1,3
      do i2=1,3
	c34=e3s(i1)+e4s(i2)
	bv=fac(7)/(c34**8)-avJn(c34,ax,7)
	bv=bv+ax*avJn(c34,ax,6)
	bv=bv*a3s(i1)*a4s(i2)*fa(e3s(i1),2)*fa(e4s(i2),3)
	av3s4s=av3s4s+bv
	enddo
	enddo


	az2s=av1s1s+av2s2s+av3s3s+av4s4s+2.d0*av1s2s+2.d0*av1s3s
     $+2.d0*av1s4s+2.d0*av2s3s+2.d0*av2s4s+2.d0*av3s4s
	end


	function az3s(ax)
      implicit double precision (a-h)		 
      implicit complex*16 (o-z)  
	       
	dimension e1s(2),a1s(2),e2s(2),a2s(2)
	dimension e3s(3),a3s(3),e4s(3),a4s(3)

      data e1s /34.2726d00,40.40d00/
      data a1s /-0.1555d00,0.01678d00/

      data e2s /26.8157d00,16.4472d00/
      data a2s /-0.11584d00,0.50977d00/

      data e3s /13.9631d00,8.89765d00,6.22432d00/
      data a3s /0.27066d00,-0.5713d00,-0.68418d00/

      data e4s /3.83098d00,2.56514d00,1.73148d00/
      data a4s /-0.00538d00,-0.00082d00,0.00002d00/


      av1s1s=0.d0
      do k1=1,2
      do k2=1,2
	c11=e1s(k1)+e1s(k2)
	bv=2.d0/(c11**3)-avJn(c11,ax,2)
	bv=bv+ax*avJn(c11,ax,1)
	bv=bv*a1s(k1)*a1s(k2)*fa(e1s(k1),0)*fa(e1s(k2),0)
	av1s1s=av1s1s+bv
	enddo
	enddo

	av2s2s=0.d0
      do j1=1,2
      do j2=1,2
	c22=e2s(j1)+e2s(j2)
	bv=fac(4)/(c22**5)-avJn(c22,ax,4)
	bv=bv+ax*avJn(c22,ax,3)
	bv=bv*a2s(j1)*a2s(j2)*fa(e2s(j1),1)*fa(e2s(j2),1)
	av2s2s=av2s2s+bv
	enddo
	enddo

	av3s3s=0.d0
      do j1=1,3
      do j2=1,3
	c33=e3s(j1)+e3s(j2)
	bv=fac(6)/(c33**7)-avJn(c33,ax,6)
	bv=bv+ax*avJn(c33,ax,5)
	bv=bv*a3s(j1)*a3s(j2)*fa(e3s(j1),2)*fa(e3s(j2),2)
	av3s3s=av3s3s+bv
	enddo
	enddo

	av4s4s=0.d0
      do j1=1,3
      do j2=1,3
	c44=e4s(j1)+e4s(j2)
	bv=fac(8)/(c44**9)-avJn(c44,ax,8)
	bv=bv+ax*avJn(c44,ax,7)
	bv=bv*a4s(j1)*a4s(j2)*fa(e4s(j1),3)*fa(e4s(j2),3)
	av4s4s=av4s4s+bv
	enddo
	enddo

	av1s2s=0.d0
      do i1=1,2
      do i2=1,2
	c12=e1s(i1)+e2s(i2)
	bv=fac(3)/(c12**4)-avJn(c12,ax,3)
	bv=bv+ax*avJn(c12,ax,2)
	bv=bv*a1s(i1)*a2s(i2)*fa(e1s(i1),0)*fa(e2s(i2),1)

	av1s2s=av1s2s+bv
	enddo
	enddo

	av1s3s=0.d0
      do i1=1,2
      do i2=1,3
	c13=e1s(i1)+e3s(i2)
	bv=fac(4)/(c13**5)-avJn(c13,ax,4)
	bv=bv+ax*avJn(c13,ax,3)
	bv=bv*a1s(i1)*a3s(i2)*fa(e1s(i1),0)*fa(e3s(i2),2)
	av1s3s=av1s3s+bv
	enddo
	enddo

	av1s4s=0.d0
      do i1=1,2
      do i2=1,3
	c14=e1s(i1)+e4s(i2)
	bv=fac(5)/(c14**6)-avJn(c14,ax,5)
	bv=bv+ax*avJn(c14,ax,4)
	bv=bv*a1s(i1)*a4s(i2)*fa(e1s(i1),0)*fa(e4s(i2),3)
	av1s4s=av1s4s+bv
	enddo
	enddo

	av2s3s=0.d0
      do i1=1,2
      do i2=1,3
	c23=e2s(i1)+e3s(i2)
	bv=fac(5)/(c23**6)-avJn(c23,ax,5)
	bv=bv+ax*avJn(c23,ax,4)
	bv=bv*a2s(i1)*a3s(i2)*fa(e2s(i1),1)*fa(e3s(i2),2)
	av2s3s=av2s3s+bv
	enddo
	enddo

	av2s4s=0.d0
      do i1=1,2
      do i2=1,3
	c24=e2s(i1)+e4s(i2)
	bv=fac(6)/(c24**7)-avJn(c24,ax,6)
	bv=bv+ax*avJn(c24,ax,5)
	bv=bv*a2s(i1)*a4s(i2)*fa(e2s(i1),1)*fa(e4s(i2),3)
	av2s4s=av2s4s+bv
	enddo
	enddo

	av3s4s=0.d0
      do i1=1,3
      do i2=1,3
	c34=e3s(i1)+e4s(i2)
	bv=fac(7)/(c34**8)-avJn(c34,ax,7)
	bv=bv+ax*avJn(c34,ax,6)
	bv=bv*a3s(i1)*a4s(i2)*fa(e3s(i1),2)*fa(e4s(i2),3)
	av3s4s=av3s4s+bv
	enddo
	enddo


	az3s=av1s1s+av2s2s+av3s3s+av4s4s+2.d0*av1s2s+2.d0*av1s3s
     $+2.d0*av1s4s+2.d0*av2s3s+2.d0*av2s4s+2.d0*av3s4s
	end



	function az4s(ax)
      implicit double precision (a-h)		 
      implicit complex*16 (o-z)  
	       
	dimension e1s(2),a1s(2),e2s(2),a2s(2)
	dimension e3s(3),a3s(3),e4s(3),a4s(3)

      data e1s /34.2726d00,40.40d00/
      data a1s /-0.04719d00,0.0035d00/

      data e2s /26.8157d00,16.4472d00/
      data a2s /-0.03974d00,0.16816d00/

      data e3s /13.9631d00,8.89765d00,6.22432d00/
      data a3s /0.09647d00,-0.20412d00,-0.30194d00/

      data e4s /3.83098d00,2.56514d00,1.73148d00/
      data a4s /0.47727d00,0.55546d00,0.12492d00/

      av1s1s=0.d0
      do k1=1,2
      do k2=1,2
	c11=e1s(k1)+e1s(k2)
	bv=2.d0/(c11**3)-avJn(c11,ax,2)
	bv=bv+ax*avJn(c11,ax,1)
	bv=bv*a1s(k1)*a1s(k2)*fa(e1s(k1),0)*fa(e1s(k2),0)
	av1s1s=av1s1s+bv
	enddo
	enddo

	av2s2s=0.d0
      do j1=1,2
      do j2=1,2
	c22=e2s(j1)+e2s(j2)
	bv=fac(4)/(c22**5)-avJn(c22,ax,4)
	bv=bv+ax*avJn(c22,ax,3)
	bv=bv*a2s(j1)*a2s(j2)*fa(e2s(j1),1)*fa(e2s(j2),1)
	av2s2s=av2s2s+bv
	enddo
	enddo

	av3s3s=0.d0
      do j1=1,3
      do j2=1,3
	c33=e3s(j1)+e3s(j2)
	bv=fac(6)/(c33**7)-avJn(c33,ax,6)
	bv=bv+ax*avJn(c33,ax,5)
	bv=bv*a3s(j1)*a3s(j2)*fa(e3s(j1),2)*fa(e3s(j2),2)
	av3s3s=av3s3s+bv
	enddo
	enddo

	av4s4s=0.d0
      do j1=1,3
      do j2=1,3
	c44=e4s(j1)+e4s(j2)
	bv=fac(8)/(c44**9)-avJn(c44,ax,8)
	bv=bv+ax*avJn(c44,ax,7)
	bv=bv*a4s(j1)*a4s(j2)*fa(e4s(j1),3)*fa(e4s(j2),3)
	av4s4s=av4s4s+bv
	enddo
	enddo

	av1s2s=0.d0
      do i1=1,2
      do i2=1,2
	c12=e1s(i1)+e2s(i2)
	bv=fac(3)/(c12**4)-avJn(c12,ax,3)
	bv=bv+ax*avJn(c12,ax,2)
	bv=bv*a1s(i1)*a2s(i2)*fa(e1s(i1),0)*fa(e2s(i2),1)
	av1s2s=av1s2s+bv
	enddo
	enddo

	av1s3s=0.d0
      do i1=1,2
      do i2=1,3
	c13=e1s(i1)+e3s(i2)
	bv=fac(4)/(c13**5)-avJn(c13,ax,4)
	bv=bv+ax*avJn(c13,ax,3)
	bv=bv*a1s(i1)*a3s(i2)*fa(e1s(i1),0)*fa(e3s(i2),2)
	av1s3s=av1s3s+bv
	enddo
	enddo

	av1s4s=0.d0
      do i1=1,2
      do i2=1,3
	c14=e1s(i1)+e4s(i2)
	bv=fac(5)/(c14**6)-avJn(c14,ax,5)
	bv=bv+ax*avJn(c14,ax,4)
	bv=bv*a1s(i1)*a4s(i2)*fa(e1s(i1),0)*fa(e4s(i2),3)
	av1s4s=av1s4s+bv
	enddo
	enddo

	av2s3s=0.d0
      do i1=1,2
      do i2=1,3
	c23=e2s(i1)+e3s(i2)
	bv=fac(5)/(c23**6)-avJn(c23,ax,5)
	bv=bv+ax*avJn(c23,ax,4)
	bv=bv*a2s(i1)*a3s(i2)*fa(e2s(i1),1)*fa(e3s(i2),2)
	av2s3s=av2s3s+bv
	enddo
	enddo

	av2s4s=0.d0
      do i1=1,2
      do i2=1,3
	c24=e2s(i1)+e4s(i2)
	bv=fac(6)/(c24**7)-avJn(c24,ax,6)
	bv=bv+ax*avJn(c24,ax,5)
	bv=bv*a2s(i1)*a4s(i2)*fa(e2s(i1),1)*fa(e4s(i2),3)
	av2s4s=av2s4s+bv
	enddo
	enddo

	av3s4s=0.d0
      do i1=1,3
      do i2=1,3
	c34=e3s(i1)+e4s(i2)
	bv=fac(7)/(c34**8)-avJn(c34,ax,7)
	bv=bv+ax*avJn(c34,ax,6)
	bv=bv*a3s(i1)*a4s(i2)*fa(e3s(i1),2)*fa(e4s(i2),3)
	av3s4s=av3s4s+bv
	enddo
	enddo


	az4s=av1s1s+av2s2s+av3s3s+av4s4s+2.d0*av1s2s+2.d0*av1s3s
     $+2.d0*av1s4s+2.d0*av2s3s+2.d0*av2s4s+2.d0*av3s4s
	end


	function az2p(ax)
      implicit double precision (a-h)		 
      implicit complex*16 (o-z)  
	       
	dimension e2p(2),a2p(2)
	dimension e3p(3),a3p(3),e4p(4),a4p(4)

      data e2p /16.3363d00,24.8543d00/
      data a2p /0.74427d00,0.09648d00/

      data e3p /15.4451d00,9.18326d00,6.4205d00/
      data a3p /0.16587d00,0.04848d00,-0.0215d00/

      data e4p /5.25899d00,3.07835d00,1.9621d00,1.422d00/
      data a4p /0.00775d00,-0.00253d00,0.00152d00,-0.0006d00/




	av2p2p=0.d0
      do j1=1,2
      do j2=1,2
	c22=e2p(j1)+e2p(j2)
	bv=fac(4)/(c22**5)-avJn(c22,ax,4)
	bv=bv+ax*avJn(c22,ax,3)
	bv=bv*a2p(j1)*a2p(j2)*fa(e2p(j1),1)*fa(e2p(j2),1)
	av2p2p=av2p2p+bv
	enddo
	enddo

	av3p3p=0.d0
      do j1=1,3
      do j2=1,3
	c33=e3p(j1)+e3p(j2)
	bv=fac(6)/(c33**7)-avJn(c33,ax,6)
	bv=bv+ax*avJn(c33,ax,5)
	bv=bv*a3p(j1)*a3p(j2)*fa(e3p(j1),2)*fa(e3p(j2),2)
	av3p3p=av3p3p+bv
	enddo
	enddo

	av4p4p=0.d0
      do j1=1,4
      do j2=1,4
	c44=e4p(j1)+e4p(j2)
	bv=fac(8)/(c44**9)-avJn(c44,ax,8)
	bv=bv+ax*avJn(c44,ax,7)
	bv=bv*a4p(j1)*a4p(j2)*fa(e4p(j1),3)*fa(e4p(j2),3)
	av4p4p=av4p4p+bv
	enddo
	enddo

	av2p3p=0.d0
      do i1=1,2
      do i2=1,3
	c23=e2p(i1)+e3p(i2)
	bv=fac(5)/(c23**6)-avJn(c23,ax,5)
	bv=bv+ax*avJn(c23,ax,4)
	bv=bv*a2p(i1)*a3p(i2)*fa(e2p(i1),1)*fa(e3p(i2),2)
	av2p3p=av2p3p+bv
	enddo
	enddo

	av2p4p=0.d0
      do i1=1,2
      do i2=1,4
	c24=e2p(i1)+e4p(i2)
	bv=fac(6)/(c24**7)-avJn(c24,ax,6)
	bv=bv+ax*avJn(c24,ax,5)
	bv=bv*a2p(i1)*a4p(i2)*fa(e2p(i1),1)*fa(e4p(i2),3)
	av2s4s=av2s4s+bv
	enddo
	enddo

	av3p4p=0.d0
      do i1=1,3
      do i2=1,4
	c34=e3p(i1)+e4p(i2)
	bv=fac(7)/(c34**8)-avJn(c34,ax,7)
	bv=bv+ax*avJn(c34,ax,6)
	bv=bv*a3p(i1)*a4p(i2)*fa(e3p(i1),2)*fa(e4p(i2),3)
	av3p4p=av3p4p+bv
	enddo
	enddo


	az2p=av2p2p+av3p3p+av4p4p+2.d0*av2p3p+2.d0*av2p4p+2.d0*av3p4p
	end

	function az3p(ax)
      implicit double precision (a-h)		 
      implicit complex*16 (o-z)  
	       
	dimension e2p(2),a2p(2)
	dimension e3p(3),a3p(3),e4p(4),a4p(4)

      data e2p /16.3363d00,24.8543d00/
      data a2p /-0.316d00,-0.03474d00/

      data e3p /15.4451d00,9.18326d00,6.4205d00/
      data a3p /-0.12839d00,0.3318d00,0.73159d00/

      data e4p /5.25899d00,3.07835d00,1.9621d00,1.422d00/
      data a4p /0.0878d00,0.00399d00,-0.00065d00,0.00026d00/

	av2p2p=0.d0
      do j1=1,2
      do j2=1,2
	c22=e2p(j1)+e2p(j2)
	bv=fac(4)/(c22**5)-avJn(c22,ax,4)
	bv=bv+ax*avJn(c22,ax,3)
	bv=bv*a2p(j1)*a2p(j2)*fa(e2p(j1),1)*fa(e2p(j2),1)
	av2p2p=av2p2p+bv
	enddo
	enddo

	av3p3p=0.d0
      do j1=1,3
      do j2=1,3
	c33=e3p(j1)+e3p(j2)
	bv=fac(6)/(c33**7)-avJn(c33,ax,6)
	bv=bv+ax*avJn(c33,ax,5)
	bv=bv*a3p(j1)*a3p(j2)*fa(e3p(j1),2)*fa(e3p(j2),2)
	av3p3p=av3p3p+bv
	enddo
	enddo

	av4p4p=0.d0
      do j1=1,4
      do j2=1,4
	c44=e4p(j1)+e4p(j2)
	bv=fac(8)/(c44**9)-avJn(c44,ax,8)
	bv=bv+ax*avJn(c44,ax,7)
	bv=bv*a4p(j1)*a4p(j2)*fa(e4p(j1),3)*fa(e4p(j2),3)
	av4p4p=av4p4p+bv
	enddo
	enddo

	av2p3p=0.d0
      do i1=1,2
      do i2=1,3
	c23=e2p(i1)+e3p(i2)
	bv=fac(5)/(c23**6)-avJn(c23,ax,5)
	bv=bv+ax*avJn(c23,ax,4)
	bv=bv*a2p(i1)*a3p(i2)*fa(e2p(i1),1)*fa(e3p(i2),2)
	av2p3p=av2p3p+bv
	enddo
	enddo

	av2p4p=0.d0
      do i1=1,2
      do i2=1,4
	c24=e2p(i1)+e4p(i2)
	bv=fac(6)/(c24**7)-avJn(c24,ax,6)
	bv=bv+ax*avJn(c24,ax,5)
	bv=bv*a2p(i1)*a4p(i2)*fa(e2p(i1),1)*fa(e4p(i2),3)
	av2s4s=av2s4s+bv
	enddo
	enddo

	av3p4p=0.d0
      do i1=1,3
      do i2=1,4
	c34=e3p(i1)+e4p(i2)
	bv=fac(7)/(c34**8)-avJn(c34,ax,7)
	bv=bv+ax*avJn(c34,ax,6)
	bv=bv*a3p(i1)*a4p(i2)*fa(e3p(i1),2)*fa(e4p(i2),3)
	av3p4p=av3p4p+bv
	enddo
	enddo


	az3p=av2p2p+av3p3p+av4p4p+2.d0*av2p3p+2.d0*av2p4p+2.d0*av3p4p
	end

	function az4p(ax)
      implicit double precision (a-h)		 
      implicit complex*16 (o-z)  
	       
	dimension e2p(2),a2p(2)
	dimension e3p(3),a3p(3),e4p(4),a4p(4)

      data e2p /16.3363d00,24.8543d00/
      data a2p /0.09018d00,0.00779d00/

      data e3p /15.4451d00,9.18326d00,6.4205d00/
      data a3p /0.03172d00,-0.08481d00,-0.26625d00/

      data e4p /5.25899d00,3.07835d00,1.9621d00,1.422d00/
      data a4p /0.0349d00,0.54322d00,0.40037d00,0.16488d00/

	av2p2p=0.d0
      do j1=1,2
      do j2=1,2
	c22=e2p(j1)+e2p(j2)
	bv=fac(4)/(c22**5)-avJn(c22,ax,4)
	bv=bv+ax*avJn(c22,ax,3)
	bv=bv*a2p(j1)*a2p(j2)*fa(e2p(j1),1)*fa(e2p(j2),1)
	av2p2p=av2p2p+bv
	enddo
	enddo

	av3p3p=0.d0
      do j1=1,3
      do j2=1,3
	c33=e3p(j1)+e3p(j2)
	bv=fac(6)/(c33**7)-avJn(c33,ax,6)
	bv=bv+ax*avJn(c33,ax,5)
	bv=bv*a3p(j1)*a3p(j2)*fa(e3p(j1),2)*fa(e3p(j2),2)
	av3p3p=av3p3p+bv
	enddo
	enddo

	av4p4p=0.d0
      do j1=1,4
      do j2=1,4
	c44=e4p(j1)+e4p(j2)
	bv=fac(8)/(c44**9)-avJn(c44,ax,8)
	bv=bv+ax*avJn(c44,ax,7)
	bv=bv*a4p(j1)*a4p(j2)*fa(e4p(j1),3)*fa(e4p(j2),3)
	av4p4p=av4p4p+bv
	enddo
	enddo

	av2p3p=0.d0
      do i1=1,2
      do i2=1,3
	c23=e2p(i1)+e3p(i2)
	bv=fac(5)/(c23**6)-avJn(c23,ax,5)
	bv=bv+ax*avJn(c23,ax,4)
	bv=bv*a2p(i1)*a3p(i2)*fa(e2p(i1),1)*fa(e3p(i2),2)
	av2p3p=av2p3p+bv
	enddo
	enddo

	av2p4p=0.d0
      do i1=1,2
      do i2=1,4
	c24=e2p(i1)+e4p(i2)
	bv=fac(6)/(c24**7)-avJn(c24,ax,6)
	bv=bv+ax*avJn(c24,ax,5)
	bv=bv*a2p(i1)*a4p(i2)*fa(e2p(i1),1)*fa(e4p(i2),3)
	av2s4s=av2s4s+bv
	enddo
	enddo

	av3p4p=0.d0
      do i1=1,3
      do i2=1,4
	c34=e3p(i1)+e4p(i2)
	bv=fac(7)/(c34**8)-avJn(c34,ax,7)
	bv=bv+ax*avJn(c34,ax,6)
	bv=bv*a3p(i1)*a4p(i2)*fa(e3p(i1),2)*fa(e4p(i2),3)
	av3p4p=av3p4p+bv
	enddo
	enddo


	az4p=av2p2p+av3p3p+av4p4p+2.d0*av2p3p+2.d0*av2p4p+2.d0*av3p4p
	end


	function az3d(ax)
      implicit double precision (a-h)		 
      implicit complex*16 (o-z)  
	       
	dimension e3d(5),a3d(5)

      data e3d /5.19134d00,3.34983d00,7.68953d00,10.3162d00,16.7379d00/
      data a3d /0.47749d00,0.10325d00,0.28079d00,0.20824d00,0.03159d00/


	av3d3d=0.d0
      do j1=1,5
      do j2=1,5
	c33=e3d(j1)+e3d(j2)
	bv=fac(6)/(c33**7)-avJn(c33,ax,6)
	bv=bv+ax*avJn(c33,ax,5)
	bv=bv*a3d(j1)*a3d(j2)*fa(e3d(j1),2)*fa(e3d(j2),2)
	av3d3d=av3d3d+bv
	enddo
	enddo

	az3d=av3d3d
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function qkr4p(az,akx,aky,akz,akex,akey,akez,am)
 	implicit double precision (a-h)
      implicit complex*16 (o-z)
	dimension al(9),ap(9)
	data ap/0.09018d00,0.00779d0,0.003172,-0.08481d00,-0.26625d00,
     $0.0349d00,0.54322d00,0.40037d00,0.16488d00/
	data al/16.3363d00,24.8543d00,15.4451d00,9.18326d00,6.4205d00,
     $5.25899d00,3.07835d00,1.9621d00,1.422d00/
	s1=(0.0d00,0.0d00)
	s2=(0.0d00,0.0d00)
	s3=(0.0d00,0.0d00)
	do i=1,2
	ta2p=ap(i)
     $*tcp(al(i),az,akx,aky,akz,akex,akey,akez,2,am)
 	 s1=s1+ta2p
	enddo
	do j=3,5
	ta3p=ap(j)
     $*tcp(al(j),az,akx,aky,akz,akex,akey,akez,3,am)
 	s2=s2+ta3p
	enddo 
	do k=6,9
	ta4p=ap(k)
     $*tcp(al(k),az,akx,aky,akz,akex,akey,akez,4,am)
 	s3=s3+ta4p
	enddo

	qkr4p=s1+s2+s3
	return
	end
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

	function fn(al,ax0,n)
	implicit double precision (a-h)
	implicit complex*16(o-z)
	dimension f(1000)
 	asom=0.d0
	do i=1,n
	f(i)=((al*ax0)**(i))/fac(i)
	asom=asom+f(i)
	enddo
	fn=asom*dexp(-al*ax0)*fac(n)/(al**(n+1))
	return
	end
cc fonction factorielle
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
	subroutine qf1(azr,akex,akey,akez,nlg,nlg1,ax,atheta
	1,aphi,sq,qxx) 

      implicit double precision (a-h)
      implicit complex*16 (o-z)
	dimension ax(nlg),azr(nlg)
	1,alpr(nlg),zar(nlg),szr(nlg)
	1,atheta(nlg1),aphi(nlg1),axx(nlg,nlg1,nlg1)
	1,axy(nlg,nlg1,nlg1),axz(nlg,nlg1),akerv(nlg,nlg1,nlg1)
     1,akerm(nlg),zc(nlg,nlg1,nlg1)
     1,qxx(nlg,nlg1,nlg1)

	xi=(0.d00,1.d00)
	api=dacos(-1.d0)

	ake=dsqrt(akex**2+akey**2+akez**2)

	do i=1,nlg
	do j=1,nlg1
	do k=1,nlg1
	axx(i,j,k)=ax(i)*dsin(atheta(j))*dcos(aphi(k))
	axy(i,j,k)=ax(i)*dsin(atheta(j))*dsin(aphi(k))
	axz(i,j)=ax(i)*dcos(atheta(j))

	akerm(i)=ake*ax(i)
	akerv(i,j,k)=akex*axx(i,j,k)+akey*axy(i,j,k)+akez*axz(i,j)
	qkie=xi*akerv(i,j,k)

	alpr(i)=-azr(i)/ake
	zar(i)=xi*alpr(i)
      zb=1.d0+xi*0.d0
	zc(i,j,k)=-xi*(akerm(i)+akerv(i,j,k))

	szr(i)=1.d0-xi*alpr(i)

      az1=2.d00
      alp1=-az1/ake
      za1=xi*alp1



      qx1=cdexp(qkie)*dexp(-api*alpr(i)/2.d0)
	1*uf(zar(i),zb,zc(i,j,k))*uc(szr(i))/sq
     	
	qx2=cdexp(qkie)*uf(za1,zb,zc(i,j,k))


      qxx(i,j,k)=qx1-qx2     
      enddo
	enddo
	enddo
	return
	end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ccccccc Cette fonction permet de calculer l'integrale(0--r0)de "(Q(Z(r))/exp(iK.r)/Qnlm)-(Q(z=1)/exp(iK.r)/Qnlm)
      function tdN(fx,az,akx,aky,akz,ake,dw1,dw2,dw3,nlg,nlg1
	1,l,am,ax,atheta,aphi,qxfr)
      implicit double precision (a-h)		 
      implicit complex*16 (o-z)
	parameter (imax=100)    
      dimension dw1(nlg),dw2(nlg1),dw3(nlg1),ax(nlg),atheta(nlg1)
	1,aphi(nlg1),qki(nlg,nlg1,nlg1),qxfr(nlg,nlg1,nlg1)
     1,fx(nlg)

	api=dacos(-1.d0)
	xi=(0.d00,1.d00)


ccccccc on integre /theta et phi ccccccccc
	alp1=az/ake
c	za=xi*alp1
c	zb=1.d0+xi*0.d0
c	sz=1.d0-xi*alp1


	qss=(0.d0,0.d0)
	do i=1,nlg
	do j=1,nlg1
	do k=1,nlg1

c	ax(i)=dxr(i)

	afk=aphi(k)
	atk=dcos(atheta(j))

ccccccc on integre /r ccccccccc
ccccccc on divise l'intervale de r cccccccccc
	ss=(0.d0,0.d0)


	axx=ax(i)*dsin(atheta(j))*dcos(aphi(k))
	axy=ax(i)*dsin(atheta(j))*dsin(aphi(k))
	axz=ax(i)*dcos(atheta(j))


	akrv=akx*axx+aky*axy+akz*axz
	qki(i,j,k)=xi*akrv

      ylmm=ylm(l,am,atk,afk)

c	qxf1=uf(za,zb,zc(i,j,k))

	qss=qss+cdexp(qki(i,j,k))*fx(i)*(ax(i)**2)*ylmm
	1*dsin(atheta(j))*dconjg(qxfr(i,j,k))*dw1(i)*dw2(j)*dw3(k)
 

	enddo
	enddo
	enddo
	tdN=qss/dsqrt((2.d0*api)**3)
	return
	end
