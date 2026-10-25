C .  TEST 3 : partial-wave Z(r) correction (PWDL+JTAB+CCOEF+PWQ3) vs a
C .  brute-force 3D integral of conj(phi_Z(r)-phi_2) R4p Y1m exp(iq.r)
      program t
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      parameter (lmx=30)
      double precision rgr(400),wgr(400),r4g(400),xx(400),ww(400)
      double precision dxt(80),dwt(80),dxf(80),dwf(80),ZKR,R4P
      complex*16 dl(400,0:lmx),jt(0:4000,0:lmx,2)
      complex*16 cc(-1:1,0:lmx,2,-lmx-1:lmx+1),res(-1:1),bf(-1:1)
      api=dacos(-1.d0) ; xi=(0.d0,1.d0)
      call KRWF(1)
      az=2.d0 ; dz=8.d0 ; ndr=8 ; nlr=12 ; lmax=16
      ak=dsqrt(7.d0*2.d0/27.21165d0)
      bkx=ak*dsin(1.1d0) ; bky=0.d0 ; bkz=ak*dcos(1.1d0)
      nr=0
      do id=1,ndr
        call GAULEG((id-1)*dz/ndr,id*dz/ndr,xx,ww,nlr,1d-15)
        do i=1,nlr
          nr=nr+1 ; rgr(nr)=xx(i) ; wgr(nr)=ww(i) ; r4g(nr)=R4P(xx(i))
        enddo
      enddo
      qmax=12.d0 ; nq=2400 ; hq=qmax/nq
      call PWDL(ak,az,nr,rgr,lmax,dl)
      call JTAB(nr,rgr,wgr,r4g,dl,lmax,nq,hq,jt)
      call CCOEF(bkx,bky,bkz,lmax,cc)
      call GAULEG(0.d0,api,dxt,dwt,48,1d-15)
      call GAULEG(0.d0,2.d0*api,dxf,dwf,48,1d-15)
      alp=az/ak
      sq=dexp(api*alp/2.d0)*uc(1.d0+xi*alp)
      do it=1,4
        bqx=0.9d0*dsin(0.8d0*it) ; bqy=0.4d0*dcos(1.3d0*it)
        bqz=1.1d0*dcos(0.5d0*it)
        if (it.eq.4) then
          bqx=2.5d0 ; bqy=-1.0d0 ; bqz=1.5d0
        endif
        res=(0.d0,0.d0)
        call PWQ3(bqx,bqy,bqz,lmax,cc,jt,nq,hq,res)
        bf=(0.d0,0.d0)
        do i=1,nr
          ar=rgr(i) ; alr=ZKR(ar)/ak
          qgr=dexp(api*alr/2.d0)*uc(1.d0+xi*alr)/sq
          do j=1,48
          do k=1,48
            ax=ar*dsin(dxt(j))*dcos(dxf(k))
            ay=ar*dsin(dxt(j))*dsin(dxf(k)) ; bz=ar*dcos(dxt(j))
            akr=bkx*ax+bky*ay+bkz*bz
            zc=-xi*(ak*ar+akr)
            qd=dconjg(cdexp(xi*akr)*(qgr*uf(-xi*alr,(1.d0,0.d0),zc)
     $        -uf(-xi*alp,(1.d0,0.d0),zc)))
            qw=cdexp(xi*(bqx*ax+bqy*ay+bqz*bz))*r4g(i)*ar*ar*qd
     $        *dsin(dxt(j))*wgr(i)*dwt(j)*dwf(k)/dsqrt((2.d0*api)**3)
            do m=-1,1
              bf(m)=bf(m)+qw*ylm(1,dble(m),dcos(dxt(j)),dxf(k))
            enddo
          enddo
          enddo
        enddo
        do m=-1,1
          write(*,'(A,I2,A,I2,2(2E13.5,2X),A,E9.2)') ' q',it,' m=',m
     $     ,res(m),bf(m),' rel.diff ',cdabs(res(m)-bf(m))/cdabs(bf(m))
        enddo
      enddo
      end
