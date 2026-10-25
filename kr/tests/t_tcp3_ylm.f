C .  TEST 1 : TCP3 = tcp of the original code (n=2,3,4 ; m=-1,0,1)
C .  TEST 2 : YALL = ylm of the original code (L <= 12, all M)
      program t
      implicit double precision (a-h)
      implicit complex*16 (o-z)
      complex*16 res(-1:1),y(0:31,-31:31)
      emax=0.d0
      do it=1,20
        bqx=dsin(1.3d0*it) ; bqy=dcos(0.7d0*it) ; bqz=0.5d0*dsin(2.1d0*it)
        bkx=0.7d0*dcos(it*0.9d0) ; bky=0.3d0*dsin(it*1.7d0)
        bkz=0.6d0*dcos(it*2.3d0)
        do n=2,4
          al=1.4d0+0.37d0*n+0.11d0*it
          call TCP3(al,2.d0,bqx,bqy,bqz,bkx,bky,bkz,n,res)
          do m=-1,1
            r0=tcp(al,2.d0,bqx,bqy,bqz,bkx,bky,bkz,n,dble(m))
            emax=max(emax,cdabs(res(m)-r0)/cdabs(r0))
          enddo
        enddo
      enddo
      write(*,'(A,E10.3)') ' TCP3 vs tcp  : max relative difference ',emax
      emax=0.d0
      do it=1,10
        ct=dcos(0.31d0*it) ; ph=0.77d0*it
        call YALL(ct,ph,12,y)
        do l=0,12
        do m=-l,l
          emax=max(emax,cdabs(y(l,m)-ylm(l,dble(m),ct,ph)))
        enddo
        enddo
      enddo
      write(*,'(A,E10.3)') ' YALL vs ylm  : max absolute difference ',emax
      end
