      PROGRAM DRV
      IMPLICIT DOUBLEPRECISION (a-h)
      IMPLICIT COMPLEX*16 (o-z)
      PARAMETER (norbx=40,nstox=400,nshx=4)
      DOUBLEPRECISION ZSCR,rr
      CHARACTER*8 alab(norbx)
      CHARACTER*16 aname
      DOUBLE PRECISION aqsh(nshx),arsh(nshx),aocn(norbx),aoci(norbx)
     $,esz(nstox),csz(nstox)
      INTEGER lorb(norbx),morb(norbx),iono(norbx),ist(norbx),nst(norbx)
     $,nsn(nstox)
      do iorb=1,3
      CALL TARGET(iorb,aname,acen,Eio,nsh,aqsh,arsh,norb,alab
     $,lorb,morb,aocn,aoci,iono,ist,nst,nsn,esz,csz,norbx,nstox,nshx)
      do ir=0,1000
        rr=0.01d0*ir+1d-6
        write(10+iorb,*) rr,ZSCR(rr,2,acen,nsh,aqsh,arsh,norb,aocn,aoci
     $  ,ist,nst,nsn,esz,csz)
      enddo
      enddo
      end
