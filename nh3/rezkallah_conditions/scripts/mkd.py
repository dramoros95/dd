# DDCS variant: theta_s read from ets.in ; mkd.py name iorb Ee
import sys
name,iorb,Ee=sys.argv[1:4]
src=open('/home/user/dd/nh3/nh3_opt_claude.for').read()
rep=[('      igeom=1\n','      igeom=2\n'),('      iorb=1\n',f'      iorb={iorb}\n'),
 ('      Ete=api/2.d0\n','      Ete=120.d0*arad\n'),
 ('      ets=6.d00*arad',"      open(9,file='ets.in',status='old')\n      read(9,*) etsdg\n      close(9)\n      ets=etsdg*arad"),
 ('      Ee=74.d00 ; Ei=500.d00 ; Esfix=0.d00',f'      Ee={Ee}d0 ; Ei=500.d0 ; Esfix=0.d0'),
 ('      nh=30\n','      nh=4\n'),('cp(6)=4.d0','cp(6)=dmax1(4.d0,1.02d0*ad+3.d0)')]
for o,n in rep:
    assert src.count(o)==1,o; src=src.replace(o,n)
open(name+'.for','w').write(src)
