# make a variant of nh3_opt_claude.for: mk.py name iorb igeom ets_deg Ete_deg Ee Ei Esfix nh [nleg md]
import sys
a=sys.argv; src=open('/home/user/dd/nh3/nh3_opt_claude.for').read()
name,iorb,igeom,ets,ete,Ee,Ei,Esf,nh=a[1:10]
rep=[('      igeom=1\n',f'      igeom={igeom}\n'),('      iorb=1\n',f'      iorb={iorb}\n'),
 ('      Ete=api/2.d0\n',f'      Ete={ete}d0*arad\n'),('      ets=6.d00*arad',f'      ets={ets}d0*arad'),
 ('      Ee=74.d00 ; Ei=500.d00 ; Esfix=0.d00',f'      Ee={Ee}d0 ; Ei={Ei}d0 ; Esfix={Esf}d0'),
 ('      nh=30\n',f'      nh={nh}\n'),('cp(6)=4.d0','cp(6)=dmax1(4.d0,1.02d0*ad+3.d0)')]
if len(a)>10: rep.append(('PARAMETER (nleg=10,nd=4,md=1,',f'PARAMETER (nleg={a[10]},nd=4,md={a[11]},'))
for o,n in rep:
    assert src.count(o)==1,(o,src.count(o)); src=src.replace(o,n)
open(name+'.for','w').write(src)
