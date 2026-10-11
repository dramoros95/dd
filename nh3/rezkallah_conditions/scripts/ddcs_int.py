import numpy as np, glob, os
from scipy.integrate import simpson
c=2/27.2114; AU=0.529177e-8**2/27.2114          # a.u. -> cm2/(eV sr2)
IP={1:10.85,2:16.4,3:27.0}; nm={1:'3a1',2:'1e',3:'2a1'}; w={1:1,2:2,3:1}
out={}
for E in [5,20,60,150]:
    tot=0; parts={}
    for o in (1,2,3):
        d=f'r_dd_o{o}_E{E}'
        if not os.path.isdir(d): break
        ki=np.sqrt(500*c); ks=np.sqrt((500-E-IP[o])*c)
        K=np.exp(np.linspace(np.log(ki-ks),np.log(ki+ks),17)); g=[]
        for j in range(1,18):
            f=glob.glob(f'{d}/n{j}/nh3_*_phie_*_{E}ev.dat')
            if not f: break
            a=np.loadtxt(f[0]); 
            if a.ndim<2 or len(a)<5: break
            g.append(2*simpson(a[:,1],x=np.radians(a[:,0])))
        if len(g)<17: break
        dd=simpson(np.array(g)*K*K/(ki*ks),x=np.log(K))*w[o]*AU
        parts[nm[o]]=dd; tot+=dd
    if len(parts)==3: out[E]=(tot,parts)
for E,(t,p) in out.items(): print(E,'%.3e'%t,{k:'%.2e'%v for k,v in p.items()})
np.savetxt('ddcs_model.dat',[[E,t,p['3a1'],p['1e'],p['2a1']] for E,(t,p) in out.items()],fmt='%.4e',
 header='Ee(eV) DDCS_total 3a1 1e(x2) 2a1  [cm2/(eV sr)], Ei=500 eV, theta_e=120 deg, BBK3CW-Z')
