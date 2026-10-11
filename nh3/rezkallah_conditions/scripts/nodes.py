# theta_s nodes (deg), equally spaced in ln K : nodes.py Ee IP
import sys,numpy as np
Ee,IP=map(float,sys.argv[1:3]); c=2/27.2114
ki=np.sqrt(500*c); ks=np.sqrt((500-Ee-IP)*c)
K=np.exp(np.linspace(np.log(ki-ks),np.log(ki+ks),17))
ct=np.clip((ki*ki+ks*ks-K*K)/(2*ki*ks),-1,1); th=np.degrees(np.arccos(ct)); th[0]=max(th[0],1e-3); th[-1]=min(th[-1],179.999)
print(' '.join('%.6f'%t for t in th))
