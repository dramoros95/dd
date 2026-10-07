import numpy as np, mpmath as mp
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss, Legendre
d=np.loadtxt('fort.11'); rz,zz=d[:,0],d[:,1]
Z=lambda r: np.interp(r,rz,zz) if r<10 else 1.0
k=np.sqrt(2*74/27.2114); eta=-1/k
cf=lambda l,x: float(mp.coulombf(l,eta,x)); cg=lambda l,x: float(mp.coulombg(l,eta,x))
def dw(l,R=40.):
    r0=1e-4
    f=lambda r,y:[y[1],(l*(l+1)/r**2-2*Z(r)/r-k*k)*y[0]]
    s=solve_ivp(f,[r0,R],[r0**(l+1),(l+1)*r0**l],rtol=1e-10,atol=1e-14,dense_output=True,max_step=0.02)
    u,up=s.y[0,-1],s.y[1,-1]; x=k*R; h=1e-5
    F,G=cf(l,x),cg(l,x); Fp=k*(cf(l,x+h)-cf(l,x-h))/(2*h); Gp=k*(cg(l,x+h)-cg(l,x-h))/(2*h)
    W=F*Gp-G*Fp; a=(u*Gp-up*G)/W; b=(up*F-u*Fp)/W
    return np.arctan2(b,a), (lambda r,A=np.hypot(a,b): s.sol(r)[0]/A)
sig=lambda l: float(mp.arg(mp.gamma(l+1+1j*eta)))
x,w=leggauss(40)
def cwz(l,r):
    a=Z(r)/k; P=Legendre.basis(l); v=0
    for xi,wi in zip(x,w):
        v+=wi*complex(mp.exp(mp.pi*a/2)*mp.gamma(1-1j*a)*mp.exp(1j*k*r*xi)*mp.hyp1f1(1j*a,1,1j*k*r*(1-xi)))*P(xi)
    return v/2/(1j**l)
print('k=%.3f au. Short-range phase shifts (beyond Coulomb Z=1):'%k)
D={}
for l in range(9):
    D[l]=dw(l); print(' l=%d delta=%6.3f'%(l,D[l][0]))
print('|A_l(r)| and phase: exact DW / 3CWZ variable-charge wave / Coulomb Z=1')
for l in range(4):
    dl,u=D[l]
    for r in [0.3,0.6,1.0,1.5,2.5]:
        ex=np.exp(1j*(sig(l)+dl))*u(r)/(k*r); c1=np.exp(1j*sig(l))*cf(l,k*r)/(k*r); cz=cwz(l,r)
        print(' l=%d r=%.1f  DW %.4f(%5.2f)  3CWZ %.4f(%5.2f)  Z1 %.4f(%5.2f)'%(l,r,abs(ex),np.angle(ex),abs(cz),np.angle(cz),abs(c1),np.angle(c1)))
for t in np.radians([30,60,90,120,150,180]):
    fc=-eta/(2*k*np.sin(t/2)**2)*np.exp(-1j*eta*np.log(np.sin(t/2)**2)+2j*sig(0))
    fs=sum((2*l+1)*np.exp(2j*sig(l))*(np.exp(2j*D[l][0])-1)/(2j*k)*Legendre.basis(l)(np.cos(t)) for l in D)
    print(' theta=%3d  elastic |f|^2 / Rutherford(Z=1) = %.1f'%(round(np.degrees(t)),abs(fc+fs)**2/abs(fc)**2))
