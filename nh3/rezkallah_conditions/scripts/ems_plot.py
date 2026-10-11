import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, sys
D='../rezdata/'; c=2/27.2114
IP={'3a1':10.85,'1e':16.4,'2a1':27.0}; run={'3a1':'ems1','1e':'ems2','2a1':'ems3'}
fine=len(sys.argv)>1
fig,axs=plt.subplots(1,3,figsize=(11,3.9))
for ax,o in zip(axs,['1e','2a1','3a1']):
    ke=np.sqrt(600*c); ki=np.sqrt((1200+IP[o])*c); th=np.radians(45)
    m=np.loadtxt(f'r_{run[o]}/nh3_{o}_phie_45_600ev.dat'); Ph=np.radians(m[:,0])
    q=np.sqrt((2*ke*np.cos(th)-ki)**2+4*ke**2*np.sin(th)**2*np.sin(Ph/2)**2)
    y=m[:,1]*1e3*(2 if o=='1e' else 1)          # 1e: both components (1ex+1ey)
    e=np.loadtxt(D+f'EMS_Bawagan_NH3_{o}.dat'); r=np.loadtxt(D+f'EMS_Rezkallah1CW_NH3_{o}.dat')
    ax.plot(r[:,0],r[:,1],color='#8a8985',lw=1.6,ls='--',label='Rezkallah 2011 (PWBA, as plotted)')
    ax.plot(q,y,color='#2a78d6',lw=2,label='this code (BBK3CW-Z), absolute')
    ax.plot(q,y*e[:,1].max()/y.max() if False else y*1.75,color='#2a78d6',lw=1.2,ls=':',label='this code ×1.75 (shape)')
    ax.plot(e[:,0],e[:,1],'s',ms=4,color='#eb6834',label='Bawagan et al. 1988 (EMS)')
    print(o,'model peak %.3f at q=%.2f ; paper curve peak %.3f at q=%.2f ; exp peak %.2f at q=%.2f'%(y.max(),q[y.argmax()],r[:,1].max(),r[r[:,1].argmax(),0],e[:,1].max(),e[e[:,1].argmax(),0]))
    if fine:
        f=np.loadtxt(f'r_emsf1/nh3_3a1_phie_45_600ev.dat')
        if o=='3a1': ax.plot(q,f[:,1]*1e3,color='#1baf7a',lw=1.2,ls=':',label='this code, finer grids')
    ax.set_xlim(0,2.6); ax.set_title(f'NH$_3$ {o}',loc='left',fontsize=10); ax.set_xlabel('q (a.u.)')
    for s in ('top','right'): ax.spines[s].set_visible(False)
    ax.grid(True,color='#e5e5e5',lw=0.8)
axs[0].set_ylabel('TDCS (10$^{-3}$ a.u.)')
h,l=axs[2].get_legend_handles_labels(); fig.legend(h,l,loc='upper center',ncol=4,frameon=False,fontsize=8)
fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig('ems_compare.png',dpi=150)
