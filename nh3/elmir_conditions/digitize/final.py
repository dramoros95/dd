import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
dig,outd,thd=sys.argv[1],sys.argv[2],sys.argv[3]
panels={'1':('2a1','2a1','Ei601.0eV_IP27.0eV',(1e-17,1e-15)),'2':('1e','1e','Ei590.4eV_IP16.4eV',None),'3':('3a1','3a1','Ei584.85eV_IP10.85eV',None)}
lim={'1':(2.0e-17,3.0e-16),'2':(3.0e-17,4.4e-16),'3':(5.0e-17,5.5e-16)}
fig,axs=plt.subplots(3,1,figsize=(8,10),sharex=True)
rep=[]
for ax,(k,(orb,lab,th,_)) in zip(axs,panels.items()):
    d=np.loadtxt(f'{dig}/raw_{k}.dat')
    lo,hi=lim[k]
    flt={'1':(1.5e-17,3.2e-16),'2':(2e-17,4.5e-16),'3':(4e-17,6e-16)}[k]
    d=d[(d[:,0]>=0)&(d[:,0]<=360)&(d[:,1]>flt[0])&(d[:,1]<flt[1])]
    ang=d[:,0]; dev=0.0
    o=np.argsort(ang); ang=ang[o]; y=d[o,1]
    fn=f'EXP_ElMir_NH3_{orb}_Es500eV_Ee74eV_thetas-6deg.dat'
    np.savetxt(f'{outd}/{fn}',np.c_[ang,y],fmt=['%8.2f','%.4e'])
    t=np.loadtxt(f'{outd}/NH3_{orb}_Es500eV_Ee74eV_thetas-6deg_{th}.dat')
    bx=(ang<=150); be=y[bx].max(); bt=t[t[:,0]<=150,1].max()
    s=be/bt; ts=t[:,1]*s
    # experimental and theoretical binary/recoil
    rx=(ang>=180)&(ang<=330); re=y[rx].max(); rt=ts[(t[:,0]>=180)&(t[:,0]<=330)].max()
    rep.append((orb,len(ang),ang[bx][y[bx].argmax()],t[t[:,0]<=150,0][t[t[:,0]<=150,1].argmax()],
                ang[rx][y[rx].argmax()],be/re,bt*s/rt))
    ax.semilogy(t[:,0],ts,color='#2a78d6',lw=2,label='BBK3CW-Z (this code), scaled to exp. binary peak')
    ax.semilogy(ang,y,'D',ms=5,color='#eb6834',mec='#fcfcfb',mew=0.8,label='Experiment, El Mir et al. (digitized)')
    ax.set_ylim(min(ts.min(),y.min())/1.4,max(ts.max(),y.max())*1.6); ax.set_xlim(0,360)
    ax.set_title(f'{lab}   (Es = 500 eV, θs = −6°, Ee = 74 eV)',loc='left',fontsize=11,color='#0b0b0b')
    ax.grid(True,color='#e5e5e5',lw=0.8); ax.set_axisbelow(True)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
    for sp in ('left','bottom'): ax.spines[sp].set_color('#8a8985')
    ax.tick_params(colors='#52514e'); ax.set_ylabel('TDCS (arb. units)',color='#52514e')
h,l=axs[0].get_legend_handles_labels(); fig.legend(h,l,loc='upper center',ncol=2,frameon=False,fontsize=9,bbox_to_anchor=(0.5,1.0))
axs[-1].set_xlabel('Ejection angle θe (deg)',color='#52514e'); axs[-1].set_xticks(range(0,361,30))
fig.patch.set_facecolor('#fcfcfb')
for a in axs: a.set_facecolor('#fcfcfb')
fig.tight_layout(rect=(0,0,1,0.97)); fig.savefig(f'{outd}/compare_NH3_Ee74eV_vs_ElMir.png',dpi=150)
for r in rep: print('%-4s n=%d  binary exp %.0f theo %.0f | recoil exp max at %.0f | B/R exp %.2f theo %.2f'%r)
