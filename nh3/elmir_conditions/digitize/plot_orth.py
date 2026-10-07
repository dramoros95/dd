import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
D='/home/user/dd/nh3/elmir_conditions/'
O='/tmp/claude-0/-home-user-dd/f5edd380-5c99-5625-bceb-482d3cf9840e/scratchpad/orth/'
cases=[('2a1','o3','Ei601.0eV_IP27.0eV'),('1e','o2','Ei590.4eV_IP16.4eV'),('3a1','o1','Ei584.85eV_IP10.85eV')]
def br(a): return a[a[:,0]<=150,1].max()/a[a[:,0]>=180,1].max()
fig,axs=plt.subplots(3,1,figsize=(8,10),sharex=True)
for ax,(o,c,th) in zip(axs,cases):
    e=np.loadtxt(D+f'EXP_ElMir_NH3_{o}_Es500eV_Ee74eV_thetas-6deg.dat')
    base=D+f'NH3_{o}_Es500eV_Ee74eV_thetas-6deg_{th}'
    cur={'BBK3CW-Z':np.loadtxt(base+'.dat'),
         'ORTH':np.loadtxt(O+f'r_{c}_sr0/nh3_{o}_coplanar_6_74ev.dat'),
         'ORTH+SR':np.loadtxt(O+f'r_{c}_sr1/nh3_{o}_coplanar_6_74ev.dat')}
    np.savetxt(base+'_ORTH.dat',cur['ORTH']); np.savetxt(base+'_ORTH_SR.dat',cur['ORTH+SR'])
    print(o,'exp B/R %.2f'%br(e),' '.join('%s %.2f'%(k,br(v)) for k,v in cur.items()))
    eb=e[e[:,0]<=150,1].max(); allv=[e[:,1]]
    for (k,v),col,ls in zip(cur.items(),['#2a78d6','#1baf7a','#8c5bd6'],['--','-','-.']):
        s=eb/v[v[:,0]<=150,1].max(); ax.semilogy(v[:,0],v[:,1]*s,color=col,lw=2,ls=ls,label=k); allv.append(v[:,1]*s)
    ax.semilogy(e[:,0],e[:,1],'D',ms=5,color='#eb6834',mec='#fcfcfb',mew=0.8,label='Experiment (El Mir et al.)')
    allv=np.concatenate(allv); ax.set_ylim(allv.min()/1.4,allv.max()*1.6); ax.set_xlim(0,360)
    ax.set_title(f'{o}  (Es=500 eV, θs=−6°, Ee=74 eV), scaled to exp. binary max',loc='left',fontsize=10)
    ax.grid(True,color='#e5e5e5',lw=0.8)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
    ax.set_ylabel('TDCS (arb. units)')
h,l=axs[0].get_legend_handles_labels(); fig.legend(h,l,loc='upper center',ncol=4,frameon=False,fontsize=8.5)
axs[-1].set_xlabel('Ejection angle θe (deg)'); axs[-1].set_xticks(range(0,361,30))
fig.tight_layout(rect=(0,0,1,0.97)); fig.savefig(D+'compare_NH3_Ee74eV_ORTH_vs_ElMir.png',dpi=150)
