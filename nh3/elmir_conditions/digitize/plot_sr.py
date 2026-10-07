import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
D='/home/user/dd/nh3/elmir_conditions/'
cases=[('2a1','Ei601.0eV_IP27.0eV','sr/r_sr_o3/nh3_2a1_coplanar_6_74ev.dat'),
       ('1e','Ei590.4eV_IP16.4eV','sr/r_sr_o2/nh3_1e_coplanar_6_74ev.dat'),
       ('3a1','Ei584.85eV_IP10.85eV','sr/r_sr_o1/nh3_3a1_coplanar_6_74ev.dat')]
fig,axs=plt.subplots(3,1,figsize=(8,10),sharex=True)
for ax,(o,th,fsr) in zip(axs,cases):
    e=np.loadtxt(D+f'EXP_ElMir_NH3_{o}_Es500eV_Ee74eV_thetas-6deg.dat')
    t0=np.loadtxt(D+f'NH3_{o}_Es500eV_Ee74eV_thetas-6deg_{th}.dat'); t1=np.loadtxt(fsr)
    eb=e[e[:,0]<=150,1].max()
    s0=eb/t0[t0[:,0]<=150,1].max(); s1=eb/t1[t1[:,0]<=150,1].max()
    ax.semilogy(t0[:,0],t0[:,1]*s0,color='#2a78d6',lw=2,ls='--',label='BBK3CW-Z (isr=0)')
    ax.semilogy(t1[:,0],t1[:,1]*s1,color='#1baf7a',lw=2,label='BBK3CW-Z + short-range potential (isr=1)')
    ax.semilogy(e[:,0],e[:,1],'D',ms=5,color='#eb6834',mec='#fcfcfb',mew=0.8,label='Experiment, El Mir et al. (digitized)')
    allv=np.r_[t0[:,1]*s0,t1[:,1]*s1,e[:,1]]
    ax.set_ylim(allv.min()/1.4,allv.max()*1.6); ax.set_xlim(0,360)
    ax.set_title(f'{o}   (Es = 500 eV, θs = −6°, Ee = 74 eV), each curve scaled to exp. binary max',loc='left',fontsize=10,color='#0b0b0b')
    ax.grid(True,color='#e5e5e5',lw=0.8); ax.set_axisbelow(True)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
    for sp in ('left','bottom'): ax.spines[sp].set_color('#8a8985')
    ax.tick_params(colors='#52514e'); ax.set_ylabel('TDCS (arb. units)',color='#52514e'); ax.set_facecolor('#fcfcfb')
h,l=axs[0].get_legend_handles_labels(); fig.legend(h,l,loc='upper center',ncol=3,frameon=False,fontsize=8.5,bbox_to_anchor=(0.5,1.0))
axs[-1].set_xlabel('Ejection angle θe (deg)',color='#52514e'); axs[-1].set_xticks(range(0,361,30))
fig.patch.set_facecolor('#fcfcfb'); fig.tight_layout(rect=(0,0,1,0.97))
fig.savefig(D+'compare_NH3_Ee74eV_SR_vs_ElMir.png',dpi=150)
