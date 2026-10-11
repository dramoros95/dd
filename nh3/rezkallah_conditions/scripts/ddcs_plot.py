import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
e=np.loadtxt('../rezdata/DDCS_Opal_NH3_Ei500_the120.dat'); r=np.loadtxt('../rezdata/DDCS_Rezkallah1CW_NH3_Ei500_the120.dat'); m=np.loadtxt('ddcs_model.dat',ndmin=2)
fig,ax=plt.subplots(figsize=(6.4,4.4))
ax.semilogy(r[:,0],r[:,1],color='#8a8985',ls='--',lw=1.6,label='Rezkallah 2011, 1CW (as plotted)')
ax.semilogy(m[:,0],m[:,1],'-o',color='#2a78d6',lw=2,ms=5,label='this code (BBK3CW-Z), total')
for j,(lab,col) in enumerate([('3a$_1$','#1baf7a'),('1e','#8c5bd6'),('2a$_1$','#c2410c')]):
    ax.semilogy(m[:,0],m[:,2+j],':',color=col,lw=1.3,label=f'  {lab} part')
ax.semilogy(e[:,0],e[:,1],'s',ms=4.5,color='#eb6834',label='Opal et al. 1972')
ax.set_xlabel('Ejection energy (eV)'); ax.set_ylabel('DDCS (cm$^2$ eV$^{-1}$ sr$^{-1}$)'); ax.set_xlim(0,205)
ax.set_title('NH$_3$, E$_i$ = 500 eV, θ$_e$ = 120°',loc='left',fontsize=10)
for s in ('top','right'): ax.spines[s].set_visible(False)
ax.grid(True,color='#e5e5e5',lw=0.8); ax.legend(frameon=False,fontsize=7.5)
fig.tight_layout(); fig.savefig('ddcs_compare.png',dpi=150)
