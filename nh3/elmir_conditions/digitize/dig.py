import sys, numpy as np
from PIL import Image
from scipy import ndimage as nd
# calibration: x0 (0 deg) and px per 50 deg from grid cols ; log-y from ticks
def fit(px,val):
    A=np.vstack([px,np.ones(len(px))]).T; s,c=np.linalg.lstsq(A,val,rcond=None)[0]; r=val-(A@[s,c]); return s,c,np.abs(r).max()
cal={
 '1':dict(gx=[178,261,344,675],gd=[50,100,150,550],  # placeholder fixed below
          ty=[24,62.5,127.5,137.5,148.5,161,176,193,214,241,279.5],
          tv=[3e-16,2e-16,1e-16,9e-17,8e-17,7e-17,6e-17,5e-17,4e-17,3e-17,2e-17]),
 '2':dict(gx=[193,276,359,690],
          ty=[36.5,66.5,110,183.5,195,207.5,221,238,257,281,311.5],
          tv=[4e-16,3e-16,2e-16,1e-16,9e-17,8e-17,7e-17,6e-17,5e-17,4e-17,3e-17]),
 '3':dict(gx=[192,275,358,689],
          ty=[15.5,43,78.5,128,212.5,225.5,240,256.5,275.5,297.5],
          tv=[5e-16,4e-16,3e-16,2e-16,1e-16,9e-17,8e-17,7e-17,6e-17,5e-17]),
}
gd=[50,100,150,350]
out={}
for k in ('1','2','3'):
    c=cal[k]
    sx,cx,ex=fit(np.array(c['gx'],float),np.array(gd,float))
    sy,cy,ey=fit(np.array(c['ty']),np.log10(c['tv']))
    a=np.asarray(Image.open(sys.argv[1]+f'/{k}.png').convert('RGB')).astype(int)
    m=(abs(a[:,:,0]-67)<40)&(a[:,:,1]<50)&(abs(a[:,:,2]-83)<40)
    lab,n=nd.label(m)
    sl=nd.find_objects(lab); areas=nd.sum(m,lab,range(1,n+1))
    med=np.median(areas[areas>20])
    pts=[]
    for i,s in enumerate(sl):
        ar=areas[i]
        if ar<0.3*med: continue
        ys,xs=np.nonzero(lab[s]==i+1); ys=ys+s[0].start; xs=xs+s[1].start
        nm=int(round(ar/med))
        if nm<=1:
            pts.append((xs.mean(),ys.mean(),1))
        else:
            # split merged markers with k-means along points
            from scipy.cluster.vq import kmeans2
            P=np.vstack([xs,ys]).T.astype(float)
            cen,_=kmeans2(P,nm,minit='++',seed=1)
            for cc in cen: pts.append((cc[0],cc[1],nm))
    pts.sort()
    res=[(sx*x+cx,10**(sy*y+cy),f) for x,y,f in pts]
    out[k]=res
    print(f'panel {k}: x-cal resid {ex:.2f} deg, y-cal resid {ey:.4f} dec, marker area {med:.0f}px, {len(res)} points, merged groups {sum(1 for p in pts if p[2]>1)}')
    with open(sys.argv[2]+f'/raw_{k}.dat','w') as fo:
        for x,y,f in res: fo.write(f'{x:9.2f} {y:.4e} {f}\n')
