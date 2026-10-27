import sys, numpy as np, scipy.io
import verify_spectral_graphs as V
from src.spectral_graphs import build_graph
cfgs = {
 "baseline(S=4,full)": ("baseline", dict(num_segments=4, k_lap=1000, order_filt=5)),
 "baseline(S=1,full)": ("baseline", dict(num_segments=1, k_lap=1000, order_filt=5)),
 "baseline(S=4,k10)":  ("baseline", dict(num_segments=4, k_lap=10, order_filt=5)),
}
for k,knn in [(10,"mutual"),(10,"union"),(30,"mutual"),(1000,"none")]:
    for h in [0,1,2,4]:
        cfgs[f"abslog h{h} k{k}{knn[0]}"]=("abs",dict(feature="log",h=h,k=k,knn=knn))
cfgs["abslog h2 k30m unsigned"]=("abs",dict(feature="log",h=2,k=30,knn="mutual",signed=False))
cfgs["abslog h2 k1000 unsigned"]=("abs",dict(feature="log",h=2,k=1000,knn="none",signed=False))
cfgs["abslin h2 k30m"]=("abs",dict(feature="abs",h=2,k=30,knn="mutual"))
cfgs["coh sep1 k10"]=("coherence",dict(min_sep=1))
cfgs["coh sep2 k10"]=("coherence",dict(min_sep=2))
cfgs["coh sep2 k30"]=("coherence",dict(min_sep=2,k=30))
cfgs["coh sep2 k30 unsigned"]=("coherence",dict(min_sep=2,k=30,signed=False))
for image in sys.argv[1:]:
    d=scipy.io.loadmat(f'dataset/{image}/g0.1_ps0.05_pt0.05_pd0.01.mat', simplify_cells=True)
    C=d['HSI_clean'].astype(float); Y=d['HSI_noisy'].astype(float); Ng=np.asarray(d['deg']['gaussian_noise'],float); Yg=C+Ng
    dg=V.gt_delta(C)
    print(image, 'tikhonov_g=%.2f'%V.tikhonov_reference(Yg,C)[0])
    print(f"{'config':28s} corrW  | GT: Rc   sel  psnr_g w_sf sgn | NOISY: Rc   sel  psnr_g w_sf sgn iso")
    for nm,(b,kw) in cfgs.items():
        gg=build_graph(b,C,**kw); gn=build_graph(b,Y,**kw)
        c=V.consistency(gg,gn); out=[]
        for g in (gg,gn):
            sm=V.smoothness(g,C,Y,Ng); it=V.intent(g,dg); pg=V.proxy_denoise(g,Yg,C)[0]
            out.append((sm['R_clean'],sm['sel'],pg,it['w_sf'],it['sign_ok'],it['iso_steep']))
        a,bb=out
        print(f"{nm:28s} {c['corr_W']:5.2f}  | {a[0]:.3f} {a[1]:5.2f} {a[2]:6.2f} {a[3]:.2f} {a[4]:.2f} | {bb[0]:.3f} {bb[1]:5.2f} {bb[2]:6.2f} {bb[3]:.2f} {bb[4]:.2f} {bb[5]:.2f}", flush=True)
