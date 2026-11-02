# Canonical GreedyMini order on (k+2i)-mers, applied to one window, in the terms of the Methods (canonical string C(x), its
# representative, projection Phi, rank_rho): the w overlapping (k+2i)-mers u = xyz of an L-window v, each ranked by
# rank_rho(C(Phi(y))) of its center y (primary key) and, between equal ranks, by its canonical representative (secondary key,
# an injective key of C(u)); the minimal (k+2i)-mer is selected. rho is a real GreedyMini order: the canonical order on binary
# 6-mers for w = 5 (orders/GreedyMini/w5_k6_canonical.gm, the t-mer order of the Kraken2 mod-minimizer), here k = 6, i = 1.
# python3 paper-export/method_figures/greedymini_order.py [--out DIR] [--vector-out DIR]  (needs rsvg-convert)
import argparse, os, tempfile
_ap = argparse.ArgumentParser(description="Draw a figure of the Methods of the paper (PDF, and optionally its SVG source).")
_ap.add_argument("--out", default=".", help="folder of the PDF (default: the current folder)")
_ap.add_argument("--vector-out", help="also keep the SVG source in this folder")
_args = _ap.parse_args()
OUT = os.path.abspath(_args.out)
SVG_DIR = os.path.abspath(_args.vector_out) if _args.vector_out else tempfile.mkdtemp()
os.makedirs(OUT, exist_ok=True); os.makedirs(SVG_DIR, exist_ok=True)
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import struct
import subprocess
GM=os.path.join(REPO, "orders", "GreedyMini", "w5_k6_canonical.gm")
d=open(GM,"rb").read(); (N,)=struct.unpack_from("<Q",d,0); RAW=struct.unpack_from(f"<{N}Q",d,8)
Phi=lambda s:''.join('0' if c in 'AT' else '1' for c in s)
rc=lambda s:s[::-1].translate(str.maketrans('ACGT','TGCA'))
rank=lambda b:None if RAW[int(b,2)]>=N else RAW[int(b,2)]   # strand-symmetric: rank(b) = rank(b reversed)
k=6; i=1; n=k+2*i; w=5; L=w+n-1
v="GTAGATCAAGAG"; assert len(v)==L
U=[v[p:p+n] for p in range(w)]; Y=[u[i:i+k] for u in U]; B=[Phi(y) for y in Y]; R=[rank(b) for b in B]
assert all(rank(b)==rank(b[::-1]) for b in B) and None not in R
REP=[min(u,rc(u)) for u in U]
sel=min(range(w),key=lambda p:(R[p],REP[p])); tie=[p for p in range(w) if R[p]==R[sel]]
MONO="DejaVu Sans Mono, Courier New, monospace"; SANS="DejaVu Sans, Helvetica, Arial, sans-serif"
# colours shared with the other method figures: neutral bases, the projection in two colours (A,T and 0: blue; C,G and 1: orange),
# ranks in purple (darker: smaller), the selected k-mer in teal green
GREEN="#1b9e77"; NEUTRAL="#f3f3f6"; NEDGE="#a3a3ad"; P0="#d6e4f5"; P1="#fbdcc4"; YEDGE="#555555"; TIE="#3f007d"
it=lambda t:f'<tspan font-style="italic">{t}</tspan>'
sub=lambda t:f'<tspan dy="5" font-size="13">{t}</tspan><tspan dy="-5"></tspan>'
CS=25; RH=40
el=[]
def rect(x,y_,w_,h,fill="none",stroke="#000",sw=1,rx=3,dash=None):
    dd=f' stroke-dasharray="{dash}"' if dash else ""
    el.append(f'<rect x="{x}" y="{y_}" width="{w_}" height="{h}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" rx="{rx}"{dd}/>')
def text(x,y_,t,size=19,anchor="middle",weight="normal",fill="#000",family=SANS):
    el.append(f'<text x="{x}" y="{y_+size*0.35}" font-size="{size}" text-anchor="{anchor}" font-weight="{weight}" fill="{fill}" font-family="{family}">{t}</text>')
def arrow(x1,x2,y_,label=""):
    el.append(f'<polyline points="{x1},{y_} {x2-9},{y_}" fill="none" stroke="#444" stroke-width="2.4"/>')
    el.append(f'<path d="M{x2},{y_} l-11,-6 l0,12 z" fill="#444"/>')
    if label: text((x1+x2)/2,y_-13,label,15,weight="bold",fill="#222")
def lcells(x,y_,seq,sz=22,proj=None):
    """DNA letters; the ones in proj (positions) are coloured by their projection, the others neutral."""
    for q,ch in enumerate(seq):
        on=proj is None or q in proj
        rect(x+q*CS,y_,sz,sz,(P1 if ch in "CG" else P0) if on else NEUTRAL,"#666" if on else NEDGE,1,3); text(x+q*CS+sz/2,y_+sz/2,ch,16,weight="bold",family=MONO)
def bcells(x,y_,b,sz=22):
    for q,ch in enumerate(b):
        rect(x+q*CS,y_,sz,sz,P1 if ch=="1" else P0,"#666",1,3); text(x+q*CS+sz/2,y_+sz/2,ch,16,weight="bold",family=MONO)
def pair(x,y_,a,b):
    text(x,y_+11,"(",22,"middle",fill="#333"); bcells(x+10,y_,a)
    text(x+14+len(a)*CS,y_+13,",",22,"middle",fill="#333"); bcells(x+22+len(a)*CS,y_,b); text(x+24+2*len(a)*CS,y_+11,")",22,"middle",fill="#333")
    return x+30+2*len(a)*CS
X0=40; X1=X0+L*CS+50; X2=X1+2*k*CS+30+60; X3=X2+70
W=X3+n*CS+20
yh=16; yv=yh+28                      # headers; the window v
text(X1+k*CS+15,yh,f"{it('C')}(Φ({it('y')}))",16,weight="bold",fill="#444")
text(X2+8,yh,f"rank{sub('π')}",16,weight="bold",fill="#444")
text(X3+n*CS/2,yh,f"{it('C')}({it('u')})",16,weight="bold",fill="#444")
text(X0-12,yv+11,it("v"),19,"end",fill="#444"); lcells(X0,yv,v,proj=())
rect(X0+sel*CS-5,yv-5,n*CS+7,32,"none",GREEN,3,5)                   # the selected 8-mer, in the window
y1=yv+44                             # labels x, y, z over the first 8-mer
xs=X0
for lab,a,b in (("x",0,i),("y",i,i+k),("z",i+k,n)):
    text(xs+(a+b)*CS/2-1.5,y1,it(lab),17,fill=YEDGE if lab=="y" else "#444",weight="bold")
for p in range(w):
    yy=y1+16+p*RH; ym=yy+11; xu=X0+p*CS
    text(X0-12,ym,f"{it('u')}{sub(p)}",17,"end",fill="#444")
    lcells(xu,yy,U[p],proj=range(i,i+k))                         # only the center y is projected
    rect(xu+i*CS-3,yy-3,k*CS+2,28,"none",YEDGE,2,4)              # the center y
    arrow(X0+L*CS+8,X1-8,ym,"Φ" if p==0 else "")
    xe=pair(X1,yy,B[p],B[p][::-1])
    if p in tie and len(tie)>1: rect(X1-8,yy-5,xe-X1+12,32,"none",TIE,2,5,"6,4")   # the tied canonical binary k-mers
    arrow(xe+4,X2-10,ym)
    rect(X2-6,yy,28,22,"#3f007d" if R[p]==R[sel] else "#bcbddc","none",0,4)
    text(X2+8,ym,str(R[p]),16,weight="bold",fill="#fff" if R[p]==R[sel] else "#2a2340")
    arrow(X2+26,X3-8,ym)
    lcells(X3,yy,REP[p],proj=())
    if p==sel: rect(X3-3,yy-3,n*CS+3,28,"none",GREEN,3,4)
if len(tie)>1:                                                       # the tied minimal ranks, each in a dashed frame
    for p in tie: rect(X2-10,y1+16+p*RH-4,36,30,"none",TIE,2,5,"6,4")
H=int(y1+16+w*RH)
svg=f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}"><rect width="100%" height="100%" fill="#fff"/>'+"".join(el)+'</svg>'
open(os.path.join(SVG_DIR, "greedymini_order.svg"),"w").write(svg)
subprocess.check_call(["rsvg-convert","-f","pdf",os.path.join(SVG_DIR, "greedymini_order.svg"),"-o",os.path.join(OUT, "greedymini_order.pdf")])
print(W,H,"ranks",R,"selected",sel,"tie",tie,REP)
