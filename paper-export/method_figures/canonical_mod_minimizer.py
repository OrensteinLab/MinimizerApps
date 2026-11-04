# The canonical mod-minimizer (Algorithm alg:mod) on a window v and its reverse complement (canonical_mod_minimizer), in the
# notation of the Proposition and its proof: the t-mer x at position m in [0, T), T = w+k-t, its offset m mod w, the selected
# offset i in v and w-1-i in v^rc; k = 9, t = 4, w = 5, L = 13. Same layout as the center-rule figure (center_rule.py):
# v on top, v^rc mirrored below, one offset track between them; the center offset (w-1)/2 on a dashed gray line; the
# selection of v in green, of v^rc in orange.
# t-mer ranks coloured by value (darkest: minimal); solid frames: selected, dashed: tied but not selected.
# python3 paper-export/method_figures/canonical_mod_minimizer.py [--out DIR] [--vector-out DIR]   (needs rsvg-convert)
import argparse, os, tempfile
_ap = argparse.ArgumentParser(description="Draw a figure of the Methods of the paper (PDF, and optionally its SVG source).")
_ap.add_argument("--out", default=".", help="folder of the PDF (default: the current folder)")
_ap.add_argument("--vector-out", help="also keep the SVG source in this folder")
_args = _ap.parse_args()
OUT = os.path.abspath(_args.out)
SVG_DIR = os.path.abspath(_args.vector_out) if _args.vector_out else tempfile.mkdtemp()
os.makedirs(OUT, exist_ok=True); os.makedirs(SVG_DIR, exist_ok=True)
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import subprocess,hashlib,re
comp={'A':'T','T':'A','C':'G','G':'C'}
rc=lambda s:''.join(comp[c] for c in reversed(s))
canon=lambda s:min(s,rc(s))
k=9; t=4; w=5; L=w+k-1; T=w+k-t
assert t%w==k%w and T%w==0
v="TTCACTGCGTGAA"   # minimal t-mers TCAC (m = 1) and its reverse complement GTGA (m = 8) tie at offsets 1 and 3
assert len(v)==L
vr=rc(v)
# canonical order on t-mers: canonical t-mers ordered by a hash; the ranks shown are dense over the canonical t-mers of v
_cls=sorted({canon(v[m:m+t]) for m in range(T)},key=lambda c:hashlib.md5(c.encode()).hexdigest())
def rank(x): return _cls.index(canon(x))
def mix(a,b,f): return "#"+"".join(f"{round(int(a[i:i+2],16)*(1-f)+int(b[i:i+2],16)*f):02x}" for i in (1,3,5))
PAL=[mix("#3f007d","#f0eff7",q/(len(_cls)-1)) for q in range(len(_cls))]   # rank 0 (darkest) .. the largest
SANS="DejaVu Sans, Helvetica, Arial, sans-serif"; MONO="DejaVu Sans Mono, Courier New, monospace"
GREEN="#1b9e77"; ORANGE="#d95f02"   # v: teal green, v^rc: orange (a colour-blind-safe pair, apart from the purple ranks)
BASE_FILL="#f3f3f6"; BASE_EDGE="#a3a3ad"   # neutral bases: the colours are kept for the selections
BW=40; X0=150; W=X0+L*BW+14; BH=30; RH=28
el=[]
def rect(x,y,w_,h,fill="none",stroke="#000",sw=1,rx=3,dash=None):
    d=f' stroke-dasharray="{dash}"' if dash else ""
    el.append(f'<rect x="{x}" y="{y}" width="{w_}" height="{h}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" rx="{rx}"{d}/>')
def text(x,y,t_,size=19,anchor="middle",weight="normal",fill="#000",family=SANS):
    el.append(f'<text x="{x}" y="{y+size*0.35}" font-size="{size}" text-anchor="{anchor}" font-weight="{weight}" fill="{fill}" font-family="{family}">{re.sub(r"&(?!(gt|lt|amp);)","&amp;",t_)}</text>')
it=lambda s:f'<tspan font-style="italic">{s}</tspan>'
RCSUP='<tspan dy="-8" font-size="14">rc</tspan><tspan dy="8"></tspan>'
X=lambda c:X0+c*BW
def select(seq):
    """Algorithm alg:mod: (selected offset i, offsets P, minimal positions m, ranks)."""
    r=[rank(seq[m:m+t]) for m in range(T)]; M=[m for m in range(T) if r[m]==min(r)]
    P=sorted({m%w for m in M}); d=lambda o:abs(2*o-(w-1)); b=min(map(d,P)); c=[o for o in P if d(o)==b]
    return (c[0] if len(c)==1 or seq<=rc(seq) else c[-1]),P,M,r
def label(y,s,size=17): text(X0-12,y,s,size,"end",fill="#444")
def frame(x,y,w_,h,col,sw,rx,dash=None):
    rect(x,y,w_,h,"none","#fff",sw+3,rx)          # a white halo, so the frame reads on the dark rank cells
    rect(x,y,w_,h,"none",col,sw,rx,dash)
def bases(y0,seq,name,i,P,col):
    label(y0+BH/2,name,21)
    for q,ch in enumerate(seq):
        rect(X(q),y0,BW,BH,BASE_FILL,BASE_EDGE,1,3); text(X(q)+BW/2,y0+BH/2,ch,20,weight="bold",family=MONO)
    for o in P:
        if o!=i: rect(X(o)-6,y0-8,k*BW+12,BH+16,"none",col,2,7,"6,4")   # the k-mer at a tied offset, dashed
    rect(X(i)-2,y0-3,k*BW+4,BH+6,"none",col,3,5)                         # the selected k-mer
def positions(y):
    for m in range(T): text(X(m)+BW/2,y,str(m),16,fill="#666")
    label(y,it("m"),18)
def ranks(y,r,M,i,col):
    for m in range(T):
        rect(X(m)+3,y,BW-6,RH,PAL[r[m]],"none",0,4)
        text(X(m)+BW/2,y+RH/2,str(r[m]),16,weight="bold",fill="#fff" if r[m]<len(PAL)/2 else "#2a2340")
    for m in M:   # the minimal t-mers: solid for the one that gives the selected offset, dashed otherwise
        if m%w==i: frame(X(m)+3,y-2,BW-6,RH+4,col,3,4)
        else: frame(X(m)+3,y-2,BW-6,RH+4,col,2.4,4,"5,3")
    label(y+RH/2,f"{it('t')}-mer rank")
i0,P0,M0,r0=select(v); ir,Pr,Mr,rr=select(vr)
assert i0+ir==w-1 and P0==Pr and all(r0[m]==rr[T-1-m] for m in range(T))
# v: bases, m, t-mer ranks; one offset track m mod w (shared: both strands have the offsets P); v^rc mirrored below
yb0=14; bases(yb0,v,it("v"),i0,P0,GREEN)
ym0=yb0+BH+20; positions(ym0)
yt0=ym0+14; ranks(yt0,r0,M0,i0,GREEN)
yo=yt0+RH+26                                                        # the offset track (center line of the numbers)
ytr=yo+26
label(yo,f"{it('m')} mod {it('w')}")
for o in range(w):
    xo=X(o)+BW/2
    col=GREEN if o==i0 else ORANGE if o==ir else None
    if col: el.append(f'<circle cx="{xo}" cy="{yo}" r="14" fill="{mix(col,"#ffffff",0.85)}" stroke="{col}" stroke-width="2.6"/>')
    if 2*o==w-1:                                                    # the center (w-1)/2: a dashed gray line through the track
        for y1,y2 in ((yt0+RH+4,yo-13),(yo+13,ytr-4)):
            el.append(f'<line x1="{xo}" y1="{y1}" x2="{xo}" y2="{y2}" stroke="#777" stroke-width="1.6" stroke-dasharray="3,3"/>')
    text(xo,yo,str(o),17,weight="bold" if o in P0 else "normal",fill=col or ("#222" if o in P0 else "#888"))
ranks(ytr,rr,Mr,ir,ORANGE)
ymr=ytr+RH+16; positions(ymr)
ybr=ymr+24; bases(ybr,vr,f"{it('v')}{RCSUP}",ir,Pr,ORANGE)
H=int(ybr+BH+16)
svg=f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}"><rect width="100%" height="100%" fill="#fff"/>'+"".join(el)+'</svg>'
name="canonical_mod_minimizer"
open(os.path.join(SVG_DIR, name + ".svg"),"w").write(svg)
subprocess.check_call(["rsvg-convert","-f","pdf",os.path.join(SVG_DIR, name + ".svg"),"-o",os.path.join(OUT, name + ".pdf")])
print(name,W,H,"i =",i0,"i_rc =",ir,"P =",P0,"M =",M0,"v<=v^rc:",v<=vr)
