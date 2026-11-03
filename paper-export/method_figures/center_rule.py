# The center rule (Algorithm alg:center, Proposition p:center) on a minimizer with a canonical order on k-mers, for
# a window v0 and its reverse complement (center_rule): the positions j and w-1-j hold the same canonical k-mer and so the same
# rank (dashed lines), and the center rule selects mirrored positions; k = 8, w = 6, L = 13 (odd).
# The ranks are coloured by value (darkest: minimal); the dashed line is the center (w-1)/2; the rank of the selected k-mer is framed.
# python3 paper-export/method_figures/center_rule.py [--out DIR] [--vector-out DIR]   (needs rsvg-convert)
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
k=8; w=6; L=w+k-1
assert L%2==1
S="AGCTCAGCTCAGCG"   # v0 = S[0:13]: the minimal canonical k-mer AGCTCAGC occurs twice, at j = 0 and j = 5, equally close to the center
# (with k = 8, w = 6 a k-mer and its reverse complement cannot sit at mirrored positions: their overlap would be an odd-length reverse palindrome)
assert len(S)==L+1
# canonical order on k-mers: canonical k-mers ordered by a hash; the ranks shown are dense over the canonical k-mers of the window
_cls=sorted({canon(S[q:q+k]) for q in range(w)},key=lambda c:hashlib.md5(c.encode()).hexdigest())
def rank(x): return _cls.index(canon(x))
PAL=["#3f007d","#5b3a9a","#7a68b5","#9e9ac8","#c3c1df","#dddcee","#f0eff7"]   # rank 0 (darkest) .. the largest; purple, apart from the base colours
SANS="DejaVu Sans, Helvetica, Arial, sans-serif"; MONO="DejaVu Sans Mono, Courier New, monospace"
GREEN="#1b9e77"; ORANGE="#d95f02"; CENTER="#777"   # v: teal green, v^rc: orange; as in the canonical mod-minimizer figure
BASE_FILL="#f3f3f6"; BASE_EDGE="#a3a3ad"   # neutral bases, as in the canonical mod-minimizer figure
BW=40; X0=120; RIGHT=60; W=X0+(L+1)*BW+RIGHT; BH=30; RH=28
el=[]
def rect(x,y,w_,h,fill="none",stroke="#000",sw=1,rx=3,dash=None):
    d=f' stroke-dasharray="{dash}"' if dash else ""
    el.append(f'<rect x="{x}" y="{y}" width="{w_}" height="{h}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" rx="{rx}"{d}/>')
def text(x,y,t_,size=19,anchor="middle",weight="normal",fill="#000",family=SANS): el.append(f'<text x="{x}" y="{y+size*0.35}" font-size="{size}" text-anchor="{anchor}" font-weight="{weight}" fill="{fill}" font-family="{family}">{re.sub(r"&(?!(gt|lt|amp);)","&amp;",t_)}</text>')
it=lambda s:f'<tspan font-style="italic">{s}</tspan>'
RCSUP='<tspan dy="-8" font-size="14">rc</tspan><tspan dy="8"></tspan>'
v0,v1=S[:L],S[1:]
def window(y0,start,seq,name,note,labels=False,center=True,frame=False):
    ranks=[rank(seq[j:j+k]) for j in range(w)]; P=[j for j,r in enumerate(ranks) if r==min(ranks)]
    d=lambda j:abs(2*j-(w-1)); b=min(d(j) for j in P); close=[j for j in P if d(j)==b]
    sel=close[0] if len(close)==1 or seq<=rc(seq) else close[-1]   # CenterRule(P, v)
    X=lambda c:X0+(start+c)*BW
    text(X0-12,y0+BH/2,name,21,"end")
    for q,ch in enumerate(seq):
        rect(X(q),y0,BW,BH,BASE_FILL,BASE_EDGE,1,3); text(X(q)+BW/2,y0+BH/2,ch,20,weight="bold",family=MONO)
    yt=y0+BH+10
    for j in range(w):
        r=ranks[j]; rect(X(j)+3,yt,BW-6,RH,PAL[r],"none",0,4)
        text(X(j)+BW/2,yt+RH/2,str(r),16,weight="bold",fill="#fff" if r<3 else "#2a2340")
    yj=yt+RH+12
    for j in range(w): text(X(j)+BW/2,yj,str(j),16,fill="#666")
    if labels: text(X0-12,yt+RH/2,"rank",17,"end",fill="#444"); text(X0-12,yj,it("j"),18,"end",fill="#444")
    xc=X((w-1)/2)+BW/2                                       # the center (w-1)/2
    if center: el.append(f'<line x1="{xc}" y1="{yt-4}" x2="{xc}" y2="{yj+10}" stroke="{RED}" stroke-width="2" stroke-dasharray="4,3"/>')
    rect(X(sel)+3,yt-2,BW-6,RH+4,"none",GREEN,3,4)          # the selected k-mer
    if frame: rect(X(sel)-2,y0-3,k*BW+4,BH+6,"none",GREEN,3,5)   # the selected k-mer, framed in the bases
    if note:                                                 # which of the tied positions is taken, under its j
        text(X(sel)+BW/2,yj+20,note,15,"middle",fill="#777")
        return start+sel,sel,yj+30
    return start+sel,sel,yj+10
def sel_of(v):
    r=[rank(v[j:j+k]) for j in range(w)]; P=[j for j in range(w) if r[j]==min(r)]
    d=lambda j:abs(2*j-(w-1)); b=min(d(j) for j in P); c=[j for j in P if d(j)==b]
    return (c[0] if len(c)==1 or v<=rc(v) else c[-1]),c,r
X=lambda c:X0+c*BW
def frame(x,y,w_,h,col,sw,rx,dash=None):
    rect(x,y,w_,h,"none","#fff",sw+3,rx)          # a white halo, so the frame reads on the dark rank cells
    rect(x,y,w_,h,"none",col,sw,rx,dash)
def bases(y0,seq,name,sel,tied,col):
    text(X0-12,y0+BH/2,name,21,"end")
    for q,ch in enumerate(seq):
        rect(X(q),y0,BW,BH,BASE_FILL,BASE_EDGE,1,3); text(X(q)+BW/2,y0+BH/2,ch,20,weight="bold",family=MONO)
    for t in tied:                                             # the tied k-mers, dashed (a little larger: they overlap the selected one)
        rect(X(t)-6,y0-8,k*BW+12,BH+16,"none",col,2,7,"6,4")
    rect(X(sel)-2,y0-3,k*BW+4,BH+6,"none",col,3,5)              # the selected k-mer, framed in the bases
def ranks(yt,r,sel,tied,col):
    for j in range(w):
        rect(X(j)+3,yt,BW-6,RH,PAL[r[j]],"none",0,4)
        text(X(j)+BW/2,yt+RH/2,str(r[j]),16,weight="bold",fill="#fff" if r[j]<3 else "#2a2340")
    for t in tied: frame(X(t)+3,yt-2,BW-6,RH+4,col,2.4,4,"5,3")   # the rank of a tied k-mer, dashed
    frame(X(sel)+3,yt-2,BW-6,RH+4,col,3,4)                       # the rank of the selected k-mer
def positions(yj):
    for j in range(w): text(X(j)+BW/2,yj,str(j),16,fill="#666")
xc=X((w-1)/2)+BW/2                                             # the center (w-1)/2
vr=rc(v0)
j0,c0,r0=sel_of(v0); jr,cr,rr=sel_of(vr)
assert j0+jr==w-1 and all(r0[j]==rr[w-1-j] for j in range(w))
V0=it("v"); VR=f"{it('v')}{RCSUP}"
# v0: bases, positions j, ranks; v0^rc mirrored below it: ranks, positions j, bases -- the two rank rows face each other
t0=[j for j in c0 if j!=j0]; tr=[j for j in cr if j!=jr]   # the tied, not selected positions
yb0=14; bases(yb0,v0,V0,j0,t0,GREEN)
yj0=yb0+BH+20; positions(yj0); text(X0-12,yj0,it("j"),18,"end",fill="#444")
yt0=yj0+14; ranks(yt0,r0,j0,t0,GREEN); text(X0-12,yt0+RH/2,"rank",17,"end",fill="#444")
ytr=yt0+RH+32; ranks(ytr,rr,jr,tr,ORANGE); text(X0-12,ytr+RH/2,"rank",17,"end",fill="#444")
yjr=ytr+RH+16; positions(yjr); text(X0-12,yjr,it("j"),18,"end",fill="#444")
ybr=yjr+24; bases(ybr,vr,VR,jr,tr,ORANGE)
# the center of both windows
el.append(f'<line x1="{xc}" y1="{yt0-4}" x2="{xc}" y2="{ytr+RH+4}" stroke="{CENTER}" stroke-width="1.6" stroke-dasharray="3,3"/>')
# the mirrored positions j and w-1-j, from rank to rank (equal ranks): the selected pair in green, and the pair 1, w-2
for j in (j0, 1, w-2):
    x1,x2=X(j)+BW/2,X(w-1-j)+BW/2
    sel=j==j0
    el.append(f'<line x1="{x1}" y1="{yt0+RH+(4 if sel else 2)}" x2="{x2}" y2="{ytr-(4 if sel else 2)}" stroke="{"#444" if sel else "#aaa"}" stroke-width="{2.6 if sel else 1.4}" stroke-dasharray="{"7,4" if sel else "4,4"}"/>')
y=ybr+BH+6
H=int(y+10)
svg=f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}"><rect width="100%" height="100%" fill="#fff"/>'+"".join(el)+'</svg>'
open(os.path.join(SVG_DIR, "center_rule.svg"),"w").write(svg)
subprocess.check_call(["rsvg-convert","-f","pdf",os.path.join(SVG_DIR, "center_rule.svg"),"-o",os.path.join(OUT, "center_rule.pdf")])
print("center_rule",W,H,j0,jr)
