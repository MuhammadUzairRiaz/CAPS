import math, random, itertools, sys
sys.argv=[""]
from screen_row14 import simulate, lab, dE, hex2rgb, to_lin
KINDS=["Protanopia","Deuteranopia","Tritanopia"]
def lum(h):
    r,g,b=[to_lin(c) for c in hex2rgb(h)]; return .2126*r+.7152*g+.0722*b
def cr(a,b):
    la,lb=sorted([lum(a),lum(b)],reverse=True); return (la+.05)/(lb+.05)
def lab2hex(L,a,b):
    fy=(L+16)/116; fx=fy+a/500; fz=fy-b/200
    fi=lambda t: t**3 if t**3>0.008856 else (t-16/116)/7.787
    X,Y,Z=fi(fx)*.95047,fi(fy),fi(fz)*1.08883
    r= 3.2406*X-1.5372*Y-0.4986*Z; g=-0.9689*X+1.8758*Y+0.0415*Z; bb=0.0557*X-0.2040*Y+1.0570*Z
    if min(r,g,bb)<-1e-4 or max(r,g,bb)>1+1e-4: return None
    s=lambda c: 12.92*c if c<=0.0031308 else 1.055*max(c,0)**(1/2.4)-0.055
    return "#%02X%02X%02X"%tuple(round(min(1,max(0,s(c)))*255) for c in (r,g,bb))
def minpair(cols):
    m=(99,None)
    for k in ["Normal"]+KINDS:
        cs=[c if k=="Normal" else simulate(c,k) for c in cols]
        for i,j in itertools.combinations(range(len(cs)),2):
            d=dE(cs[i],cs[j])
            if d<m[0]: m=(d,(k,i,j))
    return m
def cost(cols,orig,fixed,T=14.5):
    c=0
    for k in ["Normal"]+KINDS:
        cs=[x if k=="Normal" else simulate(x,k) for x in cols]
        for i,j in itertools.combinations(range(len(cs)),2):
            d=dE(cs[i],cs[j]); 
            if d<T: c+=(T-d)**2*10
    for x,o in zip(cols,orig): c+=dE(x,o)**2*0.02
    for x in cols:
        if cr(x,"#0F1113")<3: c+=200
    return c
def opt(orig,fixed=(),iters=6000,seed=1):
    random.seed(seed); cols=list(orig); best=cost(cols,orig,fixed)
    for it in range(iters):
        i=random.randrange(len(cols))
        if i in fixed: continue
        L,a,b=lab(cols[i]); s=6*(1-it/iters)+1
        h=lab2hex(L+random.gauss(0,s),a+random.gauss(0,s*1.5),b+random.gauss(0,s*1.5))
        if not h: continue
        tr=cols[:]; tr[i]=h; c=cost(tr,orig,fixed)
        if c<best: best,cols=c,tr
    return cols
if __name__=="__main__":
    E=["#8E959C","#E9ECEF","#E5534B","#4C7BD9","#D6A45E","#E3C74A","#9B7BD6","#57B26A"]  # C H O N Si S Na Cl
    CH=["#F0A83C","#6CC4D8","#E07A5F","#9B7BD6","#7CC784","#D6A45E","#E9ECEF","#4C7BD9"]
    ST=["#7CC784","#F0A83C","#FF7B72","#6CC4D8"]
    for name,p,fx in (("E",E,(1,)),("CH",CH,(0,1,6)),("ST",ST,(1,3))):
        before=minpair(p); o=opt(p,fx); after=minpair(o)
        print(name,"before %.1f"%before[0],before[1],"after %.1f"%after[0],after[1])
        for a,b in zip(p,o): print("  ",a,"->",b,"dE %.1f"%dE(a,b),"cr %.1f"%cr(b,"#0F1113"))
