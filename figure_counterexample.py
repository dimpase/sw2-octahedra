"""Draw the certified polyhedron and the two forced tiles at its cusp at infinity."""
import argparse
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from counterexample import INF, cusp_set, lorentz, outward
from verify import lor


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    points=cusp_set()
    vs=list(map(lorentz,points))
    ns=outward(vs)
    blue,red="#186594","#a63a17"
    fig,(a,b)=plt.subplots(1,2,figsize=(8.3,4.4),layout="constrained")
    for i,j in combinations(range(1,len(vs)),2):
        if sum(lor(n,vs[i])==lor(n,vs[j])==0 for n in ns)==2:
            (x,y),(u,v)=points[i],points[j]
            a.plot([x,u],[y,v],color="#859199",lw=1,zorder=1)
    a.scatter(*zip(*(pt for pt in points if pt is not INF)),color="#253843",s=12,zorder=3)
    for start,end,color,label in (((2,2),(2,0),blue,"$e_2$"),((1,0),(F(1,2),F(1,2)),red,"$e_1$")):
        a.plot([start[0],end[0]],[start[1],end[1]],color=color,lw=4,zorder=2)
        mx,my=(float(start[0]+end[0])/2,float(start[1]+end[1])/2)
        a.text(mx+.13,my+.13,label,color=color,fontweight="bold")
    a.set_title("Finite cusps and projected edges",fontsize=13)
    a.set(xlim=(-2.35,2.45),ylim=(-2.4,2.4),xticks=range(-2,3),yticks=range(-2,3))
    for size,color in ((2,blue),(1,red)):
        b.add_patch(Rectangle((0,0),size,size,facecolor=color,alpha=.17))
        b.add_patch(Rectangle((0,0),size,size,fill=False,edgecolor=color,lw=2.1))
    b.plot(.5,.5,"o",color="#302d27",ms=4)
    b.text(.28,.73,"$O_1$",color=red,fontsize=14)
    b.text(1.28,1.38,"$O_2$",color=blue,fontsize=14)
    b.set_title("Forced tiles on a horosphere at ∞",fontsize=13)
    b.set(xlim=(-.18,2.18),ylim=(-.18,2.18),xticks=[0,1,2],yticks=[0,1,2])
    for ax in (a,b):
        ax.set_aspect("equal")
        ax.set_xlabel("Re $z$")
        ax.set_ylabel("Im $z$")
        ax.spines[["top","right"]].set_visible(False)
    fig.savefig(args.output,dpi=190)
    plt.close(fig)


if __name__=="__main__":
    main()
