"""Reconstruct the symmetric 30-cusp example from seven planar points, then check its SW2 chambers and two forced tiles."""
import argparse
from fractions import Fraction as F
import gzip
from itertools import combinations
import json
from math import lcm
from pathlib import Path

from verify import (CORNERS, IDENTITY, NORMALS, act, check_boundary, cross, evaluate, lor,
                    primitive, require, supporting_planes)


INF = None
ZERO = (F(0),F(0))
GALLERY_PROBE = (7,3,2,1)


def cusp_set():
    generators=[(0,0),(1,0),(2,0),(1,1),(2,2),(F(1,2),F(1,2)),(F(7,5),F(1,5))]
    points=set()
    for a,b in generators:
        for x,y in ((a,b),(b,a)):
            points.update((F(s*x),F(t*y)) for s in (1,-1) for t in (1,-1))
    return [INF,*sorted(points)]


def lorentz(z, height=0):
    if z is INF:
        return (1,0,0,1)
    # The following translation puts the polyhedron in the saved SW2 frame.
    x,y=z[0]+3,z[1]+2
    k=x*x+y*y+height*height
    v=(k+1,2*x,2*y,k-1)
    scale=lcm(*(f.denominator for f in v))
    return primitive(tuple(int(scale*f) for f in v))


def outward(points):
    return {n if max(lor(n,p) for p in points)==0 else tuple(-x for x in n)
            for n in supporting_planes(points).values()}


def chambers(cusps):
    walls=outward(cusps)
    norms={lor(n,n) for n in walls}
    require(norms in ({1},{2}),"Expected SW2 walls of one norm")
    case="all3" if norms=={1} else "no3"
    inside=tuple(map(sum,zip(*cusps)))
    path=[]
    while True:
        face=next((i for i,n in enumerate(NORMALS,1) if lor(n,inside)>0),None)
        if face is None:
            break
        inside=act(cross(IDENTITY,face),inside)
        path.append(str(face))
    translation=evaluate("".join(reversed(path)))
    target={act(translation,n) for n in walls}
    reached={IDENTITY:""}
    queue=[IDENTITY]
    for chamber in queue:
        require(all(lor(n,act(chamber,v))<=0 for n in target for v in CORNERS),
                "An SW2 chamber crosses a proposed facet")
        for i,n in enumerate(NORMALS,1):
            if act(chamber,n) in target:
                continue
            neighbor=cross(chamber,i)
            if neighbor not in reached:
                require(all(lor(a,act(neighbor,GALLERY_PROBE))<0 for a in target),
                        "The hull has a facet outside the SW2 tessellation")
                reached[neighbor]=reached[chamber]+str(i)
                queue.append(neighbor)
    record=dict(case=case,index=len(queue),words=sorted(reached.values(),key=lambda w:(len(w),w)))
    vertices,faces,edges,_=check_boundary(record)
    require(set(vertices)=={act(translation,v) for v in cusps},"Chamber union differs from the specified hull")
    return record,faces,edges,"".join(reversed(path))


def divide(z,w):
    a,b=z
    c,d=w
    norm=c*c+d*d
    return (a*c+b*d)/norm,(b*c-a*d)/norm


def opposites_chart(z,a,b):
    """Send a,b to zero and infinity by an exact fractional linear transformation."""
    if z==b:
        return INF
    if z==a:
        return ZERO
    if a is INF:
        return divide((1,0),(z[0]-b[0],z[1]-b[1]))
    if b is INF:
        return z[0]-a[0],z[1]-a[1]
    if z is INF:
        return F(1),F(0)
    return divide((z[0]-a[0],z[1]-a[1]),(z[0]-b[0],z[1]-b[1]))


def regular_octahedra(points):
    found=set()
    for a,b in combinations(points,2):
        inverse={opposites_chart(z,a,b):z for z in points}
        require(len(inverse)==len(points),"Noninjective change of coordinates")
        for z in inverse:
            if z is INF or z==ZERO:
                continue
            x,y=z
            square=((x,y),(-y,x),(-x,-y),(y,-x))
            if all(s in inverse for s in square):
                found.add(frozenset([a,b]+[inverse[s] for s in square]))
    return found


def gap_index_input(record):
    cells={evaluate(w) for w in record["words"]}
    mirrors={}
    for w in record["words"]:
        g=evaluate(w)
        for j,n in enumerate(NORMALS,1):
            if cross(g,j) not in cells:
                mirrors.setdefault(act(g,n),[int(s) for s in w+str(j)+w[::-1]])
    return """# Reflection-subgroup index from the SW2 Coxeter presentation.
(function()
local F,s,W,r,words,H;
F:=FreeGroup(4);
s:=GeneratorsOfGroup(F);
W:=F/Concatenation(List(s,x->x^2),List([[1,2,2],[1,3,2],[1,4,3],[2,3,4],[2,4,2],[3,4,4]],
    e->(s[e[1]]*s[e[2]])^e[3]));
r:=GeneratorsOfGroup(W);
words:="""+json.dumps(list(mirrors.values()))+""";
H:=Subgroup(W,List(words,w->Product(w,i->r[i])));
Print("Reflection-subgroup index: ",Index(W,H),"\\n");
end)();
"""


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-chambers",type=Path,help="save the reconstructed SW2 words as a TOML or TOML.gz search record")
    parser.add_argument("--write-gap",type=Path,help="write a GAP reflection-subgroup index calculation")
    args=parser.parse_args()
    points=cusp_set()
    vs=list(map(lorentz,points))
    record,faces,edges,word=chambers(vs)
    index=record["index"]
    print(f"Chamber flood: {index} SW2 tetrahedra; {len(vs)} ideal vertices, {len(edges)} right-angled edges, {len(faces)} facets.")
    if args.write_chambers:
        data=(f'bound = {index}\n[visited]\nreconstruction = {index}\n[[polyhedra]]\ncase = "{record["case"]}"\nindex = {index}\n'
              + "words = "+json.dumps(record["words"])+"\n").encode()
        args.write_chambers.write_bytes(gzip.compress(data,mtime=0) if args.write_chambers.suffix==".gz" else data)
    if args.write_gap:
        args.write_gap.write_text(gap_index_input(record))
    transform=evaluate(word)
    transformed=sorted(act(transform,v) for v in vs)
    lookup={z:transformed.index(act(transform,v)) for z,v in zip(points,vs)}
    placements=regular_octahedra(points)
    forced=[]
    for length,edge in ((2,((2,2),(2,0))),(1,((1,0),(F(1,2),F(1,2))))):
        require(frozenset(lookup[z] for z in edge) in edges,"Proposed witness edge is not a boundary edge")
        options=[p for p in placements if set(edge)<=p]
        expected=frozenset((INF,(0,0),(length,0),(0,length),(length,length),(F(length,2),F(length,2))))
        require(options==[expected],"Boundary edge does not force the proposed octahedron")
        forced.append(expected)
    point=lorentz((F(1,2),F(1,2)),height=3)
    require(all(lor(n,point)<0 for tile in forced for n in outward(list(map(lorentz,tile)))),
            "Strict interior overlap failed")
    tile_ids=[{lookup[z] for z in t} for t in placements]
    unique=set()
    for edge in edges:
        options=[i for i,t in enumerate(tile_ids) if edge.issubset(t)]
        if len(options)==1:
            unique.add(options[0])
    print(f"All {len(placements)} regular ideal octahedron placements tested. Boundary edges force {len(unique)} distinct tiles.")
    print("Two forced tiles both contain (1/2, 1/2, height 3) strictly in their interiors in centered upper-half-space coordinates.")
    print("No tiling by regular ideal octahedra exists.")


if __name__=="__main__":
    main()
