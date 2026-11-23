"""Enumerate ideal right-angled SW2 polyhedra in a prescribed cusp rectangle and mirror curvature range."""
import argparse
from collections import defaultdict
from functools import cache
import gzip
from itertools import combinations
import json
from pathlib import Path
import random
import time

from verify import NORMALS, add, bracket, cusp_coord, homogeneous, lor, mul, primitive, require


INFINITY = (1, 0, 0, 1)


def bits(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length()-1
        mask ^= low


@cache
def simple_mirror(n):
    """Explicit reflection descent certifying that n is an SW2 wall."""
    probe = (7, 3, 2, 1)
    if lor(n, probe) > 0:
        n = tuple(-x for x in n)
    while n not in NORMALS:
        height = -lor(n, probe)
        for simple in NORMALS:
            factor = 2 * lor(n, simple) // lor(simple, simple)
            if factor > 0:
                other = tuple(a-factor*b for a,b in zip(n,simple))
                if 0 < -lor(other,probe) < height:
                    n = other
                    break
        else:
            return False
    return True


def cusp_triangles(v,incident,q):
    require(len(incident)==4,"Expected a rectangular cusp link")
    lengths=[]
    for a,b in combinations(incident,2):
        if lor(a,b)==-q:
            side=tuple(x+y for x,y in zip(a,b))
            m,rem=divmod(side[0],v[0])
            require(m>0 and rem==0 and side==tuple(m*x for x in v),"Invalid opposite cusp walls")
            lengths.append(m)
    require(len(lengths)==2,"Expected two pairs of opposite cusp walls")
    return 4*lengths[0]*lengths[1]//q


def chamber_index(cusps,walls):
    """Count SW2 cuspidal triangles, one per chamber, from the four walls at each cusp."""
    q=lor(walls[0],walls[0])
    require(q in (1,2) and all(lor(n,n)==q for n in walls),"Expected walls of a single norm")
    return sum(cusp_triangles(v,[n for n in walls if lor(v,n)==0],q) for v in cusps)


def edge_octahedra(a,b,c):
    """Both regular octahedra with a,b an edge and a,b,c a face, in Riemann-sphere coordinates."""
    p,q=bracket(a,c),bracket(c,b)
    bn,bd=mul(p,b[:2]),mul(p,(b[2],0))
    an,ad=mul(q,a[:2]),mul(q,(a[2],0))
    def image(z,factor):
        return homogeneous(add(mul(z,bn),mul((factor,0),an)),add(mul(z,bd),mul((factor,0),ad)))
    for sign in (-1,1):
        yield (a,b,c,image((0,sign),1),image((1,sign),1),image((1,sign),2))


class CuspSearch:
    def __init__(self, case, width, height, bend, origin=(0,0)):
        require(case in ("all3","no3"), "Unknown wall case")
        require(min(width,height,bend)>0, "Expected positive rectangle dimensions and curvature bound")
        self.case, self.q = case, (1 if case=="all3" else 2)
        self.left,self.bottom = origin
        self.right,self.top = self.left+width,self.bottom+height
        self.bend = bend
        self.walls,self.circles = self.make_walls()
        self.edges = {}
        self.orthogonal = [0]*len(self.walls)
        self.compatible = [0]*len(self.walls)
        self.prepare()

    def cusp(self, a, b, d):
        if self.case == "no3":
            a,b,d = a+b,a-b,2*d
        s = a*a+b*b
        return primitive((s+d*d,2*a*d,2*b*d,s-d*d))

    def wall(self, a, b, d):
        # In chart coordinates the center is (a/d,b/d) and the radius is q/d.
        nx,ny = a,b
        if self.case=="no3":
            if (a+b)%2:
                return None
            nx,ny = (a+b)//2,(a-b)//2
        numer=d*d+nx*nx+ny*ny-self.q
        if numer%(2*d):
            return None
        t=numer//(2*d)
        return (t,nx,ny,t-d)

    def make_walls(self):
        axes=((1,0),(0,1)) if self.case=="all3" else ((1,1),(1,-1))
        walls=[]
        for (a,b),lo,hi in zip(axes,(self.left,self.bottom),(self.right,self.top)):
            walls.extend([(-lo,-a,-b,-lo),(hi,a,b,hi)])
        require(all(simple_mirror(n) for n in walls),"Rectangle side outside SW2")
        corners=[self.cusp(x,y,1) for x in (self.left,self.right) for y in (self.bottom,self.top)]
        sides=list(walls)
        circles=[None]*4
        for d in range(1,self.bend+1):
            for a in range(self.left*d,self.right*d+1):
                for b in range(self.bottom*d,self.top*d+1):
                    n=self.wall(a,b,d)
                    if n is None:
                        continue
                    if not all(z==0 or z<=-self.q for z in (lor(n,m) for m in sides)):
                        continue
                    if any(lor(n,v)>0 for v in corners):
                        continue
                    if simple_mirror(n):
                        walls.append(n)
                        circles.append((a,b,d))
        return walls,circles

    def endpoints(self, i, j):
        # i<j, and the two supporting planes are perpendicular.
        if j<4:
            return (INFINITY,self.cusp((self.left,self.right)[i],(self.bottom,self.top)[j-2],1))
        a,b,d=self.circles[j]
        if i<4:
            dx,dy=(0,self.q) if i<2 else (self.q,0)
            return self.cusp(a+dx,b+dy,d),self.cusp(a-dx,b-dy,d)
        x,y,e=self.circles[i]
        u,v=x*d-a*e,y*d-b*e
        denom=d*d+e*e
        cx,cy=a*denom+e*u,b*denom+e*v
        return self.cusp(cx-d*v,cy+d*u,d*denom),self.cusp(cx+d*v,cy-d*u,d*denom)

    def prepare(self):
        gaps=defaultdict(int)
        incident={}
        permitted={}
        for j,m in enumerate(self.walls):
            for i,n in enumerate(self.walls[:j]):
                p=lor(n,m)
                if p==0 or p<=-self.q:
                    self.compatible[i] |= 1<<j
                    self.compatible[j] |= 1<<i
                if p==-self.q:
                    gaps[i] |= 1<<j
                    gaps[j] |= 1<<i
                if p!=0:
                    continue
                ends=self.endpoints(i,j)
                require(all(lor(v,v)==lor(n,v)==lor(m,v)==0 for v in ends) and ends[0]!=ends[1],
                        "Incorrect ideal endpoints")
                self.orthogonal[i] |= 1<<j
                self.orthogonal[j] |= 1<<i
                self.edges[i,j]=ends
        for ends in self.edges.values():
            for v in ends:
                if v not in incident:
                    zero,inside=0,0
                    for k,n in enumerate(self.walls):
                        p=lor(v,n)
                        if p==0:
                            zero |= 1<<k
                        if p<=0:
                            inside |= 1<<k
                    incident[v],permitted[v]=zero,inside
        self.constraints={}
        for (i,j), (v,w) in self.edges.items():
            needed=set()
            for end in (v,w):
                needed.add(incident[end] & gaps[i] & self.orthogonal[j])
                needed.add(incident[end] & gaps[j] & self.orthogonal[i])
            self.constraints[i,j]=permitted[v]&permitted[w],needed
        self.sites=tuple(incident)
        self.site_id={v:i for i,v in enumerate(self.sites)}
        self.sphere=tuple(map(cusp_coord,self.sites))
        self.sphere_id={v:i for i,v in enumerate(self.sphere)}
        self.wall_sites=[set() for n in self.walls]
        for (i,j),ends in self.edges.items():
            ids=[self.site_id[v] for v in ends]
            self.wall_sites[i].update(ids)
            self.wall_sites[j].update(ids)
        self.edge_tiles={}
        self.cusp_counts={}

    def solutions(self, seed):
        rng=random.Random(seed)
        def insert(k,selected,available,needs):
            available &= self.compatible[k]
            pending=needs.copy()
            for i in bits(selected & self.orthogonal[k]):
                keep,extra=self.constraints[min(i,k),max(i,k)]
                if selected & ~keep:
                    return None
                available &= keep
                pending.update(extra)
            selected |= 1<<k
            pending={mask for mask in pending if not mask&selected}
            return selected,available,pending

        def visit(selected,available,needs):
            if not needs:
                yield tuple(bits(selected))
                return
            options=min((mask&available for mask in needs),key=lambda mask:(mask.bit_count(),mask))
            candidates=list(bits(options))
            rng.shuffle(candidates)
            for k in candidates:
                state=insert(k,selected,available,needs)
                if state is not None:
                    yield from visit(*state)
                available &= ~(1<<k)

        state=0,(1<<len(self.walls))-1,set()
        for k in range(4):
            state=insert(k,*state)
            if state is None:
                return
        yield from visit(*state)

    def record(self, chosen):
        walls=[self.walls[i] for i in chosen]
        cusps=sorted({v for pair in combinations(chosen,2) for v in self.edges.get(pair,())})
        require(all(sum(lor(v,n)==0 for n in walls)==4 and all(lor(v,n)<=0 for n in walls) for v in cusps),
                "Incomplete cusp link")
        require(all(sum(lor(v,n)==0 for v in cusps)>=3 for n in walls),"Degenerate facet")
        require(len(walls)==len(cusps)+2,"Incorrect Euler characteristic")
        return dict(cusps=cusps,walls=walls,index=chamber_index(cusps,walls))

    def candidates_at_edge(self,pair):
        if pair not in self.edge_tiles:
            i,j=pair
            a,b=(self.sphere[self.site_id[v]] for v in self.edges[pair])
            tiles=set()
            for k in min((self.wall_sites[i],self.wall_sites[j]),key=len):
                c=self.sphere[k]
                if c==a or c==b:
                    continue
                for points in edge_octahedra(a,b,c):
                    if all(p in self.sphere_id for p in points):
                        tiles.add(sum(1 << self.sphere_id[p] for p in points))
            self.edge_tiles[pair]=tuple(sorted(tiles))
        return self.edge_tiles[pair]

    def screen(self, chosen):
        past=0
        pairs=[]
        corners=defaultdict(int)
        for j in chosen:
            for i in bits(self.orthogonal[j]&past):
                pair=(i,j)
                pairs.append(pair)
                for v in self.edges[pair]:
                    corners[self.site_id[v]] |= (1<<i)|(1<<j)
            past |= 1<<j
        index=0
        mask=0
        for k,incident in corners.items():
            key=k,incident
            if key not in self.cusp_counts:
                self.cusp_counts[key]=cusp_triangles(self.sites[k],[self.walls[i] for i in bits(incident)],self.q)
            index+=self.cusp_counts[key]
            mask |= 1<<k
        if index%48:
            return index,"volume"
        forced=set()
        alternatives=[]
        for pair in pairs:
            options=[t for t in self.candidates_at_edge(pair) if t&mask==t]
            if not options:
                return index,"uncovered edge"
            if len(options)==1:
                forced.add(options[0])
            else:
                alternatives.append(options)
        if len(forced)>index//48 or any(sum(t in forced for t in options)>1 for options in alternatives):
            return index,"incompatible forced tiles"
        return index,None


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",choices=("all3","no3"))
    parser.add_argument("width",type=int)
    parser.add_argument("height",type=int)
    parser.add_argument("bend",type=int,help="maximum n_t-n_z of a hemispherical wall")
    parser.add_argument("output",type=Path)
    parser.add_argument("--origin",nargs=2,type=int,default=(0,0))
    parser.add_argument("--limit",type=int,help="stop after this many embedded polyhedra")
    parser.add_argument("--seed",type=int,default=0)
    parser.add_argument("--screen",action="store_true",help="count indices and save cases failing the volume or boundary-edge tests")
    args=parser.parse_args()
    require(args.limit is None or args.limit>0,"Expected a positive sample size")
    start=time.monotonic()
    search=CuspSearch(args.case,args.width,args.height,args.bend,args.origin)
    count=0
    indices=defaultdict(int)
    obstructions=defaultdict(int)
    complete=True
    seen=set()
    opener=gzip.open if args.output.suffix==".gz" else open
    with opener(args.output,"wt") as out:
        for chosen in search.solutions(args.seed):
            require(chosen not in seen,"Duplicate configuration")
            seen.add(chosen)
            if args.screen:
                index,failure=search.screen(chosen)
            else:
                record=search.record(chosen)
                index=record["index"]
                failure="volume" if index%48 else None
            indices[index]+=1
            if failure is not None:
                obstructions[failure]+=1
                if args.screen:
                    record=search.record(chosen)
                print(f"OBSTRUCTION at index {index}: {failure}",flush=True)
            if not args.screen or failure is not None:
                print(json.dumps(record,separators=(",",":")),file=out)
            count+=1
            if count==args.limit:
                complete=False
                break
    print(json.dumps(dict(case=args.case,width=args.width,height=args.height,origin=args.origin,bend=args.bend,seed=args.seed,
                          limit=args.limit,walls=len(search.walls),configurations=count,exhausted=complete,
                          indices=dict(sorted(indices.items())),screen=args.screen,obstructions=dict(obstructions),
                          seconds=round(time.monotonic()-start,3))))


if __name__=="__main__":
    main()
