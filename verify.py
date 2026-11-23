"""Independently replay SW2 chamber words, check the solid, and test all whole-octahedron tilings."""
import argparse
from collections import Counter, defaultdict
from functools import cache
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path
import gzip
import hashlib
import json
import tomllib


IDENTITY = tuple(tuple(int(i == j) for i in range(4)) for j in range(4))
NORMALS = ((0, -1, 1, 0), (1, 1, 1, 1), (0, 0, 0, -1), (0, 0, -1, 1))
CORNERS = ((1, 1, 0, 0), (1, 0, 0, 0), (3, 1, 1, 1), (2, 1, 1, 0))


def require(test, message):
    if not test:
        raise ValueError(message)


@cache
def cross(g, face):
    t, x, y, z = g
    if face == 1:
        return t, y, x, z
    if face == 4:
        return t, x, z, y
    if face == 3:
        return t, x, y, tuple(-a for a in z)
    if face != 2:
        raise ValueError("Unknown face")
    n = tuple(map(sum, zip(*g)))
    return tuple(tuple(a+b if j == 0 else a-b for a,b in zip(col,n)) for j,col in enumerate(g))


word_frames={"":IDENTITY}


def evaluate(word):
    length=len(word)
    while word[:length] not in word_frames:
        length-=1
    frame=word_frames[word[:length]]
    for i in range(length,len(word)):
        frame=cross(frame,int(word[i]))
        word_frames[word[:i+1]]=frame
    return frame


def act(g, v):
    return tuple(sum(a*b for a,b in zip(row,v)) for row in zip(*g))


def dot(a, b):
    return sum(x*y for x,y in zip(a,b))


def lor(a, b):
    return dot(a,b)-2*a[0]*b[0]


def primitive(x):
    return tuple(a//gcd(*x) for a in x)


def det3(a,b,c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))


@cache
def through(a,b,c):
    # Euclidean covector of a projective plane. The zero covector means collinear.
    return tuple((-1)**i * det3(*(v[:i]+v[i+1:] for v in (a,b,c))) for i in range(4))


def check_boundary(record):
    words = record["words"]
    cells = {evaluate(word) for word in words}
    require(len(words) == len(cells) == record["index"], "Repeated chambers or incorrect index")
    reached = {IDENTITY}
    queue = [IDENTITY]
    points, walls, labels = set(), set(), set()
    require(IDENTITY in cells, "Missing base chamber")
    for g in queue:
        points.update(act(g,p) for p in CORNERS)
        for i in range(1,5):
            h = cross(g,i)
            if h not in cells:
                walls.add(act(g,NORMALS[i-1]))
                labels.add(i)
            elif h not in reached:
                queue.append(h)
                reached.add(h)
    require(cells == reached, "Disconnected chamber union")
    require(record["case"] in ("raw","raw-even") or
            (record["case"] == "all3" and labels == {3}) or
            (record["case"] == "no3" and 3 not in labels), "Incorrect boundary labels")
    facets_by_point = defaultdict(set)
    for n in walls:
        for p in points:
            side = lor(n,p)
            require(side <= 0, "Boundary face fails its global supporting-halfspace test")
            if side == 0:
                facets_by_point[p].add(n)
    cusps = sorted(p for p,fs in facets_by_point.items() if len(fs)>=3 and
                   any(any(through(*tri)) for tri in combinations(fs,3)))
    require(all(lor(v,v)==0 for v in cusps), "Nonideal vertex")
    faces = {n:frozenset(i for i,v in enumerate(cusps) if n in facets_by_point[v]) for n in walls}
    require(all(len(f)>=3 for f in faces.values()), "Degenerate facet")
    edges = set()
    angles = {}
    for a,b in combinations(walls,2):
        edge = faces[a] & faces[b]
        if len(edge)>=2:
            require(len(edge)==2,"Nonpolygonal edge")
            d=lor(a,b)
            if d==0:
                angles[edge]=2
            else:
                require(record["case"]=="raw-even" and d<0 and 2*d*d==lor(a,a)*lor(b,b),
                        "Boundary angle is outside requested family")
                angles[edge]=4
            edges.add(edge)
    require(len(cusps)-len(edges)+len(walls)==2, "Boundary is not a polyhedral sphere")
    for i,v in enumerate(cusps):
        orders=[m for e,m in angles.items() if i in e]
        require(len(orders)==len(facets_by_point[v]) and sum(Fraction(1,m) for m in orders)==len(orders)-2,
                "Cusp link fails the Euclidean polygon angle sum")
    return cusps, list(faces.values()),edges,angles


def incidence_code(faces,edges,angles):
    # Flags and their three adjacency involutions. The least rooted code is
    # a complete combinatorial isomorphism invariant, including reflections.
    flags = [(v,e,f) for f in faces for e in edges if e<=f for v in e]
    groups = defaultdict(list)
    for i,flag in enumerate(flags):
        for j in range(3):
            groups[j,flag[:j]+flag[j+1:]].append(i)
    adj = [[] for _ in flags]
    for i,flag in enumerate(flags):
        for j in range(3):
            pair=groups[j,flag[:j]+flag[j+1:]]
            require(len(pair)==2,"Invalid flag adjacency")
            adj[i].append(pair[0] if pair[1]==i else pair[1])
    best=None
    for root in range(len(flags)):
        ids={root:0}
        todo=[root]
        code=[]
        for i in todo:
            code.append(angles[flags[i][1]])
            for j in adj[i]:
                if j not in ids:
                    ids[j]=len(todo)
                    todo.append(j)
                code.append(ids[j])
        require(len(todo)==len(flags),"Disconnected boundary flags")
        if best is None or code<best:
            best=code
    return tuple(best)


# A cusp on the Riemann sphere is (x+iy)/q. Use [1:0] for infinity.
def cusp_coord(v):
    t,x,y,z=v
    return primitive((x,y,t-z)) if t!=z else (1,0,0)


def mul(a,b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def add(a,b):
    return a[0]+b[0],a[1]+b[1]


def homogeneous(a,b):
    if b==(0,0):
        require(a!=(0,0),"Undefined projective point")
        return (1,0,0)
    return primitive((*mul(a,(b[0],-b[1])), b[0]**2+b[1]**2))


def bracket(a,b):
    return a[0]*b[2]-b[0]*a[2], a[1]*b[2]-b[1]*a[2]


def octahedra(vertices):
    coords=[cusp_coord(v) for v in vertices]
    ids={z:i for i,z in enumerate(coords)}
    require(len(ids)==len(vertices),"Repeated cusp")
    candidates=set()
    for i,j in combinations(range(len(coords)),2):
        a,b=coords[i],coords[j]
        for k,c in enumerate(coords):
            if k in (i,j):
                continue
            # Map 0,infinity,1 to a,b,c. Compute the images of -1,i,-i.
            p,q=bracket(a,c),bracket(c,b)
            bn,bd=mul(p,b[:2]),mul(p,(b[2],0))
            an,ad=mul(q,a[:2]),mul(q,(a[2],0))
            tile=[i,j,k]
            for z in ((-1,0),(0,1),(0,-1)):
                w=homogeneous(add(mul(z,bn),an),add(mul(z,bd),ad))
                if w not in ids:
                    break
                tile.append(ids[w])
            if len(tile)==6:
                candidates.add(tuple(sorted(tile)))
    return sorted(candidates)


def separator(a,b):
    # A nonzero weakly separating covector exists iff one is orthogonal to
    # three independent vertex vectors. The cone of separators is pointed.
    for triple in combinations(sorted(set(a)|set(b)),3):
        n=through(*triple)
        if not any(n):
            continue
        left=[dot(n,v) for v in a]
        right=[dot(n,v) for v in b]
        if max(left)<=0<=min(right):
            return primitive(n)
        if max(right)<=0<=min(left):
            return primitive(tuple(-x for x in n))
    return None


def tile_test(vertices,index):
    if index%48:
        return "volume obstruction", [], [], None
    placements=octahedra(vertices)
    count=index//48
    @cache
    def sep(i,j):
        return separator([vertices[k] for k in placements[i]],[vertices[k] for k in placements[j]])
    def select(options,needed):
        if needed==0:
            return []
        for n,i in enumerate(options):
            if len(options)-n<needed:
                break
            rest=[j for j in options[n+1:] if sep(i,j) is not None] if needed>1 else []
            chosen=select(rest,needed-1)
            if chosen is not None:
                return [i,*chosen]
        return None
    chosen=select(list(range(len(placements))),count)
    if chosen is None:
        return "no octahedral tiling", [], [], len(placements)
    return ("tiled",[placements[i] for i in chosen],
            [sep(i,j) for i,j in combinations(chosen,2)],len(placements))


def supporting_planes(vs):
    """Facets of a full-dimensional convex hull, using exact projective triples."""
    facets={}
    for tri in combinations(vs,3):
        cov=through(*tri)
        if not any(cov):
            continue
        vals=[dot(cov,v) for v in vs]
        if min(vals)<0<max(vals) or all(x==0 for x in vals):
            continue
        f=frozenset(i for i,val in enumerate(vals) if val==0)
        cov=primitive(cov)
        facets[f]=(-cov[0],*cov[1:])
    require(len(facets)>=4,"The cusp hull is not three-dimensional")
    return facets


def klein_volume(vs,faces):
    """Six times Euclidean volume after normalising to t=1, invariant under subdivision."""
    center=tuple(map(sum,zip(*vs)))
    volume=Fraction(0)
    for face in faces:
        if len(face)==3:
            polygon=sorted(face)
        else:
            near={i:[] for i in face}
            for other in faces:
                common=face & other
                if len(common)==2:
                    a,b=common
                    near[a].append(b)
                    near[b].append(a)
            require(all(len(ns)==2 for ns in near.values()),"Degenerate polygon")
            polygon=[min(face)]
            polygon.append(min(near[polygon[0]]))
            while True:
                a,b=polygon[-2:]
                c=next(i for i in near[b] if i!=a)
                if c==polygon[0]:
                    break
                require(c not in polygon,"Invalid polygon cycle")
                polygon.append(c)
            require(len(polygon)==len(face),"Disconnected polygon")
        a=vs[polygon[0]]
        for i,j in zip(polygon[1:],polygon[2:]):
            b,c=vs[i],vs[j]
            volume+=Fraction(abs(dot(center,through(a,b,c))),center[0]*a[0]*b[0]*c[0])
    return volume


def verify_witness(cusps,record,faces):
    # Independent of the Mobius placement search. Construct each tile's
    # supporting planes directly from triples, then check its right angles.
    for tile in record["tiles"]:
        require(len(tile)==6 and all(type(i) is int and 0<=i<len(cusps) for i in tile),
                "Expected six cusp indices in the input range")
    tiles=[[cusps[i] for i in tile] for tile in record["tiles"]]
    if "index" in record:
        require(len(tiles)*48==record["index"],"SW2 volume mismatch")
    tile_volume=Fraction(0)
    for vs in tiles:
        require(len(set(vs))==6 and all(lor(v,v)==0 for v in vs),"Invalid tile vertices")
        facets=supporting_planes(vs)
        require(len(facets)==8 and all(len(f)==3 for f in facets),"Tile is not an octahedron")
        for a,b in combinations(facets,2):
            if len(a & b)==2:
                require(lor(facets[a],facets[b])==0,"Tile is not right angled")
        tile_volume+=klein_volume(vs,facets)
    require(tile_volume==klein_volume(cusps,faces),"Klein volumes do not cover the cusp hull")
    pairs=list(combinations(tiles,2))
    require(len(pairs)==len(record["separators"]),"Missing separator")
    for (left,right),n in zip(pairs,record["separators"]):
        require(len(n)==4 and all(type(x) is int for x in n),"Expected four integer plane coordinates")
        require(any(n) and all(dot(n,v)<=0 for v in left) and all(dot(n,v)>=0 for v in right),
                "Tile interiors overlap or separator is invalid")


def read_file(path):
    opener=gzip.open if path.suffix==".gz" else open
    with opener(path,"rb") as f:
        return f.read()


def load_input(path):
    return tomllib.loads(read_file(path).decode())


def ambient_code(words):
    """Least chamber-root code, up to SW2 itself (all four panel labels fixed)."""
    cells={evaluate(w) for w in words}
    masks={g:sum(1<<j for j in range(4) if cross(g,j+1) in cells) for g in cells}
    limit=min(masks.values())
    best=None
    for root in cells:
        if masks[root]!=limit:
            continue
        seen={root}
        todo=[root]
        code=[]
        for g in todo:
            code.append(masks[g])
            for i in range(4):
                if masks[g] & (1<<i):
                    h=cross(g,i+1)
                    if h not in seen:
                        seen.add(h)
                        todo.append(h)
        candidate=bytes(code)
        if best is None or candidate<best:
            best=candidate
    return best


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("search",nargs="+",type=Path)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write",type=Path,help="find tilings and save witnesses")
    mode.add_argument("--check",type=Path,help="replay previously saved witnesses")
    args=parser.parse_args()
    saved=json.loads(read_file(args.check)) if args.check else None
    results=[]
    shapes={}
    counts=Counter()
    for path in args.search:
        for group in load_input(path)["polyhedra"]:
            cusps,faces,edges,angles=check_boundary(group)
            code=incidence_code(faces,edges,angles)
            tag=hashlib.sha256(repr(code).encode()).hexdigest()
            shapes.setdefault(group["index"],set()).add(code)
            signature=hashlib.sha256("\n".join(sorted(group["words"])).encode()).hexdigest()
            if saved is None:
                status,tiles,planes,nplace=(tile_test(cusps,group["index"]) if 4 not in angles.values()
                                            else ("angle extension",[],[],None))
                result=dict(index=group["index"],case=group["case"],cusps=cusps,shape=tag,words_sha256=signature,
                            status=status,tiles=tiles,separators=planes,possible_octahedra=nplace)
            else:
                result=saved[len(results)]
                require((result["index"],result["case"],result["shape"],result["words_sha256"])==
                        (group["index"],group["case"],tag,signature) and result["cusps"]==list(map(list,cusps)),
                        "Witness does not match the input chambers")
            if result["status"]=="tiled":
                verify_witness(cusps,result,faces)
            elif result["status"]=="angle extension":
                require(4 in angles.values(),"Right-angled candidate classified as an angle extension")
            else:
                require(tile_test(cusps,group["index"])[0]==result["status"],"Incorrect negative result")
            counts[group["index"],result["status"]]+=1
            results.append(result)
    if saved is not None:
        require(len(results)==len(saved),"Extra witnesses")
    else:
        encoded=(json.dumps(results,separators=(",",":"))+"\n").encode()
        args.write.write_bytes(gzip.compress(encoded,mtime=0) if args.write.suffix==".gz" else encoded)
    for index,distinct in sorted(shapes.items()):
        noun="isometry type" if len(distinct)==1 else "isometry types"
        print(f"SW2 index {index}: {len(distinct)} {noun}; "
              + ", ".join(f"{n} {status}" for (i,status),n in sorted(counts.items()) if i==index))


if __name__=="__main__":
    main()
