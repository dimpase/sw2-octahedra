"""Check cusp-search output and decide the octahedral tiling question once per isometry type."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import cache
import gzip
import hashlib
from itertools import combinations
import json
from pathlib import Path

from cusp_search import bits, chamber_index, simple_mirror
from verify import (dot, incidence_code, klein_volume, lor, octahedra, primitive, require, separator,
                    supporting_planes, through, verify_witness)


def hull_faces(points):
    """Incremental exact 3D convex hull in homogeneous Klein coordinates."""
    a,b,c=points[:3]
    fourth=next((i for i in range(3,len(points)) if dot(through(a,b,c),points[i])!=0),None)
    require(fourth is not None,"Cusp hull is flat")
    initial=(0,1,2,fourth)
    center=tuple(map(sum,zip(*(points[i] for i in initial))))
    def oriented(tri):
        cov=through(*(points[i] for i in tri))
        require(dot(cov,center)!=0,"Degenerate hull triangle")
        if dot(cov,center)>0:
            cov=tuple(-x for x in cov)
        return primitive(cov)
    triangles={t:oriented(t) for t in combinations(initial,3)}
    for i,point in enumerate(points):
        if i in initial:
            continue
        visible=[tri for tri,cov in triangles.items() if dot(cov,point)>0]
        require(visible,"Repeated or nonextreme cusp")
        horizon=set()
        for tri in visible:
            del triangles[tri]
            horizon.symmetric_difference_update(combinations(tri,2))
        for a,b in horizon:
            tri=tuple(sorted((a,b,i)))
            triangles[tri]=oriented(tri)
    return {frozenset(i for i,v in enumerate(points) if dot(n,v)==0) for n in set(triangles.values())}


def boundary(cusps,walls):
    require(len(set(cusps))==len(cusps) and all(v[0]>0 and lor(v,v)==0 for v in cusps),"Invalid cusps")
    faces=[]
    for n in walls:
        require(simple_mirror(n),"Non-SW2 facet")
        vals=[lor(n,v) for v in cusps]
        require(max(vals)==0,"Nonsupporting facet")
        f=frozenset(i for i,value in enumerate(vals) if value==0)
        require(len(f)>=3,"Degenerate facet")
        faces.append(f)
    require(len(set(faces))==len(faces),"Repeated facets")
    edges={}
    for i,j in combinations(range(len(walls)),2):
        common=faces[i]&faces[j]
        if len(common)>=2:
            require(len(common)==2 and lor(walls[i],walls[j])==0,"Non-right edge")
            edges[common]=2
    require(len(cusps)-len(edges)+len(faces)==2,"Incorrect boundary Euler characteristic")
    require(all(sum(i in e for e in edges)==4 for i in range(len(cusps))),"Nonrectangular cusp")
    code=incidence_code(faces,edges,edges)
    return faces,edges,hashlib.sha256(repr(code).encode()).hexdigest()


def tile(cusps,faces,edges,index):
    if index%48:
        return dict(status="volume obstruction",index=index)
    placements=octahedra(cusps)
    vertices=[[cusps[i] for i in indices] for indices in placements]
    facets=[supporting_planes(vs) for vs in vertices]
    volumes=[klein_volume(vs,f) for vs,f in zip(vertices,facets)]
    total=klein_volume(cusps,faces)
    covectors=[]
    for vs,fs in zip(vertices,facets):
        oriented=[]
        for face,n in fs.items():
            cov=(-n[0],*n[1:])
            if dot(cov,next(v for i,v in enumerate(vs) if i not in face))>0:
                cov=tuple(-x for x in cov)
            oriented.append(cov)
        covectors.append(oriented)

    @cache
    def separate(i,j):
        for cov in covectors[i]:
            if all(dot(cov,v)>=0 for v in vertices[j]):
                return cov
        for cov in covectors[j]:
            if all(dot(cov,v)>=0 for v in vertices[i]):
                return tuple(-x for x in cov)
        return separator(vertices[i],vertices[j])

    # Each target edge is a full ideal geodesic and has angle pi/2.
    # A tiling contains exactly one octahedral edge on it.
    around=[sum(1<<k for k,p in enumerate(placements) if e.issubset(p)) for e in edges]
    def solve(options,remaining,required):
        if remaining==0:
            return [] if not required else None
        if not options:
            return None
        if required:
            choices=min((m&options for m in required),key=lambda m:(m.bit_count(),m))
        else:
            if sum((volumes[k] for k in bits(options)),Fraction(0))<remaining:
                return None
            choices=options
        for k in bits(choices):
            options &= ~(1<<k)
            if volumes[k]>remaining:
                continue
            rest=sum(1<<j for j in bits(options) if separate(min(k,j),max(k,j)) is not None)
            ans=solve(rest,remaining-volumes[k],[m for m in required if not m&(1<<k)])
            if ans is not None:
                return [k]+ans
        return None
    answer=solve((1<<len(placements))-1,total,around)
    result=dict(status="no octahedral tiling" if answer is None else "tiled",possible_octahedra=len(placements),index=index)
    if answer is not None:
        answer.sort()
        result.update(tiles=[placements[i] for i in answer],separators=[separate(i,j) for i,j in combinations(answer,2)])
        verify_witness(cusps,result,faces)
    return result


def records(path):
    opener=gzip.open if path.suffix==".gz" else open
    with opener(path,"rt") as source:
        for line in source:
            yield json.loads(line)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",nargs="+",type=Path)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write",type=Path)
    group.add_argument("--check",type=Path)
    args=parser.parse_args()
    known={}
    if args.check:
        for r in records(args.check):
            require(r["shape"] not in known,"Repeated certificate")
            known[r["shape"]]=r
    seen=set()
    sizes=Counter()
    checked=0
    output=None
    if args.write:
        opener=gzip.open if args.write.suffix==".gz" else open
        output=opener(args.write,"wt")
    try:
        for path in args.input:
            for record in records(path):
                checked+=1
                cusps=list(map(tuple,record["cusps"]))
                walls=list(map(tuple,record["walls"]))
                faces,edges,shape=boundary(cusps,walls)
                if shape in seen:
                    continue
                seen.add(shape)
                # Rebuild the full hull from the point set independently.
                require(hull_faces(cusps)==set(faces),"Facet list differs from cusp hull")
                index=chamber_index(cusps,walls)
                require(record.get("index",index)==index,"Incorrect cusp index")
                if args.check:
                    result=known[shape]
                    require(result["cusps"]==record["cusps"],"Different certificate representative")
                    require(result["index"]==index,"Incorrect certificate index")
                    if result["status"]=="tiled":
                        verify_witness(cusps,result,faces)
                    else:
                        require(tile(cusps,faces,edges,index)["status"]==result["status"],"Wrong negative result")
                else:
                    result=tile(cusps,faces,edges,index)
                    result.update(shape=shape,cusps=cusps)
                    print(json.dumps(result,separators=(",",":")),file=output,flush=True)
                if result["status"]!="tiled":
                    print(f"UNTILEABLE {shape} in {path} at record {checked}",flush=True)
                sizes[48*len(result["tiles"]) if result["status"]=="tiled" else "untileable"]+=1
                if len(seen)%100==0:
                    print(f"{checked} inputs, {len(seen)} isometry types checked",flush=True)
                if len(seen)%64==0:
                    through.cache_clear()
    finally:
        if output:
            output.close()
    if args.check:
        require(seen==set(known),"Extra or missing certificates")
    print(f"{checked} inputs, {len(seen)} isometry types. Counts by SW2 index: {dict(sorted(sizes.items(),key=lambda t:str(t[0]).zfill(9)))}")


if __name__=="__main__":
    main()
