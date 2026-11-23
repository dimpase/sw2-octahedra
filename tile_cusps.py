"""Decide whether the convex hull of given rational ideal points is tiled by regular ideal octahedra."""
import argparse
from fractions import Fraction
from functools import cache
from itertools import combinations
import json
from pathlib import Path
import sys

from verify import klein_volume, lor, octahedra, primitive, read_file, require, separator, supporting_planes, verify_witness


def decide(points):
    require(all(len(p)==4 and all(type(x) is int for x in p) and p[0]>0 and lor(p,p)==0
                for p in points),"Expected future null integer vectors [t,x,y,z]")
    vs=[primitive(tuple(p)) for p in points]
    require(len(vs)==len(set(vs)),"Duplicate projective cusp")
    faces=supporting_planes(vs)
    total=klein_volume(vs,faces)
    candidates=octahedra(vs)
    shapes=[[vs[i] for i in c] for c in candidates]
    weights=[klein_volume(s,supporting_planes(s)) for s in shapes]

    @cache
    def separate(i,j):
        return separator(shapes[i],shapes[j])

    def cover(options,remaining):
        if remaining==0:
            return []
        if sum((weights[i] for i in options),Fraction(0))<remaining:
            return None
        for pos,i in enumerate(options):
            if weights[i]>remaining:
                continue
            rest=[j for j in options[pos+1:] if separate(i,j) is not None]
            solution=cover(rest,remaining-weights[i])
            if solution is not None:
                return [i,*solution]
        return None

    chosen=cover(list(range(len(candidates))),total)
    result=dict(status="no octahedral tiling" if chosen is None else "tiled",cusps=vs,
                possible_octahedra=len(candidates),six_times_klein_volume=str(total))
    if chosen is not None:
        result["tiles"]=[candidates[i] for i in chosen]
        result["separators"]=[separate(i,j) for i,j in combinations(chosen,2)]
        verify_witness(vs,result,faces)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cusps",nargs="?",type=Path,help="JSON array of [t,x,y,z] vectors, default stdin")
    parser.add_argument("--census",type=Path,nargs="+",help="retest one cusp hull per isometry type in the given witness files")
    args=parser.parse_args()
    if args.census:
        if args.cusps:
            parser.error("Supply a cusp array or a census")
        census={}
        for path in args.census:
            for record in json.loads(read_file(path)):
                census.setdefault(record["shape"],record)
        tiled=0
        for record in census.values():
            result=decide(record["cusps"])
            expected=record["status"]=="tiled"
            require((result["status"]=="tiled")==expected,"Cusp-only result disagrees with the census")
            tiled+=expected
        print(f"Cusp-only tests on {len(census)} isometry types: {tiled} tiled, {len(census)-tiled} untileable.")
    else:
        points=json.loads(args.cusps.read_text() if args.cusps else sys.stdin.read())
        print(json.dumps(decide(points),indent=2))
