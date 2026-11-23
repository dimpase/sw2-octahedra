"""Compare every chamber embedding in a raw SW2 run with the two reduced searches."""
import argparse
from pathlib import Path
from verify import ambient_code, load_input, require


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw",type=Path)
    parser.add_argument("reduced",type=Path,nargs="+")
    args=parser.parse_args()
    raw=load_input(args.raw)
    require(set(raw["visited"])=={"raw"},"Expected a raw right-angled run")
    bound=raw["bound"]
    runs=[load_input(path) for path in args.reduced]
    require(all(run["bound"]>=bound for run in runs),"Reduced bound is too small")
    require(set().union(*(r["visited"] for r in runs))=={"all3","no3"},"Both boundary cases are required")
    def classes(records):
        return {ambient_code(record["words"]) for record in records if record["index"]<=bound}
    first=classes(raw["polyhedra"])
    second=classes(p for run in runs for p in run["polyhedra"])
    require(first==second,"The two enumerations found different embeddings")
    print(f"Raw SW2 and both reduced searches agree through {bound} chambers, "
          f"including all {len(first)} embeddings up to SW2.")


if __name__=="__main__":
    main()
