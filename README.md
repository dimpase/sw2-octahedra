# Octahedral tilings in the SW2 tessellation

A **528-chamber counterexample** to Conjecture 7.11 of Allcock, Devlin, Felikson, Kontorovich and Whitehead has 30 ideal vertices and 60 right-angled edges. Two boundary edges force octahedra whose interiors overlap. The construction uses the symmetries of a square and seven rational points in the plane. [Construction and proof](COUNTEREXAMPLE.md) · [Arithmetic Polyhedra, §7.4](https://arxiv.org/html/2609.05349v1#S7.SS4)

```sh
python3 counterexample.py
```

## Smallest sizes

The exhaustive search found **29 isometry types** with **at most 252 chambers**. All 29 admit tilings by whole regular ideal octahedra. Each has a tiling certificate.

| SW2 chambers | Whole octahedra per polyhedron | Isometry types | Representatives checked |
| ---: | ---: | ---: | ---: |
| 48 | 1 | 1 | 5 |
| 96 | 2 | 1 | 18 |
| 144 | 3 | 2 | 35 |
| 192 | 4 | 7 | 164 |
| 240 | 5 | 18 | 537 |

The table includes every index that occurred. The two boundary cases in Lemma 7.9 use different fundamental units. The last column counts distinct chamber unions up to their respective unit reflection groups. The isometry column compares their complete face incidences. An enumeration directly in the SW2 tetrahedron recovered the same 16 embeddings up to SW2 through index 144.

In the case where **all boundary faces have type 3**, the exhaustive search reached index **288**, finding 77 isometry types at that index. Each admits a tiling by six octahedra. In this case an orbit-stabilizer argument shows that every candidate, at any index, has an even number of quarter-octahedra. In the other boundary case a run through **252** visited 8,711,718 convex configurations in 18 minutes 20 seconds, with a peak resident set of 10.2 GiB. The first possible index without an exhaustive result in both cases is **264**. [Proof and enumeration method](METHOD.md)

## Run a new search

Run these commands from this directory with Julia 1.10+ and Python 3.11+.

```sh
julia --startup-file=no search.jl 96 example.toml
python3 verify.py example.toml --write example-tilings.json
```

The first command enumerates new chamber unions. The second independently reconstructs their boundary geometry, finds regular octahedra from their ideal vertices, and certifies disjoint interiors and full coverage. For this example the output gives the octahedron and its double across a face. The programs use only the standard libraries of Julia and Python. A fresh index-240 search took about ten minutes on the test host.

## Test a polyhedron given by its cusps

Supply the ideal vertices as a JSON array of homogeneous integer vectors `[t,x,y,z]`, with positive t and t²=x²+y²+z². The program decides whether their convex hull admits a tiling by regular ideal octahedra. It works independently of the SW2 search, allows non-face-to-face tilings, and returns the tiles as zero-based indices in the input list. This example uses the six coordinate directions.

```sh
python3 tile_cusps.py <<'JSON'
[[1,1,0,0],[1,-1,0,0],[1,0,1,0],[1,0,-1,0],[1,0,0,1],[1,0,0,-1]]
JSON
```

This gives one tile. Adding `[5,3,4,0]` gives a larger hull containing that octahedron. The program finds the one possible placement and proves it leaves part of the hull uncovered. [Placement and volume tests](METHOD.md#what-a-positive-tiling-certificate-proves)

## Replay the results

```sh
python3 verify.py results/right-angled-240.toml.gz --check results/right-angled-240.json.gz
python3 verify.py results/no3-252.toml.gz --check results/no3-252.json.gz
python3 verify.py results/all3-288.toml.gz --check results/all3-288.json.gz
python3 compare_runs.py results/raw-144.toml.gz results/right-angled-240.toml.gz
python3 test.py
```

The files in `results/` hold the chamber words and certificates. The search takes its data from the four mirror generators.

## Larger examples

The cusp search fixes the rectangle at infinity and enumerates compatible hemispherical mirrors up to a specified curvature. It can find examples with overlapping forced tiles. This command reproduces the 563 configurations that include the symmetric counterexample.

```sh
python3 cusp_search.py all3 4 4 4 rectangles.jsonl.gz --origin 1 0
python3 check_cusps.py rectangles.jsonl.gz --write rectangle-tilings.jsonl.gz
```

For a faster obstruction test add `--screen` to the first command. It writes just the configurations that fail a necessary tiling condition. The combined search tested 5,422 marked configurations with **2,735 distinct isometry types** and reached 2,304 chambers. It found **three untileable types**, at 528, 672 and 720 chambers. The other 2,732 types have tiling certificates. The search bounds apply to the specified cusp rectangles and curvatures. [Search method](METHOD.md#search-from-a-cusp) · [Run parameters and measured results](results/larger-runs.json)

```sh
python3 check_cusps.py results/larger-cusps.jsonl.gz --check results/larger-certificates.jsonl.gz
python3 verify.py results/other-counterexamples.toml.gz --check results/other-counterexamples.json.gz
```

## Allowing angles of π/4

```sh
julia --startup-file=no search.jl 144 angles.toml raw-even
python3 verify.py angles.toml --write angle-results.json
```

Within SW2, the search found **46 isometry types** with angles drawn from {π/2, π/4} at index at most 144. Four are right-angled and 42 have an angle of π/4. The smallest new example is an ideal tetrahedron made of 12 SW2 chambers. Its cusp triangles have angles π/2, π/4, π/4. These 42 examples fall outside the hypotheses of Conjecture 7.11. The broader classification across ambient reflection groups is among the questions in [§8 of the paper](https://arxiv.org/html/2609.05349v1#S8).

This command repeats the tiling test using only the cusp coordinates for one representative of each geometric type.

```sh
python3 tile_cusps.py --census results/right-angled-240.json.gz results/all3-288.json.gz results/even-144.json.gz
```

## GAP and HAP

From this directory, run `gap -q -b coxeter_checks.g` in the shell or `Read("coxeter_checks.g");` at the GAP prompt. The script loads HAP for the Davis complex. If HAP is unavailable, it reports the skipped check and still runs the subgroup and cup-product checks.

With GAP.jl, the same script can be run from Julia.

```julia
using GAP
GAP.Globals.Read(GAP.GapObj(abspath("coxeter_checks.g")))
```

GAP obtains subgroup indices 12, 48 and 96 for the quarter-octahedron, whole octahedron and face double. HAP's Davis complex for SW2 has cell-orbit counts [14,43,48,18].

Doubling the square antiprism across a square face gives an arithmetic reflection group. The hexagonal antiprism gives a nonarithmetic one. Their mod-2 Betti numbers agree in all degrees. GAP found **six** degree-one classes for which multiplication from H¹ to H² has rank five in the first group and **zero** in the second. The counts run over all 16,384 degree-one classes, so the distinction is independent of a choice of basis. [Arithmeticity of the antiprisms](https://arxiv.org/html/2609.05349v1#S1.SS5) · [Cup-product calculation](METHOD.md#degree-one-cup-products)
