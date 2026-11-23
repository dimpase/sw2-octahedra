# Exact enumeration in the SW2 tessellation

This experiment concerns Conjecture 7.11 of Allcock, Devlin, Felikson, Kontorovich and Whitehead. An SW2 polyhedron is a convex union of finitely many chambers in the reflection tessellation by the tetrahedron SW2. The conjecture concerns such polyhedra with all vertices ideal and all dihedral angles equal to π/2. A whole-octahedron tiling may have nonmatching faces. [Definitions and Proposition 7.10](https://arxiv.org/html/2609.05349v1#S7.SS4)

## How the counterexample was found

The exhaustive census used a standalone Julia program to grow convex unions of SW2 chambers. The four face reflections already specify the tessellation and the neighbour across each face. The program uses Julia’s standard library. We did not use ArithmeticVoronoi.jl or compute perfect forms for this search. [Chamber enumeration](search.jl)

The 528-chamber example came from a second search written in Python. We put one ideal vertex at infinity, fixed its four vertical walls and enumerated compatible hemispherical mirrors up to a curvature bound. For the rectangle [1,5] × [0,4] and maximum curvature 4 in the all-type-3 case, this gave 563 configurations. In one of them the boundary edges required 12 distinct regular octahedra although the volume allowed only 11. Examining those forced tiles gave the shorter proof using two overlapping octahedra. Both the search and the geometric verification use the Python standard library. [Cusp search](cusp_search.py) · [Proof and reproduction commands](COUNTEREXAMPLE.md)

GAP independently checked the chamber count. Given the SW2 Coxeter presentation and the 32 facet-reflection words from the geometry, GAP’s `Index(W,H)` returned 528. HAP’s contribution was a separate calculation of the Davis complex of the ambient SW2 group, with 14, 43, 48 and 18 cell orbits in dimensions 0 through 3. The cup-product comparison for the two antiprism groups used ordinary GAP and linear algebra over the field of two elements. [GAP index check](results/counterexample-528.g) · [Davis complex and cup-product calculations](coxeter_checks.g)

## Exact coordinates

Use future projective vectors (t,x,y,z) with Lorentz form −t²+x²+y²+z². The open Klein ball consists of negative lines, and its ideal boundary consists of null lines. The base chamber lies in x ≥ y ≥ z ≥ 0 and x+y+z ≤ t. In the paper's face order the outward normals are

    (0,-1,1,0)  (1,1,1,1)  (0,0,0,-1)  (0,0,-1,1).

Reflections 1 and 4 swap x,y and y,z respectively. Reflection 3 negates z. Reflection 2 sends v to v − ⟨v,n₂⟩n₂. The noncommuting orders are m₁₄=3 and m₄₃=m₃₂=4. The four vertices, ordered by omitted face, are

    (1,1,0,0)  (1,0,0,0)  (3,1,1,1)  (2,1,1,0).

The first vertex is ideal. The full orbit under reflections 1, 3, 4 about the second vertex gives the regular ideal octahedron, of volume 48 SW2 chambers. Julia uses integer matrices with checked coordinate arithmetic. The checker reconstructs these matrices independently by column swaps, sign changes and additions using Python integers.

## Search units

Lemma 7.9 says that the faces of an ideal right-angled SW2 polyhedron either all have type 3 or none has type 3. This gives two exhaustive runs.

| Run | Mirrors across which every candidate necessarily continues | Fundamental unit | SW2 chambers per unit |
| --- | --- | --- | ---: |
| `all3` | 1, 2, 4 | Right-angled triangular bipyramid | 12 |
| `no3` | 3 | Double of the tetrahedron across face 3 | 2 |
| `raw` or `raw-even` | None imposed | SW2 tetrahedron | 1 |

For each reduced run the forced generators generate a finite standard parabolic U. Sending the other generators to the identity defines a retraction W → U. Every Coxeter exponent between a retained and a deleted generator is even. The union of the U-images of SW2 is a fundamental chamber for the kernel. Its boundary reflections are conjugates under U of the remaining generators. Consequently the units form a face-to-face reflection tessellation, and every sought boundary mirror occurs in that tessellation. A search representative is taken up to translation by this unit reflection group. Different labelled embeddings can give isometric polyhedra. The raw run uses W itself and makes no appeal to Lemma 7.9. [SW2 geometry](https://arxiv.org/html/2609.05349v1#S7.SS4)

## Why the bounded search is exhaustive

1. Translate a target to contain the base unit. A chamber union with convex boundary is the intersection of the mirror halfspaces that contain it. To take the smallest such union containing a finite connected set, open precisely the mirrors crossed between adjacent units of that set. Flood from the root through these mirrors. The component is the intersection of the halfspaces containing the root for every unopened mirror. A generic segment from the root to an interior point of that intersection crosses only open mirrors. This proves both convexity and minimality. Stop when the number of units exceeds the requested bound.
2. At a finite unit vertex inspect the full finite tessellation of the spherical link. In an ideal right-angled target this point can lie in the interior, in a face, or in an edge. The selected incident units must therefore fill the whole sphere, a hemisphere, or a quadrant bounded by perpendicular mirrors. The implementation precomputes all such sectors in both orientations. For the optional angle extension it also includes sectors of angle π/4, using the exact sign and squared Lorentz inner product. The program checks that every non-right edge of the base unit has a finite endpoint. An edge with two ideal endpoints and right-angled unit sectors has only right or flat convex boundary angles.
3. At a defective link, branch on every inclusion-minimal admissible sector containing the current incident units and take mirror closure again. For any target containing the current state, at least one such branch lies within the target. A lower bound adds mandatory numbers of new units in disjoint regions of different vertex stars. In each star it counts only units outside all earlier charged regions and takes the minimum over feasible local choices. This bound can only discard supersets larger than the chamber limit.
4. When all the links pass, record the solid **and continue the search across each distinct boundary mirror**. Every proper convex chamber superset contains a face-neighbour across one of them. The forced closure across a given mirror is independent of where that mirror was first crossed. Thus the recursion reaches both minimal targets and their admissible supersets. Every branch strictly increases the count.
5. Duplicate elimination uses the full chamber adjacency with generator labels. From each possible root, breadth-first traversal records the incident labels at each visited unit in fixed generator order. Given the root and the Coxeter group, this reconstructs the chamber subset by multiplying words. The least rooted code therefore identifies exactly the translates under the unit reflection group. Choosing any one defective link remains exhaustive up to this equivalence. The raw cross-check compares the actual chamber sets up to W at the smaller bound.

## What a positive tiling certificate proves

The independent checker expands the recorded SW2 words, finds neighbours, tests connectedness, and verifies every exposed face is a globally supporting plane. This last test certifies convexity. A segment from an interior chamber to a point outside the union must leave across an exposed face, using a generic segment and closure for boundary cases. It then checks every extreme vertex is null, every ridge has the requested angle, every horospherical vertex link has Euclidean angle sum, and the boundary has Euler characteristic two. Isometry types are counted from the full vertex–edge–face incidence, with dihedral orders on edges. The least rooted flag code allows orientation reversal. Finite-volume hyperbolic Coxeter rigidity identifies these labelled combinatorial types with isometry types. [Rigidity and Coxeter presentations](https://arxiv.org/html/2609.05349v1#S2.SS1)

A candidate of index d can be tiled by regular ideal octahedra only if 48 divides d. For a candidate passing that test there is a **finite, complete placement search**, even when octahedra may be non-face-to-face.

* All six vertices of any such tile must be cusps of the candidate. In the Klein ball a point on the strictly convex unit sphere is an extreme point of any subset of the closed ball that contains it.
* Identify the sphere at infinity with the Riemann sphere. A regular ideal octahedron has vertices 0, ∞, 1, −1, i, −i. For each pair of possible opposite cusps and each choice of a third cusp, there is a unique Möbius map with these images of 0, ∞, 1. Evaluate it at the other three vertices using Gaussian integer pairs in homogeneous coordinates. Interchanging the poles or conjugating the standard octahedron preserves its vertex set. The checker keeps every six-point placement lying in the cusp set. Its convex hull is contained in the candidate.
* Search for d/48 placements with pairwise disjoint interiors. For two full-dimensional convex polytopes in projective Klein coordinates, a plane that weakly separates their vertices certifies disjoint interiors. To find one, testing planes through triples of vertices from the two tiles suffices. The cone of separating covectors is pointed because the vertex vectors span four dimensions. Any nonzero such cone has an extreme ray with three independent active constraints.
* A certificate lists the six cusp indices for each tile and an integer separating covector for every pair. A separate replay checks eight triangular supporting facets on each tile, all right dihedral angles and all separating inequalities. It checks both the SW2 chamber count and exact Euclidean volume in the Klein model by summing rational tetrahedral determinants. A finite closed union inside the candidate with this full volume covers its interior and hence its closure. This applies even when neighbouring octahedra meet along only part of a face.

The standalone cusp solver uses exact Klein volumes as weights in the disjoint-placement search. Consequently it can decide the tiling question for the convex hull of a supplied finite set of rational cusps without first computing a hyperbolic volume or identifying an ambient group.

The census proves the bounded statement. A 528-chamber example violates the unbounded conjecture. [Construction and two-edge proof](COUNTEREXAMPLE.md)

## Search from a cusp

Fix an ideal vertex at infinity in upper half-space. Four vertical facets give an integral rectangle. In the all-type-3 case use the Cartesian coordinates (x,y) in the boundary plane. In the case with no type-3 facets use (x+y,x−y). Write q for the Lorentz norm of the primitive facet normals, respectively 1 and 2. In these rectangle coordinates a hemispherical facet with primitive normal n has radius q/d, where d = nₜ−n_z is a positive integer. We enumerate integral normals up to a chosen d, retaining walls proved conjugate to a simple SW2 mirror by explicit reflections.

A possible hemisphere center must lie in the rectangle. To see this, take an edge of its projected polygon. For a vertical neighbouring facet the chord lies on a diameter and the center lies on the line of the chord. For another hemisphere the two circles are orthogonal. On the first circle, the arc exterior to the second disk lies on the same side of their common chord as the first center. These conditions, one per edge, put the center in the closed polygon.

The search extends the four selected vertical facets by completing cusp links. Two distinct facet normals in a right-angled chamber have Lorentz product zero or at most −q. Indeed their reflections in the chamber Coxeter system respectively commute or generate an infinite dihedral group, with outward normals in the chamber's simple root system. If two selected normals a,b are perpendicular, their walls meet in an edge. Orthogonal projection of an interior timelike vector onto the two walls increases it by positive multiples of a and b and preserves the strict inequalities for the other facets. Because all vertices are ideal, the two ideal endpoints of the full edge must be retained. At each endpoint the Euclidean link is a rectangle. Thus each incident side has an opposite parallel side, tangent to it at that cusp and perpendicular to the adjoining side. The program branches over all compatible choices in its finite mirror list. It rejects walls which cut off endpoints already required. A finished polyhedron has a connected edge graph, so following the edges and completing these rectangles reaches every facet. Runs with a configuration limit give samples whenever that limit interrupts the enumeration.

The SW2 index can also be computed from the cusp rectangles. For a primitive future null vector v at a cusp, the two pairs of opposite primitive outward normals sum to mv and nv, respectively, with positive integers m and n. An integral SW2 isometry takes v to (1,0,0,1). In these coordinates the rectangle has Euclidean side lengths m/√q and n/√q. The SW2 horospherical triangle has area 1/4, as seen from its vertices (0,0), (1,0), (1/2,1/2) up to symmetry. Every SW2 tetrahedron has exactly one ideal vertex. Consequently the total chamber count is the sum of 4mn/q over the cusps. This agreed with the direct chamber count for all 106 previously known isometry types through 240 and in the all-type-3 extension through 288.

The edge screen tests whether each boundary edge can occur in an octahedron with vertices among the target cusps, and whether the resulting forced tiles meet the volume and edge-count constraints. For a full tiling decision the checker rebuilds the hull, enumerates every cusp-compatible octahedron, and searches for disjoint tiles with the exact total Klein volume. It also uses the fact that each boundary edge requires exactly one tile. The counterexample has a separate reconstruction and a proof using just two edges. [Proof](COUNTEREXAMPLE.md)

## An index restriction at every size

In the all-type-3 boundary case let L be the reflection group of the quarter-octahedral triangular bipyramid and H the reflection group of the candidate. Then [L:H] is its number of quarters. Let S be the stabilizer in L of a finite apex of the bipyramid. The three walls there are mutually perpendicular, so S is an elementary abelian group of order eight.

Let S act by right multiplication on the right cosets Hg of H in L. The stabilizer of Hg is S ∩ g⁻¹Hg, conjugate to a subgroup of the H-stabilizer of the point g(apex). Because H has an ideal right-angled fundamental polyhedron, its point stabilizers **in hyperbolic space** have orders 1, 2 or 4, according as the point lies in a chamber interior, a face, or an edge. Thus each S-orbit has 2, 4 or 8 cosets. Their sum [L:H] is even. In SW2 units this proves **24 divides the chamber count** throughout the all-type-3 case. The whole-octahedron conjecture would imply divisibility by 48 and an actual geometric tiling.

## Degree-one cup products

For a right-angled Coxeter group the mod-2 cohomology is the face ring of the clique complex of its commuting graph. In degree two, a basis consists of vertex squares and products indexed by edges. In multiplication by a class supported on A, the square coordinates force every kernel element to vanish on A. The edge coordinates then force it to vanish on every neighbour of A. Thus the rank equals the size of the closed neighbourhood of A. The GAP script counts this rank for every class and checks the calculation against literal matrices over GF(2) for the zero class, every singleton and every pair, and a class of weight five. [Davis–Januszkiewicz, Theorem 4.11](https://people.math.osu.edu/davis.12/old_papers/DJ_toric.dmj.pdf)

Both face-adjacency graphs have 14 vertices, 24 edges and no triangles. In positive degree k the monomials supported on one vertex give 14 basis elements, and the monomials supported on an edge give 24(k−1) more. This proves the equality of all the mod-2 Betti numbers used in the comparison.
