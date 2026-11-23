# Independent GAP checks for SW2 and the two antiprisms.
#(function()
#local m, relators, F, W, r, walls, kind, b, outside, doubled, counts, C, M;
m := [[1,2,2],[1,3,2],[1,4,3],[2,3,4],[2,4,2],[3,4,4]];
relators := function(gens,labels)
    return Concatenation(List(gens,x->x^2),
        List(Filtered(m,e->e[1] in labels and e[2] in labels),
             e->(gens[Position(labels,e[1])]*gens[Position(labels,e[2])])^e[3]));
end;
F := FreeGroup(4);
W := F / relators(GeneratorsOfGroup(F),[1..4]);
r := GeneratorsOfGroup(W);
walls := function(k)
    local labels, F3, f, P, p, commuting, central, reps, lifts;
    labels := Difference([1..4],[k]);
    F3 := FreeGroup(3);
    f := GeneratorsOfGroup(F3);
    P := F3 / relators(f,labels);
    p := GeneratorsOfGroup(P);
    commuting := Filtered(labels,i->ForAny(m,e->Set(e{[1,2]})=Set([i,k]) and e[3]=2));
    central := Subgroup(P,p{List(commuting,i->Position(labels,i))});
    reps := AsList(RightTransversal(P,central));
    lifts := List(reps,x->MappedWord(UnderlyingElement(x),f,r{labels}));
    return List(lifts,x->r[k]^x);
end;
counts := [];
for kind in [3,2] do
    b := walls(kind);
    Add(counts,Index(W,Subgroup(W,b)));
od;
outside := b{[2..Length(b)]};
doubled := Concatenation(outside,List(outside,x->x^b[1]));
Add(counts,Index(W,Subgroup(W,doubled)));
if counts<>[12,48,96] then Error("Unexpected reflection-subgroup indices"); fi;
Print("Quarter, octahedron, face double: ",counts," SW2 chambers\n");
if LoadPackage("HAP") = true then
    M := [[1,2,2,3],[2,1,4,2],[2,4,1,4],[3,2,4,1]];
    # This function body is parsed before HAP supplies DavisComplex.
    C := ValueGlobal("DavisComplex")(M);
    counts := List([0..3],C!.dimension);
    if counts<>[14,43,48,18] then Error("Unexpected HAP Davis cell counts"); fi;
    Print("SW2 Davis complex, cell orbits in degrees 0..3: ",counts,"\n");
else
    Print("Skipping the SW2 Davis complex check because GAP could not load HAP.\n");
fi;
#end)();

# The mod-2 face-ring presentation is Davis--Januszkiewicz, Theorem 4.11.
# Count the ranks of multiplication H^1 -> H^2 by every degree-one class.
(function()
local faces, squareDouble, profiles, betti, hist, fs, f, edges, near, A, closed, i,
      rank, e, k, one, mat;
faces := function(n)
    local result,i,j;
    result := [[1..n],[n+1..2*n]];
    for i in [1..n] do
        j := (i mod n)+1;
        Add(result,Set([i,j,n+i]));
        Add(result,Set([n+i,n+j,j]));
    od;
    return result;
end;
squareDouble := function()
    local base, result, f, reflected;
    base := faces(4);
    result := [];
    for f in base{[2..Length(base)]} do
        reflected := List(f,function(x) if x>4 then return x+4; fi; return x; end);
        if Size(Intersection(f,base[1]))=2 then
            Add(result,Union(f,reflected));
        else
            Append(result,[f,reflected]);
        fi;
    od;
    return result;
end;
profiles := [];
betti := [];
one := One(GF(2));
for fs in [squareDouble(),faces(6)] do
    f := Length(fs);
    edges := Filtered(Combinations([1..f],2),e->Length(Intersection(fs[e[1]],fs[e[2]]))=2);
    if f<>14 or Length(edges)<>24 then Error("Unexpected dual graph"); fi;
    near := List([1..f],i->Union([i],Flat(Filtered(edges,e->i in e))));
    if ForAny(Combinations([1..f],3),A->ForAll(Combinations(A,2),e->e in edges)) then
        Error("Unexpected triangle in the face-adjacency graph");
    fi;
    # Monomials in each positive degree have support on one vertex or one edge.
    Add(betti,List([1..4],degree->f+(degree-1)*Length(edges)));
    hist := List([0..f],k->0);
    for A in Combinations([1..f]) do
        closed := [];
        for i in A do UniteSet(closed,near[i]); od;
        rank := Length(closed);
        hist[rank+1] := hist[rank+1]+1;
        # Check the closed-neighbourhood formula against literal GF(2) matrices.
        if Length(A)<=2 or A=[1,3,5,8,9] then
            mat := NullMat(f,f+Length(edges),GF(2));
            for i in A do mat[i][i] := one; od;
            for k in [1..Length(edges)] do
                e := edges[k];
                if e[1] in A then mat[e[2]][f+k] := one; fi;
                if e[2] in A then mat[e[1]][f+k] := one; fi;
            od;
            if RankMat(mat)<>rank then Error("Cup-product rank mismatch"); fi;
        fi;
    od;
    Add(profiles,hist);
od;
if profiles[1]=profiles[2] then Error("Cup-rank distributions agree"); fi;
if betti[1]<>betti[2] then Error("Ordinary mod-2 ranks differ"); fi;
Print("Doubled square and hexagonal antiprism, common mod-2 Betti numbers in degrees 1..4: ",betti[1],"\n");
Print("Number of degree-one classes with cup rank five: ",List(profiles,x->x[6]),"\n");
Print("Number of classes giving injective H^1 -> H^2: ",List(profiles,x->x[Length(x)]),"\n");
Print("Rank distributions (degrees of rank 0..14): ",profiles,"\n");
end)();
