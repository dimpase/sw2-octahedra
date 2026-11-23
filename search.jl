# Julia 1.10+, standard library only. Run with --help for the command line.
using TOML
using Base.Checked: checked_add, checked_sub, checked_mul, checked_neg

const V = NTuple{4,Int}
const Frame = NTuple{4,V}  # columns in coordinates (time,x,y,z)
const I4 = ntuple(j -> ntuple(i -> Int(i == j),4),4)
const NORMALS = ((0,-1,1,0), (1,1,1,1), (0,0,0,-1), (0,0,-1,1))
const CORNERS = ((1,1,0,0), (1,0,0,0), (3,1,1,1), (2,1,1,0))

plus(a::V,b::V) = map(checked_add,a,b)
minus(a::V,b::V) = map(checked_sub,a,b)
scale(k::Int,a::V) = map(x -> checked_mul(k,x),a)
lor(a::V,b::V) = checked_add(checked_sub(checked_mul(a[2],b[2]),checked_mul(a[1],b[1])),
                              checked_add(checked_mul(a[3],b[3]),checked_mul(a[4],b[4])))
plane(a::V) = a[findfirst(!iszero,a)] > 0 ? a : map(checked_neg,a)
function act(g::Frame, a::V)
    plus(plus(scale(a[1],g[1]),scale(a[2],g[2])),plus(scale(a[3],g[3]),scale(a[4],g[4])))
end
compose(g::Frame,h::Frame) = map(a -> act(g,a),h)
function reflection(n::V)
    d = lor(n,n)
    d in (1,2) || error("Unexpected root length $d")
    return map(e -> minus(e,scale(checked_mul(2÷d,lor(n,e)),n)),I4)
end
const R = map(reflection,NORMALS)

# A unit is the orbit of the simplex under the indicated finite parabolic group.
# A word is in the four face labels of the paper, in multiplication order.
function parabolic(labels)
    words = Dict{Frame,String}(I4 => "")
    todo = [I4]
    for g in todo, j in labels
        h = compose(g,R[j])
        if !haskey(words,h)
            words[h] = words[g] * string(j)
            push!(todo,h)
        end
    end
    return todo,words
end

struct Unit
    normals::Vector{V}
    mirrors::Vector{Frame}
    words::Vector{String}
    corners::Vector{V}
    finite::Vector{Int}
    orbit::Vector{String}
end

function makeunit(labels)
    sub,words = parabolic(labels)
    boundary = Dict{V,String}()
    allpoints = Set{V}()
    for g in sub
        union!(allpoints,act(g,p) for p in CORNERS)
        for j in setdiff(1:4,labels)
            boundary[act(g,NORMALS[j])] = words[g] * string(j) * reverse(words[g])
        end
    end
    ns = sort!(collect(keys(boundary)))
    points = sort!([p for p in allpoints if count(n -> lor(n,p)==0,ns)>=3])
    for n in ns
        all(p -> lor(n,p)<=0,allpoints) || error("Unit has a nonsupporting facet")
    end
    for (i,n) in enumerate(ns), m in ns[i+1:end]
        # A non-right dihedral angle will be detected at a finite endpoint.
        common = [p for p in points if lor(n,p)==lor(m,p)==0]
        if length(common)>=2 && all(p -> lor(p,p)==0,common)
            @assert lor(n,m)==0
        end
    end
    Unit(ns,reflection.(ns),[boundary[n] for n in ns],points,
         findall(p -> lor(p,p)<0,points),[words[g] for g in sub])
end

struct Tile
    frame::Frame
    walls::Vector{Int}
    neighbors::Vector{Int}
    finitepoints::Vector{V}
end
mutable struct Tiling
    unit::Unit
    ids::Dict{Frame,Int}
    planeids::Dict{V,Int}
    tiles::Vector{Tile}
end
function tile!(t::Tiling,g::Frame)
    get!(t.ids,g) do
        walls = [get!(() -> length(t.planeids)+1,t.planeids,plane(act(g,n))) for n in t.unit.normals]
        push!(t.tiles,Tile(g,walls,zeros(Int,length(walls)),
                          [act(g,t.unit.corners[i]) for i in t.unit.finite]))
        length(t.tiles)
    end
end
function neighbor!(t::Tiling,id::Int,f::Int)
    a = t.tiles[id]
    if a.neighbors[f]==0
        b = tile!(t,compose(a.frame,t.unit.mirrors[f]))
        a.neighbors[f] = b
        t.tiles[b].neighbors[f] = id
    end
    return a.neighbors[f]
end

# Each finite vertex link is a spherical tessellation. List the chambers in
# every hemisphere and every allowed wedge. A point of an ideal polyhedron can lie
# only in one of these sectors, or in the interior (the whole sphere).
struct Star
    frames::Vector{Frame}
    sectors::Vector{UInt64}
end
function makestar(u::Unit,v::V; quarters=false)
    anchor = reduce(plus,u.corners)
    localfaces = findall(n -> lor(n,v)==0,u.normals)
    frames = [I4]
    known = Set(frames)
    through = Set{V}()
    for g in frames, j in localfaces
        push!(through,plane(act(g,u.normals[j])))
        h = compose(g,u.mirrors[j])
        if !(h in known)
            push!(known,h)
            push!(frames,h)
        end
    end
    @assert length(frames)<64
    full = (UInt64(1)<<length(frames))-UInt64(1)
    pts = [act(g,anchor) for g in frames]
    ns = sort!(collect(through))
    halves = UInt64[]
    for n in ns
        @assert all(p -> lor(n,p)!=0,pts)
        mask = sum(UInt64(1)<<(i-1) for (i,p) in enumerate(pts) if lor(n,p)<0)
        push!(halves,mask)
    end
    # A later search state can approach this same vertex from a different
    # chamber. Cache the entire link, with both sides of every mirror.
    sectors = Set([full; halves; [full ⊻ s for s in halves]])
    for (i,n) in enumerate(ns), j in 1:i-1
        d = lor(n,ns[j])
        right = d==0
        diagonal = quarters && checked_mul(2,checked_mul(d,d))==checked_mul(lor(n,n),lor(ns[j],ns[j]))
        if right || diagonal
            for (signa,a) in ((1,halves[i]),(-1,full ⊻ halves[i])), (signb,b) in ((1,halves[j]),(-1,full ⊻ halves[j]))
                # Outward normals make angle 3pi/4 for an interior pi/4 sector.
                right || signa*signb*d<0 || continue
                push!(sectors,a & b)
            end
        end
    end
    Star(frames,sort!(collect(sectors);by=x->(count_ones(x),x)))
end
struct PlacedStar
    ids::Vector{Int}
    sectors::Vector{UInt64}
end
function placed!(t,templates,cache,id,k)
    a = t.tiles[id]
    get!(cache,a.finitepoints[k]) do
        star = templates[k]
        PlacedStar([tile!(t,compose(a.frame,g)) for g in star.frames],star.sectors)
    end
end
function needed(star::PlacedStar,cells)
    mask = UInt64(0)
    for (i,id) in enumerate(star.ids)
        id in cells && (mask |= UInt64(1)<<(i-1))
    end
    mask in star.sectors && return UInt64[]
    choices = UInt64[]
    for sector in star.sectors
        if (sector & mask)==mask && !any(s -> (s & sector)==s,choices)
            push!(choices,sector)
        end
    end
    @assert !isempty(choices)
    return [s & ~mask for s in choices]
end

# Return the component of the root after opening the given global mirrors.
# This is the convex hull in the Coxeter tessellation of the connected
# inputs whose internal mirrors were opened. Stop immediately past the bound.
function closehull(t,start,opened,firststeps,bound)
    cells = copy(start)
    todo = Int[]
    for id in firststeps
        if !(id in cells)
            push!(cells,id)
            push!(todo,id)
        end
    end
    length(cells)>bound && return nothing
    for id in todo
        a = t.tiles[id]
        for (j,w) in enumerate(a.walls)
            w in opened || continue
            b = neighbor!(t,id,j)
            if !(b in cells)
                push!(cells,b)
                length(cells)>bound && return nothing
                push!(todo,b)
            end
        end
    end
    return cells
end
function boundary(t,cells,opened)
    out = Dict{Int,Tuple{Int,Int}}()
    for id in cells, (j,w) in enumerate(t.tiles[id].walls)
        w in opened || (out[w]=(id,j))
    end
    return out
end

# All roots realizing the lexicographically least edge-coloured neighbourhood
# yield the same BFS code precisely for two translates of one union. Facet
# labels are invariant under the unit reflection group (left multiplication).
function shapecode(t,cells)
    masks = Dict(id=>sum((n in cells ? 1<<(j-1) : 0) for (j,n) in enumerate(t.tiles[id].neighbors)) for id in cells)
    least = minimum(values(masks))
    best = UInt8[]
    for (root,m) in masks
        m==least || continue
        todo = [root]
        visited = Set(todo)
        code = UInt8[]
        equal_prefix = !isempty(best)
        for id in todo
            push!(code,masks[id])
            if equal_prefix
                next_value = best[length(code)]
                code[end]>next_value && break
                equal_prefix = code[end]==next_value
            end
            for n in t.tiles[id].neighbors
                if n in cells && !(n in visited)
                    push!(visited,n)
                    push!(todo,n)
                end
            end
        end
        length(code)==length(cells) && (isempty(best) || isless(code,best)) && (best=code)
    end
    return best
end

function choose_defect(t,templates,starcache,cells,bound)
    visited = Set{V}()
    best = nothing
    bestscore = (typemax(Int),0)
    charged = Set{Int}()
    lowerbound = 0
    slack = bound-length(cells)
    for id in cells, (k,v) in enumerate(t.tiles[id].finitepoints)
        v in visited && continue
        push!(visited,v)
        star = placed!(t,templates,starcache,id,k)
        additions = needed(star,cells)
        isempty(additions) && continue
        keep = filter(mask -> length(cells)+count_ones(mask)<=bound,additions)
        isempty(keep) && return Vector{Vector{Int}}()
        # Each surviving sector requires at least this many tiles outside the
        # regions already counted at earlier vertices. These charges are disjoint.
        offered = reduce(|,keep)
        uncharged = offered
        bits = offered
        while !iszero(bits)
            i = trailing_zeros(bits)+1
            star.ids[i] in charged && (uncharged &= ~(UInt64(1)<<(i-1)))
            bits &= bits-UInt64(1)
        end
        extra = minimum(mask -> count_ones(mask & uncharged),keep)
        if extra>0
            lowerbound += extra
            lowerbound>slack && return Vector{Vector{Int}}()
            while !iszero(uncharged)
                push!(charged,star.ids[trailing_zeros(uncharged)+1])
                uncharged &= uncharged-UInt64(1)
            end
        end
        score = (length(keep),-minimum(count_ones,keep))
        if score<bestscore
            best = [[star.ids[i] for i in eachindex(star.ids) if !iszero((mask>>(i-1))&1)] for mask in keep]
            bestscore=score
        end
    end
    return best
end

function internal!(t,opened,cells,added)
    allcells = union(cells,added)
    for id in added, j in eachindex(t.tiles[id].walls)
        neighbor!(t,id,j) in allcells && push!(opened,t.tiles[id].walls[j])
    end
end
function chamberwords(t,cells)
    todo=[1]
    paths=Dict(1=>"")
    for id in todo, (j,n) in enumerate(t.tiles[id].neighbors)
        if n in cells && !haskey(paths,n)
            paths[n]=paths[id]*t.unit.words[j]
            push!(todo,n)
        end
    end
    @assert length(paths)==length(cells)
    sort!([path*word for path in values(paths) for word in t.unit.orbit];by=x->(length(x),x))
end

function search(kind::String,budget::Int)
    labels = kind=="all3" ? [1,2,4] : kind=="no3" ? [3] : Int[]
    u=makeunit(labels)
    bound=budget÷length(u.orbit)
    t=Tiling(u,Dict{Frame,Int}(),Dict{V,Int}(),Tile[])
    tile!(t,I4)
    templates=[makestar(u,u.corners[i];quarters=kind=="raw-even") for i in u.finite]
    starcache=Dict{V,PlacedStar}()
    seen=Set{Vector{UInt8}}()
    results=Dict{String,Any}[]
    function visit(cells,opened)
        # Fill the neighbours on internal mirrors for a complete adjacency code.
        for id in cells, (j,w) in enumerate(t.tiles[id].walls)
            if w in opened
                @assert neighbor!(t,id,j) in cells
            end
        end
        code=shapecode(t,cells)
        code in seen && return
        push!(seen,code)
        if length(seen)%1_000_000==0
            println(stderr,kind," through ",budget," chambers: ",length(seen),
                    " convex configurations visited, ",length(results)," ideal representatives found")
        end
        choices=choose_defect(t,templates,starcache,cells,bound)
        if choices===nothing
            push!(results,Dict("case"=>kind,"index"=>length(cells)*length(u.orbit),"words"=>chamberwords(t,cells)))
            choices=[[neighbor!(t,id,j)] for (id,j) in values(boundary(t,cells,opened))]
        end
        for added in choices
            extended=copy(opened)
            internal!(t,extended,cells,added)
            # A newly crossed mirror must open every old facet in that plane.
            seeds=copy(added)
            for id in cells, (j,w) in enumerate(t.tiles[id].walls)
                w in extended && !(w in opened) && push!(seeds,neighbor!(t,id,j))
            end
            closed=closehull(t,cells,extended,seeds,bound)
            if closed!==nothing
                @assert length(closed)>length(cells)
                visit(closed,extended)
            end
        end
    end
    bound>=1 && visit(Set([1]),Set{Int}())
    sort!(results;by=r->(r["index"],r["words"]))
    println(kind,": ",length(seen)," convex configurations visited; ideal representatives by SW2 index: ",
            [(i,count(r->r["index"]==i,results)) for i in sort!(unique(Int[r["index"] for r in results]))])
    return results,length(seen)
end

function main(args)
    if length(args) ∉ (2,3) || "--help" in args
        println("julia search.jl MAX_SW2_CHAMBERS OUTPUT.toml [both|all3|no3|raw|raw-even]")
        return
    end
    budget=parse(Int,args[1])
    budget>0 || error("The chamber bound must be positive")
    mode=get(args,3,"both")
    mode in ("both","all3","no3","raw","raw-even") || error("Unknown wall case")
    results=Dict{String,Any}[]
    visits=Dict{String,Int}()
    for kind in (mode=="both" ? ["all3","no3"] : [mode])
        records,visits[kind]=search(kind,budget)
        append!(results,records)
    end
    mktemp(dirname(abspath(args[2]))) do temporary,io
        TOML.print(io,Dict("bound"=>budget,"visited"=>visits,"polyhedra"=>results);sorted=true)
        close(io)
        mv(temporary,args[2];force=true)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    main(ARGS)
end
