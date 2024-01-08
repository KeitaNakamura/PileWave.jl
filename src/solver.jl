struct FEMCondition{F}
    # input
    gravity::Float64
    t_stop::Float64
    dt_cr::Float64
    input_load::F
    # output
    num_data::Int
    outdir::String
    paraview::Bool
    histinds::Vector{Int}
    show_progress::Bool
    # Newmark-beta method
    β::Float64
    γ::Float64
end

function FEMCondition(file::TOMLFile, grids::Vector{<: Grid})
    # input
    Input = file.Input
    gravity = Input.gravity
    t_stop = Input.t_stop
    ## critical time step
    dx_min = minimum(p->p.length/p.num_elements, file.Pile)
    c = maximum(p->√(p.youngs_modulus/p.density), file.Pile) # sound speed
    dt_cr = file.Advanced.CFL * (dx_min / c)
    ## input_load
    if Input.input_load isa Function
        input_load = Input.input_load
    else # Input.input_load isa String
        path = Input.input_load
        path = isabspath(path) ? path : joinpath(dirname(file.name), path)
        @assert isfile(path)
        data = CSV.File(path; comment = "#")
        input_load = linear_interpolation(data.time, data.load; extrapolation_bc=0.0)
    end

    # output
    Output = file.Output
    num_data = Output.num_data
    ## output directory
    outdir = Output.directory
    outdir = isabspath(outdir) ? outdir : joinpath(dirname(file.name), outdir)
    ## paraview
    paraview = Output.paraview
    ## histinds
    histpts = Output.history_points
    histinds = map(histpts) do pt
        # find closest node index
        grid = first(grids) # any grid is ok because we only want to use `get_allnodes`
        nodes = get_allnodes(grid)
        value, index = findmin(x->abs((only(x)-pt) - only(first(nodes))), nodes)
        index
    end
    ## progress
    show_progress = Output.show_progress

    # parameters for Newmark-beta method
    NewmarkBeta = file.Advanced.NewmarkBeta
    β = NewmarkBeta.beta
    γ = NewmarkBeta.gamma

    FEMCondition(gravity, t_stop, dt_cr, input_load,
                 num_data, outdir, paraview, histinds, show_progress,
                 β, γ)
end

function generate_grids(shape::Femto.Line, Pile::Vector{TOMLPile}, embedded_depth::Real)
    pile_length = sum(pile->pile.length, Pile)
    ndivs = Femto.num_nodes(shape) - 1

    nodes = Vec{1,Float64}[Vec(embedded_depth-pile_length)]
    nodeinds_set = UnitRange{Int}[]
    stop = 1
    for pile in Pile
        start = stop
        n = pile.num_elements
        dx = pile.length / n
        for _ in 1:n, _ in 1:ndivs
            push!(nodes, nodes[end] .+ dx/ndivs)
            stop += 1
        end
        push!(nodeinds_set, start:stop)
    end

    domains = map(nodeinds_set) do nodeinds
        conns = Femto.generate_connectivities(shape, nodeinds)
        Femto.DomainInfo(shape, conns, collect(nodeinds))
    end
    generate_gridset(nodes, domains)
end

function eachlayer(f, layers::Vector)
    # downward is positive
    top = 0.0
    for layer in layers
        bottom = top + layer.thickness
        ret = f(layer, top, bottom)
        ret !== nothing && return ret
        top = bottom
    end
end
function generate_elementstate(grids::Vector{<: Grid}, layers::Vector, piles::Vector{TOMLPile})
    # shaft state
    ests = map(grid->Femto.generate_elementstate(ElementState, grid), grids)
    for (pile, grid, est) in zip(piles, grids, ests)
        for i in eachindex(est)
            set_elementstate!(LazyRow(est, i), pile)
        end
        coords = collect(interpolate(Vf(), grid, reinterpret(Float64, get_allnodes(grid))))
        @assert size(coords) == size(est)
        for i in eachindex(coords, est)
            z = only(coords[i])
            z < 0 && continue
            eachlayer(layers) do layer, top, bottom
                if top ≤ z ≤ bottom
                    set_elementstate!(LazyRow(est, i), layer.shaft)
                    return true
                end
                nothing
            end
        end
    end

    # bottom state
    depth, btm = findmax(only, get_allnodes(last(grids)))
    layer_bottom = eachlayer(layers) do layer, top, bottom
        (top≤depth≤bottom || depth≈bottom) ? layer : nothing
    end
    estbtm = create_elementstatebottom(piles[end], layer_bottom.bottom, btm)

    ests, estbtm
end

#########
# setup #
#########

function setup(file::TOMLFile)
    grids = generate_grids(file.Advanced.shape,
                           file.Pile,
                           file.Input.embedded_depth)
    ests, estbtm = generate_elementstate(grids,
                                         file.SoilLayer,
                                         file.Pile)
    femcond = FEMCondition(file, grids)
    femcond, grids, ests, estbtm
end

#########
# solve #
#########

struct Solution
    t::Vector{Float64}
    z::Vector{Float64}
    u::Matrix{Float64}
    v::Matrix{Float64}
    a::Matrix{Float64}
    f::Matrix{Float64}
    Z::Matrix{Float64}
end

function Solution(num_nodes::Int, num_timestamps::Int)
    t = Vector{Float64}(undef, num_timestamps)
    z = Vector{Float64}(undef, num_nodes)
    u = Matrix{Float64}(undef, num_nodes, num_timestamps)
    v = Matrix{Float64}(undef, num_nodes, num_timestamps)
    a = Matrix{Float64}(undef, num_nodes, num_timestamps)
    f = Matrix{Float64}(undef, num_nodes, num_timestamps)
    Z = Matrix{Float64}(undef, num_nodes, num_timestamps)
    Solution(t, z, u, v, a, f, Z)
end

solve(path::String; return_solution::Bool=false) = solve(read_inputfile(path); return_solution)
solve(dict::Dict{String, Any}; return_solution::Bool=false) = solve(read_input(dict); return_solution)

function solve(file::TOMLFile; return_solution::Bool = false)
    femcond, grids, ests, estbtm = setup(file)
    solve(femcond, grids, ests, estbtm; return_solution)
end 

function solve(
        femcond::FEMCondition,
        grids::Vector{<: Grid},
        ests::Vector{<: StructArray{<: ElementState}},
        estbtm::ElementStateBottom;
        return_solution::Bool,
    )

    g = femcond.gravity
    β = femcond.β
    γ = femcond.γ
    field = ScalarField()

    grid_entire = reduce(grids) do grid1, grid2
        conns = [grid1.connectivities; grid2.connectivities]
        inds = [grid1.nodeindices; grid2.nodeindices]
        Grid(grid1.nodes, grid1.shape, conns, inds)
    end

    ndofs = num_dofs(field, grid_entire)

    K = mapreduce(+, grids, ests) do grid, est
        integrate(field, grid) do i, w, u
            E = est.E[i]
            A = est.A[i]
            E*A*∇(w)⋅∇(u)
        end
    end

    M = Diagonal(mapreduce(+, grids, ests) do grid, est
        integrate(field, grid) do i, w
            ρ = est.ρ[i]
            A = est.A[i]
            w*ρ*A
        end
    end)

    fγ = mapreduce(+, grids, ests) do grid, est
        integrate(field, grid) do i, w
            ρ = est.ρ[i]
            A = est.A[i]
            w*ρ*g*A
        end
    end

    f = zeros(ndofs)
    C_tan = fill!(copy(K), 0)
    K_tan = fill!(copy(K), 0)

    a = zeros(ndofs)
    v = zeros(ndofs)
    u = zeros(ndofs)
    uₙ = zeros(ndofs)
    Δu = zeros(ndofs)

    dirichlet = falses(ndofs)
    if estbtm.σ̄ᵤ == Inf
        dirichlet[estbtm.index] = true
    else
        # do nothing
    end

    # setup outputs
    isdir(femcond.outdir) && rm(femcond.outdir; recursive=true, force=true)
    if femcond.paraview
        mkpath(joinpath(femcond.outdir, "paraview"))
        pvdfile = joinpath(femcond.outdir, "paraview", "fepile1d")
        closepvd(openpvd(pvdfile))
    end
    if !isempty(femcond.histinds)
        mkpath(joinpath(femcond.outdir, "history"))
        for i in eachindex(femcond.histinds)
            index = femcond.histinds[i]
            open(joinpath(femcond.outdir, "history", "history_$i.csv"), "w") do io
                nodes = get_allnodes(grid_entire)
                z = only(nodes[index] - nodes[1])
                write(io, "# data at z = $z\n")
                write(io, "time,displacement,velocity,acceleration,force,force_down,force_up\n")
            end
        end
    end

    timestamps = LinRange(0, femcond.t_stop, round(Int, femcond.t_stop/femcond.dt_cr))
    Δt = step(timestamps)

    savepoints = collect(LinRange(0, femcond.t_stop, femcond.num_data))
    savecounts = 1

    if return_solution
        primarynodes = get_allnodes(grid_entire, 1)
        nprimarynodes = length(primarynodes)
        sol = Solution(nprimarynodes, femcond.num_data)
        sol.z .= map(z->only(z-first(primarynodes)), primarynodes)
    end

    if femcond.show_progress
        prog = ProgressMeter.Progress(length(timestamps))
    end
    for (step, t) in enumerate(timestamps)
        f .= fγ
        f[begin] += femcond.input_load(t)

        uₙ .= u
        foreach(elementstate_startup!, ests)
        elementstate_startup!(estbtm)

        # predictor
        ũ = u + Δt*v + Δt^2*(1-2β)*a/2
        ṽ = v + (1-γ)*Δt*a
        @. u = ũ
        @. v = ṽ

        # corrector
        function R!(ψ, u)
            @. a = (u - ũ) / (Δt^2*β)
            @. v = ṽ + γ*Δt*a
            @. Δu = u - uₙ

            p = K * u

            for (grid, est) in zip(grids, ests)
                Δu˜ = interpolate(field, grid, Δu)
                v˜ = interpolate(field, grid, v)
                assemble!(p, grid, Δu˜, v˜, est)
            end
            assemble!(p, Δu, v, estbtm)

            ψ .= M*a + p - f
        end
        function J!(K★, u)
            C_tan .= 0
            K_tan .= K

            for (grid, est) in zip(grids, ests)
                Δu˜ = interpolate(field, grid, Δu)
                v˜ = interpolate(field, grid, v)
                assemble!(C_tan, K_tan, grid, Δu˜, v˜, est)
            end
            assemble!(C_tan, K_tan, Δu, v, estbtm)

            K★ .= M/(Δt^2*β) + γ*C_tan/(Δt*β) + K_tan
        end
        nlsolve!(R!, J!, u, dirichlet; symmetric=true, f_tol=1e-8, dx_tol=1e-12)

        if t ≥ first(savepoints)
            popfirst!(savepoints)

            V = zeros(ndofs)
            FV = zeros(ndofs)
            ZV = zeros(ndofs) # projection for impedance `Z`
            for (grid, est) in zip(grids, ests)
                ū = interpolate(field, grid, u)
                integrate!((i,w)->one(w), V, field, grid)
                integrate!(FV, field, grid) do i, w
                    E = est.E[i]
                    A = est.A[i]
                    σ = -E * ∇(ū[i]) # compression is positive
                    only(A * σ)
                end
                integrate!(ZV, field, grid) do i, w
                    E = est.E[i]
                    A = est.A[i]
                    ρ = est.ρ[i]
                    c = √(E/ρ)
                    E*A/c
                end
            end
            Fᵢ = FV ./ V
            Zᵢ = ZV ./ V

            if return_solution
                sol.t[savecounts] = t
                sol.u[:,savecounts] .= view(u, 1:nprimarynodes)
                sol.v[:,savecounts] .= view(v, 1:nprimarynodes)
                sol.a[:,savecounts] .= view(a, 1:nprimarynodes)
                sol.f[:,savecounts] .= view(Fᵢ, 1:nprimarynodes)
                sol.Z[:,savecounts] .= view(Zᵢ, 1:nprimarynodes)
            end

            if femcond.paraview
                openpvd(pvdfile; append=true) do pvd
                    openvtk(joinpath(femcond.outdir, "paraview", "fepile1d_$savecounts"), grid_entire) do vtk
                        vtk["Displacement"] = u
                        vtk["Force"] = Fᵢ
                        pvd[t] = vtk
                    end
                end
            end

            if !isempty(femcond.histinds)
                for i in eachindex(femcond.histinds)
                    index = femcond.histinds[i]
                    open(joinpath(femcond.outdir, "history", "history_$i.csv"), "a") do io
                        displacement = u[index]
                        velocity = v[index]
                        acceleration = v[index]
                        force = Fᵢ[index]
                        impedance = Zᵢ[index]
                        force_down = (force + impedance*velocity) / 2
                        force_up = (force - impedance*velocity) / 2
                        write(io, "$t,$displacement,$velocity,$acceleration,$force,$force_down,$force_up\n")
                    end
                end
            end

            savecounts += 1
        end

        if femcond.show_progress
            ProgressMeter.update!(prog, step)
        end
    end

    return_solution ? sol : nothing
end
