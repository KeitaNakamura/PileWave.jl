struct FEMCondition{F}
    # input
    gravity::Float64
    t_stop::Float64
    dt_cr::Float64
    inputload::F
    # output
    num_data::Int
    outdir::String
    histinds::Vector{Int}
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
    ## inputload
    if Input.inputload isa Function
        inputload = Input.inputload
    else # Input.inputload isa String
        path = Input.inputload
        path = isabspath(path) ? path : joinpath(dirname(file.name), path)
        @assert isfile(path)
        data = CSV.File(path; comment = "#")
        inputload = linear_interpolation(data.time, data.load; extrapolation_bc=0.0)
    end

    # output
    Output = file.Output
    num_data = Output.num_data
    ## output directory
    outdir = Output.directory
    outdir = isabspath(outdir) ? outdir : joinpath(dirname(file.name), outdir)
    ## histinds
    histpts = Output.history_points
    histinds = map(histpts) do pt
        # find closest node index
        grid = first(grids) # any grid is ok because we only want to use `get_allnodes`
        nodes = get_allnodes(grid)
        value, index = findmin(x->abs((only(x)-pt) - only(first(nodes))), nodes)
        index
    end

    # parameters for Newmark-beta method
    NewmarkBeta = file.Advanced.NewmarkBeta
    β = NewmarkBeta.beta
    γ = NewmarkBeta.gamma

    FEMCondition(gravity, t_stop, dt_cr, inputload, num_data, outdir, histinds, β, γ)
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
function generate_elementstate(::Type{Model}, grids::Vector{<: Grid}, layers::Vector, piles::Vector{TOMLPile}) where {Model}
    # shaft state
    ests = map(grid->Femto.generate_elementstate(get_elementstate_type(Model), grid), grids)
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
    estbtm = create_elementstatebottom(Model, piles[end], layer_bottom.bottom, btm)

    ests, estbtm
end

#########
# setup #
#########

function setup(file::TOMLFile)
    grids = generate_grids(file.Advanced.shape,
                           file.Pile,
                           file.Input.embedded_depth)
    ests, estbtm = generate_elementstate(file.Input.soilmodel,
                                         grids,
                                         file.SoilLayer,
                                         file.Pile)
    femcond = FEMCondition(file, grids)
    femcond, grids, ests, estbtm
end

#########
# solve #
#########

solve(path::String; return_u::Bool=false) = solve(read_inputfile(path); return_u)
solve(dict::Dict{String, Any}; return_u::Bool=false) = solve(read_input(dict); return_u)

function solve(file::TOMLFile; return_u::Bool = false)
    femcond, grids, ests, estbtm = setup(file)
    solve(femcond, grids, ests, estbtm; return_u)
end 

function solve(
        femcond::FEMCondition,
        grids::Vector{<: Grid},
        ests::Vector{<: StructArray{<: ElementState}},
        estbtm::ElementStateBottom;
        return_u::Bool,
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
    mkpath(joinpath(femcond.outdir, "paraview"))
    mkpath(joinpath(femcond.outdir, "history"))
    ## paraview
    pvdfile = joinpath(femcond.outdir, "paraview", "fepile1d")
    openpvd(identity, pvdfile)
    ## history
    for i in eachindex(femcond.histinds)
        index = femcond.histinds[i]
        open(joinpath(femcond.outdir, "history", "history_$i.csv"), "w") do io
            z = only(get_allnodes(grid_entire)[index])
            write(io, "# data at z = $z\n")
            write(io, "time,displacement,velocity,acceleration,force,force_down,force_up\n")
        end
    end

    timestamps = LinRange(0, femcond.t_stop, round(Int, femcond.t_stop/femcond.dt_cr))
    Δt = step(timestamps)

    @showprogress for (step, t) in enumerate(timestamps)
        f .= fγ
        f[begin] += femcond.inputload(t)

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
        nlsolve!(R!, J!, u, dirichlet; symmetric=true, f_tol=1e-8, x_tol=1e-12)

        if step == 1 || step % max(1, length(timestamps)÷femcond.num_data) == 0
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
            openpvd(pvdfile; append=true) do pvd
                openvtk(joinpath(femcond.outdir, "paraview", "fepile1d_$step"), grid_entire) do vtk
                    vtk["Displacement"] = u
                    vtk["Force"] = Fᵢ
                    pvd[t] = vtk
                end
            end
            # history
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
    end

    return_u ? u : nothing
end
