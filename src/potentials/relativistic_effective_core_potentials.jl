struct RelativisticEffectiveCorePotential{T} <: AbstractEffectiveCorePotential{T}
    name::String
    gst_config::Configuration{Orbital{Int}}
    Q::Int
    V₋::Vector{GaussianExpansion{T}} # For negative κ
    V₊::Vector{GaussianExpansion{T}} # For positive κ
    reference::String
end

charge(pp::RelativisticEffectiveCorePotential) = num_electrons(pp.gst_config)
ground_state(pp::RelativisticEffectiveCorePotential) = pp.gst_config

Base.show(io::IO, pp::RelativisticEffectiveCorePotential{T}) where {T} =
    write(io, "Relativistic effective core potential for $(pp.name) ($(pp.gst_config)), Z = $(charge(pp))")

function κrange(pp::RelativisticEffectiveCorePotential)
    κ₋max = length(pp.V₋)
    κ₊max = length(pp.V₊)

    sort(vcat(-1:-1:-κ₋max,1:κ₊max), lt=(a,b)->abs(a)<abs(b) || abs(a)==abs(b) && a<0)
end

function Base.show(io::IO, ::MIME"text/plain", pp::RelativisticEffectiveCorePotential)
    show(io, pp)
    println(io)
    κ₋max = length(pp.V₋)
    κ₊max = length(pp.V₊)
    print(io, "Long-range Q = $(pp.Q), κ₋ ∈ -1:-1:$(-κ₋max)")
    κ₊max > 0 && print(io, ", κ₊ ∈ 1:$(κ₊max)")
    println(io, "\nData from \"$(pp.reference)\"")

    κs = κrange(pp)

    data = map(κs) do κ
        data = κ > 0 ? pp.V₊[κ] : pp.V₋[-κ]
        nk = length(data)
        ℓ = AtomicLevels.kappa_to_ℓ(κ)
        j = AtomicLevels.kappa_to_j(κ)
        hcat([spectroscopic_label(ℓ);repeat([""],nk-1)],
             [j;repeat([""],nk-1)],
             1:nk, data.n, data.B, data.β)
    end |> v -> vcat(v...)

    headers = ["ℓ", "j", "k", "n", "B", "β"]

    pretty_table(io, data, headers)
end

function (pp::RelativisticEffectiveCorePotential{T})(orb::RelativisticOrbital, r::AbstractVector{T}) where T
    V = -pp.Q./r

    κ = orb.κ

    data = if κ < 0
        -κ > length(pp.V₋) && return V
        pp.V₋[-κ]
    else
        κ > length(pp.V₊) && return V
        pp.V₊[κ]
    end
    V += data(r)

    V
end

function spin_orbit_potential(pp::RelativisticEffectiveCorePotential{T}, r::AbstractVector{T},
                              a::SpinOrbital{<:Orbital}, b::SpinOrbital{<:Orbital}) where T
    # This returns the (off-)diagonal, spin–orbit part of the
    # relativistic ECP.
    s = half(1)
    ℓ,mℓa = jmⱼ(a)
    mas  = spin(a)
    ℓb,mℓb = jmⱼ(b)
    mbs = spin(b)

    V = zeros(T, length(r))

    ℓ == ℓb || return V

    κs = κrange(pp)

    for j = [ℓ+s, ℓ-s]
        κ = AtomicLevels.ℓj_to_kappa(ℓ,j)
        j ≥ abs(mℓa+mas) && κ ∈ κs || continue

        coeff = clebschgordan(
            ℓ,mℓa,s,mas,j
        )*clebschgordan(
            ℓ,mℓb,s,mbs,j
        )
        V += coeff*((κ > 0 ? pp.V₊[κ] : pp.V₋[-κ])(r))
    end

    V
end

function spin_orbit_potential(pp::RelativisticEffectiveCorePotential{T}, r::AbstractVector{T},
                              a::SpinOrbital{<:RelativisticOrbital}, b::SpinOrbital{<:RelativisticOrbital}) where T
    # This returns the (off-)diagonal, spin–orbit part of the
    # relativistic ECP.
    κ,κb = a.orb.κ,b.orb.κ

    V = zeros(T, length(r))

    a.orb.ℓ == b.orb.ℓ && a.orb.j == b.orb.j && a.m == b.m || return V

    κs = κrange(pp)

    if κ ∈ κs
        V += (κ > 0 ? pp.V₊[κ] : pp.V₋[-κ])(r)
    end

    V
end

function (pp::RelativisticEffectiveCorePotential{T})(orb::SpinOrbital, r::AbstractVector{T}) where T
    V = -pp.Q./r

    V += spin_orbit_potential(pp, r, orb, orb)

    V
end


(pp::RelativisticEffectiveCorePotential{T})(orb::AbstractOrbital, r::T) where T = pp(orb, [r])[1]
