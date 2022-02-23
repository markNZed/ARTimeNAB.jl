# This code is derived from https://github.com/AP6YC/AdaptiveResonance.jl 
# Some of the changes for ARTime made it difficult to request changes to AdaptiveResonance
# Rather than forking, the essential feature were "compacted" into this file
# Please support the excellent AdaptiveResonance project (licensed under MIT)

module AdaptiveResonance

using Parameters    # ARTopts are parameters (@with_kw)
using Logging       # Logging utils used as main method of terminal reporting
using LinearAlgebra: norm   # Trace and norms
using Statistics: median, mean  # Medians and mean for linkage methods

# Abstract types
abstract type ARTOpts end               # ART module options
abstract type ARTModule end             # ART modules
abstract type ART <: ARTModule end      # ART (unsupervised)

mutable struct DataConfig
    setup::Bool
    dim::Integer
    dim_comp::Integer
end

function DataConfig()
    DataConfig(
        false,                      # setup
        0,                          # dim
        0                           # dim_comp
    )
end

@with_kw mutable struct opts_DVFA <: ARTOpts @deftype Float64
    # Lower-bound vigilance parameter: [0, 1]
    rho_lb = 0.0; @assert rho_lb >= 0.0 && rho_lb <= 1.0
    # Upper bound vigilance parameter: [0, 1]
    rho_ub = 0.0; @assert rho_ub >= 0.0 && rho_ub <= 1.0
end # opts_DVFA

mutable struct DVFA <: ART
    # Get parameters
    opts::opts_DVFA
    config::DataConfig
    # Working variables
    labels::Vector{Integer}
    W::AbstractArray{Float64, 2}
    Wx::AbstractArray{Float64, 2}
    M::Vector{Float64}
    Me::Vector{Float64}
    A::Vector{Float64}
    Ae::Vector{Float64}
    map::Vector{Integer}
    bmu::Vector{Integer}
    n_categories::Integer
    n_clusters::Integer
end

function DVFA()
    opts = opts_DVFA()
    DVFA(opts)
end # DVFA()

function DVFA(opts::opts_DVFA)
    DVFA(
        opts,                           # opts
        DataConfig(),                   # config
        Array{Integer}(undef, 0),       # labels
        Array{Float64}(undef, 0, 0),    # W
        Array{Float64}(undef, 0, 0),    # Wx
        Array{Float64}(undef, 0),       # M
        Array{Float64}(undef, 0),       # Me
        Array{Float64}(undef, 0),       # A
        Array{Float64}(undef, 0),       # Ae
        Array{Integer}(undef, 0),       # map
        Array{Integer}(undef, 0),       # bmu
        0,                              # n_categories
        0,                              # n_clusters
    )
end 

function train!(art::DVFA, x; learning::Bool=true)
    # Data information and setup
    if ndims(x) > 1
        n_samples = size(x)[2]
    else
        n_samples = 1
    end
    
    x = vcat(x, 1 .- x) # complement code
    if n_samples == 1
        y_hat = zero(Integer)
    else
        y_hat = zeros(Integer, n_samples)
    end
    # Initialization
    if isempty(art.W)
        # Set the first label as either 1 or the first provided label
        local_label = 1
        # Add the local label to the output vector
        if n_samples == 1
            y_hat = local_label
        else
            y_hat[1] = local_label
        end
        # Create a new category and cluster
        art.W = ones(art.config.dim_comp, 1)
        art.Wx = zeros(art.config.dim_comp, 1)
        art.n_categories = 1
        art.n_clusters = 1
        push!(art.labels, local_label)
        # Skip the first training entry
        push!(art.A, 0.0)
        push!(art.Ae, 0.0)
        push!(art.bmu, 1)
        skip_first = true
    else
        skip_first = false
    end
    for i in 1:n_samples
        # Skip the first sample if we just initialized
        (i == 1 && skip_first) && continue
        # Grab the sample slice
        sample = x[:, i]
        max_bmu = 0.0
        max_bmue = 0.0
        bmu_with_max = -1
        # Compute the activation and match for all categories
        activation_match!(art, sample)
        # Sort activation function values in descending order
        index = sortperm(art.M, rev=true)
        # Default to mismatch
        mismatch_flag = true
        label = -1
        # Loop over all categories
        for j = 1:art.n_categories
            # Best matching unit, order does not matter
            bmu = index[j]
            # Vigilance test upper bound
            if art.M[bmu] > max_bmu
                bmu_with_max = bmu
                max_bmue = art.Me[bmu]
                max_bmu = art.M[bmu]
            end
            if !learning
                # no learning
            elseif art.M[bmu] >= art.opts.rho_ub
                # learn with fast commit
                art.W[:, bmu] = art.Wx[:, bmu]
                # Update sample label for output`
                if art.M[bmu] >= max_bmu
                    label = art.labels[bmu]
                end
                mismatch_flag = false
            # Vigilance test lower bound
            elseif art.M[bmu] >= art.opts.rho_lb && mismatch_flag
                if art.M[bmu] >= max_bmu
                    label = art.labels[bmu]
                end
                push!(art.labels, label)
                # Fast commit the sample, same as per mismatch
                art.W = hcat(art.W, sample)
                art.Wx = hcat(art.W, zeros(art.config.dim_comp, 1))
                art.n_categories += 1
                # No mismatch
                mismatch_flag = false
                break
            else
                break
            end
        end
        push!(art.A, max_bmu)
        push!(art.Ae, max_bmue)
        push!(art.bmu, bmu_with_max)
        # If there was no resonant category, make a new one
        if mismatch_flag && learning
            label = -1
            push!(art.labels, art.n_clusters + 1)
            # Fast commit the sample
            art.W = hcat(art.W, sample)
            art.Wx = hcat(art.W, zeros(art.config.dim_comp, 1))
            # Increment the number of categories and clusters
            art.n_categories += 1
            art.n_clusters += 1
        end

        if n_samples == 1
            y_hat = label
        else
            y_hat[i] = label
        end
    end
    return y_hat
end 

function activation_match!(art::DVFA, x::Vector)
    art.M = zeros(art.n_categories)
    art.Me = zeros(art.n_categories)
    ei = length(x)÷2
    for jx = 1:art.n_categories
        W = art.W[:, jx]
        em = minimum([x W], dims = 2)
        art.Wx[:, jx] = em #stored because this can be reused in the learning (fast commit)
        numerator = norm(em, 1)
        nW = norm(W, 1)
        if nW == 0 nW = 0.001 end
        feature_similarity = numerator/nW 
        eme = [em[ei];em[end]]
        We = [W[ei];W[end]]
        neme = norm(eme, 1)
        nWe = norm(We, 1)
        energy_similarity = neme / nWe
        art.M[jx] = feature_similarity^3 * energy_similarity^2
        art.Me[jx] = energy_similarity
    end 
end 

end