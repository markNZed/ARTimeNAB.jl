module ARTime

using RedefStructs
using Wavelets
using Statistics
using OnlineStats
import LinearAlgebra: norm

include("AdaptiveResonanceCompact.jl")

@kwredef mutable struct ClassifyState
    art::AdaptiveResonance.DVFA = AdaptiveResonance.DVFA()
    last_anomaly_i::Int = 0
    last_anomaly_sim::Float64 = 0.0
    last_rho_update_i::Int = 0
    mask_after_cat::Bool = false
    no_new_cat_count::Int = 0
    trend_window_f = []
    anomaly_sim_history::Vector{Float64} = []
    sim_diff_window::Vector{Float64} = []
    rho_ub_mean::Mean{Float64, EqualWeight} = OnlineStats.Mean()
    sim_window::Vector{Float64} = []
    ds_window::Vector{Float64} = []
    ds_moving_average::Float64 = 0.0
    medbin::Vector{Int} = []
    medlevel::Int = 0
    belowmed::Int = 0
    abovemed::Int = 0
    f_window::Vector{Float64} = []
    dsi::Int = 1
end

@kwredef mutable struct P
    i::Int = 1
    cs::ClassifyState = ClassifyState()
    wavelett = wavelet(WT.haar)
    datafile::String = ""
    dmin::Float64 = 0.0
    dmax::Float64 = 0.0
    dlength::Int = 0
    window::Int = 8
    probationary_period::Int = 0
    windows_per_pb = 13
    sstep::Int = 1
    discretize_chomp::Float64 = 0.075
    nlevels::Int = 80
    mask_rho_after_anomaly::Int = 0 
    trend_window::Int = 0
    initial_rho::Float64 = 0.80
end

function init(dmin, dmax, dlength, p=p)
    p.dmin = dmin
    p.dmax = dmax
    p.dlength = dlength
    probationary_period = dlength < 5000 ? Int.(floor(0.15 * dlength)) : 750
    p.probationary_period = probationary_period - mod(probationary_period, 2) # make an even number
    p.sstep = max(1, round(Int, div(p.probationary_period, p.window * p.windows_per_pb)))
    p.trend_window = floor(Int, p.probationary_period / p.sstep)
    p.mask_rho_after_anomaly = p.window * 1.5
    # initialise detector state variables
    p.cs.sim_window = ones(p.trend_window ÷ 2 + 1)
    p.cs.sim_diff_window = zeros(p.trend_window + 1)
    p.cs.ds_window = zeros(p.sstep) # downsampling window
    p.cs.medbin = zeros(Int, p.nlevels + 1) 
    p.cs.f_window = zeros(p.window)
    return true
end

function process_sample!(di, p=p)
    i = p.i
    p.cs.ds_window = [p.cs.ds_window[2:end];di]
    anomaly = 0.0
    if mod(i,p.sstep) == 0
        # Downsample
        mean = Statistics.mean(p.cs.ds_window)
        max = maximum(p.cs.ds_window)
        min = minimum(p.cs.ds_window)
        if p.cs.dsi == 1 p.cs.ds_moving_average = mean end
        p.cs.ds_moving_average = (p.cs.ds_moving_average + mean) / 2
        ds = max
        # Spike below the mean
        if abs(max - p.cs.ds_moving_average) < abs(min - p.cs.ds_moving_average)
            ds = min
        end
        if abs(ds - mean) < (0.1 * mean) # spike must be at least 10%
            ds = mean
        end
        # Normalize
        ds = ds - p.dmin
        if (p.dmax - p.dmin) != 0 ds = ds / (p.dmax - p.dmin) end
        # Discretize
        level = round(Int, ds * p.nlevels)
        ds = level / p.nlevels
        # Levelize
        p.cs.medbin[level + 1] += 1
        medpos = p.cs.dsi ÷ 2
        if p.cs.dsi == 1
            p.cs.medlevel = level
        end
        if p.cs.medlevel > level
            p.cs.belowmed += 1
        elseif p.cs.medlevel < level
            p.cs.abovemed += 1
        end
        # Not strictly a running median but close enough
        if medpos < p.cs.abovemed
            p.cs.belowmed += p.cs.medbin[p.cs.medlevel + 1]
            p.cs.medlevel += 1
            while p.cs.medbin[p.cs.medlevel+1] == 0
                p.cs.medlevel += 1
            end
            p.cs.abovemed -= p.cs.medbin[p.cs.medlevel + 1]
        elseif medpos < p.cs.belowmed
            p.cs.abovemed += p.cs.medbin[p.cs.medlevel + 1]
            p.cs.medlevel -= 1
            while p.cs.medbin[p.cs.medlevel + 1] == 0
                p.cs.medlevel -= 1
            end
            p.cs.belowmed -= p.cs.medbin[p.cs.medlevel + 1]
        end
        med = p.cs.medlevel / p.nlevels
        if Base.abs(ds - med) < p.discretize_chomp
            ds = med
        end
        # Extract features
        features = zeros(p.window * 2)
        p.cs.f_window = [p.cs.f_window[2:end];ds]
        if p.cs.dsi >= p.window
            dw = copy(p.cs.f_window)
            dw_min = minimum(dw)
            dw = dw .- dw_min
            dw_max = maximum(dw)
            if dw_max != 0 dw = dw ./ dw_max end
            fw = dwt(dw, p.wavelett)
            fw_min = minimum(fw)
            fw = (fw .- fw_min)
            fw_max = maximum(fw)
            if fw_max != 0 fw = fw ./ fw_max end
            features = [fw;p.cs.f_window]
        end
        anomaly = process_features!(features, p.cs.dsi, p)
        p.cs.dsi += 1
    end
    p.i += 1
    return anomaly
end

function process_features!(f, i, p)
    anomaly = 0.0
    if i <= p.trend_window
        # Here we could build a matrix instead of an array
        push!(p.cs.trend_window_f, f)
        # Batch process the probationary period
        if i == p.trend_window
            features_mat = hcat(p.cs.trend_window_f...)
            p.cs.art.config.dim = length(f)
            p.cs.art.config.dim_comp = 2 * p.cs.art.config.dim
            p.cs.art.config.setup = true
            rho = init_rho(features_mat[:,p.window:p.trend_window], p)
            update_rho!(rho, rho, p.cs.art)
            for (fi, ff) in enumerate(eachcol(features_mat))
                detect!(ff, fi, p)
            end
        end 
    else
        anomaly = detect!(f, i, p)
    end
    return anomaly
end

function detect!(f, i, p)
    update_rho_after_anomaly = (i - p.cs.last_anomaly_i) == p.mask_rho_after_anomaly
    update_rho_for_trend = (i - p.cs.last_rho_update_i) >= p.trend_window ÷ 2
    mask_after_anomaly = (i - p.cs.last_anomaly_i) <= p.mask_rho_after_anomaly
    if i > p.trend_window + p.mask_rho_after_anomaly && p.cs.mask_after_cat
        if p.cs.no_new_cat_count >= p.mask_rho_after_anomaly
            p.cs.mask_after_cat = false
        end
    end
    # The samples prior to a complete feature window are not used for training
    # Call train! anyway but don't learn - this keeps ART indexes and arrays aligned with input data
    if i < p.window
        AdaptiveResonance.train!(p.cs.art, f, learning=false)
        cat = -1
    else
        cat = AdaptiveResonance.train!(p.cs.art, f)
    end
    p.cs.no_new_cat_count = cat == -1 ? 0 : p.cs.no_new_cat_count + 1
    OnlineStats.fit!(p.cs.rho_ub_mean, p.cs.art.A[i]) # running mean
    p.cs.sim_window = [p.cs.sim_window[2:end];p.cs.art.A[i]]
    p.cs.sim_diff_window = [p.cs.sim_diff_window[2:end];p.cs.art.opts.rho_ub - p.cs.art.A[i]]
    # Store the smallest similarity during the masking window for each anomaly
    if (i - p.cs.last_anomaly_i) < p.mask_rho_after_anomaly && length(p.cs.anomaly_sim_history) > 0 
        if p.cs.art.A[i] < p.cs.anomaly_sim_history[end]
            p.cs.anomaly_sim_history[end] = p.cs.art.A[i]
        end
    end
    masking_anomaly = p.cs.mask_after_cat || mask_after_anomaly
    below_last_scale = mask_after_anomaly ? 0.90 : 0.70
    below_last = p.cs.art.A[i] < p.cs.last_anomaly_sim * below_last_scale
    anomaly_with_cat = cat == -1 && (!masking_anomaly || below_last)    
    if i > p.trend_window && anomaly_with_cat
        anomaly = confidence(p.cs.art.A[i], p.cs.art.Ae[i], p)
        push!(p.cs.anomaly_sim_history, p.cs.art.A[i])
        p.cs.last_anomaly_sim = p.cs.art.A[i]
        p.cs.last_anomaly_i = i
    else
        anomaly = 0.0
    end
    p.cs.mask_after_cat = cat == -1 || p.cs.mask_after_cat
    # ART could use supervised learning to improve here, but NAB does not allow this
    if i > p.trend_window && (update_rho_after_anomaly || update_rho_for_trend)
        min_sim_in_trend_window = minimum(p.cs.sim_window)
        new_rho_ub = OnlineStats.value(p.cs.rho_ub_mean) 
        new_rho_ub = min(0.97, new_rho_ub) # capping ub
        prev_rho_lb = p.cs.art.opts.rho_lb
        if prev_rho_lb <= min_sim_in_trend_window
            incr = (min_sim_in_trend_window - prev_rho_lb) * 0.19
            new_rho_lb = prev_rho_lb + incr
        else
            decr = 0.0
            if i > p.trend_window * 2
                below_rho_idxs = findall(x -> x > 0.05, p.cs.sim_diff_window)
                below_rho = p.cs.sim_diff_window[below_rho_idxs]
                below_rho = map(x -> min(0.37, x), below_rho)
                if length(below_rho) > 0 decr = mean(below_rho) end
            end
            decr = max(0.01, decr)
            new_rho_lb = prev_rho_lb - (decr / 2)
            if length(p.cs.anomaly_sim_history) > 0
                new_rho_lb = max(mean(p.cs.anomaly_sim_history), new_rho_lb)
            end
        end
        new_rho_lb = min(new_rho_ub, new_rho_lb)
        update_rho!(new_rho_lb, new_rho_ub, p.cs.art)
        p.cs.last_rho_update_i = i
        p.cs.mask_after_cat = true
    end
    return anomaly
end

function confidence(features_sim, energy_sim, p)
    features_sim = min(0.999, features_sim)
    ub = ((1 - features_sim) - (1-p.cs.art.opts.rho_ub))/(1-features_sim)
    lb = ((1 - features_sim) - (1-p.cs.art.opts.rho_lb))/(1-features_sim)
    s = (ub*0.35 + lb*0.65) + (1.0 - energy_sim)*1.5
    s = min(1.0, s)
    return round(s, digits=6) 
end

function init_rho(raw_x_optim, p)
    lengthx = length(raw_x_optim[1,:])
    raw_x_sort = raw_x_optim
    # Build a similarity matrix
    sim = ones(lengthx, lengthx)
    for i in 1:lengthx
        for j in 1:lengthx
            if j > i continue end # symmetrical so save some computation
            if i == j continue end
            sim_score = similarity(raw_x_optim[:,i], raw_x_optim[:,j])
            sim[i,j] = sim_score
            sim[j,i] = sim_score
        end
    end
    sim_sum = zeros(lengthx)
    for i in 1:lengthx
        sim_sum[i] = sum(sim[:,i])
    end
    sim_order = sortperm(sim_sum) 
    raw_x_sort = copy(raw_x_optim)
    for (i1,i2) in enumerate(sim_order)
        raw_x_sort[:,i1] = raw_x_optim[:,i2]
    end
    # Find initial rho
    art = AdaptiveResonance.DVFA()
    art.config = p.cs.art.config
    opt_rho = p.initial_rho
    update_rho!(opt_rho, opt_rho, art) 
    AdaptiveResonance.train!(art, raw_x_sort)
    return mean(art.A[p.trend_window ÷ 2:end])
end

# Close to cosine similarity
function similarity(t1, t2)
    nt1 = norm(t1)
    nt2 = norm(t2)
    if nt1 == 0.0
        s = 1 - nt2
    else
        if nt2 == 0.0
            s = 1 - nt1
        else
            s = sum(t1 .* t2) / (nt1 * nt2)
        end
    end
    return s
end

function update_rho!(rho_lb, rho_ub, art)
    art.opts.rho_lb = Float64.(rho_lb)
    art.opts.rho_ub = Float64.(rho_ub)
end

end
