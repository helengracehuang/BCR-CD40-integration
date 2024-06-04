using Interpolations;
using Distributions, Random;

# Define shapes of input signals : step, exponential, Hill functions
#-------------------------------------------------------------------------
function specifyIKKCurve(; IKK_shape = "CD40") # Deprecated
    time = [0, 15, 30, 45, 60, 90, 120, 240, 360, 600, 1440, 1440*2, 1440*3, 1440*4, 1440*5, 1440*6, 1440*10] /60;
    if (IKK_shape == "CpG250")
        values = [1, 2, 7, 13, 19, 32, 29, 27, 25, 24, 22, 20, 19, 14, 11, 5, 1] /100;
    elseif (IKK_shape == "CpG10")
        values = [1, 2, 7, 13, 19, 32, 27, 25, 24, 22, 20, 10, 7, 5, 3, 2, 1] /100;
    elseif (IKK_shape == "CpG2")
        values = [1, 2, 7, 13, 19, 32, 29, 27, 25, 24, 22, 20, 7, 3, 2, 1, 1] /100;
    elseif (IKK_shape == "LPS") # LPS45p in review model
        values = [1, 2, 7, 13, 19, 25, 7, 6, 3, 1, 1] /100;
        time = [0, 15, 30, 45, 60, 90, 120, 240, 360, 600, 2400] /60;
    elseif (IKK_shape == "IgM")
        values = [1, 5, 15, 25, 45, 65, 75, 100, 90, 75, 50, 45, 20, 15, 2*5, 1*10, 1*10, 1*10, 1*10] /100;
        time = [0, 5, 10, 15, 20, 25, 30, 60, 90, 120, 180, 210, 240, 360, 480, 600, 1440, 2160, 2880] /60;
    elseif (IKK_shape == "BAFF")
        values = [1, 2*2, 3*2, 5*2, 7*2, 8*2, 9*2, 10*2, 9*2, 8*2, 7*2, 6*1.5, 5, 4, 1, 1, 1, 1, 1] /100;
        time = [0, 5, 10, 15, 20, 25, 30, 60, 90, 120, 180, 210, 240, 360, 480, 600, 1440, 2160, 2880] /60;
    elseif (IKK_shape == "CD40") # TENTATIVE CD40 curve!!!!!!!!!!!!!!!
        values = [1, 2, 7*2, 19*2, 13*2, 7*2, 7, 6, 3, 1, 1] /100;
        time = [0, 15, 30, 45, 60, 90, 120, 240, 360, 600, 2400] /60;
    else
        print("Please input the correct IKK shape!");
        return nothing;
    end
    return LinearInterpolation(time, values, extrapolation_bc=Line());
end

function specifyNIKCurve(; NIK_shape = "CD40") # Deprecated
    if (NIK_shape == "CD40")
        values = [10, 10, 12, 25, 40, 50, 70, 80, 90, 95, 100] /100;
        time = [0, 120, 240, 360, 480, 960, 1440, 1920, 2400, 2880, 3360] /60;
    else
        print("Please input the correct IKK shape!");
        return nothing;
    end
    return LinearInterpolation(time, values, extrapolation_bc=Line());
end

# Define function to evaluate general Michaelis-Menten expression
#------------------------------------------------------------------
function michaelisMenten(; K::Float64, X::Float64)
    ans::Float64 = (X / (K + X));
    return ans;
end

# Define function to evaluate general Hill activation expression
#------------------------------------------------------------------
function hillActivation(; Ka::Float64, na::Float64, X::Float64)
    ans::Float64 = ((X^na)/((Ka^na) + (X^na)));
    return ans;
end

# Define function to evaluate general Hill repression expression
#------------------------------------------------------------------
function hillRepression(; Kr::Float64, nr::Float64, X::Float64)
    ans::Float64 = ((Kr^nr)/((Kr^nr) + (X^nr)));
    return ans;
end

# Define function to evaluate general Hill induction expression (2 params)
#------------------------------------------------------------------
function hillInduction(w1::Int64, w2::Int64, X1::Float64, X2::Float64; Kd::Float64, hill::Float64)
    ans::Float64 = (1 + (w1 * abs(X1 / Kd) ^ hill) + (w2 * abs(X2 / Kd) ^ hill)) / (1 + (abs(X1 / Kd) ^ hill) + (abs(X2 / Kd) ^ hill));
    return ans;
end

# Define function to evaluate general Hill induction expression (3 params)
#------------------------------------------------------------------
function hillInduction(w1::Int64, w2::Int64, w3::Int64, X1::Float64, X2::Float64, X3::Float64; Kd::Float64, hill::Float64)
    ans::Float64 = (1 + (w1 * abs(X1 / Kd) ^ hill) + (w2 * abs(X2 / Kd) ^ hill) + (w3 * abs(X3 / Kd) ^ hill)) / (1 + (abs(X1 / Kd) ^ hill) + (abs(X2 / Kd) ^ hill) + (abs(X3 / Kd) ^ hill));
    return ans;
end

# Define function to evaluate general Hill induction expression (4 params)
#------------------------------------------------------------------
function hillInduction(w1::Int64, w2::Int64, w3::Int64, w4::Int64, X1::Float64, X2::Float64, X3::Float64, X4::Float64; Kd::Float64, hill::Float64)
    ans::Float64 = (1 + (w1 * abs(X1 / Kd) ^ hill) + (w2 * abs(X2 / Kd) ^ hill) + (w3 * abs(X3 / Kd) ^ hill) + (w4 * abs(X4 / Kd) ^ hill)) / (1 + (abs(X1 / Kd) ^ hill) + (abs(X2 / Kd) ^ hill) + (abs(X3 / Kd) ^ hill)+ (abs(X4 / Kd) ^ hill));
    return ans;
end

# Define function to evaluate general Hill induction expression (4 params)
#------------------------------------------------------------------
function hillInduction_NFkB(w1::Int64, w2::Int64, w3::Int64, w4::Int64, X1::Float64, X2::Float64, X3::Float64, X4::Float64; Kd::Float64, hill::Float64)
    ans::Float64 = ((w1 * abs(X1 / Kd) ^ hill) + (w2 * abs(X2 / Kd) ^ hill) + (w3 * abs(X3 / Kd) ^ hill) + (w4 * abs(X4 / Kd) ^ hill)) / (1 + (abs(X1 / Kd) ^ hill) + (abs(X2 / Kd) ^ hill) + (abs(X3 / Kd) ^ hill)+ (abs(X4 / Kd) ^ hill));
    return ans;
end

# For p100: Define function to evaluate general Hill induction expression (4 params)
#------------------------------------------------------------------
function hillInduction_p100(w1::Float64, w2::Float64, w3::Float64, w4::Float64, X1::Float64, X2::Float64, X3::Float64, X4::Float64; Kd::Float64, hill::Float64)
    ans::Float64 = ((w1 * abs(X1 / Kd) ^ hill) + (w2 * abs(X2 / Kd) ^ hill) + (w3 * abs(X3 / Kd) ^ hill) + (w4 * abs(X4 / Kd) ^ hill)) / (1 + (abs(X1 / Kd) ^ hill) + (abs(X2 / Kd) ^ hill) + (abs(X3 / Kd) ^ hill)+ (abs(X4 / Kd) ^ hill));
    return ans;
end

# Rate computation functions for variable rates
#------------------------------------------------
function RateCycBDecay(concentration, Srates) # V2 in Shokhirev2015
    ans::Float64 = ((Srates[CYCB, DECAY] * concentration[CDC20P]) + (Srates[CYCB, BASALDECAY] * (1.0-concentration[CDH1])) + (Srates[CYCB, PHOSPHORYLATION] * concentration[CDH1]));
    return ans;
end
#------------------------------------------------
function RateCdh1Decay(concentration, Srates) # V4 in Shokhirev2015
    ans::Float64 = (Srates[CDH1, DEPHOSPHORYLATION] * ((0.3 * concentration[CYCA]) + (1.0 * concentration[CYCB])));
    return ans;
end
#------------------------------------------------
function RateP27Decay(concentration, Srates) # V6 in Shokhirev2015
    ans::Float64 = (Srates[P27, BASALDECAY] + (Srates[P27, PHOSPHORYLATION] * ((0.5*concentration[CYCA]) + (1.0*concentration[CYCB]) + (0.5*concentration[CYCE]))));
    return ans;
end
#------------------------------------------------
function RateCycEDecay(concentration, Srates) # V8 in Shokhirev2015
    ans::Float64 = (Srates[CYCE, DECAY] + (Srates[CYCE, PHOSPHORYLATION] * (concentration[CYCA] + concentration[CYCE] + 0.05*concentration[CYCB])/(0.15 + concentration[CYCE] + concentration[CYCEP27])));
    return ans;
end

# Define parameter distribution functions
#--------------------------------------------
function distribute_params(rates)
    # 1: initialConc ~ lognormal(log(stdConc), sqrt(log(CV^2)))
    # 2: scale (no change)
    # 3-16: rates ~ normal(stdRates, stdRates * CV) ~ stdRates * normal(1, CV)

    New_rates = similar(rates);
    @inbounds for i in 1:TOTAL_SPECIES
        # 1: initial concentrations
        if FOUNDER_CONC_CV > 0
            mu = log(rates[i, INITIALCONC]);
            sigma = sqrt(log(FOUNDER_CONC_CV*FOUNDER_CONC_CV+1));
            New_rates[i, INITIALCONC] = rand(LogNormal(mu, sigma));
        else
            New_rates[i, INITIALCONC] = rates[i, INITIALCONC];
        end

        # 2: scale
        New_rates[i, SCALE] = rates[i, SCALE];

        # 3-end: rate constants (synthesis, decay, etc.)
        if i <= 20 # receptor module is distributed differently from other modules
            New_rates[i, 3:end] .= rates[i, 3:end] .* [max(0, rand(Normal(1, FOUNDER_RECEPTOR_CV))) for i in 1:(TOTAL_PARAMS-2)];
        else
            New_rates[i, 3:end] .= rates[i, 3:end] .* [max(0, rand(Normal(1, FOUNDER_RATE_CV))) for i in 1:(TOTAL_PARAMS-2)];
        end
    end
    return New_rates;
end

function distribute_daughter_params(finalConcentration, rates, division_factor, which_daughter)
    New_rates = similar(rates)
    # 1: initial concentrations (inheritted from parent w/ asymmetric division)
    # all concentrations redistributed according to concentration CV
    if DAUGHTER_CONC_CV > 0
        @inbounds for i in 1:TOTAL_SPECIES
            mu = log(finalConcentration[i]);
            sigma = sqrt(log(DAUGHTER_CONC_CV*DAUGHTER_CONC_CV+1));
            New_rates[i, INITIALCONC] = rand(LogNormal(mu, sigma));
        end
    else
        New_rates[:, INITIALCONC] .= finalConcentration;
    end

    # asymmetric division partition cells such that if one daughter is large the other is small
    if which_daughter == 1
        @. New_rates[:, INITIALCONC] = (@view New_rates[:, INITIALCONC]) * division_factor;
    elseif which_daughter == 2
        @. New_rates[:, INITIALCONC] = (@view New_rates[:, INITIALCONC]) * (2-division_factor);
    else
        print("Please input the correct daughter index! There can only be 2 daughter cells.")
    end

    # restore CDH1 & cell machinery and mass divide by 2
    New_rates[GEN, INITIALCONC] = finalConcentration[GEN];
    New_rates[CDH1, INITIALCONC] = 1.0;
    New_rates[GM, INITIALCONC] = New_rates[GM, INITIALCONC] / 2;
    # New_rates[MASS, INITIALCONC] = New_rates[MASS, INITIALCONC] / 2; # Already done by the callback


    # empty out the contents of the nucleus (IkBs, NFkB monomers, NFkB dimers, IkB-NFkB complexes)
    @. New_rates[IKBA:IKBD, INITIALCONC] = ((@view New_rates[IKBA:IKBD, INITIALCONC]) + (@view New_rates[NIKBA:NIKBD, INITIALCONC])) / 2;
    @. New_rates[NIKBA:NIKBD, INITIALCONC] = (@view New_rates[IKBA:IKBD, INITIALCONC]);
    @. New_rates[RELA:P52, INITIALCONC] = ((@view New_rates[RELA:P52, INITIALCONC]) + (@view New_rates[NRELA:NP52, INITIALCONC])) / 2;
    @. New_rates[NRELA:NP52, INITIALCONC] = (@view New_rates[RELA:P52, INITIALCONC]);

    @. New_rates[AA:P52P52, INITIALCONC] = ((@view New_rates[AA:P52P52, INITIALCONC]) + (@view New_rates[NAA:NP52P52, INITIALCONC])) / 2;
    @. New_rates[NAA:NP52P52, INITIALCONC] = (@view New_rates[AA:P52P52, INITIALCONC]);

    @. New_rates[IKBAAA:IKBDAA, INITIALCONC] = ((@view New_rates[IKBAAA:IKBDAA, INITIALCONC]) + (@view New_rates[NIKBAAA:NIKBDAA, INITIALCONC])) / 2;
    @. New_rates[NIKBAAA:NIKBDAA, INITIALCONC] = (@view New_rates[IKBAAA:IKBDAA, INITIALCONC]);
    @. New_rates[IKBAA50:IKBDA50, INITIALCONC] = ((@view New_rates[IKBAA50:IKBDA50, INITIALCONC]) + (@view New_rates[NIKBAA50:NIKBDA50, INITIALCONC])) / 2;
    @. New_rates[NIKBAA50:NIKBDA50, INITIALCONC] = (@view New_rates[IKBAA50:IKBDA50, INITIALCONC]);
    @. New_rates[IKBAA52:IKBDA52, INITIALCONC] = ((@view New_rates[IKBAA52:IKBDA52, INITIALCONC]) + (@view New_rates[NIKBAA52:NIKBDA52, INITIALCONC])) / 2;
    @. New_rates[NIKBAA52:NIKBDA52, INITIALCONC] = (@view New_rates[IKBAA52:IKBDA52, INITIALCONC]);

    @. New_rates[IKBAB50:IKBDB50, INITIALCONC] = ((@view New_rates[IKBAB50:IKBDB50, INITIALCONC]) + (@view New_rates[NIKBAB50:NIKBDB50, INITIALCONC])) / 2;
    @. New_rates[NIKBAB50:NIKBDB50, INITIALCONC] = (@view New_rates[IKBAB50:IKBDB50, INITIALCONC]);
    @. New_rates[IKBAB52:IKBDB52, INITIALCONC] = ((@view New_rates[IKBAB52:IKBDB52, INITIALCONC]) + (@view New_rates[NIKBAB52:NIKBDB52, INITIALCONC])) / 2;
    @. New_rates[NIKBAB52:NIKBDB52, INITIALCONC] = (@view New_rates[IKBAB52:IKBDB52, INITIALCONC]);
    # @. New_rates[IKBA5050:IKBD5050, INITIALCONC] = ((@view New_rates[IKBA5050:IKBD5050, INITIALCONC]) + (@view New_rates[NIKBA5050:NIKBD5050, INITIALCONC])) / 2;
    # @. New_rates[NIKBA5050:NIKBD5050, INITIALCONC] = (@view New_rates[IKBA5050:IKBD5050, INITIALCONC]);
    @. New_rates[IKBAC50:IKBDC50, INITIALCONC] = ((@view New_rates[IKBAC50:IKBDC50, INITIALCONC]) + (@view New_rates[NIKBAC50:NIKBDC50, INITIALCONC])) / 2;
    @. New_rates[NIKBAC50:NIKBDC50, INITIALCONC] = (@view New_rates[IKBAC50:IKBDC50, INITIALCONC]);
    @. New_rates[IKBAC52:IKBDC52, INITIALCONC] = ((@view New_rates[IKBAC52:IKBDC52, INITIALCONC]) + (@view New_rates[NIKBAC52:NIKBDC52, INITIALCONC])) / 2;
    @. New_rates[NIKBAC52:NIKBDC52, INITIALCONC] = (@view New_rates[IKBAC52:IKBDC52, INITIALCONC]);

    # 2-end: scale & rate constants remain the same as parents
    New_rates[:, 2:end] .= (@view rates[:, 2:end]);

    # 3-end: rate constants (synthesis, decay, etc.) are redistributed accroding to rate CV
    if DAUGHTER_RATE_CV > 0
        @. New_rates[:, 3:end] = (@view New_rates[:, 3:end]) * [[max(0, rand(Normal(1, DAUGHTER_RATE_CV))) for i in 1:(TOTAL_PARAMS-2)] for j in 1:TOTAL_SPECIES];
    end
    return New_rates;
end
