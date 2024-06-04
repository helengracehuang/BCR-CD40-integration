# include("HelperFunctions.jl");
# Compute reaction fluxes for differentiation module
#--------------------------------------------

function computeDiffFluxes!(concentration, reactionFlux, Srates, phase)
    @inbounds begin

        #--------------------------------------------
        # MODULE 1: GC phenotype (Pax5-Bcl6) module
        #--------------------------------------------
        # B01. 0 ----> Pax5 : Basal synthesis
        reactionFlux[PAX5, BASALSYNTHESIS] = Srates[PAX5, BASALSYNTHESIS];
        # B02. Pax5 ----> 0 : Degradation
        reactionFlux[PAX5, DECAY] = Srates[PAX5, DECAY] * concentration[PAX5];
        # B03. Blimp1 ----| Pax5 : Repression
        reactionFlux[PAX5, REPRESSION] = hillRepression(Kr = 1.0, nr = 2.0, X = concentration[BLIMP1]);
        #--------------------------------------------
        # initialize Bcl6 : Synthesis
        reactionFlux[BCL6, SYNTHESIS] = 0
        # B04. 0 --Pax5--> Bcl6 : Synthesis
        reactionFlux[BCL6, SYNTHESIS] += Srates[BCL6, SYNTHESIS] * hillActivation(Ka = 1.0, na = 1.0, X = concentration[PAX5])
        # B05. nC50 ----> Bcl6 : Induced synthesis
        reactionFlux[BCL6, SYNTHESIS] += Srates[BCL6, SYNTHESIS] * hillActivation(Ka = 3.0, na = 1.0, X = concentration[NC50])
        # B06. Bcl6 ----> 0 : Degradation
        reactionFlux[BCL6, DECAY] = Srates[BCL6, DECAY] * concentration[BCL6]
        # B07. IRF4 && Blimp1 ----| Bcl6 : Repression : repressors taken together, cannot use normal Hill repression form for this, so hard-coding it here
        reactionFlux[BCL6, REPRESSION] = (1.0^2.0) / ((1.0^2.0) + (concentration[IRF4]^2.0) + (concentration[BLIMP1]^2.0))
        # B08. Bcl6 ----| Bcl6 : Autorepression
        # reactionFlux[BCL6].repression *= hillRepression(Kr = Srates[BCL6].Kr[BCL6], nr = Srates[BCL6].nr[BCL6], X = concentration[BCL6]);

        #-----------------------------------------------------
        # MODULE 2: Plasmablast phenotype (Blimp1-IRF4) module
        #-----------------------------------------------------
        # initialize Blimp1 : Synthesis
        reactionFlux[BLIMP1, SYNTHESIS] = 0
        # C01. IRF4 ----> Blimp1 : Induction by IRF4
        reactionFlux[BLIMP1, SYNTHESIS] += Srates[BLIMP1, SYNTHESIS] * hillActivation(Ka = 3.0, na = 2.0, X = concentration[IRF4]);
        # C02. nA50 ----> Blimp1 : Induction by nA50
        reactionFlux[BLIMP1, SYNTHESIS] += Srates[BLIMP1, SYNTHESIS] * hillActivation(Ka = 50.0, na = 2.0, X = concentration[NA50]);
        # C03. Blimp1 ----> 0 : Degradation
        reactionFlux[BLIMP1, DECAY] = Srates[BLIMP1, DECAY] * concentration[BLIMP1];
        # C04. Bcl6 ----| Blimp1 : Repression by Bcl6
        reactionFlux[BLIMP1, REPRESSION] = hillRepression(Kr = 1.0, nr = 2.0, X = concentration[BCL6]);
        #--------------------------------------------
        # C05. 0 ----> IRF4 : Basal expression
        reactionFlux[IRF4, BASALSYNTHESIS] = Srates[IRF4, BASALSYNTHESIS];
        # initialize IRF4 : Synthesis
        reactionFlux[IRF4, SYNTHESIS] = 0
        # C06. Blimp1 ----> IRF4 : Induction by Blimp1
        reactionFlux[IRF4, SYNTHESIS] += (1/3) * Srates[IRF4, SYNTHESIS] * hillActivation(Ka = 1.0, na = 2.0, X = concentration[BLIMP1]);
        # C07. nC50 ----> IRF4 : Induction by nC50
        reactionFlux[IRF4, SYNTHESIS] += (1/40) * Srates[IRF4, SYNTHESIS] * hillActivation(Ka = 3.0, na = 2.0, X = concentration[NC50]);
        # C08. nA50 ----> IRF4 : Induction by nA50
        reactionFlux[IRF4, SYNTHESIS] += (1/40) * Srates[IRF4, SYNTHESIS] * hillActivation(Ka = 50.0, na = 2.0, X = concentration[NA50]);
        # C09. Myc ----> IRF4 : Induction by Myc
        reactionFlux[IRF4, SYNTHESIS] += 0;
        # reactionFlux[IRF4, SYNTHESIS] += Srates[IRF4, SYNTHESIS] * hillActivation(Ka = Srates[IRF4].Ka[MYC], na = Srates[IRF4].na[MYC], X = concentration[MYC]);
        # C10. IRF4 ----> 0 : Degradation
        reactionFlux[IRF4, DECAY] = Srates[IRF4, DECAY] * concentration[IRF4];
        # C11. IRF-4 ----> IRF4 : Autoinduction
        #reactionFlux[IRF4, SYNTHESIS] += (0.7) * Srates[IRF4, SYNTHESIS] * hillActivation(Ka = Srates[IRF4].Ka[IRF4], na = Srates[IRF4].na[IRF4], X = concentration[IRF4]);

    end
    nothing
end

function computeDiffNettFluxes!(nettFlux, concentration, reactionFlux, Srates, phase)
    #------------------------------------------------
    # Compute differentiation reaction fluxes
    computeDiffFluxes!(concentration, reactionFlux, Srates, phase);

    # Compute differentiation net fluxes
    @inbounds begin
        # 1 : Pax5
        nettFlux[PAX5] = (reactionFlux[PAX5, BASALSYNTHESIS] * reactionFlux[PAX5, REPRESSION]) - reactionFlux[PAX5, DECAY];
        # 2 : Bcl6
        nettFlux[BCL6] = (reactionFlux[BCL6, SYNTHESIS] * reactionFlux[BCL6, REPRESSION]) - reactionFlux[BCL6, DECAY];
        # 3 : Blimp1
        nettFlux[BLIMP1] = (reactionFlux[BLIMP1, SYNTHESIS] * reactionFlux[BLIMP1, REPRESSION]) - reactionFlux[BLIMP1, DECAY];
        # 4 : IRF4
        nettFlux[IRF4] = reactionFlux[IRF4, BASALSYNTHESIS] + reactionFlux[IRF4, SYNTHESIS] - reactionFlux[IRF4, DECAY];

    end
    nothing
end
