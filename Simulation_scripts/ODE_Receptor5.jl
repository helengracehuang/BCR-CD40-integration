# successor of ODE_Receptor3.jl, compatible with main3.jl
# For sequential stimulation, make BCR-induced IKK stop at 1,3,5,8hrs but not AICD effects

# include("HelperFunctions.jl");
# BCR module do not include CBM complex in this version
# Compute reaction fluxes for CD40 module
#--------------------------------------------
function computeBCRFluxes!(concentration, reactionFlux, Srates, phase, time)
    @inbounds begin
        # MODULE 1 : BCR antigen-receptor interaction
        # 1 : Antigen basal decay
        reactionFlux[ANTIGEN, BASALDECAY] = Srates[ANTIGEN, BASALDECAY] * concentration[ANTIGEN];
        # 2 : BCR basal synthesis
        reactionFlux[BCR, BASALSYNTHESIS] = Srates[BCR, BASALSYNTHESIS];
        # 3 : BCR basal decay
        reactionFlux[BCR, BASALDECAY] = Srates[BCR, BASALDECAY] * concentration[BCR];
        # 4 : Antigen-BCR association
        reactionFlux[ABCR, ASSOCIATION] = Srates[ABCR, ASSOCIATION] * concentration[ANTIGEN] * concentration[BCR];
        # 5 : Antigen-BCR dissociation
        reactionFlux[ABCR, DISSOCIATION] = Srates[ABCR, DISSOCIATION] * concentration[ABCR];
        # 6 : Antigen-BCR basal decay
        reactionFlux[ABCR, BASALDECAY] = Srates[ABCR, BASALDECAY] * concentration[ABCR];

        # MODULE 2 : CARMA1-Bcl10-MALT1 complex dynamics
        # 7 : CBM activation (induced by ABCR)
        reactionFlux[CBM, ACTIVATION] = Srates[CBM, ACTIVATION] * concentration[ABCR] * concentration[CBM];
        # if phase == 2
        #     if BCR_END > 0 && time >= BCR_END
        #         reactionFlux[CBM, ACTIVATION] = 0.0;
        #     end
        # end
        # 8 : CBM activation (induced by activated IKK)
        reactionFlux[ACBM, SYNTHESIS] = Srates[ACBM, SYNTHESIS] * concentration[IKK] * concentration[CBM];
        # 9 : CBM deactivation
        reactionFlux[ACBM, DEACTIVATION] = Srates[ACBM, DEACTIVATION] * concentration[ACBM];
        # 10 : CBM inhibition (induced by activated IKK)
        reactionFlux[ICBM, SYNTHESIS] = Srates[ICBM, SYNTHESIS] * concentration[IKK] * concentration[ACBM];
        # 11 : CBM renewal
        reactionFlux[CBM, BASALSYNTHESIS] = Srates[CBM, BASALSYNTHESIS] * concentration[ICBM];

    end
    nothing
end

function computeCD40Fluxes!(concentration, reactionFlux, Srates, phase, time)
    @inbounds begin
        # MODULE 1 : CD40 ligand-receptor interaction
        # 1 : CD40L basal decay
        reactionFlux[CD40L, BASALDECAY] = Srates[CD40L, BASALDECAY] * concentration[CD40L];
        # 2 : CD40R basal synthesis
        reactionFlux[CD40R, BASALSYNTHESIS] = Srates[CD40R, BASALSYNTHESIS];
        # 3 : CD40R basal decay
        reactionFlux[CD40R, BASALDECAY] = Srates[CD40R, BASALDECAY] * concentration[CD40R];
        # 4 : CD40L-R association
        reactionFlux[CD40LR, ASSOCIATION] = Srates[CD40LR, ASSOCIATION] * concentration[CD40L] * concentration[CD40R];
        # 5 : CD40L-R dissociation
        reactionFlux[CD40LR, DISSOCIATION] = Srates[CD40LR, DISSOCIATION] * concentration[CD40LR];
        # 6 : CD40L-R basal decay
        reactionFlux[CD40LR, BASALDECAY] = Srates[CD40LR, BASALDECAY] * concentration[CD40LR];

        # MODULE 2 : TRAF6 dynamics
        # 7 : TRAF6 activation (induced by CD40L-R complex)
        reactionFlux[TRAF6, ACTIVATION] = Srates[TRAF6, ACTIVATION] * concentration[CD40LR] * concentration[TRAF6_OFF];
        # 8 : TRAF6 deactivation
        reactionFlux[TRAF6, DEACTIVATION] = Srates[TRAF6, DEACTIVATION] * concentration[TRAF6];
        # 9 : TAK1 activation (induced by TRAF6)
        reactionFlux[TAK1, ACTIVATION] = Srates[TAK1, ACTIVATION] * concentration[TRAF6] * concentration[TAK1];

        # MODULE 3 : Noncanonical pathway activation
        # 10 : TRAF3-cIAP complex basal synthesis
        reactionFlux[TRAF3, BASALSYNTHESIS] = Srates[TRAF3, BASALSYNTHESIS];
        # 11 : TRAF3-cIAP complex basal decay
        reactionFlux[TRAF3, BASALDECAY] = Srates[TRAF3, BASALDECAY] * concentration[TRAF3];
        # 12 : TRAF3-cIAP complex decay (induced by CD40L-R complex)
        reactionFlux[TRAF3, DECAY] = Srates[TRAF3, DECAY] * concentration[CD40LR] * concentration[TRAF3];

    end
    nothing
end

function computeKinaseFluxes!(concentration, reactionFlux, Srates, phase, time)
    @inbounds begin
        # MODULE 1 : TAK1 dynamics
        # 1 : TAK1 activation (induced by ABCR-complex)
        reactionFlux[ATAK1, ACTIVATION] = Srates[ATAK1, ACTIVATION] * concentration[ACBM] * michaelisMenten(K = 5.52, X = concentration[TAK1]);
        # 2 : TAK1 activation (induced by IKK2)
        reactionFlux[ATAK1, SYNTHESIS] = Srates[ATAK1, SYNTHESIS] * concentration[IKK2] * michaelisMenten(K = 0.102, X = concentration[TAK1]);
        # 3 : TAK1 activation (induced by IKK3) The naming is strange because there are so many TAK1 activation equations!
        reactionFlux[TAK1, DECAY] = Srates[TAK1, DECAY] * concentration[IKK3] * michaelisMenten(K = 0.102, X = concentration[TAK1]);
        # 4 : TAK1 basal activation
        reactionFlux[ATAK1, BASALSYNTHESIS] = Srates[ATAK1, BASALSYNTHESIS] * concentration[TAK1];
        # 5 : TAK1 deactivation
        reactionFlux[ATAK1, DEACTIVATION] = Srates[ATAK1, DEACTIVATION] * michaelisMenten(K = 143.0, X = concentration[ATAK1]);

        # MODULE 2 : Canonical pathway activation (IKK dynamics)
        # 6 : IKK2 activation (induced by ATAK1)
        reactionFlux[IKK2, ACTIVATION] = Srates[IKK2, ACTIVATION] * concentration[ATAK1] * michaelisMenten(K = 0.986*140-10, X = concentration[IKK_OFF]);
        # 7 : IKK2 deactivation
        reactionFlux[IKK2, DEACTIVATION] = Srates[IKK2, DEACTIVATION] * michaelisMenten(K = 0.202*140, X = concentration[IKK2]);
        # 8 : IKK3 basal activation
        reactionFlux[IKK3, BASALSYNTHESIS] = Srates[IKK3, BASALSYNTHESIS] * michaelisMenten(K = 3.56*140-10, X = concentration[IKK2]);
        # 9 : IKK3 activation (self induction)
        reactionFlux[IKK3, ACTIVATION] = Srates[IKK3, ACTIVATION] * concentration[IKK3] * michaelisMenten(K = 1.56*140-10, X = concentration[IKK2]);
        # 10 : IKK3 deactivation
        reactionFlux[IKK3, DEACTIVATION] = Srates[IKK3, DEACTIVATION] * michaelisMenten(K = 1.76*140, X = concentration[IKK3]);
        # 11 : IKK inhibition
        reactionFlux[IIKK, BASALSYNTHESIS] = Srates[IIKK, BASALSYNTHESIS] * michaelisMenten(K = 0.45*140, X = concentration[IKK3]);
        # 12 : IKK renewal
        reactionFlux[IKK_OFF, BASALSYNTHESIS] = Srates[IKK_OFF, BASALSYNTHESIS] * michaelisMenten(K = 2.6*140-10, X = concentration[IIKK]);
        # 13 : total IKK activity (IKK2 + IKK3) * IKK_modifier
        concentration[IKK] = IKK_MOD * (concentration[IKK2] + concentration[IKK3]);

        # MODULE 3 : Noncanonical pathway activation (NIK dynamics)
        # 14 : NIK basal synthesis
        reactionFlux[NIK, BASALSYNTHESIS] = Srates[NIK, BASALSYNTHESIS];
        # 15 : NIK basal decay
        reactionFlux[NIK, BASALDECAY] = Srates[NIK, BASALDECAY] * concentration[NIK];
        # 16 : NIK decay (induced by TRAF3 complex)
        reactionFlux[NIK, DECAY] = Srates[NIK, DECAY] * concentration[NIK] * hillActivation(Ka = 0.5, na = 2.0, X = concentration[TRAF3]);

    end
    nothing
end

function ReceptorNettFluxes!(nettFlux, reactionFlux)
    @inbounds begin
        # 1 : ANTIGEN
        nettFlux[ANTIGEN] = - reactionFlux[ANTIGEN, BASALDECAY] - reactionFlux[ABCR, ASSOCIATION] * SCALE_CELLULAR2MEDIA + reactionFlux[ABCR, DISSOCIATION] * SCALE_CELLULAR2MEDIA;
        # 2 : BCR
        nettFlux[BCR] = reactionFlux[BCR, BASALSYNTHESIS] - reactionFlux[BCR, BASALDECAY] - reactionFlux[ABCR, ASSOCIATION] + reactionFlux[ABCR, DISSOCIATION];
        # 3 : Antigen-BCR
        nettFlux[ABCR] = reactionFlux[ABCR, ASSOCIATION] - reactionFlux[ABCR, DISSOCIATION] - reactionFlux[ABCR, BASALDECAY];
        # 4 : inactivated CBM
        nettFlux[CBM] = - reactionFlux[CBM, ACTIVATION] - reactionFlux[ACBM, SYNTHESIS] + reactionFlux[ACBM, DEACTIVATION] + reactionFlux[CBM, BASALSYNTHESIS];
        # 5 : activated CBM
        nettFlux[ACBM] = reactionFlux[CBM, ACTIVATION] + reactionFlux[ACBM, SYNTHESIS] - reactionFlux[ACBM, DEACTIVATION] - reactionFlux[ICBM, SYNTHESIS];
        # 6 : inhibited CBM
        nettFlux[ICBM] = reactionFlux[ICBM, SYNTHESIS] - reactionFlux[CBM, BASALSYNTHESIS];

        # 1 : CD40 ligand
        nettFlux[CD40L] = - reactionFlux[CD40L, BASALDECAY] - reactionFlux[CD40LR, ASSOCIATION] * SCALE_CELLULAR2MEDIA + reactionFlux[CD40LR, DISSOCIATION] * SCALE_CELLULAR2MEDIA;
        # 2 : CD40 receptor
        nettFlux[CD40R] = reactionFlux[CD40R, BASALSYNTHESIS] - reactionFlux[CD40R, BASALDECAY] - reactionFlux[CD40LR, ASSOCIATION] + reactionFlux[CD40LR, DISSOCIATION];
        #------------------------------------------------
        # 3 : CD40L-R complex
        nettFlux[CD40LR] = reactionFlux[CD40LR, ASSOCIATION] - reactionFlux[CD40LR, DISSOCIATION] - reactionFlux[CD40LR, BASALDECAY];

        # 4 : inactivated TRAF6
        nettFlux[TRAF6_OFF] = - reactionFlux[TRAF6, ACTIVATION] - reactionFlux[TRAF6, SYNTHESIS] + reactionFlux[TRAF6, DEACTIVATION];
        # 5 : activated TRAF6
        nettFlux[TRAF6] = reactionFlux[TRAF6, ACTIVATION] + reactionFlux[TRAF6, SYNTHESIS] - reactionFlux[TRAF6, DEACTIVATION];

        # 4 : inactivated TAK1
        nettFlux[TAK1] = - reactionFlux[TAK1, ACTIVATION] - reactionFlux[ATAK1, ACTIVATION] - reactionFlux[ATAK1, SYNTHESIS] - reactionFlux[TAK1, DECAY] - reactionFlux[ATAK1, BASALSYNTHESIS] + reactionFlux[ATAK1, DEACTIVATION];;
        # 5 : activated TAK1
        nettFlux[ATAK1] = reactionFlux[TAK1, ACTIVATION] + reactionFlux[ATAK1, ACTIVATION] + reactionFlux[ATAK1, SYNTHESIS] + reactionFlux[TAK1, DECAY] + reactionFlux[ATAK1, BASALSYNTHESIS] - reactionFlux[ATAK1, DEACTIVATION];
        # 6 : IKK (IKK1)
        nettFlux[IKK_OFF] = reactionFlux[IKK_OFF, BASALSYNTHESIS] - reactionFlux[IKK2, ACTIVATION] + reactionFlux[IKK2, DEACTIVATION];
        # 7 : pIKK (IKK2)
        nettFlux[IKK2] = reactionFlux[IKK2, ACTIVATION] - reactionFlux[IKK2, DEACTIVATION] - reactionFlux[IKK3, BASALSYNTHESIS] - reactionFlux[IKK3, ACTIVATION] + reactionFlux[IKK3, DEACTIVATION];
        # 8 : ppIKK (IKK3)
        nettFlux[IKK3] = reactionFlux[IKK3, BASALSYNTHESIS] + reactionFlux[IKK3, ACTIVATION] - reactionFlux[IKK3, DEACTIVATION] - reactionFlux[IIKK, BASALSYNTHESIS];
        # 9 : pppIKK (IKK4)
        nettFlux[IIKK] = reactionFlux[IIKK, BASALSYNTHESIS] - reactionFlux[IKK_OFF, BASALSYNTHESIS];

        # 11 : TRAF3-cIAP1/2 complex
        nettFlux[TRAF3] = reactionFlux[TRAF3, BASALSYNTHESIS] - reactionFlux[TRAF3, BASALDECAY] - reactionFlux[TRAF3, DECAY];
        # 12 : NIK
        nettFlux[NIK] = reactionFlux[NIK, BASALSYNTHESIS] - reactionFlux[NIK, BASALDECAY] - reactionFlux[NIK, DECAY];

    end
    nothing
end

# time-independent ODE (pre-simulation: phase = 1)
function computeReceptorNettFluxes!(nettFlux, concentration, reactionFlux, Srates, phase)
    #------------------------------------------------
    # Compute BCR reaction fluxes
    computeBCRFluxes!(concentration, reactionFlux, Srates, phase, nothing);
    # Compute CD40 reaction fluxes
    computeCD40Fluxes!(concentration, reactionFlux, Srates, phase, nothing);
    # Compute Kinase reaction fluxes
    computeKinaseFluxes!(concentration, reactionFlux, Srates, phase, nothing);
    # Compute all net fluxes
    ReceptorNettFluxes!(nettFlux, reactionFlux);
    nothing
end

# time-dependent ODE (simulation, w/ delay: phase = 2)
function computeReceptorNettFluxes!(nettFlux, concentration, reactionFlux, Srates, phase, time)
    #------------------------------------------------
    # Compute BCR reaction fluxes
    computeBCRFluxes!(concentration, reactionFlux, Srates, phase, time);
    # Compute CD40 reaction fluxes
    computeCD40Fluxes!(concentration, reactionFlux, Srates, phase, time);
    # Compute Kinase reaction fluxes
    computeKinaseFluxes!(concentration, reactionFlux, Srates, phase, time);
    # Compute all net fluxes
    ReceptorNettFluxes!(nettFlux, reactionFlux);
    nothing
end
