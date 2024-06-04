# Hari's version
# include("HelperFunctions.jl");
# Compute reaction fluxes for proliferation (cell-cycle) module
#--------------------------------------------

function computeProlifFluxes!(concentration, reactionFlux, Srates, phase)
    @inbounds begin

        #-------------------------------------------------
        # MODULE 1: Growth and division (tMyc-Myc) module
        #-------------------------------------------------
        # D01. 0 --NFkB--> tMyc : Myc transcript induction by NFkB
        # reactionFlux[TMYC, SYNTHESIS] = Srates[TMYC, SYNTHESIS] * hillActivation(Ka = 1.0, na = 2.0, X = ((concentration[NA50] + concentration[NC50])/MYCTHR));
        # reactionFlux[TMYC, SYNTHESIS] = Srates[TMYC, SYNTHESIS] * hillActivation(Ka = 1.0, na = 2.0, X = ((0.45*concentration[NA50] + 0.45*concentration[NC50] + 0.1*concentration[IKK])/MYCTHR));
        reactionFlux[TMYC, SYNTHESIS] = Srates[TMYC, SYNTHESIS] * hillActivation(Ka = 1.0, na = 2.0, X = ((0.45*(concentration[NA50]+concentration[NA52]) + 0.45*(concentration[NC50]+concentration[NC52]) + 0.1*concentration[IKK])/MYCTHR));
        # D02. tMyc ----> 0 : mRNA degradation
        reactionFlux[TMYC, DECAY] = Srates[TMYC, DECAY] * concentration[TMYC];
        # D03. tMyc ----> Myc : Synthesis by translation
        reactionFlux[MYC, TRANSLATION] = EPS * Srates[MYC, TRANSLATION] * concentration[TMYC];
        # initialize Myc : Synthesis
        reactionFlux[MYC, SYNTHESIS] = 0;
        # D04. 0 --IRF4--> Myc : Induction by IRF4
        reactionFlux[MYC, SYNTHESIS] += 0;
        # reactionFlux[MYC, SYNTHESIS] += Srates[MYC, SYNTHESIS] * hillActivation(Ka = Srates[MYC].Ka[IRF4], na = Srates[MYC].na[IRF4], X = concentration[IRF4]);
        # D05. Myc ----> 0 : Protein degradation
        reactionFlux[MYC, DECAY] = Srates[MYC, DECAY] * concentration[MYC];

        #-----------------------------------------------------------
        # Execution layer : cell cycle module
        #-----------------------------------------------------------

        #----------------------------------------
        # MODULE 2: Cyclin D module
        #----------------------------------------
        # D01. 0 --Myc--> CycD : protein induction by Myc
        reactionFlux[CYCD, SYNTHESIS] = EPS * Srates[CYCD, SYNTHESIS] * michaelisMenten(K = 0.15, X = concentration[MYC]);
        # D02. CycD ----> 0 : protein degradation
        reactionFlux[CYCD, DECAY] = Srates[CYCD, DECAY] * concentration[CYCD];
        # D03. CycD + p27 ----> CycD:p27 : complex formation
        reactionFlux[CYCDP27, ASSOCIATION] = Srates[CYCDP27, ASSOCIATION] * concentration[CYCD] * concentration[P27];
        # D04. CycD:p27 ----> CycD + p27 : complex dissociation
        reactionFlux[CYCDP27, DISSOCIATION] = Srates[CYCDP27, DISSOCIATION] * concentration[CYCDP27];
        # D05. CycD:p27 ----> p27 : protein degradation
        reactionFlux[CYCDP27, DECAY] = Srates[CYCD, DECAY] * concentration[CYCDP27];

        #----------------------------------------
        # MODULE 3: Cyclin E module
        #----------------------------------------
        # E01. 0 --E2F--> CycE : protein synthesis induced by E2F
        reactionFlux[CYCE, SYNTHESIS] = EPS * (Srates[CYCE, BASALSYNTHESIS] + Srates[CYCE, SYNTHESIS] * michaelisMenten(K = 0.1, X = concentration[E2F]));
        # E02. CycE ----> 0 : protein degradation
        reactionFlux[CYCE, DECAY] = RateCycEDecay(concentration, Srates) * concentration[CYCE];
        # E03. CycE + p27 ----> CycE:p27 : complex formation
        reactionFlux[CYCEP27, ASSOCIATION] = Srates[CYCEP27, ASSOCIATION] * concentration[CYCE] * concentration[P27];
        # E04. CycE:p27 ----> CycE + p27 : complex dissociation
        reactionFlux[CYCEP27, DISSOCIATION] = Srates[CYCEP27, DISSOCIATION] * concentration[CYCEP27];
        # E05. CycE:p27 --CycE/A/B--> p27 : protein degradation
        reactionFlux[CYCEP27, DECAY] = RateCycEDecay(concentration, Srates) * concentration[CYCEP27];

        #----------------------------------------
        # MODULE 4: E2F-Rb module
        #----------------------------------------
        # F01. 0 --Myc--> tE2F : E2F transcript induction by Myc
        reactionFlux[TE2F, BASALSYNTHESIS] = Srates[TE2F, BASALSYNTHESIS] * michaelisMenten(K = 2.5, X = concentration[MYC]);
        # F02. 0 --Myc&&E2F--> tE2F : E2F transcript induction by synergistic Myc and E2F
        reactionFlux[TE2F, SYNTHESIS] = Srates[TE2F, SYNTHESIS] * michaelisMenten(K = 0.15, X = concentration[MYC]) * michaelisMenten(K = 0.15, X = concentration[E2F]);
        # F03. tE2F ----> 0 : E2F transcript degradation
        reactionFlux[TE2F, DECAY] = Srates[TE2F, DECAY] * concentration[TE2F];
        # F04. tE2F ----> E2F : translation
        reactionFlux[E2F, TRANSLATION] = EPS * Srates[E2F, TRANSLATION] * concentration[TE2F];
        # F05. E2F --CycA,CycB--> 0 : free E2F degradation mediated by Cyclin A,B
        reactionFlux[E2F, DECAY] = (Srates[E2F, BASALDECAY] + Srates[E2F, DECAY] * (concentration[CYCA] + concentration[CYCB])) * (concentration[E2F]);
        #----------------------------------------
        # F06. 0 ----> Rb : protein synthesis
        reactionFlux[RB, SYNTHESIS] = EPS * Srates[RB, SYNTHESIS];
        # F07. Rb ----> 0 : protein degradation
        reactionFlux[RB, DECAY] = Srates[RB, DECAY] * concentration[RB];
        # F08. Rb --CycD,CycE,CycA,CycB--> ppRb : phosphorylation mediated by all cyclins
        reactionFlux[RB, PHOSPHORYLATION] = Srates[RB, PHOSPHORYLATION] * michaelisMenten(K = 0.92, X = concentration[RB]) * ((3.3 * (concentration[CYCD] + concentration[CYCDP27])) + (5.0 * concentration[CYCE]) + (3.0 * concentration[CYCA]) + (5.0 * concentration[CYCB]));
        # F09. ppRb ----> Rb : dephosphorylation
        reactionFlux[PPRB, DEPHOSPHORYLATION] = (20.0 * concentration[PPRB]) / (1.0 + (25.0*(concentration[CYCE] + concentration[CYCA]) + 2.0*concentration[CYCB]));
        # F10. ppRb ----> 0 : phosphorylated protein degradation
        reactionFlux[PPRB, DECAY] = Srates[PPRB, DECAY] * concentration[PPRB];
        #----------------------------------------
        # F11. E2F + Rb ----> E2F:Rb : complex formation
        reactionFlux[E2FRB, ASSOCIATION] = Srates[E2FRB, ASSOCIATION] * concentration[E2F] * concentration[RB];
        # F12. E2F:Rb --CycD,CycE,CycA,CycB--> E2F + ppRb : complex dissociation due to Rb phosphorylation mediated by all cyclins
        reactionFlux[E2FRB, PHOSPHORYLATION] = Srates[RB, PHOSPHORYLATION] * michaelisMenten(K = 0.92, X = concentration[E2FRB]) * ((3.3 * (concentration[CYCD] + concentration[CYCDP27])) + (5.0 * concentration[CYCE]) + (3.0 * concentration[CYCA]) + (5.0 * concentration[CYCB]));
        # F13. E2F:Rb ----> 0 : complex degradation
        reactionFlux[E2FRB, DECAY] = Srates[E2FRB, DECAY] * concentration[E2FRB];

        #----------------------------------------
        # MODULE 5: Growth module
        #----------------------------------------
        # G01. 0 ----> GM : general machinery synthesis : Basal growth to maintain size => isolated
        reactionFlux[GM, BASALSYNTHESIS] = 0.01 * Srates[GM, SYNTHESIS];
        # G02. 0 --Myc&&Rb--> GM : general machinery synthesis : Active/Total Rb < 0.8 ? Additional non-linear growth : No additional growth
        reactionFlux[GM, SYNTHESIS] = (((concentration[RB] + concentration[E2FRB]) < (4.0*concentration[PPRB])) ? 0.0 : 0.99 * Srates[GM, SYNTHESIS] * (2.8*concentration[MASS]^2)*exp(-0.14*concentration[MASS]^2)) * hillActivation(Ka = 1.0, na = 2.0, X = (concentration[MYC]/GROWTHR));
        # G03. GM ----> 0 : general machinery degradation
        reactionFlux[GM, DECAY] = Srates[GM, DECAY] * concentration[GM];
        # G04. 0 --GM--> Mass : anabolism
        reactionFlux[MASS, SYNTHESIS] = EPS * Srates[MASS, SYNTHESIS] * concentration[GM];
        # G05. Mass ----> 0 : catabolism
        reactionFlux[MASS, DECAY] = Srates[MASS, DECAY] * concentration[MASS];

        #----------------------------------------
        # MODULE 6: Cyclin A module
        #----------------------------------------
        # H01. 0 --E2F&&Myc--> CycA : protein synthesis induced by E2F : GrowthMultiplier = hillActivation(Ka = 1, na = 2, X = (concentration[MYC]/GROWTHR))
        reactionFlux[CYCA, SYNTHESIS] = concentration[MASS] <= 0.5 ? 0.0 : (EPS * Srates[CYCA, SYNTHESIS] * concentration[E2F] * concentration[MASS]);
        # H02. CycA --Cdc20--> 0 : protein degradation
        reactionFlux[CYCA, DECAY] = Srates[CYCA, DECAY] * concentration[CDC20] * concentration[CYCA];
        # H03. CycA + p27 ----> CycA:p27 : complex formation
        reactionFlux[CYCAP27, ASSOCIATION] = Srates[CYCAP27, ASSOCIATION] * concentration[CYCA] * concentration[P27];
        # H04. CycA:p27 ----> CycA + p27 : complex dissociation
        reactionFlux[CYCAP27, DISSOCIATION] = Srates[CYCAP27, DISSOCIATION] * concentration[CYCAP27];
        # H05. CycA:p27 --CycE/A/B--> p27 : protein degradation
        reactionFlux[CYCAP27, DECAY] = Srates[CYCA, DECAY] * concentration[CDC20] * concentration[CYCAP27];

        #----------------------------------------
        # MODULE 7: p27 module
        #----------------------------------------
        # I01. 0 ----> p27 : protein synthesis
        reactionFlux[P27, SYNTHESIS] = EPS * Srates[P27, SYNTHESIS];
        # I02. p27 --CycE/A/B--> 0 : protein degradation
        # I03. CycD:p27 --CycE/A/B--> CycD : protein degradation
        # I04. CycE:p27 --CycE/A/B--> CycE : protein degradation
        # I05. CycA:p27 --CycE/A/B--> CycA : protein degradation
        reactionFlux[P27, DECAY] = RateP27Decay(concentration, Srates); # *** Write as an overall rate and multiply with appropriate concentrations as required ***

        #----------------------------------------
        # MODULE 8: Cyclin B module
        #----------------------------------------
        # J01. 0 --CycB--> CycB : protein synthesis
        reactionFlux[CYCB, SYNTHESIS] = EPS * (Srates[CYCB, BASALSYNTHESIS] + Srates[CYCB, SYNTHESIS] * hillActivation(Ka = 0.1, na = 2.0, X = concentration[CYCB]));
        # J02. CycB --Cdc20,Cdh1--> 0 : protein degradation
        reactionFlux[CYCB, DECAY] = RateCycBDecay(concentration, Srates) * concentration[CYCB];
        # J03. 0 ----> PPX : phosphatase synthesis
        reactionFlux[PPX, SYNTHESIS] = EPS * Srates[PPX, SYNTHESIS];
        # J04. PPX ----> 0 : phosphatase degradation
        reactionFlux[PPX, DECAY] = Srates[PPX, DECAY] * concentration[PPX];
        # J05. 0 --CycB--> IEP : intermediate enzyme synthesis
        reactionFlux[IEP, PHOSPHORYLATION] = Srates[IEP, PHOSPHORYLATION] * concentration[CYCB] * michaelisMenten(K = 0.01, X = (1.0 - concentration[IEP]));
        # J06. IEP --PPX--> 0 : intermediate enzyme degradation
        reactionFlux[IEP, DEPHOSPHORYLATION] = Srates[IEP, DEPHOSPHORYLATION] * concentration[PPX] * michaelisMenten(K = 0.01, X = concentration[IEP]);
        # J07. 0 --CycB--> Cdc20T : protein synthesis
        reactionFlux[CDC20, SYNTHESIS] = EPS * (Srates[CDC20, BASALSYNTHESIS] + Srates[CDC20, SYNTHESIS] * concentration[CYCB]);
        # J08. Cdc20T ----> 0 : protein degradation
        reactionFlux[CDC20, DECAY] = Srates[CDC20, DECAY] * concentration[CDC20];
        # J09. 0 --IEP--> Cdc20 : protein synthesis
        reactionFlux[CDC20, PHOSPHORYLATION] = Srates[CDC20, PHOSPHORYLATION] * concentration[IEP] * michaelisMenten(K = 0.005, X = (concentration[CDC20] - concentration[CDC20P]));
        # J10. Cdc20 --Cdc20--> 0 : protein autocatalytic degradation
        reactionFlux[CDC20P, DEPHOSPHORYLATION] = Srates[CDC20P, DEPHOSPHORYLATION] * michaelisMenten(K = 0.005, X = concentration[CDC20P]);
        # J11. Cdc20 ----> 0 : protein degradation
        reactionFlux[CDC20P, DECAY] = Srates[CDC20, DECAY] * concentration[CDC20P];
        # J12. 0 --Cdc20--> Cdh1 : protein activation
        reactionFlux[CDH1, PHOSPHORYLATION] = concentration[CDH1] > 1.0 ? 0.0 : ((Srates[CDH1, BASALSYNTHESIS] + Srates[CDH1, PHOSPHORYLATION] * concentration[CDC20P]) * michaelisMenten(K = 0.01, X = (1.0 - concentration[CDH1])));
        # J13. Cdh1 --CycE/A/B--> 0 : protein deactivation
        reactionFlux[CDH1, DEPHOSPHORYLATION] = RateCdh1Decay(concentration, Srates) * michaelisMenten(K = 0.01, X = concentration[CDH1]);

    end
    nothing
end

function computeProlifNettFluxes!(nettFlux, concentration, reactionFlux, Srates, phase)
    #------------------------------------------------
    # Compute proliferation reaction fluxes
    computeProlifFluxes!(concentration, reactionFlux, Srates, phase);

    # Compute proliferation net fluxes
    @inbounds begin
        # 1 : tMyc
        nettFlux[TMYC] = reactionFlux[TMYC, SYNTHESIS] - reactionFlux[TMYC, DECAY];
        # 2 : Myc
        nettFlux[MYC] = reactionFlux[MYC, TRANSLATION] + reactionFlux[MYC, SYNTHESIS] - reactionFlux[MYC, DECAY];
        #------------------------------------------------
        # 3 : tE2F
        nettFlux[TE2F] = reactionFlux[TE2F, BASALSYNTHESIS] + reactionFlux[TE2F, SYNTHESIS] - reactionFlux[TE2F, DECAY];
        # 4 : E2F
        nettFlux[E2F] = reactionFlux[E2F, TRANSLATION] - reactionFlux[E2F, DECAY] - reactionFlux[E2FRB, ASSOCIATION] + reactionFlux[E2FRB, PHOSPHORYLATION];
        # 5 : Rb
        nettFlux[RB] = reactionFlux[RB, SYNTHESIS] - reactionFlux[RB, DECAY] - reactionFlux[RB, PHOSPHORYLATION] + reactionFlux[PPRB, DEPHOSPHORYLATION] - reactionFlux[E2FRB, ASSOCIATION];
        # 6 : ppRB
        nettFlux[PPRB] = reactionFlux[RB, PHOSPHORYLATION] - reactionFlux[PPRB, DEPHOSPHORYLATION] - reactionFlux[PPRB, DECAY] + reactionFlux[E2FRB, PHOSPHORYLATION];
        # 7 : E2F:Rb
        nettFlux[E2FRB] = reactionFlux[E2FRB, ASSOCIATION] - reactionFlux[E2FRB, PHOSPHORYLATION] - reactionFlux[E2FRB, DECAY];
        # 8 : CycD
        nettFlux[CYCD] = reactionFlux[CYCD, TRANSLATION] + reactionFlux[CYCD, SYNTHESIS] - reactionFlux[CYCD, DECAY] - reactionFlux[CYCDP27, ASSOCIATION] + reactionFlux[CYCDP27, DISSOCIATION] + (reactionFlux[P27, DECAY] * concentration[CYCDP27]);
        # 9 : CycE
        nettFlux[CYCE] = reactionFlux[CYCE, SYNTHESIS] - reactionFlux[CYCE, DECAY] - reactionFlux[CYCEP27, ASSOCIATION] + reactionFlux[CYCEP27, DISSOCIATION] + (reactionFlux[P27, DECAY] * concentration[CYCEP27]);
        # 10 : CycA
        nettFlux[CYCA] = reactionFlux[CYCA, SYNTHESIS] - reactionFlux[CYCA, DECAY] - reactionFlux[CYCAP27, ASSOCIATION] + reactionFlux[CYCAP27, DISSOCIATION] + (reactionFlux[P27, DECAY] * concentration[CYCAP27]);
        # 11 : CycB
        nettFlux[CYCB] = reactionFlux[CYCB, SYNTHESIS] - reactionFlux[CYCB, DECAY];
        # 12 : p27
        nettFlux[P27] = reactionFlux[P27, SYNTHESIS] - (reactionFlux[P27, DECAY] * concentration[P27]) + reactionFlux[CYCDP27, DECAY] + reactionFlux[CYCEP27, DECAY] + reactionFlux[CYCAP27, DECAY] + reactionFlux[CYCDP27, DISSOCIATION] + reactionFlux[CYCEP27, DISSOCIATION] + reactionFlux[CYCAP27, DISSOCIATION] - reactionFlux[CYCDP27, ASSOCIATION] - reactionFlux[CYCEP27, ASSOCIATION] - reactionFlux[CYCAP27, ASSOCIATION];
        # 13 : CycD:p27
        nettFlux[CYCDP27] = reactionFlux[CYCDP27, ASSOCIATION] - reactionFlux[CYCDP27, DISSOCIATION] - reactionFlux[CYCDP27, DECAY] - (reactionFlux[P27, DECAY] * concentration[CYCDP27]);
        # 14 : CycE:p27
        nettFlux[CYCEP27] = reactionFlux[CYCEP27, ASSOCIATION] - reactionFlux[CYCEP27, DISSOCIATION] - reactionFlux[CYCEP27, DECAY] - (reactionFlux[P27, DECAY] * concentration[CYCEP27]);
        # 15 : CycA:p27
        nettFlux[CYCAP27] = reactionFlux[CYCAP27, ASSOCIATION] - reactionFlux[CYCAP27, DISSOCIATION] - reactionFlux[CYCAP27, DECAY] - (reactionFlux[P27, DECAY] * concentration[CYCAP27]);
        # 16 : IEP (active)
        nettFlux[IEP] = reactionFlux[IEP, PHOSPHORYLATION] - reactionFlux[IEP, DEPHOSPHORYLATION];
        # 17 : PPX
        nettFlux[PPX] = reactionFlux[PPX, SYNTHESIS] - reactionFlux[PPX, DECAY];
        # 18 : Cdc20 (total)
        nettFlux[CDC20] = reactionFlux[CDC20, SYNTHESIS] - reactionFlux[CDC20, DECAY]; # - reactionFlux[CDC20, PHOSPHORYLATION];
        # 19 : Cdc20P (active)
        nettFlux[CDC20P] = reactionFlux[CDC20, PHOSPHORYLATION] - reactionFlux[CDC20P, DEPHOSPHORYLATION] - reactionFlux[CDC20P, DECAY];
        # 20 : Cdh1 (active)
        nettFlux[CDH1] = reactionFlux[CDH1, PHOSPHORYLATION] - reactionFlux[CDH1, DEPHOSPHORYLATION];
        #------------------------------------------------
        # 21 : GM
        nettFlux[GM] = reactionFlux[GM, BASALSYNTHESIS] + reactionFlux[GM, SYNTHESIS] - reactionFlux[GM, DECAY];
        # 22 : Mass
        nettFlux[MASS] = reactionFlux[MASS, SYNTHESIS] - reactionFlux[MASS, DECAY];
        # 23 : Rb Growth Switch
        concentration[RBGS] = (concentration[RB] + concentration[E2FRB]) < (4*concentration[PPRB]) ? 0 : 1;
        # 24 : Generations : To be set by interrupt
        #------------------------------------------------

    end
    nothing
end
