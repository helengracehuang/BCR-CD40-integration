# 3rd edition of NFkB model (added non-canonical pathway according to Mitchell 2022 review paper)
# include("HelperFunctions.jl");
# Compute reaction fluxes for NFkB module
#--------------------------------------------

function computeNFkBFluxes!(concentration, reactionFlux, Srates, phase, delay, historicFlux, time)
    @inbounds begin
        # get historic values for nuclear NFkB concentrations (for delayed transcription)
        if phase == 2
            p = (Srates, reactionFlux, historicFlux);
            hist_NA50_025 = delay(p, time-0.25/CONVERSION; idxs=NA50);
            hist_NA52_025 = delay(p, time-0.25/CONVERSION; idxs=NA52);

            hist_NA50_075 = delay(p, time-0.75/CONVERSION; idxs=NA50);
            hist_NA52_075 = delay(p, time-0.75/CONVERSION; idxs=NA52);
            hist_NC50_075 = delay(p, time-0.75/CONVERSION; idxs=NC50);
            hist_NC52_075 = delay(p, time-0.75/CONVERSION; idxs=NC52);

            hist_NA50_1 = delay(p, time-1.0/CONVERSION; idxs=NA50);
            hist_NA52_1 = delay(p, time-1.0/CONVERSION; idxs=NA52);
            hist_NC50_1 = delay(p, time-1.0/CONVERSION; idxs=NC50);
            hist_NC52_1 = delay(p, time-1.0/CONVERSION; idxs=NC52);

            hist_NA50_40 = delay(p, time-4.0/CONVERSION; idxs=NA50);
            hist_NA52_40 = delay(p, time-4.0/CONVERSION; idxs=NA52);
            hist_NC50_40 = delay(p, time-4.0/CONVERSION; idxs=NC50);
            hist_NC52_40 = delay(p, time-4.0/CONVERSION; idxs=NC52);

            hist_NA50_120 = delay(p, time-12.0/CONVERSION; idxs=NA50);
            hist_NA52_120 = delay(p, time-12.0/CONVERSION; idxs=NA52);
            hist_NC50_120 = delay(p, time-12.0/CONVERSION; idxs=NC50);
            hist_NC52_120 = delay(p, time-4.0/CONVERSION; idxs=NC52);

            # delay(historicFlux, p, time-0.25/CONVERSION);
            # hist_NA50_025 = historicFlux[NA50];
            # hist_NA52_025 = historicFlux[NA52];
            # delay(historicFlux, p, time-0.75/CONVERSION);
            # hist_NA50_075 = historicFlux[NA50];
            # hist_NA52_075 = historicFlux[NA52];
            # hist_NC50_075 = historicFlux[NC50];
            # hist_NC52_075 = historicFlux[NC52];
            # delay(historicFlux, p, time-1.0/CONVERSION);
            # hist_NA50_1 = historicFlux[NA50];
            # hist_NA52_1 = historicFlux[NA52];
            # hist_NC50_1 = historicFlux[NC50];
            # hist_NC52_1 = historicFlux[NC52];
            # delay(historicFlux, p, time-4.0/CONVERSION);
            # hist_NA50_40 = historicFlux[NA50];
            # hist_NA52_40 = historicFlux[NA52];
            # hist_NC50_40 = historicFlux[NC50];
            # hist_NC52_40 = historicFlux[NC52];
            # delay(historicFlux, p, time-12.0/CONVERSION);
            # hist_NA50_120 = historicFlux[NA50];
            # hist_NA52_120 = historicFlux[NA52];
            # hist_NC50_120 = historicFlux[NC50];
            # hist_NC52_120 = historicFlux[NC52];
        end

        #--------------------------------------------
        # MODULE 1: IkB Reactions
        #--------------------------------------------
        # IkB mRNA
        # 0 --nNFkB--> tIkB : nuclear NFkB dimer-induced synthesis
        # const rate * ( 1 + sum( w*([d]/Kd)^hill)) / ( 1 + sum( ([d]/Kd)^hill))
        # tIkBa induction weights changed from (200, 200) to (100, 100)
        reactionFlux[TIKBA, BASALSYNTHESIS] = Srates[TIKBA, BASALSYNTHESIS] * hillInduction(200, 200, concentration[NA50], concentration[NA52]; Kd = 150.0, hill = 1.1);
        reactionFlux[TIKBB, BASALSYNTHESIS] = Srates[TIKBB, BASALSYNTHESIS];
        reactionFlux[TIKBE, BASALSYNTHESIS] = Srates[TIKBE, BASALSYNTHESIS] * hillInduction(5, 5, 25, 25, concentration[NA50], concentration[NA52], concentration[NC50], concentration[NC52]; Kd = 150.0, hill = 1.1);
        if phase == 1  # pre-simulation phase (no delay)
            for i in TIKBA:TIKBE
                reactionFlux[i, SYNTHESIS] = reactionFlux[i, BASALSYNTHESIS];
            end
        elseif phase == 2
            reactionFlux[TIKBB, SYNTHESIS] = reactionFlux[TIKBB, BASALSYNTHESIS]; # no induction for tIkbb
            if time > 0.25
                reactionFlux[TIKBA, SYNTHESIS] = Srates[TIKBA, BASALSYNTHESIS] * hillInduction(200, 200, hist_NA50_025, hist_NA52_025; Kd = 150.0, hill = 1.1);
            else # transcription is at basal level
                reactionFlux[TIKBA, SYNTHESIS] = Srates[TIKBA, BASALSYNTHESIS];
            end
            if time > 0.75
                reactionFlux[TIKBE, SYNTHESIS] = Srates[TIKBE, BASALSYNTHESIS] * hillInduction(5, 5, 25, 25, hist_NA50_075, hist_NA52_075, hist_NC50_075, hist_NC52_075; Kd = 150.0, hill = 1.1);
            else # transcription is at basal level
                reactionFlux[TIKBE, SYNTHESIS] = Srates[TIKBE, BASALSYNTHESIS];
            end
        end
        # tIkB ----> 0 : Basal decay
        for i in TIKBA:TIKBE
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        # IkB free in cytoplasm
        for i in IKBA:IKBD
            # tIkB ----> IkB : translation
            if i == IKBD  # IkBd is associated from p100, not translated from RNA
                reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[P100] ^ 2;
                reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            else
                reactionFlux[i, TRANSLATION] = EPS * Srates[i, TRANSLATION] * concentration[i-IKBA+TIKBA];
            end
            # IkB ----> 0 : Basal decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
            # IkB + IKK ----> IKK : IKK-mediated free IkB decay (NIK-mediated for IkBd)
            if i == IKBD
                reactionFlux[i, DECAY] = Srates[i, DECAY] * concentration[NIK] * concentration[i] / (NIK_KM * (1 + concentration[P100] / NIK_KM) + concentration[i]);
            else
                reactionFlux[i, DECAY] = concentration[IKK] * Srates[i, DECAY] * concentration[i];
            end
            # IkB ----> nIkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
        end

        # IkB free in nucleus
        for i in NIKBA:NIKBD
            if i == NIKBD  # nIkBd can be associated from np100
                reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[NP100] ^ 2;
                reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            end
            # nIkB ----> 0 : Basal decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
            # nIkB ----> IkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 2: NFkB monomer Reactions
        #--------------------------------------------
        # NFkB mRNA
        # 0 --nNFkB--> tNFkB : nuclear NFkB dimer-induced synthesis
        # const rate  +  sum( w*([d]/Kd)^hill) / ( 1 + sum( ([d]/Kd)^hill)) / 100000
        reactionFlux[TRELA, BASALSYNTHESIS] = Srates[TRELA, BASALSYNTHESIS];
        reactionFlux[TP50, BASALSYNTHESIS] = Srates[TP50, BASALSYNTHESIS] + hillInduction_NFkB(20, 20, 20, 20, concentration[NA50], concentration[NA52], concentration[NC50], concentration[NC52]; Kd = 150.0, hill = 1.0) / 100000;
        reactionFlux[TRELB, BASALSYNTHESIS] = Srates[TRELB, BASALSYNTHESIS] + hillInduction_NFkB(5, 5, 100, 100, concentration[NA50], concentration[NA52], concentration[NC50], concentration[NC52]; Kd = 150.0, hill = 1.0) / 100000;
        reactionFlux[TP100, BASALSYNTHESIS] = Srates[TP100, BASALSYNTHESIS] + hillInduction_p100(1000*P100_MOD, 1000*P100_MOD, 1500*P100_MOD, 1500*P100_MOD, concentration[NA50], concentration[NA52], concentration[NC50], concentration[NC52]; Kd = 150.0, hill = 1.0);
        reactionFlux[TCREL, BASALSYNTHESIS] = Srates[TCREL, BASALSYNTHESIS] + hillInduction_NFkB(20, 20, 20, 20, concentration[NA50], concentration[NA52], concentration[NC50], concentration[NC52]; Kd = 150.0, hill = 1.0) / 100000;
        if phase == 1  # pre-simulation phase (no delay)
            for i in TRELA:TCREL
                reactionFlux[i, SYNTHESIS] = reactionFlux[i, BASALSYNTHESIS];
            end
        elseif phase == 2  # simulation phase (with delay)
            reactionFlux[TRELA, SYNTHESIS] = reactionFlux[TRELA, BASALSYNTHESIS]; # no delay for tRelA
            if time > 1
                reactionFlux[TP50, SYNTHESIS] = Srates[TP50, BASALSYNTHESIS] + hillInduction_NFkB(20, 20, 20, 20, hist_NA50_1, hist_NA52_1, hist_NC50_1, hist_NC52_1; Kd = 150.0, hill = 1.0) / 100000;
                reactionFlux[TRELB, SYNTHESIS] = Srates[TRELB, BASALSYNTHESIS] + hillInduction_NFkB(5, 5, 100, 100, hist_NA50_1, hist_NA52_1, hist_NC50_1, hist_NC52_1; Kd = 150.0, hill = 1.0) / 100000;
            else # transcription is at basal level (with basal induction)
                reactionFlux[TP50, SYNTHESIS] = Srates[TP50, BASALSYNTHESIS];
                reactionFlux[TRELB, SYNTHESIS] = Srates[TRELB, BASALSYNTHESIS];
            end
            if time > 4
                reactionFlux[TP100, SYNTHESIS] = Srates[TP100, BASALSYNTHESIS] + hillInduction_p100(1000*P100_MOD, 1000*P100_MOD, 1500*P100_MOD, 1500*P100_MOD, hist_NA50_40, hist_NA52_40, hist_NC50_40, hist_NC52_40; Kd = 150.0, hill = 1.0);
            else # transcription is at basal level (with basal induction)
                reactionFlux[TP100, SYNTHESIS] = Srates[TP100, BASALSYNTHESIS];
            end
            if time > 12
                reactionFlux[TCREL, SYNTHESIS] = Srates[TCREL, BASALSYNTHESIS] + hillInduction_NFkB(20, 20, 20, 20, hist_NA50_120, hist_NA52_120, hist_NC50_120, hist_NC52_120; Kd = 150.0, hill = 1.0) / 100000;
            else # transcription is at basal level (with basal induction)
                reactionFlux[TCREL, SYNTHESIS] = Srates[TCREL, BASALSYNTHESIS];
            end
        end
        # Blimp1 ----| tcRel : Repression (multiplicative factor)
    #         reactionFlux[TCREL, REPRESSION] = hillRepression(Kr = 1.0, nr = 2.0, X = concentration[BLIMP1]);
        reactionFlux[TCREL, REPRESSION] = 1.0;
        reactionFlux[TCREL, SYNTHESIS] = reactionFlux[TCREL, REPRESSION] * reactionFlux[TCREL, SYNTHESIS];
        for i in TRELA:TCREL
            # tNFkB ----> 0 : Basal decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        # NFkB free in cytoplasm
        for i in RELA:P52
            # tNFkB ----> NFkB : translation
            if i == P52
                # p100 --NIK--> p52 : NIK-mediated p100 processing (S6 in Mitchell )
                reactionFlux[i, TRANSLATION] = Srates[i, TRANSLATION] * concentration[NIK] * concentration[P100] / (NIK_KM * (1 + concentration[IKBD] / NIK_KM) + concentration[P100]);
            else
                reactionFlux[i, TRANSLATION] = EPS * Srates[i, TRANSLATION] * concentration[i-(RELA-TRELA)];
            end
            # NFkB ----> 0 : Basal decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        # NFkB free in nucleus
        for i in NRELA:NP52
            # nNFkB ----> 0 : Basal decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 3: NFkB monomer <-> dimer Reactions
        #--------------------------------------------
        # NFkB + NFkB ----> NFkB:NFkB : association (cytoplasm)
        reactionFlux[AA, ASSOCIATION] = Srates[AA, ASSOCIATION] * concentration[RELA] ^ 2;
        reactionFlux[A50, ASSOCIATION] = Srates[A50, ASSOCIATION] * concentration[RELA] * concentration[P50];
        reactionFlux[A52, ASSOCIATION] = Srates[A52, ASSOCIATION] * concentration[RELA] * concentration[P52];
        reactionFlux[B50, ASSOCIATION] = Srates[B50, ASSOCIATION] * concentration[RELB] * concentration[P50];
        reactionFlux[B52, ASSOCIATION] = Srates[B52, ASSOCIATION] * concentration[RELB] * concentration[P52];
        reactionFlux[C50, ASSOCIATION] = Srates[C50, ASSOCIATION] * concentration[CREL] * concentration[P50];
        reactionFlux[C52, ASSOCIATION] = Srates[C52, ASSOCIATION] * concentration[CREL] * concentration[P52];
        reactionFlux[C100, ASSOCIATION] = Srates[C100, ASSOCIATION] * concentration[CREL] * concentration[P100];
        reactionFlux[P50P50, ASSOCIATION] = Srates[P50P50, ASSOCIATION] * concentration[P50] ^ 2;
        reactionFlux[P52P52, ASSOCIATION] = Srates[P52P52, ASSOCIATION] * concentration[P52] ^ 2;
        # nNFkB + nNFkB ----> nNFkB:NFkB : association (nucleus)
        reactionFlux[NAA, ASSOCIATION] = Srates[NAA, ASSOCIATION] * concentration[NRELA] ^ 2;
        reactionFlux[NA50, ASSOCIATION] = Srates[NA50, ASSOCIATION] * concentration[NRELA] * concentration[NP50];
        reactionFlux[NA52, ASSOCIATION] = Srates[NA52, ASSOCIATION] * concentration[NRELA] * concentration[NP52];
        reactionFlux[NB50, ASSOCIATION] = Srates[NB50, ASSOCIATION] * concentration[NRELB] * concentration[NP50];
        reactionFlux[NB52, ASSOCIATION] = Srates[NB52, ASSOCIATION] * concentration[NRELB] * concentration[NP52];
        reactionFlux[NC50, ASSOCIATION] = Srates[NC50, ASSOCIATION] * concentration[NCREL] * concentration[NP50];
        reactionFlux[NC52, ASSOCIATION] = Srates[NC52, ASSOCIATION] * concentration[NCREL] * concentration[NP52];
        reactionFlux[NC100, ASSOCIATION] = Srates[NC100, ASSOCIATION] * concentration[NCREL] * concentration[NP100];
        reactionFlux[NP50P50, ASSOCIATION] = Srates[NP50P50, ASSOCIATION] * concentration[NP50] ^ 2;
        reactionFlux[NP52P52, ASSOCIATION] = Srates[NP52P52, ASSOCIATION] * concentration[NP52] ^ 2;

        for i in AA:P52P52
            # NFkB:NFkB ----> NFkB + NFkB : dissociation (cytoplasm)
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # NFkB:NFkB ----> nNFkB:NFkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
            # NFkB:NFkB ----> 0 : Basal decay (cytoplasm)
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        for i in NAA:NP52P52
            # nNFkB:NFkB ----> nNFkB + nNFkB : dissociation (nucleus)
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # nNFkB:NFkB ----> NFkB:NFkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
            # nNFkB:NFkB ----> 0 : Basal decay (nucleus)
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 4: AA:IkB Reactions
        #--------------------------------------------
        # in cytoplasm
        for i in IKBAAA:IKBDAA
            # AA + IkB ----> AA:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[AA] * concentration[i-(IKBAAA-IKBA)];
            # AA:IkB ----> AA + IkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # AA:IkB ----> nAA:IkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
            # AA:IkB + IKK ----> AA + IKK : IKK-mediated NFkB-bound IkB decay (NIK-mediated for IkBd)
            reactionFlux[i, DECAY] = concentration[IKK] * Srates[i-(IKBAAA-IKBA), DECAY] * concentration[i];
            # AA:IkB ----> AA : IKK-independent NFkB-bound IkB decay
            reactionFlux[i, DECAY] += Srates[i, DECAY] * concentration[i];
            # AA:IkB ----> IkB : IkB-bound NFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end
        # AA:IkBd + NIK ----> AA + NIK : NIK-mediated NFkB-bound IkBd decay
        reactionFlux[IKBDAA, DECAY] = concentration[NIK] * Srates[IKBD, DECAY] * concentration[IKBDAA] / (NIK_KM * (1 + concentration[P100] / NIK_KM) + concentration[IKBDAA]) +
            Srates[IKBDAA, DECAY] * concentration[IKBDAA];

        # in nucleus
        for i in NIKBAAA:NIKBDAA
            # nAA + nIkB ----> nAA:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[NAA] * concentration[i-(NIKBAAA-NIKBA)];
            # nAA:IkB ----> nAA + nIkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # nAA:IkB ----> AA:IkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
            # nAA:IkB ----> nAA : IKK-independent nNFkB-bound nIkB decay
            reactionFlux[i, DECAY] = Srates[i, DECAY] * concentration[i];
            # nAA:IkB ----> nIkB : nIkB-bound nNFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 5: A50:IkB Reactions
        #--------------------------------------------
        # in cytoplasm
        for i in IKBAA50:IKBDA50
            # A50 + IkB ----> A50:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[A50] * concentration[i-(IKBAA50-IKBA)];
            # A50:IkB ----> A50 + IkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # A50:IkB ----> nA50:IkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
            # A50:IkB + IKK ----> A50 + IKK : IKK-mediated NFkB-bound IkB decay (NIK-mediated for IkBd)
            reactionFlux[i, DECAY] = concentration[IKK] * Srates[i-(IKBAA50-IKBA), DECAY] * concentration[i];
            # A50:IkB ----> A50 : IKK-independent NFkB-bound IkB decay
            reactionFlux[i, DECAY] += Srates[i, DECAY] * concentration[i];
            # A50:IkB ----> IkB : IkB-bound NFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end
        # A50:IkBd + NIK ----> A50 + NIK : NIK-mediated NFkB-bound IkBd decay
        reactionFlux[IKBDA50, DECAY] = concentration[NIK] * Srates[IKBD, DECAY] * concentration[IKBDA50] / (NIK_KM * (1 + concentration[P100] / NIK_KM) + concentration[IKBDA50]) +
            Srates[IKBDA50, DECAY] * concentration[IKBDA50];

        # in nucleus
        for i in NIKBAA50:NIKBDA50
            # nA50 + nIkB ----> nA50:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[NA50] * concentration[i-(NIKBAA50-NIKBA)];
            # nA50:IkB ----> nA50 + nIkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # nA50:IkB ----> A50:IkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
            # nA50:IkB ----> nA50 : IKK-independent nNFkB-bound nIkB decay
            reactionFlux[i, DECAY] = Srates[i, DECAY] * concentration[i];
            # nA50:IkB ----> nIkB : nIkB-bound nNFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 6: A52:IkB Reactions
        #--------------------------------------------
        # in cytoplasm
        for i in IKBAA52:IKBDA52
            # A52 + IkB ----> A52:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[A52] * concentration[i-(IKBAA52-IKBA)];
            # A52:IkB ----> A52 + IkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # A52:IkB ----> nA52:IkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
            # A52:IkB + IKK ----> A52 + IKK : IKK-mediated NFkB-bound IkB decay (NIK-mediated for IkBd)
            reactionFlux[i, DECAY] = concentration[IKK] * Srates[i-(IKBAA52-IKBA), DECAY] * concentration[i];
            # A52:IkB ----> A52 : IKK-independent NFkB-bound IkB decay
            reactionFlux[i, DECAY] += Srates[i, DECAY] * concentration[i];
            # A52:IkB ----> IkB : IkB-bound NFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end
        # A52:IkBd + NIK ----> A52 + NIK : NIK-mediated NFkB-bound IkBd decay
        reactionFlux[IKBDA52, DECAY] = concentration[NIK] * Srates[IKBD, DECAY] * concentration[IKBDA52] / (NIK_KM * (1 + concentration[P100] / NIK_KM) + concentration[IKBDA52]) +
            Srates[IKBDA52, DECAY] * concentration[IKBDA52];

        # in nucleus
        for i in NIKBAA52:NIKBDA52
            # nA52 + nIkB ----> nA52:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[NA52] * concentration[i-(NIKBAA52-NIKBA)];
            # nA52:IkB ----> nA52 + nIkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # nA52:IkB ----> A52:IkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
            # nA52:IkB ----> nA52 : IKK-independent nNFkB-bound nIkB decay
            reactionFlux[i, DECAY] = Srates[i, DECAY] * concentration[i];
            # nA52:IkB ----> nIkB : nIkB-bound nNFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 7: B50:IkB Reactions
        #--------------------------------------------
        # in cytoplasm
        for i in IKBAB50:IKBDB50
            # B50 + IkB ----> B50:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[B50] * concentration[i-(IKBAB50-IKBA)];
            # B50:IkB ----> B50 + IkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # B50:IkB ----> nB50:IkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
            # B50:IkB + IKK ----> B50 + IKK : IKK-mediated NFkB-bound IkB decay (NIK-mediated for IkBd)
            reactionFlux[i, DECAY] = concentration[IKK] * Srates[i-(IKBAB50-IKBA), DECAY] * concentration[i];
            # B50:IkB ----> B50 : IKK-independent NFkB-bound IkB decay
            reactionFlux[i, DECAY] += Srates[i, DECAY] * concentration[i];
            # B50:IkB ----> IkB : IkB-bound NFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end
        # B50:IkBd + NIK ----> B50 + NIK : NIK-mediated NFkB-bound IkBd decay
        reactionFlux[IKBDB50, DECAY] = concentration[NIK] * Srates[IKBD, DECAY] * concentration[IKBDB50] / (NIK_KM * (1 + concentration[P100] / NIK_KM) + concentration[IKBDB50]) +
            Srates[IKBDB50, DECAY] * concentration[IKBDB50];

        # in nucleus
        for i in NIKBAB50:NIKBDB50
            # nB50 + nIkB ----> nB50:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[NB50] * concentration[i-(NIKBAB50-NIKBA)];
            # nB50:IkB ----> nB50 + nIkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # nB50:IkB ----> B50:IkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
            # nB50:IkB ----> nB50 : IKK-independent nNFkB-bound nIkB decay
            reactionFlux[i, DECAY] = Srates[i, DECAY] * concentration[i];
            # nB50:IkB ----> nIkB : nIkB-bound nNFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 8: B52:IkBd Reactions
        #--------------------------------------------
        # in cytoplasm
        # B52 + IkB ----> B52:IkB : association
        reactionFlux[IKBDB52, ASSOCIATION] = Srates[IKBDB52, ASSOCIATION] * concentration[B52] * concentration[IKBD];
        # B52:IkB ----> B52 + IkB : dissociation
        reactionFlux[IKBDB52, DISSOCIATION] = Srates[IKBDB52, DISSOCIATION] * concentration[IKBDB52];
        # B52:IkB ----> nB52:IkB : transport into nucleus
        reactionFlux[IKBDB52, TRANSPORTIN] = Srates[IKBDB52, TRANSPORTIN] * concentration[IKBDB52];
        # B52:IkB + NIK ----> B52 + NIK : NIK-mediated NFkB-bound IkBd decay
        # B52:IkB ----> B52 : NIK-independent NFkB-bound IkB decay
        reactionFlux[IKBDB52, DECAY] = concentration[NIK] * Srates[IKBD, DECAY] * concentration[IKBDB52] / (NIK_KM * (1 + concentration[P100] / NIK_KM) + concentration[IKBDB52]) +
            Srates[IKBDB52, DECAY] * concentration[IKBDB52];
        # B52:IkB ----> IkB : IkB-bound NFkB decay
        reactionFlux[IKBDB52, BASALDECAY] = Srates[IKBDB52, BASALDECAY] * concentration[IKBDB52];

        # in nucleus
        # nB52 + nIkB ----> nB52:IkB : association
        reactionFlux[NIKBDB52, ASSOCIATION] = Srates[NIKBDB52, ASSOCIATION] * concentration[NB52] * concentration[NIKBD];
        # nB52:IkB ----> nB52 + nIkB : dissociation
        reactionFlux[NIKBDB52, DISSOCIATION] = Srates[NIKBDB52, DISSOCIATION] * concentration[NIKBDB52];
        # nB52:IkB ----> B52:IkB : transport out of nucleus
        reactionFlux[NIKBDB52, TRANSPORTOUT] = Srates[NIKBDB52, TRANSPORTOUT] * concentration[NIKBDB52];
        # nB52:IkB ----> nB52 : NIK-independent nNFkB-bound nIkB decay
        reactionFlux[NIKBDB52, DECAY] = Srates[NIKBDB52, DECAY] * concentration[NIKBDB52];
        # nB52:IkB ----> nIkB : nIkB-bound nNFkB decay
        reactionFlux[NIKBDB52, BASALDECAY] = Srates[NIKBDB52, BASALDECAY] * concentration[NIKBDB52];


        #--------------------------------------------
        # MODULE 9: C50:IkB Reactions
        #--------------------------------------------
        # in cytoplasm
        for i in IKBAC50:IKBDC50
            # C50 + IkB ----> C50:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[C50] * concentration[i-(IKBAC50-IKBA)];
            # C50:IkB ----> C50 + IkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # C50:IkB ----> nC50:IkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
            # C50:IkB + IKK ----> C50 + IKK : IKK-mediated NFkB-bound IkB decay (NIK-mediated for IkBd)
            reactionFlux[i, DECAY] = concentration[IKK] * Srates[i-(IKBAC50-IKBA), DECAY] * concentration[i];
            # C50:IkB ----> C50 : IKK-independent NFkB-bound IkB decay
            reactionFlux[i, DECAY] += Srates[i, DECAY] * concentration[i];
            # C50:IkB ----> IkB : IkB-bound NFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end
        # C50:IkBd + NIK ----> C50 + NIK : NIK-mediated NFkB-bound IkBd decay
        reactionFlux[IKBDC50, DECAY] = concentration[NIK] * Srates[IKBD, DECAY] * concentration[IKBDC50] / (NIK_KM * (1 + concentration[P100] / NIK_KM) + concentration[IKBDC50]) +
            Srates[IKBDC50, DECAY] * concentration[IKBDC50];

        # in nucleus
        for i in NIKBAC50:NIKBDC50
            # nC50 + nIkB ----> nC50:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[NC50] * concentration[i-(NIKBAC50-NIKBA)];
            # nC50:IkB ----> nC50 + nIkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # nC50:IkB ----> C50:IkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
            # nC50:IkB ----> nC50 : IKK-independent nNFkB-bound nIkB decay
            reactionFlux[i, DECAY] = Srates[i, DECAY] * concentration[i];
            # nC50:IkB ----> nIkB : nIkB-bound nNFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 10: C52:IkB Reactions
        #--------------------------------------------
        # in cytoplasm
        for i in IKBAC52:IKBDC52
            # C52 + IkB ----> C52:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[C52] * concentration[i-(IKBAC52-IKBA)];
            # C52:IkB ----> C52 + IkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # C52:IkB ----> nC52:IkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
            # C52:IkB + IKK ----> C52 + IKK : IKK-mediated NFkB-bound IkB decay (NIK-mediated for IkBd)
            reactionFlux[i, DECAY] = concentration[IKK] * Srates[i-(IKBAC52-IKBA), DECAY] * concentration[i];
            # C52:IkB ----> C52 : IKK-independent NFkB-bound IkB decay
            reactionFlux[i, DECAY] += Srates[i, DECAY] * concentration[i];
            # C52:IkB ----> IkB : IkB-bound NFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end
        # C52:IkBd + NIK ----> C52 + NIK : NIK-mediated NFkB-bound IkBd decay
        reactionFlux[IKBDC52, DECAY] = concentration[NIK] * Srates[IKBD, DECAY] * concentration[IKBDC52] / (NIK_KM * (1 + concentration[P100] / NIK_KM) + concentration[IKBDC52]) +
            Srates[IKBDC52, DECAY] * concentration[IKBDC52];

        # in nucleus
        for i in NIKBAC52:NIKBDC52
            # nC52 + nIkB ----> nC52:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[NC52] * concentration[i-(NIKBAC52-NIKBA)];
            # nC52:IkB ----> nC52 + nIkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # nC52:IkB ----> C52:IkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
            # nC52:IkB ----> nC52 : IKK-independent nNFkB-bound nIkB decay
            reactionFlux[i, DECAY] = Srates[i, DECAY] * concentration[i];
            # nC52:IkB ----> nIkB : nIkB-bound nNFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end
    end
    nothing
end

function NFkBNettFluxes!(nettFlux, reactionFlux)
    @inbounds begin
        #------------------------------------------------
        # MODULE 1
        #--------------------------------------------
        # 2-5 : IkB mRNA
        for i in TIKBA:TIKBE
            nettFlux[i] = reactionFlux[i, SYNTHESIS] - reactionFlux[i, BASALDECAY];
        end
        # 6-9 : IkB free in cytoplasm
        for i in IKBA:IKBD
            if i == IKBD
                nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION];
            else
                nettFlux[i] = reactionFlux[i, TRANSLATION];
            end
            nettFlux[i] += - reactionFlux[i, BASALDECAY] - reactionFlux[i, DECAY] -
            reactionFlux[i, TRANSPORTIN] + reactionFlux[i-IKBA+NIKBA, TRANSPORTOUT] -
            reactionFlux[i-IKBA+IKBAAA, ASSOCIATION] + reactionFlux[i-IKBA+IKBAAA, DISSOCIATION] + reactionFlux[i-IKBA+IKBAAA, BASALDECAY] -
            reactionFlux[i-IKBA+IKBAA50, ASSOCIATION] + reactionFlux[i-IKBA+IKBAA50, DISSOCIATION] + reactionFlux[i-IKBA+IKBAA50, BASALDECAY] -
            reactionFlux[i-IKBA+IKBAA52, ASSOCIATION] + reactionFlux[i-IKBA+IKBAA52, DISSOCIATION] + reactionFlux[i-IKBA+IKBAA52, BASALDECAY] -
            reactionFlux[i-IKBA+IKBAB50, ASSOCIATION] + reactionFlux[i-IKBA+IKBAB50, DISSOCIATION] + reactionFlux[i-IKBA+IKBAB50, BASALDECAY] -
            reactionFlux[i-IKBA+IKBAB52, ASSOCIATION] + reactionFlux[i-IKBA+IKBAB52, DISSOCIATION] + reactionFlux[i-IKBA+IKBAB52, BASALDECAY] -
            reactionFlux[i-IKBA+IKBAC50, ASSOCIATION] + reactionFlux[i-IKBA+IKBAC50, DISSOCIATION] + reactionFlux[i-IKBA+IKBAC50, BASALDECAY] -
            reactionFlux[i-IKBA+IKBAC52, ASSOCIATION] + reactionFlux[i-IKBA+IKBAC52, DISSOCIATION] + reactionFlux[i-IKBA+IKBAC52, BASALDECAY];
        end

        # 10-13 : IkB free in nucleus
        for i in NIKBA:NIKBD
            nettFlux[i] = - reactionFlux[i, BASALDECAY] +
            reactionFlux[i-NIKBA+IKBA, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] -
            reactionFlux[i-NIKBA+NIKBAAA, ASSOCIATION] + reactionFlux[i-NIKBA+NIKBAAA, DISSOCIATION] + reactionFlux[i-NIKBA+NIKBAAA, BASALDECAY] -
            reactionFlux[i-NIKBA+NIKBAA50, ASSOCIATION] + reactionFlux[i-NIKBA+NIKBAA50, DISSOCIATION] + reactionFlux[i-NIKBA+NIKBAA50, BASALDECAY] -
            reactionFlux[i-NIKBA+NIKBAA52, ASSOCIATION] + reactionFlux[i-NIKBA+NIKBAA52, DISSOCIATION] + reactionFlux[i-NIKBA+NIKBAA52, BASALDECAY] -
            reactionFlux[i-NIKBA+NIKBAB50, ASSOCIATION] + reactionFlux[i-NIKBA+NIKBAB50, DISSOCIATION] + reactionFlux[i-NIKBA+NIKBAB50, BASALDECAY] -
            reactionFlux[i-NIKBA+NIKBAB52, ASSOCIATION] + reactionFlux[i-NIKBA+NIKBAB52, DISSOCIATION] + reactionFlux[i-NIKBA+NIKBAB52, BASALDECAY] -
            reactionFlux[i-NIKBA+NIKBAC50, ASSOCIATION] + reactionFlux[i-NIKBA+NIKBAC50, DISSOCIATION] + reactionFlux[i-NIKBA+NIKBAC50, BASALDECAY] -
            reactionFlux[i-NIKBA+NIKBAC52, ASSOCIATION] + reactionFlux[i-NIKBA+NIKBAC52, DISSOCIATION] + reactionFlux[i-NIKBA+NIKBAC52, BASALDECAY];
            if i == NIKBD
                nettFlux[i] += reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION];
            end
        end

        #------------------------------------------------
        # MODULE 2
        #--------------------------------------------
        # NFkB mRNA
        for i in TRELA:TCREL
            nettFlux[i] = reactionFlux[i, SYNTHESIS] - reactionFlux[i, BASALDECAY];
        end
        # NFkB free in cytoplasm
        for i in RELA:P52
            nettFlux[i] = reactionFlux[i, TRANSLATION] - reactionFlux[i, BASALDECAY];
        end
        nettFlux[P100] += - reactionFlux[P52, TRANSLATION] - 2 * reactionFlux[IKBD, ASSOCIATION] + 2 * reactionFlux[IKBD, DISSOCIATION]; # NIK-mediated p100 processing

        nettFlux[RELA] += - 2 * reactionFlux[AA, ASSOCIATION] + 2 * reactionFlux[AA, DISSOCIATION];
        nettFlux[RELA] += - reactionFlux[A50, ASSOCIATION] + reactionFlux[A50, DISSOCIATION];
        nettFlux[RELA] += - reactionFlux[A52, ASSOCIATION] + reactionFlux[A52, DISSOCIATION];
        nettFlux[RELB] += - reactionFlux[B50, ASSOCIATION] + reactionFlux[B50, DISSOCIATION];
        nettFlux[RELB] += - reactionFlux[B52, ASSOCIATION] + reactionFlux[B52, DISSOCIATION];
        nettFlux[P50] += - 2 * reactionFlux[P50P50, ASSOCIATION] + 2 * reactionFlux[P50P50, DISSOCIATION];
        nettFlux[P50] += - reactionFlux[A50, ASSOCIATION] + reactionFlux[A50, DISSOCIATION];
        nettFlux[P50] += - reactionFlux[B50, ASSOCIATION] + reactionFlux[B50, DISSOCIATION];
        nettFlux[P50] += - reactionFlux[C50, ASSOCIATION] + reactionFlux[C50, DISSOCIATION];
        nettFlux[P100] += - reactionFlux[C100, ASSOCIATION] + reactionFlux[C100, DISSOCIATION];
        nettFlux[P52] += - 2 * reactionFlux[P52P52, ASSOCIATION] + 2 * reactionFlux[P52P52, DISSOCIATION];
        nettFlux[P52] += - reactionFlux[A52, ASSOCIATION] + reactionFlux[A52, DISSOCIATION];
        nettFlux[P52] += - reactionFlux[B52, ASSOCIATION] + reactionFlux[B52, DISSOCIATION];
        nettFlux[P52] += - reactionFlux[C52, ASSOCIATION] + reactionFlux[C52, DISSOCIATION];
        nettFlux[CREL] += - reactionFlux[C50, ASSOCIATION] + reactionFlux[C50, DISSOCIATION];
        nettFlux[CREL] += - reactionFlux[C52, ASSOCIATION] + reactionFlux[C52, DISSOCIATION];
        nettFlux[CREL] += - reactionFlux[C100, ASSOCIATION] + reactionFlux[C100, DISSOCIATION];
        # NFkB free in nucleus
        for i in NRELA:NP52
            nettFlux[i] = - reactionFlux[i, BASALDECAY];
        end
        nettFlux[NP100] += - 2 * reactionFlux[NIKBD, ASSOCIATION] + 2 * reactionFlux[NIKBD, DISSOCIATION]; # NIK-mediated p100 processing

        nettFlux[NRELA] += - 2 * reactionFlux[NAA, ASSOCIATION] + 2 * reactionFlux[NAA, DISSOCIATION];
        nettFlux[NRELA] += - reactionFlux[NA50, ASSOCIATION] + reactionFlux[NA50, DISSOCIATION];
        nettFlux[NRELA] += - reactionFlux[NA52, ASSOCIATION] + reactionFlux[NA52, DISSOCIATION];
        nettFlux[NRELB] += - reactionFlux[NB50, ASSOCIATION] + reactionFlux[NB50, DISSOCIATION];
        nettFlux[NRELB] += - reactionFlux[NB52, ASSOCIATION] + reactionFlux[NB52, DISSOCIATION];
        nettFlux[NP50] += - 2 * reactionFlux[NP50P50, ASSOCIATION] + 2 * reactionFlux[NP50P50, DISSOCIATION];
        nettFlux[NP50] += - reactionFlux[NA50, ASSOCIATION] + reactionFlux[NA50, DISSOCIATION];
        nettFlux[NP50] += - reactionFlux[NB50, ASSOCIATION] + reactionFlux[NB50, DISSOCIATION];
        nettFlux[NP50] += - reactionFlux[NC50, ASSOCIATION] + reactionFlux[NC50, DISSOCIATION];
        nettFlux[NP100] += - reactionFlux[NC100, ASSOCIATION] + reactionFlux[NC100, DISSOCIATION];
        nettFlux[NP52] += - 2 * reactionFlux[NP52P52, ASSOCIATION] + 2 * reactionFlux[NP52P52, DISSOCIATION];
        nettFlux[NP52] += - reactionFlux[NA52, ASSOCIATION] + reactionFlux[NA52, DISSOCIATION];
        nettFlux[NP52] += - reactionFlux[NB52, ASSOCIATION] + reactionFlux[NB52, DISSOCIATION];
        nettFlux[NP52] += - reactionFlux[NC52, ASSOCIATION] + reactionFlux[NC52, DISSOCIATION];
        nettFlux[NCREL] += - reactionFlux[NC50, ASSOCIATION] + reactionFlux[NC50, DISSOCIATION];
        nettFlux[NCREL] += - reactionFlux[NC52, ASSOCIATION] + reactionFlux[NC52, DISSOCIATION];
        nettFlux[NCREL] += - reactionFlux[NC100, ASSOCIATION] + reactionFlux[NC100, DISSOCIATION];

        #------------------------------------------------
        # MODULE 3
        #--------------------------------------------
        # NFkB free in cytoplasm
        for i in AA:P52P52
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] - reactionFlux[i, TRANSPORTIN] + reactionFlux[i-AA+NAA, TRANSPORTOUT] - reactionFlux[i, BASALDECAY];
        end
        # NFkB free in nucleus
        for i in NAA:NP52P52
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] + reactionFlux[i-NAA+AA, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] - reactionFlux[i, BASALDECAY];
        end

        #------------------------------------------------
        # MODULE 4
        #--------------------------------------------
        for i in IKBAAA:IKBDAA
            # AA free in cytoplasm
            nettFlux[AA] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # AA:IkB in cytoplasm
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] - reactionFlux[i, TRANSPORTIN] + reactionFlux[i-IKBAAA+NIKBAAA, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        for i in NIKBAAA:NIKBDAA
            # AA free in nucleus
            nettFlux[NAA] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # AA:IkB in nucleus
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] + reactionFlux[i-NIKBAAA+IKBAAA, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        #------------------------------------------------
        # MODULE 5
        #--------------------------------------------
        for i in IKBAA50:IKBDA50
            # A50 free in cytoplasm
            nettFlux[A50] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # A50:IkB in cytoplasm
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] - reactionFlux[i, TRANSPORTIN] + reactionFlux[i-IKBAA50+NIKBAA50, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        for i in NIKBAA50:NIKBDA50
            # A50 free in nucleus
            nettFlux[NA50] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # A50:IkB in nucleus
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] + reactionFlux[i-NIKBAA50+IKBAA50, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        #------------------------------------------------
        # MODULE 6
        #--------------------------------------------
        for i in IKBAA52:IKBDA52
            # A52 free in cytoplasm
            nettFlux[A52] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # A52:IkB in cytoplasm
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] - reactionFlux[i, TRANSPORTIN] + reactionFlux[i-IKBAA52+NIKBAA52, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        for i in NIKBAA52:NIKBDA52
            # A52 free in nucleus
            nettFlux[NA52] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # A52:IkB in nucleus
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] + reactionFlux[i-NIKBAA52+IKBAA52, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        #------------------------------------------------
        # MODULE 7
        #--------------------------------------------
        for i in IKBAB50:IKBDB50
            # B50 free in cytoplasm
            nettFlux[B50] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # B50:IkB in cytoplasm
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] - reactionFlux[i, TRANSPORTIN] + reactionFlux[i-IKBAB50+NIKBAB50, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        for i in NIKBAB50:NIKBDB50
            # B50 free in nucleus
            nettFlux[NB50] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # B50:IkB in nucleus
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] + reactionFlux[i-NIKBAB50+IKBAB50, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        #------------------------------------------------
        # MODULE 8
        #--------------------------------------------
        for i in IKBAB52:IKBDB52
            # B52 free in cytoplasm
            nettFlux[B52] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # B52:IkB in cytoplasm
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] - reactionFlux[i, TRANSPORTIN] + reactionFlux[i-IKBAB52+NIKBAB52, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        for i in NIKBAB52:NIKBDB52
            # B52 free in nucleus
            nettFlux[NB52] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # B52:IkB in nucleus
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] + reactionFlux[i-NIKBAB52+IKBAB52, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        #------------------------------------------------
        # MODULE 9
        #--------------------------------------------
        for i in IKBAC50:IKBDC50
            # C50 free in cytoplasm
            nettFlux[C50] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # C50:IkB in cytoplasm
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] - reactionFlux[i, TRANSPORTIN] + reactionFlux[i-IKBAC50+NIKBAC50, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        for i in NIKBAC50:NIKBDC50
            # C50 free in nucleus
            nettFlux[NC50] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # C50:IkB in nucleus
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] + reactionFlux[i-NIKBAC50+IKBAC50, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        #------------------------------------------------
        # MODULE 10
        #--------------------------------------------
        for i in IKBAC52:IKBDC52
            # C52 free in cytoplasm
            nettFlux[C52] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # C52:IkB in cytoplasm
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] - reactionFlux[i, TRANSPORTIN] + reactionFlux[i-IKBAC52+NIKBAC52, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        for i in NIKBAC52:NIKBDC52
            # C52 free in nucleus
            nettFlux[NC52] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # C52:IkB in nucleus
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] + reactionFlux[i-NIKBAC52+IKBAC52, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

    end
    nothing
end

# time-independent ODE (pre-simulation: phase = 1)
function computeNFkBNettFluxes!(nettFlux, concentration, reactionFlux, Srates, phase)
    #------------------------------------------------
    # Compute NFkB reaction fluxes
    computeNFkBFluxes!(concentration, reactionFlux, Srates, phase, nothing, nothing, nothing);
    # Compute NFkB net fluxes
    NFkBNettFluxes!(nettFlux, reactionFlux);
    nothing
end

# time-dependent ODE (simulation, w/ delay: phase = 2)
function computeNFkBNettFluxes!(nettFlux, concentration, reactionFlux, Srates, phase, delay, historicFlux, time)
    #------------------------------------------------
    # Compute NFkB reaction fluxes
    computeNFkBFluxes!(concentration, reactionFlux, Srates, phase, delay, historicFlux, time);
    # Compute NFkB net fluxes
    NFkBNettFluxes!(nettFlux, reactionFlux);
    nothing
end
