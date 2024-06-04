# 2nd edition of Apoptosis model (added utilities to solve ODE separately from NFkB model)
# eq. 17+1 is BCR-mediated C8 processing, also change the inputCurves assignment
# eq. A03+1 is BCR-mediated tBcl2 decay
# eq. 41+1 is BCR-mediated Bcl2 decay
# include("HelperFunctions.jl");
# Compute reaction fluxes for apoptosis module
#--------------------------------------------

function computeApoptosisFluxes!(concentration, reactionFlux, Srates, phase)
    @inbounds begin
        #--------------------------------------------
        # MODULE 1: Bcl2 mRNA transcript
        #--------------------------------------------
        # A01. 0 --nNFkB--> tBcl2 activity : Activation (((kA50*[A50n]+kC50*[C50n]+kIKK*[IKK])/Kd)^expn/(k1+((kA50*[A50n]+kC50*[C50n]+kIKK*[IKK])/Kd)^expn)*(1-k2)+k2 ----> IKK activity represents other transcription factor
        # reactionFlux[TBCL2, ACTIVATION] = hillActivation(Ka = 1.0, na = 2.0, X = (0.3*concentration[NA50] + 0.6*concentration[NC50] + 0.1*concentration[IKK])/BCL2THR) * (1-Srates[TBCL2, ACTIVATION]) + Srates[TBCL2, ACTIVATION];
        reactionFlux[TBCL2, ACTIVATION] = hillActivation(Ka = 1.0, na = 2.0, X = (0.3*(concentration[NA50]+concentration[NA52]) + 0.6*(concentration[NC50]+concentration[NC52]) + 0.1*concentration[IKK])/BCL2THR) * (1-Srates[TBCL2, ACTIVATION]) + Srates[TBCL2, ACTIVATION];
        # A02. tBcl2 Act ----> tBcl2 : Synthesis
        reactionFlux[TBCL2, SYNTHESIS] = Srates[TBCL2, SYNTHESIS] * reactionFlux[TBCL2, ACTIVATION];
        # A03. tBcl2 ----> 0 : Basal decay
        reactionFlux[TBCL2, BASALDECAY] = Srates[TBCL2, BASALDECAY] * concentration[TBCL2];
        # A03+1. tBcl2 + ABCR ----> ABCR : Basal decay
        reactionFlux[TBCL2, DECAY] = Srates[TBCL2, DECAY] * concentration[TBCL2] * concentration[ABCR];

        #--------------------------------------------
        # MODULE 2: Receptor & DISC system
        #--------------------------------------------
        # 01. 0 ----> L : Basal synthesis
        reactionFlux[L, BASALSYNTHESIS] = EPS * Srates[L, BASALSYNTHESIS];
        # 02. L ----> 0 : Basal decay
        reactionFlux[L, BASALDECAY] = Srates[L, BASALDECAY] * concentration[L];
        # 03. 0 ----> R : Basal synthesis
        reactionFlux[R, BASALSYNTHESIS] = EPS * Srates[R, BASALSYNTHESIS];
        # 04. R ----> 0 : Basal decay
        reactionFlux[R, BASALDECAY] = Srates[R, BASALDECAY] * concentration[R];
        # 05. L + R ----> L:R : complex formation
        reactionFlux[LR, ASSOCIATION] = Srates[LR, ASSOCIATION] * concentration[L] * concentration[R];
        # 06. L:R ----> L + R : complex dissociation
        reactionFlux[LR, DISSOCIATION] = Srates[LR, DISSOCIATION] * concentration[LR];
        # 07. L:R ----> DISC : Activation
        reactionFlux[DISC, ACTIVATION] = Srates[DISC, ACTIVATION] * concentration[LR];

        # 08. 0 ----> flip : Basal synthesis
        reactionFlux[FLIP, BASALSYNTHESIS] = EPS * Srates[FLIP, BASALSYNTHESIS];
        # 09. flip ----> 0 : Basal decay
        reactionFlux[FLIP, BASALDECAY] = Srates[FLIP, BASALDECAY] * concentration[FLIP];
        # 10. flip + DISC ----> flip:DISC : complex formation
        reactionFlux[FLIPDISC, ASSOCIATION] = Srates[FLIPDISC, ASSOCIATION] * concentration[FLIP] * concentration[DISC];
        # 11. flip:DISC ----> flip + DISC : complex dissociation
        reactionFlux[FLIPDISC, DISSOCIATION] = Srates[FLIPDISC, DISSOCIATION] * concentration[FLIPDISC];
        # 12. flip:DISC ----> 0 : Basal decay
        reactionFlux[FLIPDISC, BASALDECAY] = Srates[FLIPDISC, BASALDECAY] * concentration[FLIPDISC];

        #--------------------------------------------
        # MODULE 3: CASPASE 8 module (initiator caspase)
        #--------------------------------------------
        # 13. 0 ----> pC8 : Basal synthesis
        reactionFlux[PC8, BASALSYNTHESIS] = EPS * Srates[PC8, BASALSYNTHESIS];
        # 14. pC8 ----> 0 : Basal decay
        reactionFlux[PC8, BASALDECAY] = Srates[PC8, BASALDECAY] * concentration[PC8];
        # 15. pC8 + DISC ----> DISC:pC8 : complex formation
        reactionFlux[DISCPC8, ASSOCIATION] = Srates[DISCPC8, ASSOCIATION] * concentration[PC8] * concentration[DISC];
        # 16. DISC:pC8 ----> pC8 + DISC : complex dissociation
        reactionFlux[DISCPC8, DISSOCIATION] = Srates[DISCPC8, DISSOCIATION] * concentration[DISCPC8];
        # 17. DISC:pC8 ----> C8 + DISC : activation
        reactionFlux[C8, ACTIVATION] = Srates[C8, ACTIVATION] * concentration[DISCPC8];
        # 17+1. pC8 + BCR ----> C8 + BCR : BCR-mediated pC8 processing
        reactionFlux[C8, BASALSYNTHESIS] = Srates[C8, BASALSYNTHESIS] * concentration[PC8] * concentration[ABCR];

        # 18. 0 ----> Bar : Basal synthesis
        reactionFlux[BAR, BASALSYNTHESIS] = EPS * Srates[BAR, BASALSYNTHESIS];
        # 19. Bar ----> 0 : Basal decay
        reactionFlux[BAR, BASALDECAY] = Srates[BAR, BASALDECAY] * concentration[BAR];
        # 20. Bar + C8 ----> Bar:C8 : complex formation
        reactionFlux[BARC8, ASSOCIATION] = Srates[BARC8, ASSOCIATION] * concentration[C8] * concentration[BAR];
        # 21. Bar:C8 ----> Bar + C8 : complex dissociation
        reactionFlux[BARC8, DISSOCIATION] = Srates[BARC8, DISSOCIATION] * concentration[BARC8];
        # 22. Bar:C8 ----> 0 : Basal decay
        reactionFlux[BARC8, BASALDECAY] = Srates[BARC8, BASALDECAY] * concentration[BARC8];

        #------------------------------------------------
        # MODULE 4: Before Bax enter mitochondria
        #------------------------------------------------
        # 23. 0 ----> Bid : Basal synthesis
        reactionFlux[BID, BASALSYNTHESIS] = EPS * Srates[BID, BASALSYNTHESIS];
        # 24. Bid ----> 0 : Basal decay
        reactionFlux[BID, BASALDECAY] = Srates[BID, BASALDECAY] * concentration[BID];
        # 25. C8 + Bid ----> C8:Bid : complex formation
        reactionFlux[C8BID, ASSOCIATION] = Srates[C8BID, ASSOCIATION] * concentration[C8] * concentration[BID];
        # 26. C8:Bid ----> C8 + Bid : complex dissociation
        reactionFlux[C8BID, DISSOCIATION] = Srates[C8BID, DISSOCIATION] * concentration[C8BID];
        # 27. C8:Bid ----> C8 + tBid : activation
        reactionFlux[TBID, ACTIVATION] = Srates[TBID, ACTIVATION] * concentration[C8BID];

        # 28. 0 ----> cBcl2 : Basal synthesis
        reactionFlux[CBCL2, BASALSYNTHESIS] = EPS * Srates[CBCL2, BASALSYNTHESIS];
        # 29. cBcl2 ----> 0 : Basal decay
        reactionFlux[CBCL2, BASALDECAY] = Srates[CBCL2, BASALDECAY] * concentration[CBCL2];
        # 30. cBcl2 + tBid ----> cBcl2:tBid : complex formation
        reactionFlux[CBCL2TBID, ASSOCIATION] = Srates[CBCL2TBID, ASSOCIATION] * concentration[CBCL2] * concentration[TBID];
        # 31. cBcl2:tBid ----> cBcl2 + tBid : complex dissociation
        reactionFlux[CBCL2TBID, DISSOCIATION] = Srates[CBCL2TBID, DISSOCIATION] * concentration[CBCL2TBID];
        # 32. cBcl2:tBid ----> 0 : Basal decay
        reactionFlux[CBCL2TBID, BASALDECAY] = Srates[CBCL2TBID, BASALDECAY] * concentration[CBCL2TBID];

        #---------------------------------------------
        # MODULE 5: Oligomerization of Bax
        #---------------------------------------------
        # 33. 0 ----> Bax : Basal synthesis
        reactionFlux[BAX, BASALSYNTHESIS] = EPS * Srates[BAX, BASALSYNTHESIS];
        # 34. Bax ----> 0 : Basal decay
        reactionFlux[BAX, BASALDECAY] = Srates[BAX, BASALDECAY] * concentration[BAX];
        # 35. tBid + Bax ----> tBid:Bax : complex formation
        reactionFlux[TBIDBAX, ASSOCIATION] = Srates[TBIDBAX, ASSOCIATION] * concentration[TBID] * concentration[BAX];
        # 36. tBid:Bax ----> tBid + Bax : complex dissociation
        reactionFlux[TBIDBAX, DISSOCIATION] = Srates[TBIDBAX, DISSOCIATION] * concentration[TBIDBAX];
        # 37. tBid:Bax ----> tBid + aBax : activation
        reactionFlux[ABAX, ACTIVATION] = Srates[ABAX, ACTIVATION] * concentration[TBIDBAX];

        # 38. aBax ----> mBax : transport into mitochondria
        reactionFlux[ABAX, TRANSPORTIN] = Srates[ABAX, TRANSPORTIN] * concentration[ABAX];
        # 39. mBax ----> aBax : transport out of mitochondria
        reactionFlux[MBAX, TRANSPORTOUT] = Srates[MBAX, TRANSPORTOUT] * concentration[MBAX];

        # 40. tBcl2 ----> Bcl2 : translation
        if phase == 1 # assume that during equilibration the cell is receiving a survival signal
            reactionFlux[BCL2, TRANSLATION] = EPS * BASAL_TBCL2;
        else
            reactionFlux[BCL2, TRANSLATION] = EPS * Srates[BCL2, TRANSLATION] * concentration[TBCL2];
        end

        # 41. Bcl2 ----> 0 : Basal decay
        reactionFlux[BCL2, BASALDECAY] = Srates[BCL2, BASALDECAY] * concentration[BCL2];
        # 41+1. Bcl2 + ABCR ----> ABCR : Decay
        reactionFlux[BCL2, DECAY] = Srates[BCL2, DECAY] * concentration[BCL2] * concentration[ABCR];
        # 42. mBax + Bcl2 ----> mBax:Bcl2 : complex formation
        reactionFlux[MBAXBCL2, ASSOCIATION] = Srates[MBAXBCL2, ASSOCIATION] * concentration[MBAX] * concentration[BCL2];
        # 43. mBax:Bcl2 ----> mBax + Bcl2 : complex dissociation
        reactionFlux[MBAXBCL2, DISSOCIATION] = Srates[MBAXBCL2, DISSOCIATION] * concentration[MBAXBCL2];
        # 44. mBax:Bcl2 ----> 0 : Basal decay
        reactionFlux[MBAXBCL2, BASALDECAY] = Srates[MBAXBCL2, BASALDECAY] * concentration[MBAXBCL2];

        # 45. mBax + mBax ----> Bax2 : complex formation
        reactionFlux[BAX2, ASSOCIATION] = Srates[BAX2, ASSOCIATION] * concentration[MBAX] ^ 2;
        # 46. Bax2 ----> mBax + mBax : complex dissociation
        reactionFlux[BAX2, DISSOCIATION] = Srates[BAX2, DISSOCIATION] * concentration[BAX2];
        # 47. Bax2 + Bcl2 ----> Bax2:Bcl2 : complex formation
        reactionFlux[BAX2BCL2, ASSOCIATION] = Srates[BAX2BCL2, ASSOCIATION] * concentration[BAX2] * concentration[BCL2];
        # 48. Bax2:Bcl2 ----> Bax2 + Bcl2 : complex dissociation
        reactionFlux[BAX2BCL2, DISSOCIATION] = Srates[BAX2BCL2, DISSOCIATION] * concentration[BAX2BCL2];
        # 49. Bax2:Bcl2 ----> 0 : Basal decay
        reactionFlux[BAX2BCL2, BASALDECAY] = Srates[BAX2BCL2, BASALDECAY] * concentration[BAX2BCL2];

        # 50. Bax2 + Bax2 ----> Bax4 : complex formation
        reactionFlux[BAX4, ASSOCIATION] = Srates[BAX4, ASSOCIATION] * concentration[BAX2] ^ 2;
        # 51. Bax4 ----> Bax2 + Bax2 : complex dissociation
        reactionFlux[BAX4, DISSOCIATION] = Srates[BAX4, DISSOCIATION] * concentration[BAX4];
        # 52. Bax4 + Bcl2 ----> Bax4:Bcl2 : complex formation
        reactionFlux[BAX4BCL2, ASSOCIATION] = Srates[BAX4BCL2, ASSOCIATION] * concentration[BAX4] * concentration[BCL2];
        # 53. Bax4:Bcl2 ----> Bax4 + Bcl2 : complex dissociation
        reactionFlux[BAX4BCL2, DISSOCIATION] = Srates[BAX4BCL2, DISSOCIATION] * concentration[BAX4BCL2];
        # 54. Bax4:Bcl2 ----> 0 : Basal decay
        reactionFlux[BAX4BCL2, BASALDECAY] = Srates[BAX4BCL2, BASALDECAY] * concentration[BAX4BCL2];

        #----------------------------------------
        # MODULE 6: MOMP (pore-forming and transporting)
        #----------------------------------------
        # 55. Bax4 + Mito ----> Bax4:Mito : complex formation
        reactionFlux[BAX4MITO, ASSOCIATION] = Srates[BAX4MITO, ASSOCIATION] * concentration[BAX4] * concentration[MITO];
        # 56. Bax4:Mito ----> Bax4 + Mito : complex dissociation
        reactionFlux[BAX4MITO, DISSOCIATION] = Srates[BAX4MITO, DISSOCIATION] * concentration[BAX4MITO];
        # 57. Bax4:Mito ----> aMito : activation
        reactionFlux[AMITO, ACTIVATION] = Srates[AMITO, ACTIVATION] * concentration[BAX4MITO];

        # 58. 0 ----> mCytoC : Basal synthesis
        reactionFlux[MCYTOC, BASALSYNTHESIS] = EPS * Srates[MCYTOC, BASALSYNTHESIS];
        # 59. mCytoC ----> 0 : Basal decay
        reactionFlux[MCYTOC, BASALDECAY] = Srates[MCYTOC, BASALDECAY] * concentration[MCYTOC];
        # 60. aMito + mCytoC ----> aMito:mCytoC : complex formation
        reactionFlux[AMITOMCYTOC, ASSOCIATION] = Srates[AMITOMCYTOC, ASSOCIATION] * concentration[AMITO] * concentration[MCYTOC];
        # 61. aMito:mCytoC ----> aMito + mCytoC : complex dissociation
        reactionFlux[AMITOMCYTOC, DISSOCIATION] = Srates[AMITOMCYTOC, DISSOCIATION] * concentration[AMITOMCYTOC];
        # 62. aMito:mCytoC ----> aMito + aCytoC : activation
        reactionFlux[ACYTOC, ACTIVATION] = Srates[ACYTOC, ACTIVATION] * concentration[AMITOMCYTOC];

        # 63. 0 ----> mSmac : Basal synthesis
        reactionFlux[MSMAC, BASALSYNTHESIS] = EPS * Srates[MSMAC, BASALSYNTHESIS];
        # 64. mSmac ----> 0 : Basal decay
        reactionFlux[MSMAC, BASALDECAY] = Srates[MSMAC, BASALDECAY] * concentration[MSMAC];
        # 65. aMito + mSmac ----> aMito:mSmac : complex formation
        reactionFlux[AMITOMSMAC, ASSOCIATION] = Srates[AMITOMSMAC, ASSOCIATION] * concentration[AMITO] * concentration[MSMAC];
        # 66. aMito:mSmac ----> aMito + mSmac : complex dissociation
        reactionFlux[AMITOMSMAC, DISSOCIATION] = Srates[AMITOMSMAC, DISSOCIATION] * concentration[AMITOMSMAC];
        # 67. aMito:mSmac ----> aMito + aSmac : activation
        reactionFlux[ASMAC, ACTIVATION] = Srates[ASMAC, ACTIVATION] * concentration[AMITOMSMAC];

        # 68. aCytoC ----> cCytoC : transport out of mitochondria
        reactionFlux[ACYTOC, TRANSPORTOUT] = Srates[ACYTOC, TRANSPORTOUT] * concentration[ACYTOC];
        # 69. cCytoC ----> aCytoC : transport into mitochondria
        reactionFlux[CCYTOC, TRANSPORTIN] = Srates[CCYTOC, TRANSPORTIN] * concentration[CCYTOC];

        # 70. aSmac ----> cSmac : transport out of mitochondria
        reactionFlux[ASMAC, TRANSPORTOUT] = Srates[ASMAC, TRANSPORTOUT] * concentration[ASMAC];
        # 71. cSmac ----> aSmac : transport into mitochondria
        reactionFlux[CSMAC, TRANSPORTIN] = Srates[CSMAC, TRANSPORTIN] * concentration[CSMAC];

        # 72. aMito ----> Mito : deactivation
        reactionFlux[AMITO, DEACTIVATION] = Srates[AMITO, DEACTIVATION] * concentration[AMITO];

        #----------------------------------------
        # MODULE 7: Feedforward pathway 1 (XIAP)
        #----------------------------------------
        # 73. 0 ----> XIAP : Basal synthesis
        reactionFlux[XIAP, BASALSYNTHESIS] = EPS * Srates[XIAP, BASALSYNTHESIS];
        # 74. XIAP ----> 0 : Basal decay
        reactionFlux[XIAP, BASALDECAY] = Srates[XIAP, BASALDECAY] * concentration[XIAP];

        # 75. cSmac ----> 0 : Basal decay
        reactionFlux[CSMAC, BASALDECAY] = Srates[CSMAC, BASALDECAY] * concentration[CSMAC];

        # 76. cSmac + XIAP ----> cSmac:XIAP : complex formation
        reactionFlux[CSMACXIAP, ASSOCIATION] = Srates[CSMACXIAP, ASSOCIATION] * concentration[CSMAC] * concentration[XIAP];
        # 77. cSmac:XIAP ----> cSmac + XIAP : complex dissociation
        reactionFlux[CSMACXIAP, DISSOCIATION] = Srates[CSMACXIAP, DISSOCIATION] * concentration[CSMACXIAP];
        # 78. cSmac:XIAP ----> 0 : Basal decay
        reactionFlux[CSMACXIAP, BASALDECAY] = Srates[CSMACXIAP, BASALDECAY] * concentration[CSMACXIAP];

        #----------------------------------------
        # MODULE 8: Feedforward pathway 2 (Apoptosome)
        #----------------------------------------
        # 79. cCytoC ----> 0 : Basal decay
        reactionFlux[CCYTOC, BASALDECAY] = Srates[CCYTOC, BASALDECAY] * concentration[CCYTOC];

        # 80. Apaf + cCytoC ----> Apaf:cCytoC : complex formation
        reactionFlux[APAFCCYTOC, ASSOCIATION] = Srates[APAFCCYTOC, ASSOCIATION] * concentration[APAF] * concentration[CCYTOC];
        # 81. Apaf:cCytoC ----> Apaf + cCytoC : complex dissociation
        reactionFlux[APAFCCYTOC, DISSOCIATION] = Srates[APAFCCYTOC, DISSOCIATION] * concentration[APAFCCYTOC];
        # 82. Apaf:cCytoC ----> aApaf + cCytoC : activation
        reactionFlux[AAPAF, ACTIVATION] = Srates[AAPAF, ACTIVATION] * concentration[APAFCCYTOC];

        # 83. aApaf + pC9 ----> Apop : complex formation
        reactionFlux[APOP, ASSOCIATION] = Srates[APOP, ASSOCIATION] * concentration[AAPAF] * concentration[PC9];
        # 84. Apop ----> aApaf + pC9 : complex dissociation
        reactionFlux[APOP, DISSOCIATION] = Srates[APOP, DISSOCIATION] * concentration[APOP];

        # 85. aApaf ----> Apaf : deactivation
        reactionFlux[AAPAF, DEACTIVATION] = Srates[AAPAF, DEACTIVATION] * concentration[AAPAF];

        # 86. Apop + XIAP ----> Apop:XIAP : complex formation
        reactionFlux[APOPXIAP, ASSOCIATION] = Srates[APOPXIAP, ASSOCIATION] * concentration[APOP] * concentration[XIAP];
        # 87. Apop:XIAP ----> Apop + XIAP : complex dissociation
        reactionFlux[APOPXIAP, DISSOCIATION] = Srates[APOPXIAP, DISSOCIATION] * concentration[APOPXIAP];
        # 88. Apop:XIAP ----> Apop : deactivation
        reactionFlux[APOP, DEACTIVATION] = Srates[APOP, DEACTIVATION] * concentration[APOPXIAP];

        #----------------------------------------
        # MODULE 9: CASPASE 3 module (effector caspase)
        #----------------------------------------
        # 89. 0 ----> pC3 : Basal synthesis
        reactionFlux[PC3, BASALSYNTHESIS] = EPS * Srates[PC3, BASALSYNTHESIS];
        # 90. pC3 ----> 0 : Basal decay
        reactionFlux[PC3, BASALDECAY] = Srates[PC3, BASALDECAY] * concentration[PC3];
        # 91. C8 + pC3 ----> C8:pC3 : complex formation
        reactionFlux[C8PC3, ASSOCIATION] = Srates[C8PC3, ASSOCIATION] * concentration[C8] * concentration[PC3];
        # 92. C8:pC3 ----> C8 + pC3 : complex dissociation
        reactionFlux[C8PC3, DISSOCIATION] = Srates[C8PC3, DISSOCIATION] * concentration[C8PC3];
        # 93. C8:pC3 ----> C8 + C3 : activation
        reactionFlux[C3, ACTIVATION] = Srates[C3, ACTIVATION] * concentration[C8PC3];

        # 94. XIAP + C3 ----> XIAP:C3 : complex formation
        reactionFlux[XIAPC3, ASSOCIATION] = Srates[XIAPC3, ASSOCIATION] * concentration[XIAP] * concentration[C3];
        # 95. XIAP:C3 ----> XIAP + C3 : complex dissociation
        reactionFlux[XIAPC3, DISSOCIATION] = Srates[XIAPC3, DISSOCIATION] * concentration[XIAPC3];
        # 96. XIAP:C3 ----> XIAP + uC3 : deactivation
        reactionFlux[C3, DEACTIVATION] = Srates[C3, DEACTIVATION] * concentration[XIAPC3];
        # 97. uC3 ----> 0 : Basal decay
        reactionFlux[UC3, BASALDECAY] = Srates[UC3, BASALDECAY] * concentration[UC3];

        # 98. Apop + pC3 ----> Apop:pC3 : complex formation
        reactionFlux[APOPPC3, ASSOCIATION] = Srates[APOPPC3, ASSOCIATION] * concentration[APOP] * concentration[PC3];
        # 99. Apop:pC3 ----> Apop + pC3 : complex dissociation
        reactionFlux[APOPPC3, DISSOCIATION] = Srates[APOPPC3, DISSOCIATION] * concentration[APOPPC3];
        # 100. Apop:pC3 ----> Apop + C3 : synthesis
        reactionFlux[C3, SYNTHESIS] = Srates[C3, SYNTHESIS] * concentration[APOPPC3];

        #----------------------------------------
        # MODULE 10: CASPASE 6 feedback module
        #----------------------------------------
        # 101. 0 ----> pC6 : Basal synthesis
        reactionFlux[PC6, BASALSYNTHESIS] = EPS * Srates[PC6, BASALSYNTHESIS];
        # 102. pC6 ----> 0 : Basal decay
        reactionFlux[PC6, BASALDECAY] = Srates[PC6, BASALDECAY] * concentration[PC6];
        # 103. C3 + pC6 ----> C3:pC6 : complex formation
        reactionFlux[C3PC6, ASSOCIATION] = Srates[C3PC6, ASSOCIATION] * concentration[C3] * concentration[PC6];
        # 104. C3:pC6 ----> C3 + pC6 : complex dissociation
        reactionFlux[C3PC6, DISSOCIATION] = Srates[C3PC6, DISSOCIATION] * concentration[C3PC6];
        # 105. C3:pC6 ----> C3 + C6 : activation
        reactionFlux[C6, ACTIVATION] = Srates[C6, ACTIVATION] * concentration[C3PC6];

        # 106. C6 ----> 0 : Basal decay
        reactionFlux[C6, BASALDECAY] = Srates[C6, BASALDECAY] * concentration[C6];

        # 107. C6 + pC8 ----> C6:pC8 : complex formation
        reactionFlux[C6PC8, ASSOCIATION] = Srates[C6PC8, ASSOCIATION] * concentration[C6] * concentration[PC8];
        # 108. C6:pC8 ----> C6 + pC8 : complex dissociation
        reactionFlux[C6PC8, DISSOCIATION] = Srates[C6PC8, DISSOCIATION] * concentration[C6PC8];
        # 109. C6:pC8 ----> C6 + C8 : synthesis
        reactionFlux[C8, SYNTHESIS] = Srates[C8, SYNTHESIS] * concentration[C6PC8];

        #----------------------------------------
        # MODULE 11: Cell death module
        #----------------------------------------
        # 110. 0 ----> PARP : Basal synthesis
        reactionFlux[PARP, BASALSYNTHESIS] = EPS * Srates[PARP, BASALSYNTHESIS];
        # 111. PARP ----> 0 : Basal decay
        reactionFlux[PARP, BASALDECAY] = Srates[PARP, BASALDECAY] * concentration[PARP];

        # 112. C3 + PARP ----> C3:PARP : complex formation
        reactionFlux[C3PARP, ASSOCIATION] = Srates[C3PARP, ASSOCIATION] * concentration[C3] * concentration[PARP];
        # 113. C3:PARP ----> C3 + PARP : complex dissociation
        reactionFlux[C3PARP, DISSOCIATION] = Srates[C3PARP, DISSOCIATION] * concentration[C3PARP];
        # 114. C3:PARP ----> C3 + tPARP : deactivation
        reactionFlux[PARP, DEACTIVATION] = Srates[PARP, DEACTIVATION] * concentration[C3PARP];

        # 115. tPARP ----> 0 : Basal decay
        reactionFlux[CPARP, BASALDECAY] = Srates[CPARP, BASALDECAY] * concentration[CPARP];

    end
    nothing
end

function ApoptosisNettFluxes!(nettFlux, reactionFlux)

    # Compute apoptosis net fluxes
    @inbounds begin
        #------------------------------------------------
        # MODULE 1
        #--------------------------------------------
        # 1 : tBcl2
        nettFlux[TBCL2] = reactionFlux[TBCL2, SYNTHESIS] - reactionFlux[TBCL2, BASALDECAY] - reactionFlux[TBCL2, DECAY];

        #------------------------------------------------
        # MODULE 2
        #--------------------------------------------
        # 2 : L
        #     nettFlux[L] = reactionFlux[L, BASALSYNTHESIS] - reactionFlux[L, BASALDECAY] - reactionFlux[LR, ASSOCIATION] + reactionFlux[LR, DISSOCIATION];
        nettFlux[L] = 0;
        # 3 : R
        nettFlux[R] = reactionFlux[R, BASALSYNTHESIS] - reactionFlux[R, BASALDECAY] - reactionFlux[LR, ASSOCIATION] + reactionFlux[LR, DISSOCIATION];
        # 4 : L-R
        nettFlux[LR] = reactionFlux[LR, ASSOCIATION] - reactionFlux[LR, DISSOCIATION] - reactionFlux[DISC, ACTIVATION];
        # 5 : DISC
        nettFlux[DISC] = reactionFlux[DISC, ACTIVATION] - reactionFlux[FLIPDISC, ASSOCIATION] + reactionFlux[FLIPDISC, DISSOCIATION] - reactionFlux[DISCPC8, ASSOCIATION] + reactionFlux[DISCPC8, DISSOCIATION] + reactionFlux[C8, ACTIVATION];
        # 6 : flip
        nettFlux[FLIP] = reactionFlux[FLIP, BASALSYNTHESIS] - reactionFlux[FLIP, BASALDECAY] - reactionFlux[FLIPDISC, ASSOCIATION] + reactionFlux[FLIPDISC, DISSOCIATION];
        # 7 : flip:DISC
        nettFlux[FLIPDISC] = reactionFlux[FLIPDISC, ASSOCIATION] - reactionFlux[FLIPDISC, DISSOCIATION] - reactionFlux[FLIPDISC, BASALDECAY];

        #------------------------------------------------
        # MODULE 3
        #--------------------------------------------
        # 8 : pC8
        nettFlux[PC8] = reactionFlux[PC8, BASALSYNTHESIS] - reactionFlux[C8, BASALSYNTHESIS] - reactionFlux[PC8, BASALDECAY] - reactionFlux[DISCPC8, ASSOCIATION] + reactionFlux[DISCPC8, DISSOCIATION] - reactionFlux[C6PC8, ASSOCIATION] + reactionFlux[C6PC8, DISSOCIATION];
        # 9 : DISC:pC8
        nettFlux[DISCPC8] = reactionFlux[DISCPC8, ASSOCIATION] - reactionFlux[DISCPC8, DISSOCIATION] - reactionFlux[C8, ACTIVATION];
        # 10 : C8
        nettFlux[C8] = reactionFlux[C8, ACTIVATION] + reactionFlux[C8, BASALSYNTHESIS] - reactionFlux[BARC8, ASSOCIATION] + reactionFlux[BARC8, DISSOCIATION] - reactionFlux[C8BID, ASSOCIATION] + reactionFlux[C8BID, DISSOCIATION] + reactionFlux[TBID, ACTIVATION] - reactionFlux[C8PC3, ASSOCIATION] + reactionFlux[C8PC3, DISSOCIATION] + reactionFlux[C3, ACTIVATION] + reactionFlux[C8, SYNTHESIS];
        # 11 : Bar
        nettFlux[BAR] = reactionFlux[BAR, BASALSYNTHESIS] - reactionFlux[BAR, BASALDECAY] - reactionFlux[BARC8, ASSOCIATION] + reactionFlux[BARC8, DISSOCIATION];
        # 12 : Bar:C8
        nettFlux[BARC8] = reactionFlux[BARC8, ASSOCIATION] - reactionFlux[BARC8, DISSOCIATION] - reactionFlux[BARC8, BASALDECAY];

        #------------------------------------------------
        # MODULE 4
        #--------------------------------------------
        # 13 : Bid
        nettFlux[BID] = reactionFlux[BID, BASALSYNTHESIS] - reactionFlux[BID, BASALDECAY] - reactionFlux[C8BID, ASSOCIATION] + reactionFlux[C8BID, DISSOCIATION];
        # 14 : C8:Bid
        nettFlux[C8BID] = reactionFlux[C8BID, ASSOCIATION] - reactionFlux[C8BID, DISSOCIATION] - reactionFlux[TBID, ACTIVATION];
        # 15 : tBid
        nettFlux[TBID] = reactionFlux[TBID, ACTIVATION] - reactionFlux[CBCL2TBID, ASSOCIATION] + reactionFlux[CBCL2TBID, DISSOCIATION] - reactionFlux[TBIDBAX, ASSOCIATION] + reactionFlux[TBIDBAX, DISSOCIATION] + reactionFlux[ABAX, ACTIVATION];
        # 16 : cBcl2
        nettFlux[CBCL2] = reactionFlux[CBCL2, BASALSYNTHESIS] - reactionFlux[CBCL2, BASALDECAY] - reactionFlux[CBCL2TBID, ASSOCIATION] + reactionFlux[CBCL2TBID, DISSOCIATION];
        # 17 : cBcl2:tBid
        nettFlux[CBCL2TBID] = reactionFlux[CBCL2TBID, ASSOCIATION] - reactionFlux[CBCL2TBID, DISSOCIATION] - reactionFlux[CBCL2TBID, BASALDECAY];

        #------------------------------------------------
        # MODULE 5
        #--------------------------------------------
        # 18 : Bax
        nettFlux[BAX] = reactionFlux[BAX, BASALSYNTHESIS] - reactionFlux[BAX, BASALDECAY] - reactionFlux[TBIDBAX, ASSOCIATION] + reactionFlux[TBIDBAX, DISSOCIATION];
        # 19 : tBid:Bax
        nettFlux[TBIDBAX] = reactionFlux[TBIDBAX, ASSOCIATION] - reactionFlux[TBIDBAX, DISSOCIATION] - reactionFlux[ABAX, ACTIVATION];
        # 20 : aBax
        nettFlux[ABAX] = reactionFlux[ABAX, ACTIVATION] - reactionFlux[ABAX, TRANSPORTIN] + reactionFlux[MBAX, TRANSPORTOUT];
        # 21 : mBax
        nettFlux[MBAX] = reactionFlux[ABAX, TRANSPORTIN] - reactionFlux[MBAX, TRANSPORTOUT] - reactionFlux[MBAXBCL2, ASSOCIATION] * SCALE_MITO2CELLULAR + reactionFlux[MBAXBCL2, DISSOCIATION] - 2 * reactionFlux[BAX2, ASSOCIATION] * SCALE_MITO2CELLULAR + 2 * reactionFlux[BAX2, DISSOCIATION];
        # 22 : Bcl2
        nettFlux[BCL2] = reactionFlux[BCL2, TRANSLATION] - reactionFlux[BCL2, BASALDECAY] - reactionFlux[BCL2, DECAY] - reactionFlux[MBAXBCL2, ASSOCIATION] * SCALE_MITO2CELLULAR + reactionFlux[MBAXBCL2, DISSOCIATION] - reactionFlux[BAX2BCL2, ASSOCIATION] * SCALE_MITO2CELLULAR + reactionFlux[BAX2BCL2, DISSOCIATION] - reactionFlux[BAX4BCL2, ASSOCIATION] * SCALE_MITO2CELLULAR + reactionFlux[BAX4BCL2, DISSOCIATION];
        # 23 : mBax:Bcl2
        nettFlux[MBAXBCL2] = reactionFlux[MBAXBCL2, ASSOCIATION] * SCALE_MITO2CELLULAR - reactionFlux[MBAXBCL2, DISSOCIATION] - reactionFlux[MBAXBCL2, BASALDECAY];
        # 24 : Bax2
        nettFlux[BAX2] = reactionFlux[BAX2, ASSOCIATION] * SCALE_MITO2CELLULAR - reactionFlux[BAX2, DISSOCIATION] - reactionFlux[BAX2BCL2, ASSOCIATION] * SCALE_MITO2CELLULAR + reactionFlux[BAX2BCL2, DISSOCIATION] - 2 * reactionFlux[BAX4, ASSOCIATION] * SCALE_MITO2CELLULAR + 2 * reactionFlux[BAX4, DISSOCIATION];
        # 25 : Bax2:Bcl2
        nettFlux[BAX2BCL2] = reactionFlux[BAX2BCL2, ASSOCIATION] * SCALE_MITO2CELLULAR - reactionFlux[BAX2BCL2, DISSOCIATION] - reactionFlux[BAX2BCL2, BASALDECAY];
        # 26 : Bax4
        nettFlux[BAX4] = reactionFlux[BAX4, ASSOCIATION] * SCALE_MITO2CELLULAR - reactionFlux[BAX4, DISSOCIATION] - reactionFlux[BAX4BCL2, ASSOCIATION] * SCALE_MITO2CELLULAR + reactionFlux[BAX4BCL2, DISSOCIATION] - reactionFlux[BAX4MITO, ASSOCIATION] * SCALE_MITO2CELLULAR + reactionFlux[BAX4MITO, DISSOCIATION];
        # 27 : Bax4:Bcl2
        nettFlux[BAX4BCL2] = reactionFlux[BAX4BCL2, ASSOCIATION] * SCALE_MITO2CELLULAR - reactionFlux[BAX4BCL2, DISSOCIATION] - reactionFlux[BAX4BCL2, BASALDECAY];

        #------------------------------------------------
        # MODULE 6
        #--------------------------------------------
        # 28 : Mito
        nettFlux[MITO] = - reactionFlux[BAX4MITO, ASSOCIATION] * SCALE_MITO2CELLULAR + reactionFlux[BAX4MITO, DISSOCIATION] + reactionFlux[AMITO, DEACTIVATION];
        # 29 : Bax4:Mito
        nettFlux[BAX4MITO] = reactionFlux[BAX4MITO, ASSOCIATION] * SCALE_MITO2CELLULAR - reactionFlux[BAX4MITO, DISSOCIATION] - reactionFlux[AMITO, ACTIVATION];
        # 30 : aMito
        nettFlux[AMITO] = reactionFlux[AMITO, ACTIVATION] - reactionFlux[AMITOMCYTOC, ASSOCIATION] * SCALE_MITO2CELLULAR + reactionFlux[AMITOMCYTOC, DISSOCIATION] + reactionFlux[ACYTOC, ACTIVATION] - reactionFlux[AMITOMSMAC, ASSOCIATION] * SCALE_MITO2CELLULAR + reactionFlux[AMITOMSMAC, DISSOCIATION] + reactionFlux[ASMAC, ACTIVATION] - reactionFlux[AMITO, DEACTIVATION];
        # 31 : mCytoC
        nettFlux[MCYTOC] = reactionFlux[MCYTOC, BASALSYNTHESIS] - reactionFlux[MCYTOC, BASALDECAY] - reactionFlux[AMITOMCYTOC, ASSOCIATION] * SCALE_MITO2CELLULAR + reactionFlux[AMITOMCYTOC, DISSOCIATION];
        # 32 : aMito:mCytoC
        nettFlux[AMITOMCYTOC] = reactionFlux[AMITOMCYTOC, ASSOCIATION] * SCALE_MITO2CELLULAR - reactionFlux[AMITOMCYTOC, DISSOCIATION] - reactionFlux[ACYTOC, ACTIVATION];
        # 33 : aCytoC
        nettFlux[ACYTOC] = reactionFlux[ACYTOC, ACTIVATION] - reactionFlux[ACYTOC, TRANSPORTOUT] + reactionFlux[CCYTOC, TRANSPORTIN];
        # 34 : mSmac
        nettFlux[MSMAC] = reactionFlux[MSMAC, BASALSYNTHESIS] - reactionFlux[MSMAC, BASALDECAY] - reactionFlux[AMITOMSMAC, ASSOCIATION] * SCALE_MITO2CELLULAR + reactionFlux[AMITOMSMAC, DISSOCIATION];
        # 35 : aMito:mSmac
        nettFlux[AMITOMSMAC] = reactionFlux[AMITOMSMAC, ASSOCIATION] * SCALE_MITO2CELLULAR - reactionFlux[AMITOMSMAC, DISSOCIATION] - reactionFlux[ASMAC, ACTIVATION];
        # 36 : aSmac
        nettFlux[ASMAC] = reactionFlux[ASMAC, ACTIVATION] - reactionFlux[ASMAC, TRANSPORTOUT] + reactionFlux[CSMAC, TRANSPORTIN];

        #------------------------------------------------
        # MODULE 7
        #--------------------------------------------
        # 37 : XIAP
        nettFlux[XIAP] = reactionFlux[XIAP, BASALSYNTHESIS] - reactionFlux[XIAP, BASALDECAY] - reactionFlux[CSMACXIAP, ASSOCIATION] + reactionFlux[CSMACXIAP, DISSOCIATION] - reactionFlux[APOPXIAP, ASSOCIATION] + reactionFlux[APOPXIAP, DISSOCIATION] - reactionFlux[XIAPC3, ASSOCIATION] + reactionFlux[XIAPC3, DISSOCIATION] + reactionFlux[C3, DEACTIVATION];
        # 38 : cSmac
        nettFlux[CSMAC] = reactionFlux[ASMAC, TRANSPORTOUT] - reactionFlux[CSMAC, TRANSPORTIN] - reactionFlux[CSMAC, BASALDECAY] - reactionFlux[CSMACXIAP, ASSOCIATION] + reactionFlux[CSMACXIAP, DISSOCIATION];
        # 39 : cSmac:XIAP
        nettFlux[CSMACXIAP] = reactionFlux[CSMACXIAP, ASSOCIATION] - reactionFlux[CSMACXIAP, DISSOCIATION] - reactionFlux[CSMACXIAP, BASALDECAY];

        #------------------------------------------------
        # MODULE 8
        #--------------------------------------------
        # 40 : cCytoC
        nettFlux[CCYTOC] = reactionFlux[ACYTOC, TRANSPORTOUT] - reactionFlux[CCYTOC, TRANSPORTIN] - reactionFlux[CCYTOC, BASALDECAY] - reactionFlux[APAFCCYTOC, ASSOCIATION] + reactionFlux[APAFCCYTOC, DISSOCIATION] + reactionFlux[AAPAF, ACTIVATION];
        # 41 : Apaf
        nettFlux[APAF] = - reactionFlux[APAFCCYTOC, ASSOCIATION] + reactionFlux[APAFCCYTOC, DISSOCIATION] + reactionFlux[AAPAF, DEACTIVATION];
        # 42 : Apaf:cCytoC
        nettFlux[APAFCCYTOC] = reactionFlux[APAFCCYTOC, ASSOCIATION] - reactionFlux[APAFCCYTOC, DISSOCIATION] - reactionFlux[AAPAF, ACTIVATION];
        # 43 : aApaf
        nettFlux[AAPAF] = reactionFlux[AAPAF, ACTIVATION] - reactionFlux[APOP, ASSOCIATION] + reactionFlux[APOP, DISSOCIATION] - reactionFlux[AAPAF, DEACTIVATION];
        # 44 : pC9
        nettFlux[PC9] = - reactionFlux[APOP, ASSOCIATION] + reactionFlux[APOP, DISSOCIATION];
        # 45 : Apop
        nettFlux[APOP] = reactionFlux[APOP, ASSOCIATION] - reactionFlux[APOP, DISSOCIATION] - reactionFlux[APOPXIAP, ASSOCIATION] + reactionFlux[APOPXIAP, DISSOCIATION] + reactionFlux[APOP, DEACTIVATION] - reactionFlux[APOPPC3, ASSOCIATION] + reactionFlux[APOPPC3, DISSOCIATION] + reactionFlux[C3, SYNTHESIS];
        # 46 : Apop:XIAP
        nettFlux[APOPXIAP] = reactionFlux[APOPXIAP, ASSOCIATION] - reactionFlux[APOPXIAP, DISSOCIATION] - reactionFlux[APOP, DEACTIVATION];

        #------------------------------------------------
        # MODULE 9
        #--------------------------------------------
        # 47 : pC3
        nettFlux[PC3] = reactionFlux[PC3, BASALSYNTHESIS] - reactionFlux[PC3, BASALDECAY] - reactionFlux[C8PC3, ASSOCIATION] + reactionFlux[C8PC3, DISSOCIATION] - reactionFlux[APOPPC3, ASSOCIATION] + reactionFlux[APOPPC3, DISSOCIATION];
        # 48 : C8:pC3
        nettFlux[C8PC3] = reactionFlux[C8PC3, ASSOCIATION] - reactionFlux[C8PC3, DISSOCIATION] - reactionFlux[C3, ACTIVATION];
        # 49 : C3
        nettFlux[C3] = reactionFlux[C3, ACTIVATION] - reactionFlux[XIAPC3, ASSOCIATION] + reactionFlux[XIAPC3, DISSOCIATION] + reactionFlux[C3, SYNTHESIS] - reactionFlux[C3PC6, ASSOCIATION] + reactionFlux[C3PC6, DISSOCIATION] + reactionFlux[C6, ACTIVATION] - reactionFlux[C3PARP, ASSOCIATION] + reactionFlux[C3PARP, DISSOCIATION] + reactionFlux[PARP, DEACTIVATION]
        # 50 : XIAP:C3
        nettFlux[XIAPC3] = reactionFlux[XIAPC3, ASSOCIATION] - reactionFlux[XIAPC3, DISSOCIATION] - reactionFlux[C3, DEACTIVATION]
        # 51 : Apop:pC3
        nettFlux[APOPPC3] = reactionFlux[APOPPC3, ASSOCIATION] - reactionFlux[APOPPC3, DISSOCIATION] - reactionFlux[C3, SYNTHESIS]
        # 52 : uC3
        nettFlux[UC3] = reactionFlux[C3, DEACTIVATION] - reactionFlux[UC3, BASALDECAY]

        #------------------------------------------------
        # MODULE 10
        #--------------------------------------------
        # 53 : pC6
        nettFlux[PC6] = reactionFlux[PC6, BASALSYNTHESIS] - reactionFlux[PC6, BASALDECAY] - reactionFlux[C3PC6, ASSOCIATION] + reactionFlux[C3PC6, DISSOCIATION]
        # 54 : C3:pC6
        nettFlux[C3PC6] = reactionFlux[C3PC6, ASSOCIATION] - reactionFlux[C3PC6, DISSOCIATION] - reactionFlux[C6, ACTIVATION]
        # 55 : C6
        nettFlux[C6] = reactionFlux[C6, ACTIVATION] - reactionFlux[C6, BASALDECAY] - reactionFlux[C6PC8, ASSOCIATION] + reactionFlux[C6PC8, DISSOCIATION] + reactionFlux[C8, SYNTHESIS]
        # 56 : C6:pC8
        nettFlux[C6PC8] = reactionFlux[C6PC8, ASSOCIATION] - reactionFlux[C6PC8, DISSOCIATION] - reactionFlux[C8, SYNTHESIS]

        #------------------------------------------------
        # MODULE 11
        #--------------------------------------------
        # 57 : PARP
        nettFlux[PARP] = reactionFlux[PARP, BASALSYNTHESIS] - reactionFlux[PARP, BASALDECAY] - reactionFlux[C3PARP, ASSOCIATION] + reactionFlux[C3PARP, DISSOCIATION]
        # 58 : C3:PARP
        nettFlux[C3PARP] = reactionFlux[C3PARP, ASSOCIATION] - reactionFlux[C3PARP, DISSOCIATION] - reactionFlux[PARP, DEACTIVATION]
        # 59 : tPARP
        nettFlux[CPARP] = reactionFlux[PARP, DEACTIVATION] - reactionFlux[CPARP, BASALDECAY]
        #------------------------------------------------

    end
    nothing
end

# time-independent ODE (pre-simulation: phase = 1)
function computeApoptosisNettFluxes!(nettFlux, concentration, reactionFlux, Srates, phase)
    #------------------------------------------------
    # Compute Apoptosis reaction fluxes
    computeApoptosisFluxes!(concentration, reactionFlux, Srates, phase);
    # Compute Apoptosis net fluxes
    ApoptosisNettFluxes!(nettFlux, reactionFlux);
    nothing
end

# time-dependent ODE (simulation, phase = 2)
function computeApoptosisNettFluxes!(nettFlux, concentration, reactionFlux, Srates, phase, time, birthday, inputCurves)
    #------------------------------------------------
    if phase == 2
        concentration[ABCR] = inputCurves(time + birthday, idxs=ABCR);
        concentration[NA50] = inputCurves(time + birthday, idxs=NA50);
        concentration[NA52] = inputCurves(time + birthday, idxs=NA52);
        concentration[NC50] = inputCurves(time + birthday, idxs=NC50);
        concentration[NC52] = inputCurves(time + birthday, idxs=NC52);
        concentration[IKK] = inputCurves(time + birthday, idxs=IKK);
    end
    # Compute Apoptosis reaction fluxes
    computeApoptosisFluxes!(concentration, reactionFlux, Srates, phase);
    # Compute Apoptosis net fluxes
    ApoptosisNettFluxes!(nettFlux, reactionFlux);
    nothing
end
