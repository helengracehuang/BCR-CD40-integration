# 3rd edition of reaction rate parameters (tuned params for BCR-CBM pathway so BCR can activate NFkB) -- Works with ODE_Receptor3.jl
include("ConstantParams2.jl");
# Input values of rate constants for all species
# Rate constants are in hr-1 (converted from min-1 in Shokhirev model)
# Initial values are in [# of molecules] for apoptosis module and [nM] otherwise
#------------------------------------------------

# Function to set parameter values for any given species
function setParams!(rates_entry; initialConc=0, scale=1, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, phosphorylation=0, dephosphorylation=0, activation=0, deactivation=0, transport_in=0, transport_out=0)
    rates_entry[INITIALCONC] = initialConc
    rates_entry[SCALE] = scale
    rates_entry[BASALSYNTHESIS] = basalsynthesis * scale
    rates_entry[SYNTHESIS] = synthesis * scale
    rates_entry[TRANSLATION] = translation * scale
    rates_entry[BASALDECAY] = basaldecay * scale
    rates_entry[DECAY] = decay * scale
    rates_entry[REPRESSION] = 0.0
    rates_entry[ASSOCIATION] = association * scale
    rates_entry[DISSOCIATION] = dissociation * scale
    rates_entry[PHOSPHORYLATION] = phosphorylation * scale
    rates_entry[DEPHOSPHORYLATION] = dephosphorylation * scale
    rates_entry[ACTIVATION] = activation * scale
    rates_entry[DEACTIVATION] = deactivation * scale
    rates_entry[TRANSPORTIN] = transport_in * scale
    rates_entry[TRANSPORTOUT] = transport_out * scale
    nothing
end

function setAllRates!(rates)
    #################################
    ######## Receptor Model #########
    #################################
    # MODULE 1: B-cell Receptor species
    # 1 : Antigen
    setParams!((@view rates[ANTIGEN, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.05, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 2 : B-cell Receptor
    setParams!((@view rates[BCR, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=4.93, synthesis=0, translation=0, basaldecay=1.43, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 3 : Antigen-BCR (half-life = 2hrs)
    setParams!((@view rates[ABCR, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.35, decay=0, association=66.0, dissociation=1.26, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 4 : inactivated CBM-complex
    setParams!((@view rates[CBM, 1:end]), initialConc=100, scale=CONVERSION, basalsynthesis=6.84/100, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=66.0/100, deactivation=0, transport_in=0, transport_out=0);
    # 5 : activated CBM-complex
    setParams!((@view rates[ACBM, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=0, deactivation=1.26, transport_in=0, transport_out=0);
    # 6 : inhibited CBM-complex
    setParams!((@view rates[ICBM, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=18.06/1000, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    #--------------------------------------------
    # MODULE 2: CD40 Receptor species
    # 7 : CD40L (basaldecay: 0.115 for sCD154, 0.05 for a-CD40)
    setParams!((@view rates[CD40L, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.05, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 8 : CD40R
    setParams!((@view rates[CD40R, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=7.672, synthesis=0, translation=0, basaldecay=0.05, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 9 : CD40LR basaldecay +++
    setParams!((@view rates[CD40LR, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.17, decay=0, association=0.04, dissociation=11.3, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 10 : inactivated TRAF6
    setParams!((@view rates[TRAF6_OFF, 1:end]), initialConc=100, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 11 : activated TRAF6
    setParams!((@view rates[TRAF6, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0.5, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=0.1, deactivation=7.5, transport_in=0, transport_out=0);
    # 12 : TRAF3
    setParams!((@view rates[TRAF3, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=10, synthesis=0, translation=0, basaldecay=0.5, decay=10, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    #--------------------------------------------
    # MODULE 3: Kinase species
    # 13 : inactivated TAK1
    setParams!((@view rates[TAK1, 1:end]), initialConc=100, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=165600.0/140, association=0, dissociation=0, activation=0.006*10000, deactivation=0, transport_in=0, transport_out=0);
    # 14 : activated TAK1 (activation = 1050 for ACBM or 105000 for ABCR, basalsynthesis=24900.0/1000 for ACBM or 24900.0/100 for ABCR)
    setParams!((@view rates[ATAK1, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=24900.0/1000, synthesis=56160.0/140, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=1050*ATAK1_MOD, deactivation=258600.0, transport_in=0, transport_out=0);
    # 15 : inactivated IKK (IKK1)
    setParams!((@view rates[IKK_OFF, 1:end]), initialConc=140, scale=CONVERSION, basalsynthesis=6.84*140, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 16 : activated pIKK (IKK2)
    setParams!((@view rates[IKK2, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=57.66/100*140, deactivation=15.54*140, transport_in=0, transport_out=0);
    # 17 : activated ppIKK (IKK3)
    setParams!((@view rates[IKK3, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0.0000672*140, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=2094.0, deactivation=384.6*140, transport_in=0, transport_out=0);
    # 18 : inhibited pppIKK (IKK4)
    setParams!((@view rates[IIKK, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=18.06*140, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 19 : total activated IKK
    setParams!((@view rates[IKK, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 20 : NIK (MAP3K14)
    setParams!((@view rates[NIK, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=12*NIK_MOD, synthesis=0, translation=0, basaldecay=0.231, decay=2, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);



    #############################
    ######## NFkB Model #########
    #############################
    # 3 : tIkBa
    setParams!((@view rates[TIKBA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0.288, synthesis=0, translation=0, basaldecay=1.74, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 4 : tIkBb
    setParams!((@view rates[TIKBB, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0.072, synthesis=0, translation=0, basaldecay=0.174, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 5 : tIkBe
    setParams!((@view rates[TIKBE, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0.0043, synthesis=0, translation=0, basaldecay=0.228, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 6 : IkBa
    setParams!((@view rates[IKBA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=720.0, basaldecay=7.2, decay=0.084, association=0, dissociation=0, activation=0, deactivation=0, transport_in=3.6, transport_out=0);
    # 7 : IkBb
    setParams!((@view rates[IKBB, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=720.0, basaldecay=7.2, decay=0.027, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0.54, transport_out=0);
    # 8 : IkBe
    setParams!((@view rates[IKBE, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=720.0, basaldecay=0.72, decay=0.0204, association=0, dissociation=0, activation=0, deactivation=0, transport_in=2.7, transport_out=0);
    # 9 : IkBd
    setParams!((@view rates[IKBD, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.18, decay=36.0, association=0.72, dissociation=0.72, activation=0, deactivation=0, transport_in=2.7, transport_out=0);
    # 10 : IkBdn
    setParams!((@view rates[NIKBD, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.18, decay=0, association=0.72, dissociation=0.72, activation=0, deactivation=0, transport_in=0, transport_out=0.72);
    # 11 : IkBan
    setParams!((@view rates[NIKBA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=7.2, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0.72);
    # 12 : IkBbn
    setParams!((@view rates[NIKBB, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=7.2, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0.72);
    # 13 : IkBen
    setParams!((@view rates[NIKBE, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.72, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0.72);
    # 14 : tRelA
    setParams!((@view rates[TRELA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0.0069, synthesis=0, translation=0, basaldecay=0.174, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 15 : tp50
    setParams!((@view rates[TP50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0.002, synthesis=0, translation=0, basaldecay=0.174, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 16 : tRelB
    setParams!((@view rates[TRELB, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0.0031, synthesis=0, translation=0, basaldecay=0.174, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 17 : tP100
    setParams!((@view rates[TP100, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0.0009, synthesis=0, translation=0, basaldecay=0.0576, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 18 : tcRel
    setParams!((@view rates[TCREL, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0.0022, synthesis=0, translation=0, basaldecay=0.0576, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 19 : RelA
    setParams!((@view rates[RELA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=720.0, basaldecay=1.38, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 20 : p50
    setParams!((@view rates[P50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=720.0, basaldecay=1.38, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 21 : RelB
    setParams!((@view rates[RELB, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=720.0, basaldecay=1.38, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 22 : P100
    setParams!((@view rates[P100, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=720.0, basaldecay=1.38, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 23 : cRel
    setParams!((@view rates[CREL, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=720.0, basaldecay=1.38, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 24 : P52
    setParams!((@view rates[P52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=3.0, basaldecay=0, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 25 : RelAn
    setParams!((@view rates[NRELA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=1.38, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 26 : p50n
    setParams!((@view rates[NP50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=1.38, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 27 : RelBn
    setParams!((@view rates[NRELB, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=1.38, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 28 : P100n
    setParams!((@view rates[NP100, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=1.38, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 29 : cReln
    setParams!((@view rates[NCREL, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=1.38, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 30 : P52n
    setParams!((@view rates[NP52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=1.38, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 31 : RelA:RelA
    setParams!((@view rates[AA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0, dissociation=28.8, activation=0, deactivation=0, transport_in=324.0, transport_out=0);
    # 32 : RelA:p50
    setParams!((@view rates[A50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.114, dissociation=1.14, activation=0, deactivation=0, transport_in=324.0, transport_out=0);
    # 33 : RelA:p52
    setParams!((@view rates[A52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.0576, dissociation=0.114, activation=0, deactivation=0, transport_in=324.0, transport_out=0);
    # 34 : RelB:p50
    setParams!((@view rates[B50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.0576, dissociation=0.276, activation=0, deactivation=0, transport_in=324.0, transport_out=0);
    # 35 : RelB:p52
    setParams!((@view rates[B52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.0576, dissociation=0.84, activation=0, deactivation=0, transport_in=324.0, transport_out=0);
    # 36 : cRel:p50
    setParams!((@view rates[C50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.0576, dissociation=0.084, activation=0, deactivation=0, transport_in=324.0, transport_out=0);
    # 37 : cRel:p52
    setParams!((@view rates[C52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.0576, dissociation=0.084, activation=0, deactivation=0, transport_in=324.0, transport_out=0);
    # 38 : cRel:p100
    setParams!((@view rates[C100, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.0576, dissociation=0.084, activation=0, deactivation=0, transport_in=324.0, transport_out=0);
    # 39 : p50:p50
    setParams!((@view rates[P50P50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.108, dissociation=3.24, activation=0, deactivation=0, transport_in=324.0, transport_out=0);
    # 40 : p52:p52
    setParams!((@view rates[P52P52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.108, dissociation=3.24, activation=0, deactivation=0, transport_in=324.0, transport_out=0);
    # 41 : RelA:RelAn
    setParams!((@view rates[NAA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0, dissociation=2.88, activation=0, deactivation=0, transport_in=0, transport_out=0.288);
    # 42 : RelA:p50n
    setParams!((@view rates[NA50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.114, dissociation=0.114, activation=0, deactivation=0, transport_in=0, transport_out=0.288);
    # 43 : RelA:p52n
    setParams!((@view rates[NA52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.0576, dissociation=0.0114, activation=0, deactivation=0, transport_in=0, transport_out=0.288);
    # 44 : RelB:p50n
    setParams!((@view rates[NB50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.0576, dissociation=0.0276, activation=0, deactivation=0, transport_in=0, transport_out=0.288);
    # 45 : RelB:p52n
    setParams!((@view rates[NB52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.0576, dissociation=0.084, activation=0, deactivation=0, transport_in=0, transport_out=0.288);
    # 46 : cRel:p50n
    setParams!((@view rates[NC50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.0576, dissociation=0.0084, activation=0, deactivation=0, transport_in=0, transport_out=0.288);
    # 47 : cRel:p52n
    setParams!((@view rates[NC52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.0576, dissociation=0.0084, activation=0, deactivation=0, transport_in=0, transport_out=0.288);
    # 48 : cRel:p100n
    setParams!((@view rates[NC100, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.0576, dissociation=0.027, activation=0, deactivation=0, transport_in=0, transport_out=0.288);
    # 49 : p50:p50n
    setParams!((@view rates[NP50P50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.108, dissociation=3.24, activation=0, deactivation=0, transport_in=0, transport_out=0.288);
    # 50 : p52:p52n
    setParams!((@view rates[NP52P52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0, association=0.108, dissociation=3.24, activation=0, deactivation=0, transport_in=0, transport_out=0.288);
    # 51 : IkBa-RelA:RelA
    setParams!((@view rates[IKBAAA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0806, dissociation=0.0016, activation=0, deactivation=0, transport_in=16.56, transport_out=0);
    # 52 : IkBb-RelA:RelA
    setParams!((@view rates[IKBBAA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=18.0, dissociation=0.072, activation=0, deactivation=0, transport_in=1.656, transport_out=0);
    # 53 : IkBe-RelA:RelA
    setParams!((@view rates[IKBEAA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0128, dissociation=1.0152, activation=0, deactivation=0, transport_in=8.28, transport_out=0);
    # 54 : IkBd-RelA:RelA
    setParams!((@view rates[IKBDAA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.36, dissociation=3.6, activation=0, deactivation=0, transport_in=8.28, transport_out=0);
    # 55 : IkBa-RelA:RelAn
    setParams!((@view rates[NIKBAAA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0806, dissociation=0.0016, activation=0, deactivation=0, transport_in=0, transport_out=50.4);
    # 56 : IkBb-RelA:RelAn
    setParams!((@view rates[NIKBBAA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=18.0, dissociation=0.072, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 57 : IkBe-RelA:RelAn
    setParams!((@view rates[NIKBEAA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0128, dissociation=1.0152, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 58 : IkBd-RelA:RelAn
    setParams!((@view rates[NIKBDAA, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.36, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 59 : IkBa-RelA:p50
    setParams!((@view rates[IKBAA50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.18, dissociation=0.036, activation=0, deactivation=0, transport_in=16.56, transport_out=0);
    # 60 : IkBb-RelA:p50
    setParams!((@view rates[IKBBA50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0128, dissociation=1.0152, activation=0, deactivation=0, transport_in=1.656, transport_out=0);
    # 61 : IkBe-RelA:p50
    setParams!((@view rates[IKBEA50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0806, dissociation=0.36, activation=0, deactivation=0, transport_in=8.28, transport_out=0);
    # 62 : IkBd-RelA:p50
    setParams!((@view rates[IKBDA50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.036, dissociation=0.05, activation=0, deactivation=0, transport_in=1.656, transport_out=0);
    # 63 : IkBa-RelA:p50n
    setParams!((@view rates[NIKBAA50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.18, dissociation=0.036, activation=0, deactivation=0, transport_in=0, transport_out=50.4);
    # 64 : IkBb-RelA:p50n
    setParams!((@view rates[NIKBBA50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0128, dissociation=1.0152, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 65 : IkBe-RelA:p50n
    setParams!((@view rates[NIKBEA50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0806, dissociation=0.36, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 66 : IkBd-RelA:p50n
    setParams!((@view rates[NIKBDA50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.036, dissociation=0.05, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 67 : IkBa-RelA:p52
    setParams!((@view rates[IKBAA52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.1804, dissociation=0.036, activation=0, deactivation=0, transport_in=16.56, transport_out=0);
    # 68 : IkBb-RelA:p52
    setParams!((@view rates[IKBBA52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0128, dissociation=1.0152, activation=0, deactivation=0, transport_in=16.56, transport_out=0);
    # 69 : IkBe-RelA:p52
    setParams!((@view rates[IKBEA52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0806, dissociation=0.36, activation=0, deactivation=0, transport_in=8.28, transport_out=0);
    # 70 : IkBd-RelA:p52
    setParams!((@view rates[IKBDA52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.036, dissociation=0.0504, activation=0, deactivation=0, transport_in=1.656, transport_out=0);
    # 71 : IkBa-RelA:p52n
    setParams!((@view rates[NIKBAA52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.1804, dissociation=0.036, activation=0, deactivation=0, transport_in=0, transport_out=50.4);
    # 72 : IkBb-RelA:p52n
    setParams!((@view rates[NIKBBA52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0128, dissociation=1.0152, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 73 : IkBe-RelA:p52n
    setParams!((@view rates[NIKBEA52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0806, dissociation=0.36, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 74 : IkBd-RelA:p52n
    setParams!((@view rates[NIKBDA52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.036, dissociation=0.0504, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 75 : IkBa-RelB:p50
    setParams!((@view rates[IKBAB50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.1804, dissociation=0.036, activation=0, deactivation=0, transport_in=16.56, transport_out=0);
    # 76 : IkBb-RelB:p50
    setParams!((@view rates[IKBBB50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0, dissociation=0.0, activation=0, deactivation=0, transport_in=0.0, transport_out=0);
    # 77 : IkBe-RelB:p50
    setParams!((@view rates[IKBEB50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0806, dissociation=0.36, activation=0, deactivation=0, transport_in=8.28, transport_out=0);
    # 78 : IkBd-RelB:p50
    setParams!((@view rates[IKBDB50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0036, dissociation=0.0504, activation=0, deactivation=0, transport_in=1.656, transport_out=0);
    # 79 : IkBa-RelB:p50n
    setParams!((@view rates[NIKBAB50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.1804, dissociation=0.036, activation=0, deactivation=0, transport_in=0, transport_out=50.4);
    # 80 : IkBb-RelB:p50n
    setParams!((@view rates[NIKBBB50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0, dissociation=0.0, activation=0, deactivation=0, transport_in=0, transport_out=0.0);
    # 81 : IkBe-RelB:p50n
    setParams!((@view rates[NIKBEB50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0806, dissociation=0.36, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 82 : IkBd-RelB:p50n
    setParams!((@view rates[NIKBDB50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0036, dissociation=0.0504, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 83 : IkBa-RelB:p52
    setParams!((@view rates[IKBAB52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0, dissociation=0.0, activation=0, deactivation=0, transport_in=0.0, transport_out=0);
    # 84 : IkBb-RelB:p52
    setParams!((@view rates[IKBBB52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0, dissociation=0.0, activation=0, deactivation=0, transport_in=0.0, transport_out=0);
    # 85 : IkBe-RelB:p52
    setParams!((@view rates[IKBEB52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0, dissociation=0.0, activation=0, deactivation=0, transport_in=0.0, transport_out=0);
    # 86 : IkBd-RelB:p52
    setParams!((@view rates[IKBDB52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.036, dissociation=0.05, activation=0, deactivation=0, transport_in=1.656, transport_out=0);
    # 87 : IkBa-RelB:p52n
    setParams!((@view rates[NIKBAB52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0, dissociation=0.0, activation=0, deactivation=0, transport_in=0, transport_out=0.0);
    # 88 : IkBb-RelB:p52n
    setParams!((@view rates[NIKBBB52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0, dissociation=0.0, activation=0, deactivation=0, transport_in=0, transport_out=0.0);
    # 89 : IkBe-RelB:p52n
    setParams!((@view rates[NIKBEB52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0, dissociation=0.0, activation=0, deactivation=0, transport_in=0, transport_out=0.0);
    # 90 : IkBd-RelB:p52n
    setParams!((@view rates[NIKBDB52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.036, dissociation=0.05, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 91 : IkBa-cRel:p50
    setParams!((@view rates[IKBAC50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.1804, dissociation=0.288, activation=0, deactivation=0, transport_in=16.56, transport_out=0);
    # 92 : IkBb-cRel:p50
    setParams!((@view rates[IKBBC50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0128, dissociation=1.0152, activation=0, deactivation=0, transport_in=1.656, transport_out=0);
    # 93 : IkBe-cRel:p50
    setParams!((@view rates[IKBEC50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0806, dissociation=0.0016, activation=0, deactivation=0, transport_in=8.28, transport_out=0);
    # 94 : IkBd-cRel:p50
    setParams!((@view rates[IKBDC50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.036, dissociation=0.05, activation=0, deactivation=0, transport_in=1.656, transport_out=0);
    # 95 : IkBa-cRel:p50n
    setParams!((@view rates[NIKBAC50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.1804, dissociation=0.288, activation=0, deactivation=0, transport_in=0, transport_out=50.4);
    # 96 : IkBb-cRel:p50n
    setParams!((@view rates[NIKBBC50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0128, dissociation=1.0152, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 97 : IkBe-cRel:p50n
    setParams!((@view rates[NIKBEC50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0806, dissociation=0.0016, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 98 : IkBd-cRel:p50n
    setParams!((@view rates[NIKBDC50, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.036, dissociation=0.05, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 99 : IkBa-cRel:p52
    setParams!((@view rates[IKBAC52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.1804, dissociation=0.288, activation=0, deactivation=0, transport_in=16.56, transport_out=0);
    # 100 : IkBb-cRel:p52
    setParams!((@view rates[IKBBC52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0128, dissociation=1.0152, activation=0, deactivation=0, transport_in=16.56, transport_out=0);
    # 101 : IkBe-cRel:p52
    setParams!((@view rates[IKBEC52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0806, dissociation=0.36, activation=0, deactivation=0, transport_in=8.28, transport_out=0);
    # 102 : IkBd-cRel:p52
    setParams!((@view rates[IKBDC52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.118, dissociation=0.05, activation=0, deactivation=0, transport_in=1.656, transport_out=0);
    # 103 : IkBa-cRel:p52n
    setParams!((@view rates[NIKBAC52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.1804, dissociation=0.288, activation=0, deactivation=0, transport_in=0, transport_out=50.4);
    # 104 : IkBb-cRel:p52n
    setParams!((@view rates[NIKBBC52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0128, dissociation=1.0152, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 105 : IkBe-cRel:p52n
    setParams!((@view rates[NIKBEC52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.0806, dissociation=0.36, activation=0, deactivation=0, transport_in=0, transport_out=25.2);
    # 106 : IkBd-cRel:p52n
    setParams!((@view rates[NIKBDC52, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.0144, decay=0.0144, association=0.118, dissociation=0.05, activation=0, deactivation=0, transport_in=0, transport_out=25.2);




    ##############################
    ###### Apoptosis Model #######
    ##############################
    # MODULE 1: NFkB regulatory species :
    # 1 : tBcl2 ("activation" is the basal fractional activity Kb)
    setParams!((@view rates[TBCL2, 1:end]), initialConc=2000 * TBCL2_CONC, scale=CONVERSION, basalsynthesis=0, synthesis=480, translation=0, basaldecay=0.42 * TBCL2_MOD, decay=BCR_DEATH2, association=0, dissociation=0, activation=0.01, deactivation=0, transport_in=0, transport_out=0);
    #------------------------------------------------
    # MODULE 2: Receptor & DISC system species :
    # 2 : L (10000 in Shokhirev, 1 in Mitchell 2018??)
    setParams!((@view rates[L, 1:end]), initialConc=1, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 3 : R
    setParams!((@view rates[R, 1:end]), initialConc=1000, scale=CONVERSION, basalsynthesis=138.74634, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 4 : L-R
    setParams!((@view rates[LR, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.00144, dissociation=0.0036, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 5 : DISC
    setParams!((@view rates[DISC, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=36, deactivation=0, transport_in=0, transport_out=0);
    # 6 : flip
    setParams!((@view rates[FLIP, 1:end]), initialConc=2000, scale=CONVERSION, basalsynthesis=69.4316196, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 7 : flip:DISC
    setParams!((@view rates[FLIPDISC, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0.0036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    #--------------------------------------------
    # MODULE 3: CASPASE 8 module (initiator caspase) species :
    # 8 : pC8
    setParams!((@view rates[PC8, 1:end]), initialConc=100000, scale=CONVERSION, basalsynthesis=14011.14402, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 9 : DISC:pC8
    setParams!((@view rates[DISCPC8, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.00036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 10 : C8 (tune basalsynthesis=0.0001~0.00036)
    setParams!((@view rates[C8, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=BCR_DEATH, synthesis=3600, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=3600, deactivation=0, transport_in=0, transport_out=0);
    # 11 : Bar
    setParams!((@view rates[BAR, 1:end]), initialConc=1000, scale=CONVERSION, basalsynthesis=1270.296408, synthesis=0, translation=0, basaldecay=2.0794416, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 12 : Bar:C8
    setParams!((@view rates[BARC8, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=6.9314718, decay=0, association=0.0036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    #--------------------------------------------
    # MODULE 4: species before Bax enter mitochondria :
    # 13 : Bid
    setParams!((@view rates[BID, 1:end]), initialConc=60000, scale=CONVERSION, basalsynthesis=28625.67252, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 14 : C8:Bid
    setParams!((@view rates[C8BID, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.00036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 15 : tBid
    setParams!((@view rates[TBID, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=3600, deactivation=0, transport_in=0, transport_out=0);
    # 16 : cBcl2
    setParams!((@view rates[CBCL2, 1:end]), initialConc=20000, scale=CONVERSION, basalsynthesis=14762.72892, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 17 : cBcl2:tBid
    setParams!((@view rates[CBCL2TBID, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0.0036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    #--------------------------------------------
    # MODULE 5: Oligomerization of Bax :
    # 18 : Bax
    setParams!((@view rates[BAX, 1:end]), initialConc=80000, scale=CONVERSION, basalsynthesis=72098.4354, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 19 : tBid:Bax
    setParams!((@view rates[TBIDBAX, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.00036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 20 : aBax
    setParams!((@view rates[ABAX, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=3600, deactivation=0, transport_in=36, transport_out=0);
    # 21 : mBax
    setParams!((@view rates[MBAX, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=3600);
    # 22 : Bcl2
    setParams!((@view rates[BCL2, 1:end]), initialConc=30000 * BCL2_CONC, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=18.42, basaldecay=0.693147, decay=BCR_DEATH3, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 23 : mBax:Bcl2
    setParams!((@view rates[MBAXBCL2, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0.0036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 24 : Bax2
    setParams!((@view rates[BAX2, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.0036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 25 : Bax2:Bcl2
    setParams!((@view rates[BAX2BCL2, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0.0036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 26 : Bax4
    setParams!((@view rates[BAX4, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.0036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 27 : Bax4:Bcl2
    setParams!((@view rates[BAX4BCL2, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0.0036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    #--------------------------------------------
    # MODULE 6: MOMP (pore-forming and transporting) species :
    # 28 : Mito
    setParams!((@view rates[MITO, 1:end]), initialConc=500000, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 29 : Bax4:Mito
    setParams!((@view rates[BAX4MITO, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.0036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 30 : aMito
    setParams!((@view rates[AMITO, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=3600, deactivation=6.9314718, transport_in=0, transport_out=0);
    # 31 : mCytoC
    setParams!((@view rates[MCYTOC, 1:end]), initialConc=500000, scale=CONVERSION, basalsynthesis=346576.2318, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 32 : aMito-mCytoC
    setParams!((@view rates[AMITOMCYTOC, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.0072, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 33 : aCytoC
    setParams!((@view rates[ACYTOC, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=36000, deactivation=0, transport_in=0, transport_out=3600);
    # 34 : mSmac
    setParams!((@view rates[MSMAC, 1:end]), initialConc=100000, scale=CONVERSION, basalsynthesis=4158914.778, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 35 : aMito:mSmac
    setParams!((@view rates[AMITOMSMAC, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.0072, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 36 : aSmac
    setParams!((@view rates[ASMAC, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=36000, deactivation=0, transport_in=0, transport_out=3600);
    #--------------------------------------------
    # MODULE 7: Feedforward pathway 1 (XIAP) species :
    # 37 : XIAP
    setParams!((@view rates[XIAP, 1:end]), initialConc=100000, scale=CONVERSION, basalsynthesis=4159259.826, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 38 : cSmac
    setParams!((@view rates[CSMAC, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=36, transport_out=0);
    # 39 : cSmac:XIAP
    setParams!((@view rates[CSMACXIAP, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0.0252, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    #--------------------------------------------
    # MODULE 8: Feedforward pathway 2 (Apoptosome) species :
    # 40 : cCytoC
    setParams!((@view rates[CCYTOC, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=6.9314718, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=36, transport_out=0);
    # 41 : Apaf
    setParams!((@view rates[APAF, 1:end]), initialConc=100000, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 42 : Apaf:cCytoC
    setParams!((@view rates[APAFCCYTOC, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.0018, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 43 : aApaf
    setParams!((@view rates[AAPAF, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=3600, deactivation=6.9314718, transport_in=0, transport_out=0);
    # 44 : pC9
    setParams!((@view rates[PC9, 1:end]), initialConc=100000, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 45 : Apop
    setParams!((@view rates[APOP, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.00018, dissociation=3.6, activation=0, deactivation=0.693147, transport_in=0, transport_out=0);
    # 46 : Apop:XIAP
    setParams!((@view rates[APOPXIAP, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.0072, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    #--------------------------------------------
    # MODULE 9: CASPASE 3 module (effector caspase) species :
    # 47 : pC3
    setParams!((@view rates[PC3, 1:end]), initialConc=10000, scale=CONVERSION, basalsynthesis=416569.242, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 48 : C8:pC3
    setParams!((@view rates[C8PC3, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.00036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 49 : C3
    setParams!((@view rates[C3, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=3600, translation=0, basaldecay=0, decay=0, association=0, dissociation=0, activation=3600, deactivation=360, transport_in=0, transport_out=0);
    # 50 : XIAP:C3
    setParams!((@view rates[XIAPC3, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.0072, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 51 : Apop:pC3
    setParams!((@view rates[APOPPC3, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.000018, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 52 : uC3
    setParams!((@view rates[UC3, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    #--------------------------------------------
    # MODULE 10: CASPASE 6 feedback module species :
    # 53 : pC6
    setParams!((@view rates[PC6, 1:end]), initialConc=10000, scale=CONVERSION, basalsynthesis=416569.242, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 54 : C3:pC6
    setParams!((@view rates[C3PC6, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.00036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 55 : C6
    setParams!((@view rates[C6, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=6.9314718, decay=0, association=0, dissociation=0, activation=3600, deactivation=0, transport_in=0, transport_out=0);
    # 56 : C6:pC8
    setParams!((@view rates[C6PC8, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.00036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    #--------------------------------------------
    # MODULE 11: Dcell death species :
    # 57 : PARP
    setParams!((@view rates[PARP, 1:end]), initialConc=1000000, scale=CONVERSION, basalsynthesis=694271.958, synthesis=0, translation=0, basaldecay=0.693147, decay=0, association=0, dissociation=0, activation=0, deactivation=72000, transport_in=0, transport_out=0);
    # 58 : C3:PARP
    setParams!((@view rates[C3PARP, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=0, decay=0, association=0.0036, dissociation=3.6, activation=0, deactivation=0, transport_in=0, transport_out=0);
    # 59 : tPARP
    setParams!((@view rates[CPARP, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0, synthesis=0, translation=0, basaldecay=4.1588832, decay=0, association=0, dissociation=0, activation=0, deactivation=0, transport_in=0, transport_out=0);



    ####################################
    ###### Differentiation Model #######
    ####################################

    # Differentiation decision layer :
    # 1 : Pax-5
    setParams!((@view rates[PAX5, 1:end]), initialConc=1.972, scale=CONVERSION, basalsynthesis=2.0, synthesis=7.5, translation=0.0, basaldecay=0.0, decay=1.0, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 2 : Bcl-6
    setParams!((@view rates[BCL6, 1:end]), initialConc=2.3864, scale=CONVERSION, basalsynthesis=0.0, synthesis=7.5, translation=0.0, basaldecay=0.0, decay=1.0, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 3 : Blimp-1
    setParams!((@view rates[BLIMP1, 1:end]), initialConc=0.1191, scale=CONVERSION, basalsynthesis=0.0, synthesis=7.5, translation=0.0, basaldecay=0.0, decay=1.0, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 4 : IRF-4
    setParams!((@view rates[IRF4, 1:end]), initialConc=1.035, scale=CONVERSION, basalsynthesis=1.0, synthesis=7.5, translation=0.0, basaldecay=0.0, decay=1.0, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);



    ##################################
    ###### Proliferation Model #######
    ##################################

    # Proliferation decision layer :
    # 1 : tMyc
    setParams!((@view rates[TMYC, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=10.0 * TMYC_MOD, translation=0.0, basaldecay=0.0, decay=0.231, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0, activation=0.01, deactivation=0);
    # 2 : Myc
    setParams!((@view rates[MYC, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.0, translation=0.4, basaldecay=0.0, decay=0.2, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    #------------------------------------------------
    # Cell cycle layer :
    # 3 : tE2F
    setParams!((@view rates[TE2F, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.15 * TE2F_MOD, synthesis=4.0, translation=0.0, basaldecay=0.0, decay=0.25, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 4 : E2F
    setParams!((@view rates[E2F, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.0, translation=0.4, basaldecay=0.005, decay=1.0, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 5 : Rb
    setParams!((@view rates[RB, 1:end]), initialConc=2.983, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.18, translation=0.0, basaldecay=0.0, decay=0.06 * RB_MOD, association=0.0, dissociation=0.0, phosphorylation=18.0, dephosphorylation=0.0);
    # 6 : ppRB
    setParams!((@view rates[PPRB, 1:end]), initialConc=0.0176, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.0, translation=0.0, basaldecay=0.0, decay=0.06, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=5.0);
    # 7 : E2F:Rb
    setParams!((@view rates[E2FRB, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.0, translation=0.0, basaldecay=0.0, decay=0.03, association=10000.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 8 : CycD
    setParams!((@view rates[CYCD, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.3, translation=1.0, basaldecay=0.0, decay=5.0, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 9 : tCycD
    setParams!((@view rates[TCYCD, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.0, translation=0.0, basaldecay=0.0, decay=0.0, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 10 : CycE
    setParams!((@view rates[CYCE, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.35, translation=0.0, basaldecay=0.0, decay=0.1, association=0.0, dissociation=0.0, phosphorylation=2.0, dephosphorylation=0.0);
    # 11 : CycA
    setParams!((@view rates[CYCA, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.5, translation=0.0, basaldecay=0.0, decay=20.0, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 12 : CycB
    setParams!((@view rates[CYCB, 1:end]), initialConc=0.00508, scale=CONVERSION, basalsynthesis=0.1, synthesis=0.6 * CYCB_MOD, translation=0.0, basaldecay=0.05, decay=1.0, association=0.0, dissociation=0.0, phosphorylation=20.0, dephosphorylation=0.0);
    # 13 : p27
    setParams!((@view rates[P27, 1:end]), initialConc=1.9033, scale=CONVERSION, basalsynthesis=0.0, synthesis=20.0, translation=0.0, basaldecay=10.0, decay=0.0, association=0.0, dissociation=0.0, phosphorylation=100.0, dephosphorylation=0.0);
    # 14 : CycD:p27
    setParams!((@view rates[CYCDP27, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.0, translation=0.0, basaldecay=0.0, decay=0.0, association=1000.0, dissociation=10.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 15 : CycE:p27
    setParams!((@view rates[CYCEP27, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.0, translation=0.0, basaldecay=0.0, decay=0.0, association=1000.0, dissociation=10.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 16 : CycA:p27
    setParams!((@view rates[CYCAP27, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.0, translation=0.0, basaldecay=0.0, decay=0.0, association=1000.0, dissociation=10.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 17 : IEP
    setParams!((@view rates[IEP, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.0, translation=0.0, basaldecay=0.0, decay=0.0, association=0.0, dissociation=0.0, phosphorylation=0.7, dephosphorylation=1.8);
    # 18 : PPX
    setParams!((@view rates[PPX, 1:end]), initialConc=1.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.05, translation=0.0, basaldecay=0.0, decay=0.05, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 19 : Cdc20 (total)
    setParams!((@view rates[CDC20, 1:end]), initialConc=0.0050785, scale=CONVERSION, basalsynthesis=0.0, synthesis=1.5, translation=0.0, basaldecay=0.0, decay=1.5, association=0.0, dissociation=0.0, phosphorylation=5.0, dephosphorylation=0.0);
    # 20 : Cdc20P (active)
    setParams!((@view rates[CDC20P, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.0, translation=0.0, basaldecay=0.0, decay=0.0, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=2.5);
    # 21 : Cdh1
    setParams!((@view rates[CDH1, 1:end]), initialConc=1.0, scale=CONVERSION, basalsynthesis=7.5, synthesis=0.0, translation=0.0, basaldecay=0.0, decay=0.0, association=0.0, dissociation=0.0, phosphorylation=140.0, dephosphorylation=40.0);
    #------------------------------------------------
    # 22 : GM
    setParams!((@view rates[GM, 1:end]), initialConc=1.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.2, translation=0.0, basaldecay=0.0, decay=0.002, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 23 : Mass
    setParams!((@view rates[MASS, 1:end]), initialConc=1.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.02, translation=0.0, basaldecay=0.0, decay=0.02, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 24 : RBGS
    setParams!((@view rates[RBGS, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.0, translation=0.0, basaldecay=0.0, decay=0.0, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
    # 25 : Generation number
    setParams!((@view rates[GEN, 1:end]), initialConc=0.0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.0, translation=0.0, basaldecay=0.0, decay=0.0, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);

    # To keep track of the cell fate decision # CHANGED to 0.0 for time (Td). Need to switch to 0 for cell fate (normal)
    setParams!((@view rates[TOTAL_SPECIES, 1:end]), initialConc=0, scale=CONVERSION, basalsynthesis=0.0, synthesis=0.0, translation=0.0, basaldecay=0.0, decay=0.0, association=0.0, dissociation=0.0, phosphorylation=0.0, dephosphorylation=0.0);
#
#     # !!!!!!!!! SUBJECT TO CHANGE !!!!!!!!!!
#     # NFkB steady state as initial concentration (to avoid long pre-smulation time)
#     (@view rates[IKK:NIKBDC50, INITIALCONC]) .= [
#            0.001,
#        0.14243854850248225,
#        0.20833333333330958,
#        0.08941492524493025,
#        0.01338984282014388,
#       10.049787440438708,
#       19.663912847101077,
#       34.87908600758657,
#       68.2144371842325,
#        3.6856752411058173,
#        0.9091036284893966,
#       47.0975968339208,
#        4.802337947113564,
#        0.02499999999999703,
#        0.028256250100611162,
#        0.028256250100611162,
#        7.779565765495248,
#        4.708673448009479,
#       10.710052018346948,
#        1.0983108720214592,
#        3.5354400990669816,
#        0.1807735098829773,
#        0.013515089220393151,
#        0.46492159616093537,
#        0.016911821399738753,
#        0.2079944855847172,
#        0.22995465600894835,
#        6.384492498073419,
#       10.902549678268038,
#        2.841562523760022,
#        1.2430503691429073,
#       74.06790080834702,
#        0.10073030396476601,
#        0.14496721223916761,
#        0.4095398070034898,
#        4.99665587055513,
#        0.037054499866281905,
#        0.05542665573261573,
#       74.23613016985541,
#        0.16823820481499938,
#       48.39473939795053,
#       43.05372489412052,
#       24.495010485324777,
#        0.013442318072217905,
#       16.607131694555186,
#       19.085317776328274,
#        0.0,
#        0.0,
#        0.0,
#        0.0,
#        0.0,
#        0.0,
#        0.0,
#        0.0,
#        5.08119619519807,
#        0.07512102802002223,
#      245.81630003191387,
#       19.221710041425,
#        1.6963519672197198,
#        0.005998127363757545,
#       81.09861213430489,
#        8.5080603921139
#     ];
end
