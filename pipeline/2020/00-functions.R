get_coef_north <- function(est, subset=1, ...) {
    estn <- est[subset,]

    mu <- drop(Xage %*% estn[colnames(Xage)])
    lam1 <- exp(mu)
    lam1 <- lam1[!grepl("9", names(lam1))]
    lamCC <- lam1[grepl("CC", names(lam1))]

#    MOD <- c("ROAD", "mWell", "mSoft", "mEnSft", "mTrSft", "mSeism")
    MOD <- c("mWell", "mSoft", "mEnSft", "mTrSft", "mSeism")
    Z <- exp(estn[MOD])
    isSoft <- estn["mSoft"] != 0 & estn["mEnSft"] == 0
    #isSoft2 <- get_mid(resn)[,"Contrast"] == 3
    if (isSoft) {
        estn["mEnSft"] <- estn["mSoft"]
        estn["mTrSft"] <- estn["mSoft"]
        estn["mSeism"] <- estn["mSoft"]
    }
    pm <- c("ROAD"=1, "mWell"=0.2, "mSoft"=0.2,
        "mEnSft"=0.2, "mTrSft"=0.2, "mSeism"=0.05)
    pm <- c("ROAD"=1, "mWell"=0.1, "mSoft"=0.1,
        "mEnSft"=0.1, "mTrSft"=0.1, "mSeism"=0.025)
    for (i in MOD)
        Z[i] <- linexp(1, estn[i], pm[i])

    HFc <- c("Crop", "Industrial", "Mine", "RoughP", "Rural", "TameP", "Urban")

    FUN <- mean

    # not in HFc and not forestry!
    Xn2 <- Xn[en$DAT$mWell > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nWell <- nrow(Xn2)
    lamWell <- exp(Xn2 %*% estn[colnames(Xn2)])
    estWell <- FUN(lamWell * Z["mWell"])

    Xn2 <- Xn[en$DAT$mEnSft > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nEnSoft <- nrow(Xn2)
    lamEnSft <- exp(Xn2 %*% estn[colnames(Xn2)])
    estEnSft <- FUN(lamEnSft * Z["mEnSft"])

    ## TrSft incorporates ROAD effect as well?
    Xn2 <- Xn[en$DAT$mTrSft > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nTrSoft <- nrow(Xn2)
    lamTrSft <- exp(Xn2 %*% estn[colnames(Xn2)])
    #estTrSft <- quantile(lamTrSft * Z[,"mTrSft"] * Z[,"ROAD"], c(0.5, 0.05, 0.95))
    estTrSft <- FUN(lamTrSft * Z["mTrSft"])

    Xn2 <- Xn[en$DAT$mSeism > 0 & en$DAT$vegc %ni% HFc & en$DAT$fCC2 == 0, colnames(Xage)]
    nSeism <- nrow(Xn2)
    lamSeism <- exp(Xn2 %*% estn[colnames(Xn2)])
    estSeism <- FUN(lamSeism * Z["mSeism"])


    ## capping values
    ## early seral: young forest (R) and open/shrub
    EARLY <- c(names(lam1)[endsWith(names(lam1), "R")],
        "GrassHerb", "Shrub", "GraminoidFen", "Marsh")
    ESMAX <- max(lam1[EARLY])
    ## non vegetated HF: this should apply to wellsites besides early seral
    NONVEG <-  c("Industrial", "Rural", "Urban")
    NVMAX <- max(lam1[NONVEG], ESMAX)

    ## UrbInd
    #Xn2 <- Xn[en$DAT$vegc %ni% c("Industrial", "Mine", "Rural", "Urban"), colnames(Xage)]
    #nUI <- nrow(Xn2)
    #lamUI <- apply(exp(Xn2 %*% t(estn[,colnames(Xn2),drop=FALSE])), 2, median)
    #estUI <- quantile((nUI * lamUI + nWell * lamWell) / (nUI + nWell),
    #    c(0.5, 0.05, 0.95))

    lam1 <- c(lam1,
        Wellsites=min(NVMAX, estWell), # nonveg HF + early seral
        EnSeismic=min(ESMAX, estSeism),# early seral
        EnSoftLin=min(ESMAX, estEnSft),# early seral
        TrSoftLin=min(ESMAX, estTrSft),# early seral
        HardLin=0,
        Water=0,
        Bare=0,
        SnowIce=0,
        MineV=unname(lam1["Mine"]))
    lam1["Mine"] <- 0
    names(lam1) <- gsub("Spruce", "WhiteSpruce", names(lam1))
    names(lam1) <- gsub("Decid", "Deciduous", names(lam1))
    lam1
}

get_coef_south <- function(est, subset=1, ...) {
    ests <- est[subset,]
    LCC <- c("soilcBlowout", "soilcClaySub", "soilcCrop",
        "soilcIndustrial", "soilcMine", "soilcOther", "soilcRapidDrain",
        "soilcRoughP", "soilcRural", "soilcSandyLoam", "soilcTameP",
        "soilcThinBreak", "soilcUrban")
    LCC2 <- c("soilc2Crop", "soilc2OtherBlowThinRapid",
        "soilc2RoughP", "soilc2Rural", "soilc2TameP", "soilc2UrbInd")
    LCC1 <- c("soilc1Crop", "soilc1RoughP", "soilc1TameP", "soilc1UrbIndRur")
    if (all(ests[LCC] == 0)) {
        ## soilc2 in use
        if (all(ests[LCC1] == 0)) {
            #ests[c("Loamy", "ClaySub", "SandyLoam")] <- ests["ClaySubLoamSand"] # intercept
            ests[c("soilcBlowout", "soilcOther", "soilcRapidDrain",
                "soilcThinBreak")] <- ests["soilc2OtherBlowThinRapid"]
            ests[c("soilcIndustrial", "soilcMine", "soilcUrban")] <- ests["soilc2UrbInd"]
            ests[c("soilcCrop", "soilcRoughP", "soilcRural", "soilcTameP")] <-
              ests[c("soilc2Crop", "soilc2RoughP", "soilc2Rural", "soilc2TameP")]
        } else {
        ## soilc1 in use
#            ests[c("soilcBlowout", "soilcOther", "soilcRapidDrain", "soilcThinBreak",
#                "soilcLoamy", "soilcClaySub", "soilcSandyLoam")] <- ests["soilc1SoilNative"]
            ests[c("soilcIndustrial", "soilcMine", "soilcUrban",
                "soilcRural")] <- ests["soilc1UrbIndRur"]
            ests[c("soilcCrop", "soilcRoughP", "soulcTameP")] <-
              ests[c("soilc1Crop", "soilc1RoughP", "soulc1TameP")]
        }
    }
    muLCC <- c(ests[1], ests[1]+ests[LCC])
    names(muLCC) <- levels(es$DAT$soilc)
    lamLCC0 <- exp(muLCC) # nontreed
    o <- c("Loamy", "Blowout", "ClaySub", "RapidDrain", "SandyLoam", "ThinBreak", "Other",
        "RoughP", "TameP","Crop","Urban", "Rural", "Industrial", "Mine")
    lamLCC0 <- lamLCC0[o]

#    MOD <- c("ROAD", "mWell", "mSoft")
    MOD <- c("mWell", "mSoft")
    Z <- exp(ests[MOD])
    pm <- c("ROAD"=1, "mWell"=0.2, "mSoft"=0.2)
    for (i in MOD)
        Z[i] <- linexp(1, ests[i], pm[i])
    lamMOD <- Z
#    names(lamMOD) <- c("Road", "Well", "Soft")
    names(lamMOD) <- c("Well", "Soft")

    HFc <- c("RoughP", "TameP","Crop","Urban", "Rural", "Industrial", "Mine")

    FUN <- mean

    Xs2 <- Xs[es$DAT$mSoft > 0 & es$DAT$soilc %ni% HFc, colnames(Xs)]
    nSoft <- nrow(Xs2)
    lamSoft <- exp(Xs2 %*% ests[colnames(Xs2)])
    estSoft <- FUN(lamSoft * Z["mSoft"])

    Xs2 <- Xs[es$DAT$mWell > 0 & es$DAT$soilc %ni% HFc, colnames(Xs)]
    nWell <- nrow(Xs2)
    lamWell <- exp(Xs2 %*% ests[colnames(Xs2)])
    estWell <- FUN(lamWell * Z["mWell"])

#    Xs2 <- Xs[es$DAT$ROAD > 0 & es$DAT$soilc %ni% HFc, colnames(Xs)]
#    nROAD <- nrow(Xs2)
#    lamROAD <- exp(Xs2 %*% ests[colnames(Xs2)])
#    estROAD <- FUN(lamROAD * Z["ROAD"])

    ## capping values
    ## early seral: young forest (R) and open/shrub
    EARLY <- c("Loamy", "Blowout", "ClaySub", "RapidDrain", "SandyLoam", "ThinBreak", "Other")
    ESMAX <- max(lamLCC0[EARLY])
    ## non vegetated HF: this should apply to wellsites besides early seral
    NONVEG <-  c("Industrial", "Rural", "Urban")
    NVMAX <- max(lamLCC0[NONVEG], ESMAX)

    lam0 <- c(lamLCC0,
        Wellsites=min(NVMAX, estWell),
        EnSeismic=min(ESMAX, estSoft),
        EnSoftLin=min(ESMAX, estSoft),
        TrSoftLin=min(ESMAX, estSoft),
        HardLin=0, # survace is unsuitable as hell
        #HFor=0,
        MineV=unname(lamLCC0["Mine"]),
        Water=0)
    lam0["Mine"] <- 0
    attr(lam0, "pAspen") <- ests["pAspen"]
    lam0
}


ALLBIRDSPP <- list(
  north = c(AlderFlycatcher = "ALFL", AmericanCrow = "AMCR", AmericanGoldfinch = "AMGO",
    AmericanRedstart = "AMRE", AmericanRobin = "AMRO", AmericanThreetoedWoodpecker = "ATTW",
    BaltimoreOriole = "BAOR", BarnSwallow = "BARS", BlackAndWhiteWarbler = "BAWW",
    BlackbilledMagpie = "BBMA", BaybreastedWarbler = "BBWA", BlackbackedWoodpecker = "BBWO",
    BlackcappedChickadee = "BCCH", BrownheadedCowbird = "BHCO", BlueheadedVireo = "BHVI",
    BlackburnianWarbler = "BLBW", BlueJay = "BLJA", BlackpollWarbler = "BLPW",
    BorealChickadee = "BOCH", BrewersBlackbird = "BRBL", BrownCreeper = "BRCR",
    BlackthroatedGreenWarbler = "BTNW", CanadaWarbler = "CAWA", ClaycoloredSparrow = "CCSP",
    CedarWaxwing = "CEDW", ChippingSparrow = "CHSP", CapeMayWarbler = "CMWA",
    CommonNighthawk = "CONI", ConnecticutWarbler = "CONW", CommonRaven = "CORA",
    CommonYellowthroat = "COYE", ChestnutsidedWarbler = "CSWA", DarkeyedJunco = "DEJU",
    DownyWoodpecker = "DOWO", EasternKingbird = "EAKI", EasternPhoebe = "EAPH",
    EuropeanStarling = "EUST", EveningGrosbeak = "EVGR", FoxSparrow = "FOSP",
    GoldencrownedKinglet = "GCKI", GrayJay = "GRAJ", GrayCatbird = "GRCA",
    GreaterYellowlegs = "GRYE", HairyWoodpecker = "HAWO", HermitThrush = "HETH",
    HornedLark = "HOLA", HouseSparrow = "HOSP", HouseWren = "HOWR",
    Killdeer = "KILL", LeContesSparrow = "LCSP", LeastFlycatcher = "LEFL",
    LesserYellowlegs = "LEYE", LincolnsSparrow = "LISP", MagnoliaWarbler = "MAWA",
    MarshWren = "MAWR", MountainChickadee = "MOCH", MourningDove = "MODO",
    MourningWarbler = "MOWA", NelsonsSparrow = "NESP", NorthernFlicker = "NOFL",
    NorthernWaterthrush = "NOWA", OrangecrownedWarbler = "OCWA",
    OlivesidedFlycatcher = "OSFL", Ovenbird = "OVEN", PalmWarbler = "PAWA",
    PhiladelphiaVireo = "PHVI", PineSiskin = "PISI", PileatedWoodpecker = "PIWO",
    PurpleFinch = "PUFI", RosebreastedGrosbeak = "RBGR", RedbreastedNuthatch = "RBNU",
    RubycrownedKinglet = "RCKI", RedeyedVireo = "REVI", RockPigeon = "ROPI",
    RustyBlackbird = "RUBL", RuffedGrouse = "RUGR", RedwingedBlackbird = "RWBL",
    SavannahSparrow = "SAVS", SedgeWren = "SEWR", Sora = "SORA",
    SolitarySandpiper = "SOSA", SongSparrow = "SOSP", SpottedSandpiper = "SPSA",
    SwampSparrow = "SWSP", SwainsonsThrush = "SWTH", TennesseeWarbler = "TEWA",
    TreeSwallow = "TRES", VariedThrush = "VATH", Veery = "VEER",
    VesperSparrow = "VESP", WarblingVireo = "WAVI", WhitebreastedNuthatch = "WBNU",
    WhitecrownedSparrow = "WCSP", WesternMeadowlark = "WEME", WesternTanager = "WETA",
    WesternWoodPewee = "WEWP", WilsonsSnipe = "WISN", WilsonsWarbler = "WIWA",
    WinterWren = "WIWR", WhitethroatedSparrow = "WTSP", WhitewingedCrossbill = "WWCR",
    YellowbelliedFlycatcher = "YBFL", YellowbelliedSapsucker = "YBSA",
    YellowWarbler = "YEWA", YellowheadedBlackbird = "YHBL", YellowrumpedWarbler = "YRWA"),
  south = c(AlderFlycatcher = "ALFL", AmericanCrow = "AMCR", AmericanGoldfinch = "AMGO",
    AmericanRobin = "AMRO", BairdsSparrow = "BAIS", BankSwallow = "BANS",
    BaltimoreOriole = "BAOR", BarnSwallow = "BARS", BlackbilledMagpie = "BBMA",
    BrownheadedCowbird = "BHCO", Bobolink = "BOBO", BrewersBlackbird = "BRBL",
    BrewersSparrow = "BRSP", BrownThrasher = "BRTH", ChestnutcollaredLongspur = "CCLO",
    ClaycoloredSparrow = "CCSP", CedarWaxwing = "CEDW", ChippingSparrow = "CHSP",
    CommonGrackle = "COGR", CommonRaven = "CORA", CommonYellowthroat = "COYE",
    DownyWoodpecker = "DOWO", EasternKingbird = "EAKI", EasternPhoebe = "EAPH",
    EuropeanStarling = "EUST", GrayPartridge = "GRAP", GrasshopperSparrow = "GRSP",
    HornedLark = "HOLA", HouseSparrow = "HOSP", HouseWren = "HOWR",
    Killdeer = "KILL", LarkBunting = "LARB", LarkSparrow = "LASP",
    LongbilledCurlew = "LBCU", LeContesSparrow = "LCSP", LeastFlycatcher = "LEFL",
    LincolnsSparrow = "LISP", MarbledGodwit = "MAGO", McCownsLongspur = "MCLO",
    MourningDove = "MODO", NelsonsSparrow = "NESP", NorthernFlicker = "NOFL",
    RingneckedPheasant = "RNEP", RockPigeon = "ROPI", RedwingedBlackbird = "RWBL",
    SavannahSparrow = "SAVS", Sora = "SORA", SongSparrow = "SOSP",
    SpraguesPipit = "SPPI", SpottedSandpiper = "SPSA", SharptailedGrouse = "STGR",
    TreeSwallow = "TRES", UplandSandpiper = "UPSA", VesperSparrow = "VESP",
    WarblingVireo = "WAVI", WesternKingbird = "WEKI", WesternMeadowlark = "WEME",
    WesternWoodPewee = "WEWP", Willet = "WILL", WilsonsPhalarope = "WIPH",
    WilsonsSnipe = "WISN", WhitethroatedSparrow = "WTSP", YellowWarbler = "YEWA",
    YellowheadedBlackbird = "YHBL"))
