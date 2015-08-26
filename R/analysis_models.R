modsVeg <- list(
    "Hab"=list( # 1
        .~. + hab1,
        .~. + hab1b),
    "Age"=list( # 2
        .~. + wtAge,
        .~. + wtAge + wtAge2,
        ## here dec+mix is both as intercept
        .~. + wtAge + wtAge2 + wtAge:isCon + wtAge2:isCon,
        .~. + wtAge + wtAge2 + wtAge:isUpCon + wtAge:isBSLarch +
            wtAge2:isUpCon + wtAge2:isBSLarch,
        .~. + wtAge + wtAge2 + wtAge:isMix + wtAge:isPine + wtAge:isWSpruce + wtAge:isBSLarch +
            wtAge2:isMix + wtAge2:isPine + wtAge2:isWSpruce + wtAge2:isBSLarch,
        .~. + wtAge05,
        .~. + wtAge05 + wtAge05:isCon,
        .~. + wtAge05 + wtAge05:isUpCon + wtAge05:isBSLarch,
        .~. + wtAge05 + wtAge05:isMix + wtAge05:isPine + wtAge05:isWSpruce + wtAge05:isBSLarch,
        .~. + wtAge05 + wtAge,
        .~. + wtAge05 + wtAge + wtAge05:isCon + wtAge:isCon,
        .~. + wtAge05 + wtAge + wtAge05:isUpCon + wtAge05:isBSLarch +
            wtAge:isUpCon + wtAge:isBSLarch,
        .~. + wtAge05 + wtAge + wtAge05:isMix + wtAge05:isPine + wtAge05:isWSpruce + wtAge05:isBSLarch +
            wtAge:isMix + wtAge:isPine + wtAge:isWSpruce + wtAge:isBSLarch),
    "CC"=list( # 3
        .~. + fCC2),
    "Contrast"=list( # 4
        .~. + ROAD01,
        .~. + SoftLin_PC,
        .~. + ROAD01 + SoftLin_PC,
        .~. + ROAD01 + ROAD01:hab_lcc2,
        .~. + ROAD01 + SoftLin_PC + ROAD01:hab_lcc2,
        .~. + ROAD01 + ROAD01:hab_lcc3,
        .~. + ROAD01 + SoftLin_PC + ROAD01:hab_lcc3,
        .~. + ROAD01 + ROAD01:hab_lcc, # 5 levels
        .~. + ROAD01 + SoftLin_PC + ROAD01:hab_lcc),
    "ARU"=list( # 5
        .~. + ARU),
    "Space"=list( # 6
        .~.+ xPET,## climate only
        .~.+ xMAT,
        .~.+ xAHM,
        .~.+ xFFP,
        .~.+ xMAP + xFFP,
        .~.+ xMAP + xFFP + xMAP:xFFP,
        .~.+ xMAT + xMAP + xPET + xAHM,
        .~.+ xMAT + xMAP + xPET + xAHM + xPET:xMAP + xMAT:xAHM,
        .~.+ xMAT + xMAP,
        .~.+ xMWMT + xMCMT,
        .~.+ xAHM + xPET,
        .~.+ xlat + xlong + xlat:xlong + xPET,## linear lat-long and climate
        .~.+ xlat + xlong + xlat:xlong + xMAT,
        .~.+ xlat + xlong + xlat:xlong + xAHM,
        .~.+ xlat + xlong + xlat:xlong + xFFP,
        .~.+ xlat + xlong + xlat:xlong + xMAP + xFFP,
        .~.+ xlat + xlong + xlat:xlong + xMAP + xFFP + xMAP:xFFP,
        .~.+ xlat + xlong + xlat:xlong + xMAT + xMAP + xPET + xAHM,
        .~.+ xlat + xlong + xlat:xlong + xMAT + xMAP + xPET + xAHM + xPET:xMAP + xMAT:xAHM,
        .~.+ xlat + xlong + xlat:xlong + xMAT + xMAP,
        .~.+ xlat + xlong + xlat:xlong + xMWMT + xMCMT,
        .~.+ xlat + xlong + xlat:xlong + xAHM + xPET,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xPET,## quadratic lat-long and climate
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMAT,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xAHM,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xFFP,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMAP + xFFP,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMAP + xFFP + xMAP:xFFP,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMAT + xMAP + xPET + xAHM,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMAT + xMAP + xPET + xAHM + xPET:xMAP + xMAT:xAHM,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMAT + xMAP,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMWMT + xMCMT,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xAHM + xPET),
    "HF"=list( # 7
        .~.+ THF_KM,## HF only
        .~.+ Lin_KM + Nonlin_KM,
        .~.+ Succ_KM + Alien_KM,
        .~.+ Succ_KM + Noncult_KM + Cult_KM,
        .~.+ THF_KM + THF2_KM,
        .~.+ Lin_KM + Nonlin_KM + Nonlin2_KM,
        .~.+ Succ_KM + Alien_KM + Succ2_KM,
        .~.+ Succ_KM + Noncult_KM + Cult_KM + Succ2_KM,
        .~.+ Succ_KM + Alien_KM + Alien2_KM,
        .~.+ Succ_KM + Noncult_KM + Cult_KM + Noncult2_KM,
        .~.+ Succ_KM + Alien_KM + Succ2_KM + Alien2_KM,
        .~.+ Succ_KM + Noncult_KM + Cult_KM + Succ2_KM + Noncult2_KM),
    "Year"=list( # 8
        .~.+ YR))

modsSoil <- list(
    "Hab"=list( # 1
        .~. + soil1,
        .~. + soil1 + pAspen,
        .~. + soil1v,
        .~. + soil1v + pAspen),
    "Contrast"=list( # 2
        .~. + ROAD01,
        .~. + ROAD01 + ROAD01:hab_lcc2),
    "Space"=list( # 3
        .~.+ xPET,## climate only
        .~.+ xMAT,
        .~.+ xAHM,
        .~.+ xFFP,
        .~.+ xMAP + xFFP,
        .~.+ xMAP + xFFP + xMAP:xFFP,
        .~.+ xMAT + xMAP + xPET + xAHM,
        .~.+ xMAT + xMAP + xPET + xAHM + xPET:xMAP + xMAT:xAHM,
        .~.+ xMAT + xMAP,
        .~.+ xMWMT + xMCMT,
        .~.+ xAHM + xPET,
        .~.+ xlat + xlong + xlat:xlong + xPET,## linear lat-long and climate
        .~.+ xlat + xlong + xlat:xlong + xMAT,
        .~.+ xlat + xlong + xlat:xlong + xAHM,
        .~.+ xlat + xlong + xlat:xlong + xFFP,
        .~.+ xlat + xlong + xlat:xlong + xMAP + xFFP,
        .~.+ xlat + xlong + xlat:xlong + xMAP + xFFP + xMAP:xFFP,
        .~.+ xlat + xlong + xlat:xlong + xMAT + xMAP + xPET + xAHM,
        .~.+ xlat + xlong + xlat:xlong + xMAT + xMAP + xPET + xAHM + xPET:xMAP + xMAT:xAHM,
        .~.+ xlat + xlong + xlat:xlong + xMAT + xMAP,
        .~.+ xlat + xlong + xlat:xlong + xMWMT + xMCMT,
        .~.+ xlat + xlong + xlat:xlong + xAHM + xPET,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xPET,## quadratic lat-long and climate
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMAT,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xAHM,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xFFP,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMAP + xFFP,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMAP + xFFP + xMAP:xFFP,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMAT + xMAP + xPET + xAHM,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMAT + xMAP + xPET + xAHM + xPET:xMAP + xMAT:xAHM,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMAT + xMAP,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xMWMT + xMCMT,
        .~.+ xlat + xlong + xlat:xlong + xlat2 + xlong2 + xAHM + xPET),
    "HF"=list( # 4
        .~.+ THF_KM,## HF only
        .~.+ Lin_KM + Nonlin_KM,
        .~.+ Succ_KM + Alien_KM,
        .~.+ Succ_KM + Noncult_KM + Cult_KM,
        .~.+ THF_KM + THF2_KM,
        .~.+ Lin_KM + Nonlin_KM + Nonlin2_KM,
        .~.+ Succ_KM + Alien_KM + Succ2_KM,
        .~.+ Succ_KM + Noncult_KM + Cult_KM + Succ2_KM,
        .~.+ Succ_KM + Alien_KM + Alien2_KM,
        .~.+ Succ_KM + Noncult_KM + Cult_KM + Noncult2_KM,
        .~.+ Succ_KM + Alien_KM + Succ2_KM + Alien2_KM,
        .~.+ Succ_KM + Noncult_KM + Cult_KM + Succ2_KM + Noncult2_KM),
    "Year"=list( # 5
        .~.+ YR))