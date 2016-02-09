source("~/repos/abmianalytics/veghf/veghf-setup.R")

## compiling lookup tables -------------------------

#tab_veg <- data.frame(Label=colnames(dd1km_pred[[2]]),
#    prop_cr=colSums(dd1km_pred[[2]]) / sum(dd1km_pred[[2]]))
#write.csv(tab_veg, file=file.path(ROOT, VER, "out", "tab_veg.csv"))

ltveg <- read.csv("~/repos/abmianalytics/lookup/lookup-veg.csv")
ltsoil <- read.csv("~/repos/abmianalytics/lookup/lookup-soil.csv")

tveg <- data.frame(VEGHFAGE=colnames(dd1km_pred$veg_current))
tveg$HF <- hfgroups$HF_GROUP[match(tveg$VEGHFAGE, hfgroups$HF_GROUP)]
tveg$HF[substr(as.character(tveg$VEGHFAGE), 1, 2) == "CC"] <- "CutBlocks"
tveg$VEGAGE <- ltveg$VEGAGE[match(tveg$VEGHFAGE, ltveg$VEGAGE)]
tveg$VEGAGE <- as.character(tveg$VEGAGE)
tveg$VEGHFAGE <- as.character(tveg$VEGHFAGE)
tveg$VEGAGE[substr(tveg$VEGHFAGE, 1, 2) == "CC"] <- 
    substr(tveg$VEGHFAGE, 3, nchar(tveg$VEGHFAGE))[substr(tveg$VEGHFAGE, 1, 2) == "CC"]
tveg <- data.frame(tveg, 
    ltveg[match(tveg$VEGAGE, ltveg$VEGAGE),-1],
    hfgroups[match(tveg$HF, hfgroups$HF_GROUP),-1])
tveg$VEGAGE <- as.factor(tveg$VEGAGE)
rownames(tveg) <- tveg$VEGHFAGE
colnames(tveg)[colnames(tveg)=="Type.1"] <- "HFtype"

tsoil <- data.frame(SOILHF=colnames(dd1km_pred$soil_current))
tsoil$HF <- hfgroups$HF_GROUP[match(tsoil$SOILHF, hfgroups$HF_GROUP)]
tsoil$SOIL <- ltsoil$SOILclass[match(tsoil$SOILHF, ltsoil$SOILclass)]
tsoil <- data.frame(tsoil, 
    ltsoil[match(tsoil$SOILHF, ltsoil$SOILclass),-1],
    hfgroups[match(tsoil$SOILHF, hfgroups$HF_GROUP),-1])
rownames(tsoil) <- tsoil$SOILHF

write.csv(tveg, file="~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")
write.csv(tsoil, file="~/repos/abmianalytics/lookup/lookup-soil-hf.csv")

## compiling hab-age tables by NSR ------------------------

lt <- read.csv("~/repos/abmianalytics/lookup/lookup-veg-hf-age.csv")

load(file.path(ROOT, VER, "out/kgrid", "veg-hf_1kmgrid_fix-fire.Rdata"))
load(file.path(ROOT, VER, "out/kgrid", "kgrid_table.Rdata"))

summary(rowSums(dd1km_pred$veg_current)/10^6)
vhf <- groupSums(dd1km_pred$veg_current, 1, kgrid$NSRNAME)
vhf2 <- groupSums(dd1km_pred$veg_current, 1, kgrid$NRNAME)
veg <- groupSums(dd1km_pred$veg_reference, 1, kgrid$NSRNAME)
veg2 <- groupSums(dd1km_pred$veg_reference, 1, kgrid$NRNAME)
stopifnot(all(colnames(vhf) == rownames(lt)))

lt$AGE[lt$AGE == ""] <- NA
labs_to_keep <- is.na(lt$HF) & !is.na(lt$AGE)
cols_to_keep <- colnames(vhf)[labs_to_keep]
known_ages <- cols_to_keep[!grepl("0", cols_to_keep)]

## no CC included here
Target0 <- c("Conif0", "Decid0", "Mixwood0", "Pine0", 
    "Swamp-Conif0", "Swamp-Decid0", "Swamp-Mixwood0", "Swamp-Pine0", 
    "Wetland-BSpr0", "Wetland-Decid0", "Wetland-Larch0")
Target <- gsub("0", "", Target0)
Ages <- c("0", "R", as.character(1:9))

NSRs <- levels(kgrid$NSRNAME)
NRs <- nonDuplicated(kgrid, NSRNAME, TRUE)
NRs <- as.character(NRs[NSRs, "NRNAME"])
names(NRs) <- NSRs

ages_cr <- array(NA, c(length(Target), length(Ages), length(NSRs)),
    list(Target, Ages, NSRs))
ages_rf <- ages_cr

for (What in c("cr","rf")) {
    cat("\n-------------", What, "---------------")
    for (nsr in NSRs) {
        for (i in Target) {
            tmp <- if (What == "cr")
                vhf[nsr, paste0(i, Ages)] else veg[nsr, paste0(i, Ages)]
            if (tmp[1] > 0 && sum(tmp[-1]) <= 0) {## these are all <1%, just ignore
                cat("\nproblem:", nsr, "\t", i, "\t", 
                    round(100 * sum(tmp) / sum(vhf[nsr, cols_to_keep]), 4), 
                    ifelse(sum(vhf[nsr, known_ages]) > 0, "OK", "zero !!!"))
            } else {
                #cat("\nOK:", nsr, "\t", i)
            }
            if (tmp[1] > 0 && sum(tmp[-1]) <= 0) {
                ## we do NOT have known ages to work with
                if (sum(vhf[nsr, known_ages]) <= 0) {
                ## NO know age forest in NSR --> use NR
                    allages <- if (What == "cr")
                        vhf2[NRs[nsr], known_ages] else veg2[NRs[nsr], known_ages]
                } else {
                ## we have known age forest in NSR
                    allages <- if (What == "cr")
                        vhf[nsr, known_ages] else veg[nsr, known_ages]
                }
                dat <- data.frame(A=allages,
                    VEG=substr(names(allages), 1, nchar(names(allages))-1),
                    AGE=substr(names(allages), nchar(names(allages)),
                        nchar(names(allages))))
                mat <- as.matrix(Xtab(A ~ VEG + AGE, dat))
                mat <- mat[Target, Ages[-1]]
                tmp <- colSums(mat)
                tmp <- c("0"=0, tmp)
            } else {
                ## we have known ages to work with
                tmp[1] <- 0
            }
            #stopifnot(all(!is.na(tmp / sum(tmp))))
            tmp <- tmp / ifelse(sum(tmp) <= 0, 1, sum(tmp))
            stopifnot(all(!is.na(tmp)))
            if (What == "cr")
                ages_cr[i,,nsr] <- tmp
            if (What == "rf")
                ages_rf[i,,nsr] <- tmp
        }
    }
    cat("\n\n")
}

sum(is.na(ages_cr))
sum(is.na(ages_rf))

## we need to know availability for combining forest classes
nsr_cr <- groupSums(dd1km_pred$veg_current, 1, kgrid$NSRNAME)
nsr_rf <- groupSums(dd1km_pred$veg_reference, 1, kgrid$NSRNAME)
cn1 <- as.character(lt$VEG)
cn1[!is.na(lt$HF)] <- "HF"
cn2 <- cn1[is.na(lt$HF)]
nsr_cr <- as.matrix(groupSums(nsr_cr, 2, cn1))
nsr_rf <- as.matrix(groupSums(nsr_rf, 2, cn2))
nsr_cr <- nsr_cr[,dimnames(ages_cr)[[1]]]
nsr_rf <- nsr_rf[,dimnames(ages_cr)[[1]]]

AvgAges <- list(current=ages_cr, reference=ages_rf,
    area_cr=nsr_cr, area_rf=nsr_rf)

if (SAVE)
save(AvgAges, 
    file=file.path(ROOT, VER, "out/kgrid", "veg-hf_avgages_fix-fire.Rdata"))

