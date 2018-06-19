library(tidyverse)

# Derive generalised habitat type

# ====================================================================
# LUT

path_lut = "LUT"

# Peter's LUT
lut <- read_csv(file.path(path_lut,"lookup-veg-v6_forMetaData.csv"))
lut <- lut %>%
  rename(PreBackfill_Source=preBackfill_Source)

# Daiyuan's LUT for Cutblock not in upland
DomInEachNReg <- read_csv(file.path(path_lut,"DomInEachNReg.csv"))
upland = c('AlpineLarch', 'Decid', 'Fir', 'Mixedwood', 'Pine', 'Spruce') 

Harvest_Area <-  c("CUTBLOCK","HARVEST-AREA")

# ====================================================================
# Input and output file definition

path = "S2017_SummaryTables_Final"
pathout = str_c(path,"_habitat_SC")
dir.create(pathout)
fin =""

fout <- file.path(pathout,str_c("GHT_",fin))

# read input file
df <- read_csv(file.path( path, fin),
               col_types = cols(
  .default = col_character(),
  Pct_of_Larch = col_integer(),
  Shape_Length = col_double(),
  Shape_Area = col_double()
))

# ====================================================================
# join all tables


df2 <- df %>% 
  left_join(lut,c("Veg_Type","Moisture_Reg","PreBackfill_Source"))%>%
  left_join(DomInEachNReg,c("NSRNAME"="NSRName"))

df2 <- df2 %>%
  mutate( 
    Larch = Pct_of_Larch > 0,
    CWCS_Class2 = CWCS_Class,
    CWCS_Class = if_else(
      is.na(CWCS_Class),
      "Not Wetland",
      CWCS_Class
      ),
    Combined2 = ""
    )

#=============================================================================
# Rule 1
#============

df2 <- df2 %>% mutate(
  Combined2 = if_else(
    Larch & Combined == "TreedBog-BSpr",
    "TreedFen-Larch",
    Combined2
    # Combined
  )
)

#=============================================================================
# Rule 2
#============

cc = "TreedWetland-Mixedwood"

df2$Combined2[df2$Combined==cc & df2$CWCS_Class== "Fen"] <- "TreedFen-Mixedwood"

df2$Combined2[df2$Combined==cc & df2$CWCS_Class== "Bog"] <- "TreedBog-BSpr"

df2$Combined2[df2$Combined==cc & df2$CWCS_Class== "Swamp"] <- "TreedSwamp-Mixedwood"

df2$Combined2[df2$Combined==cc & !(df2$CWCS_Class %in%c("Fen","Bog","Swamp"))] <- "TreedFen-Mixedwood"

zo <- df2 %>% filter(Combined == cc)
#--------------------------------------------------------------------------------

#=============================================================================
# Rule 3
#============

cc = "GraminoidWetland"

df2$Combined2[df2$Combined==cc & df2$CWCS_Class== "Bog"] <- "ShrubbyBog"

df2$Combined2[df2$Combined==cc & df2$CWCS_Class== "Marsh"] <- "Marsh"

df2$Combined2[df2$Combined==cc & df2$CWCS_Class== "Fen"] <- "GraminoidFen"

df2$Combined2[df2$Combined==cc & !(df2$CWCS_Class %in% c("Fen","Bog","Marsh"))] <- "Marsh"

zo <- df2 %>% filter(Combined == cc)
#=============================================================================
# Rule 4
#============
cc = "ShrubbyWetland"

df2$Combined2[df2$Combined==cc & df2$CWCS_Class== "Fen"] <- "ShrubbyFen"

df2$Combined2[df2$Combined==cc & df2$CWCS_Class== "Bog"] <- "ShrubbyBog"

df2$Combined2[df2$Combined==cc & !(df2$CWCS_Class %in% c("Fen","Bog"))] <- "ShrubbySwamp"

zo <- df2 %>% filter(Combined == cc)
#=============================================================================
# Rule 5
#============
cc = "Muskeg"

df2$Combined2[df2$Combined==cc & df2$CWCS_Class== "Fen"] <- "GraminoidFen"

df2$Combined2[ df2$Combined==cc &  !(df2$CWCS_Class %in% c("Swamp","Fen","Bog")) & df2$PostBackfill_Source != "WBNP"] <-  "GraminoidFen"

df2$Combined2[df2$Combined==cc & df2$CWCS_Class== "Fen" & df2$PostBackfill_Source == "Phase1" ] <- "TreedFen-BSpr"

df2$Combined2[df2$Combined==cc & df2$CWCS_Class != "Fen" & df2$PostBackfill_Source == "Phase1" ] <- "TreedBog-BSpr"

df2$Combined2[df2$Combined==cc & df2$Combined2=="" ] <- "TreedFen-BSpr"

zo <- df2 %>% filter(Combined == cc)
#=============================================================================
# No more rule
#============

selection <- df2$Combined2==""
df2$Combined2[ selection ] <- df2$Combined[ selection ]

#=============================================================================
# CUtblock not on upland additional rule
#============

selection <- df2$FEATURE_TY %in% Harvest_Area &  !(df2$Combined2 %in% upland)
df2$Combined2[selection] <- df2$UplandGene[selection]

#=============================================================================
# cleaning and saving

df2$CWCS_Class <- df2$CWCS_Class2

df3 <- df2 %>%
  select(-Combined:-CWCS_Class2) %>%
  rename(Combined_ChgByCWCS = Combined2)
  
write_csv(df3,fout)
