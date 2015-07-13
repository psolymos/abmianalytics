:: switch off commad echo
@echo off
:: output folder
set out=y:\Oracle_access_2015\out\
:: switch on/off
set TEST=0
set MASTERLIST=0
set BIRDS=0
set VPLANTS=1
set MITES=0
set MOSSES=0
set LICHENS=0
set MAMMALS=0
set AQPLANTS=1
set AQINVERTS=0
set SUMMARY=0

:TEST_SUB
IF %TEST%==1 (
echo Running test (hit OK to read/write exception error if happens)
R CMD BATCH --vanilla "00test.R" "%out%OUT_test.log"
)

:MASTERLIST_SUB
IF %MASTERLIST%==1 (
echo Compiling master list for sites
R CMD BATCH --vanilla "00masterlist.R" "%out%00masterlist.log"
)

:BIRDS_SUB
IF %BIRDS%==1 (
echo Compiling species data for terrestrial birds
R CMD BATCH --vanilla "birds.R" "%out%OUT_birds.log"
)

:VPLANTS_SUB
IF %VPLANTS%==1 (
echo Compiling species data for terrestrial vascular plants
R CMD BATCH --vanilla "vplants.R" "%out%OUT_vplants.log"
)

:MITES_SUB
IF %MITES%==1 (
echo Compiling species data for terrestrial mites
R CMD BATCH --vanilla "mites.R" "%out%OUT_mites.log"
)

:MOSSES_SUB
IF %MOSSES%==1 (
echo Compiling species data for terrestrial mosses
R CMD BATCH --vanilla "mosses.R" "%out%OUT_mosses.log"
)

:LICHENS_SUB
IF %LICHENS%==1 (
echo Compiling species data for terrestrial lichens
R CMD BATCH --vanilla "lichens.R" "%out%OUT_lichens.log"
)

:MAMMALS_SUB
IF %MAMMALS%==1 (
echo Compiling species data for terrestrial mammals
R CMD BATCH --vanilla "mammals.R" "%out%OUT_mammals.log"
)

:AQPLANTS_SUB
IF %AQPLANTS%==1 (
echo Compiling species data for aquatic vascular plants
R CMD BATCH --vanilla "aqplants.R" "%out%OUT_aqplants.log"
)

:AQINVERTS_SUB
IF %AQINVERTS%==1 (
echo Compiling species data for aquatic invertebrates
R CMD BATCH --vanilla "aqinverts.R" "%out%OUT_aqinverts.log"
)

:SUMMARY_SUB
IF %SUMMARY%==1 (
echo Running summaries
R CMD BATCH --vanilla "00summary.R" "%out%00summary.log"
)

:END
:: erase output folder variable
set out=
:: Exit
exit
