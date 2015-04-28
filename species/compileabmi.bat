:: switch off commad echo
@echo off

:: test
:: echo Running test (hit OK to read/write exception error if happens)
:: R CMD BATCH --vanilla "y:\Oracle_access_2015\R\00test.R" "y:\Oracle_access_2015\out\OUT_test.log"

:: MasterList
:: echo Compiling master list for sites
:: R CMD BATCH --vanilla "y:\Oracle_access_2015\R\00masterlist.R" "y:\Oracle_access_2015\out\00masterlist.log"

:: Birds
echo Compiling species data for terrestrial birds
R CMD BATCH --vanilla "y:\Oracle_access_2015\R\birds.R" "y:\Oracle_access_2015\out\OUT_birds.log"

:: VPlants
echo Compiling species data for terrestrial vascular plants
R CMD BATCH --vanilla "y:\Oracle_access_2015\R\vplants.R" "y:\Oracle_access_2015\out\OUT_vplants.log"

:: Mites
echo Compiling species data for terrestrial mites
R CMD BATCH --vanilla "y:\Oracle_access_2015\R\mites.R" "y:\Oracle_access_2015\out\OUT_mites.log"

:: Mosses
echo Compiling species data for terrestrial mosses
R CMD BATCH --vanilla "y:\Oracle_access_2015\R\mosses.R" "y:\Oracle_access_2015\out\OUT_mosses.log"

:: Lichens
echo Compiling species data for terrestrial lichens
R CMD BATCH --vanilla "y:\Oracle_access_2015\R\lichens.R" "y:\Oracle_access_2015\out\OUT_lichens.log"

:: Mammals
:: echo Compiling species data for terrestrial mammals
:: R CMD BATCH --vanilla "y:\Oracle_access_2015\R\mammals.R" "y:\Oracle_access_2015\out\OUT_mammals.log"

:: AqVPlants
::echo Compiling species data for aquatic vascular plants
::R CMD BATCH --vanilla "y:\Oracle_access_2015\R\aqplants.R" "y:\Oracle_access_2015\out\OUT_aqplants.log"

:: AqInverts
::echo Compiling species data for aquatic invertebrates
::R CMD BATCH --vanilla "y:\Oracle_access_2015\R\aqinverts.R" "y:\Oracle_access_2015\out\OUT_aqinverts.log"

:: Summary
::echo Running summaries
::R CMD BATCH --vanilla "y:\Oracle_access_2015\R\00summary.R" "y:\Oracle_access_2015\out\00summary.log"

:: Exit
exit
