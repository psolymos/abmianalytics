:: switch off commad echo
@echo off
:: output folder
set out=y:\Oracle_access_2015\out\

:: test
:: echo Running test (hit OK to read/write exception error if happens)
:: R CMD BATCH --vanilla "00test.R" "%out%OUT_test.log"

:: MasterList
:: echo Compiling master list for sites
:: R CMD BATCH --vanilla "00masterlist.R" "%out%00masterlist.log"

:: Birds
echo Compiling species data for terrestrial birds
R CMD BATCH --vanilla "birds.R" "%out%OUT_birds.log"

:: VPlants
echo Compiling species data for terrestrial vascular plants
R CMD BATCH --vanilla "vplants.R" "%out%OUT_vplants.log"

:: Mites
echo Compiling species data for terrestrial mites
R CMD BATCH --vanilla "mites.R" "%out%OUT_mites.log"

:: Mosses
echo Compiling species data for terrestrial mosses
R CMD BATCH --vanilla "mosses.R" "%out%OUT_mosses.log"

:: Lichens
echo Compiling species data for terrestrial lichens
R CMD BATCH --vanilla "lichens.R" "%out%OUT_lichens.log"

:: Mammals
:: echo Compiling species data for terrestrial mammals
:: R CMD BATCH --vanilla "c:\Users\Peter\repos\abmianalytics\species\R\mammals.R" "%out%OUT_mammals.log"

:: AqVPlants
::echo Compiling species data for aquatic vascular plants
::R CMD BATCH --vanilla "c:\Users\Peter\repos\abmianalytics\species\R\aqplants.R" "%out%OUT_aqplants.log"

:: AqInverts
::echo Compiling species data for aquatic invertebrates
::R CMD BATCH --vanilla "c:\Users\Peter\repos\abmianalytics\species\R\aqinverts.R" "%out%OUT_aqinverts.log"

:: Summary
::echo Running summaries
::R CMD BATCH --vanilla "c:\Users\Peter\repos\abmianalytics\species\R\00summary.R" "%out%00summary.log"

:: erase output folder variable
set out=
:: Exit
exit
