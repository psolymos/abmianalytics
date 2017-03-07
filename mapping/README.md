# &species&

## Description

Predicted relative abundance for &species&
in 1 km^2 prediction grid in Alberta, Canada.

## Version

Alberta Biodiversity Monitoring Institute, species website
Version &version&, http://species.abmi.ca

## Fields in the CSV file:

* Row_Col: ID number defining raster row and column IDs (see Projection)
* Ref: normalized relative abundance under reference landscape conditions
* Curr: normalized relative abundance under current landscape conditions

Normalized relative abundance is calculated from original values
(Ref0, Curr0) as:
* Ref = round(Ref0 / max(Ref0, Curr0))
* Curr = round(Curr0 / max(Ref0, Curr0))

## Projection

The Row_Col field links the raster cells to these products:

* Geodatabase:
  http://ftp.public.abmi.ca/species.abmi.ca/gis/Grid1km_working.gdb.zip
* CSV with latitude/longitude (NAD_1983_10TM_AEP_Forest projection):
  http://ftp.public.abmi.ca/species.abmi.ca/gis/Grid1km_working.csv.zip

CRS: '+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
