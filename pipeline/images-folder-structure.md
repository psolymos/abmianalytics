# /images folder structure

Subfolders refer to the groups of species (i.e. taxa,  such as birds,  
mosses etc.).

Within subfolder,  species are organized as 
`<SpeciesID>-<type>.png`. <SpeciesID> is the unique ID for thst species 
that can be found in lookup table. <type> can take the following values:

All species with at least 3 detections:

* `useavail-north`: use-availability based on north data.
* `useavail-south`: use-availability based on south data.
* `det`: detection map (north and south combined).

All species with models:
	
* `coef-north`: coefficient figures based on north models.
* `coef-south`: coefficient figures based on south models.
* `sector-north`: sector effects in the north.
* `sector-south`: sector effects in the south.
* `map`: reference,  current,  difference maps (north and south combined).

Optionally might include:
	
* `unc`: uncertainty map.
* `spclim-north`: climate part of north models.
* `spclim-south`: climate part of south models.
