# Settings

1. number of nodes to use (12 cores per node),
2. south/north

## Loading git on jasper

`module load application/git/1.7.10.1`

## Loading R on jasper

`module load application/R/3.1.2`

## script to update the runs

```
rm ~/abmi/*
module load application/git/1.7.10.1
cd ~/repos/abmianalytics/
git pull
cd ~/abmi/
cp ~/repos/abmianalytics/wg/* ~/abmi/
```
