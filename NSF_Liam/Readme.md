Data from C2VSim
===============

Water levels
-----------

The matlab script **ExtractDataLiam.m** is used to extract the water levels.
The results are written in the **WaterLevelPLSSLiam.dat** in the folowing format:
```
CODE X1 Y1 X2 Y2 H1 H2 ... H505
```
where ```CODE``` is the field CO_MTRS

```X1 Y1``` are the coordinates of the point in EPSG:4326 
```X2 Y2``` are the coordinates of the point in EPSG:26910 (C2Vsim projection)
```H1 H2 ... H505``` are the head monthly values in feet from 9/30/1973 to 9/30/2015.

There are some points with nan values. These are located outside of the C2Vsim mesh.