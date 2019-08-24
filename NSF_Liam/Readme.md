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

Groundwater Pumping
-------
The second code snippet from the  **ExtractDataLiam.m** script is used to extract pumping.
The results are in the **PumpingPLSSLiam.dat** file. THe format of the file is identical to the previous except that the units are ft^3 and the starting date is 10/31/1973 

Agricultural and Urban Deliveries
--------------------
These deliveries in C2VSim are split into 4 categories
 - Non Ponded crops
 - Refuge
 - Rice
 - Urban

 For each category there are three files with suffix _Area.dat_, _Deliv.dat_ and _Pump.dat_ Which correspond to the area (sq ft) that the category is covering, the delivery and Pumping amounts (cu ft). 

**Important Note:**
For the groundwater pumping, Agricultural and Urban Deliveries the zero values can be due to the fact that the plss are outside the C2VSim domain, while for water levels the plss nodes outside the modeling domain get nan values.

Calculation procedure
----
We assume that the plss correspond to square mile rectangular center.

Therefore first we construct the rectangular and then we identify the intesecting elements.

![alt text][exampleAverage.png]

For each element we multiply the C2Vsim timeseries wit hthe ratio of intersected area/ element area. Then we sum all the weighted averaged timeseries

