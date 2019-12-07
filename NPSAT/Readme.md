### NPSAT
In this folder there are scripts and functions that extract and prepare the data for the NPSAT engine model based on the C2Vsim simulation.


#### Pre-process
First download the C2Vsim model from [here](https://data.cnra.ca.gov/dataset/c2vsimfg_beta2) and run the model using the default configuration. This involves three steps:

1. Run the preprocecor
2. Run the simulation
3. Run the budget

#### Geometry and property files
* Mesh input file. This is the file wit hthe initial mesh. The generation of  quadrilateral mesh is a rather complicated issue and the options for generating good quality quadrilateral meshes are quite limited ([Gmsh](http://gmsh.info/) which is free but difficult to get the desired result, [cubit](https://cubit.sandia.gov/) which results in rather good quality mesh, bit it's not free to everyone). Our approach is to print the original mesh as [obj](https://en.wikipedia.org/wiki/Wavefront_.obj_file) file and then use a standard remesher tool from 3D industry [exoside](https://exoside.com/).


