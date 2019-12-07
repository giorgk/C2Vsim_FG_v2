# First fetch some functions from another repo 
source("../../C2VsimCG/Rwrkspc/c2vsim_io.R")
source("npsat_functions.R")

# Read the nodes
XY <- c2vsim.readNodes(filename = "../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Preprocessor/C2VSimFG_Nodes.dat", ND = 30179, Nskip = 90)
# Read the mesh
MSH <- c2vsim.readMesh(filename = "../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Preprocessor/C2VSimFG_Elements.dat",
                       NE = 32537, Nskip = 142,Ncols = 6)

# Write obj file
npsat.writeMesh2obj(filename = "temp.obj",XY = cbind(XY[,-1],0),MSH = MSH[,c(-1,-6)])
# The temp.obj file was converted to quadrilateral mesh in Zbrush (THis is the C2VSIM_16K.OBJ)
# Next the quadrilateral mesh was edited in Houdini to modify the remaining triangle elements.
# The final quadrilateral only mesh is the InputMesh_Modif.obj.
QuadMesh <- npsat.readOBJmesh("InputMesh_Modif.obj")

# Write the mesh input file
npsat.Input.WriteMesh("inputfiles/c2vsimMesh.npsat", XY = QuadMesh[[1]], MSH = QuadMesh[[2]])

