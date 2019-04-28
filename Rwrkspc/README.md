The file FaceData.RData is created from the PreprocesInputs.R script and containts the following variables
- XY [30179x3] ID, X, Y The coordinates of the C2Vsim nods
- MSH [32537x6] ID, ND:1,2,3,4, S(subregion ID). The finite element mesh
- strat [30179x6] ID, GSE, L1, L2, L3, L4 The node elevations of the groundwater surface and the 4 layers
- FcLm [62719x2] The face elements. The face flows values correspond to this matrix
- faceIndex [32537x4] The index of each face in the FcLM matrix.
- faceArea [62719x4] The Area of each face. Each column corresponds to a layer
- elemArea [32537x1] The area of each mesh element

The PreprocOutputs.R script is used to extract the flow values and convert the flow volume to rates.
At the end it writes an PartTrackData.h5 file with all the required information by the particle tracking algorithm.

Last the ExtractFlowField.R extracts the vertical and horizontal flow field for selected months and prints the output to csv files
