%% Paths
c2vsim_path = '..\c2vsimfg-betapublicrelease\C2VSimFG-BETA_PublicRelease\';
gis_data = '..\gis_data\';
mat_data = '..\mat_data\';
%% Hydraulic head
C2VsimHead = readC2Vsim_headalloutput([c2vsim_path 'Results\C2VSimFG_GW_HeadAll.out']);
save([mat_data 'C2VsimHead'],'C2VsimHead');