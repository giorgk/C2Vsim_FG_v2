%% Paths
%
c2vsim_path = '..\c2vsimfg-betapublicrelease\C2VSimFG-BETA_PublicRelease\';
gis_data = '..\gis_data\';
mat_data = '..\mat_data\';
%%
c2vsim_path = '/media/giorgk/DATA/giorgk/Documents/C2Vsim_FG_v2/C2VSimFG-BETA_PublicRelease/';
mat_data = '/media/giorgk/DATA/giorgk/Documents/C2Vsim_FG_v2/mat_data/';
hou_data = '/media/giorgk/DATA/giorgk/Documents/C2Vsim_FG_v2/hou_data/';
%% Hydraulic head
C2VsimHead = readC2Vsim_headalloutput([c2vsim_path 'Results/C2VSimFG_GW_HeadAll.out']);
save([mat_data 'C2VsimHead'],'C2VsimHead');
%% GroundwaterBudgets 
GBinfo = h5info([c2vsim_path 'Results/C2VSimFG_GW_ZBudget.hdf']);
%% Create face flows
load([mat_data 'C2Vsim_Nodes.mat'], 'C2Vsim_nodes');
Xcoords = [C2Vsim_nodes.X]';
Ycoords = [C2Vsim_nodes.Y]';
faceElements = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(1).Name GBinfo.Name GBinfo.Groups(1).Datasets(17).Name])';

elementsNodes = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(1).Name GBinfo.Name GBinfo.Groups(1).Datasets(16).Name])';
%}
%%
clf
hold on
faceFlowsCoords = nan(size(faceElements,1),4);
for ii = 1:size(faceElements,1)
    is_bnd_face = false;
    if faceElements(ii,1) ~= 0
        eln1 = elementsNodes(faceElements(ii,1),:);
        if eln1(end) == 0; eln1(:,end) = [];end
    else
        is_bnd_face = true;
    end
    if faceElements(ii,2) ~= 0
        eln2 = elementsNodes(faceElements(ii,2),:);
        if eln2(end) == 0; eln2(:,end) = [];end
    else
        is_bnd_face = true;
    end
    
    if ~is_bnd_face
        c1 = [mean(Xcoords(eln1)) mean(Ycoords(eln1))]; 
        c2 = [mean(Xcoords(eln2)) mean(Ycoords(eln2))];
        faceFlowsCoords(ii,1:2) = [(c1(1) + c2(1))/2 (c1(2) + c2(2))/2];
        nrm = c2 - c1;
        faceFlowsCoords(ii,3:4) = nrm/sqrt(sum(nrm.^2));
    end
end
%% Extract X-Y Flows from all layers
for ii = 1:4
    faceflowLay{ii,1} = h5read(GBinfo.Filename, ...
        [GBinfo.Groups(ii+1).Name GBinfo.Name GBinfo.Groups(ii+1).Datasets(9).Name])';
end
%%
Nfaces = size(faceFlowsCoords,1);
frmt = '%f';
for ii = 1:9
    frmt = [frmt ' %f'];
end
frmt = [frmt ' \n'];
for ii = 1:size(faceflowLay{1,1},1)
    ii
    fid = fopen([hou_data 'faceflowlay' num2str(ii, '%04d') '.dat'],'w');
    TAB = [faceFlowsCoords(:,1) ...
           zeros(Nfaces,1) ...
           faceFlowsCoords(:,2) ...
           faceFlowsCoords(:,3) ...
           zeros(Nfaces,1) ...
           faceFlowsCoords(:,4) ...
           faceflowLay{1,1}(ii,:)' ...
           faceflowLay{2,1}(ii,:)' ...
           faceflowLay{3,1}(ii,:)' ...
           faceflowLay{4,1}(ii,:)' ];
    TAB(isnan(TAB(:,1)),:) = [];
    % X Y Z nx ny nz, lay1 lay2 lay3 lay4
    fprintf(fid, frmt, TAB');
    fclose(fid);
end
%% Deep percolation
% dee percolation - is all zeros
DeepPerc = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(5).Name])';

%% Vertical Flows
% dee percolation - is all zeros
VertFlows = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(4).Name GBinfo.Name GBinfo.Groups(4).Datasets(27).Name])';
    


