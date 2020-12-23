%% C2Vsim main Path
c2vsim_path = ['..' filesep 'c2vsimfg_beta2_publicrelease' filesep 'C2VSimFG_BETA2_PublicRelease' filesep];
%% Read constraint BC node IDs
fid = fopen([c2vsim_path 'Simulation' filesep 'Groundwater' filesep 'C2VSimFG_ConstrainedHeadBC.dat']);
C = textscan(fid,'%d %d %d %d %f %f %d %d %s %s %s',109,'HeaderLines',119);
fclose(fid);
%% Read All heads
% ALLHEADS = readC2Vsim_headalloutput([c2vsim_path filesep 'Results' filesep 'C2VSimFG_GW_HeadAll.out']);
% or
load('ALLHEADS')
%% Read node coordinates in 3310
XY_3310 = shaperead('C2VsimNodes_3310.shp');
XY_3310 = [[XY_3310.X]' [XY_3310.Y]'];
%% Read the mesh boundary nodes
mesh_bnd_nd = shaperead('C2VsimMesh_Outline_3310');
[XX, YY] = polysplit(mesh_bnd_nd.X, mesh_bnd_nd.Y);
mesh_bnd_nd = [XX{1,1}' YY{1,1}'];
%% Interpolate the head values on the boundary nodes
msh_heads = nan(size(mesh_bnd_nd,1), length(ALLHEADS));
for ii = 1:length(ALLHEADS)
    ii
    F = scatteredInterpolant(XY_3310(:,1), XY_3310(:,2), ALLHEADS{ii,2}(:,1),'linear','nearest');
    msh_heads(:,ii) = F(mesh_bnd_nd(:,1), mesh_bnd_nd(:,2));
end
%% Compute statistics
Head_stats_monthly = [mean(msh_heads,2) std(msh_heads,0,2) std(msh_heads,1,2)];
Coef_Var_monthly = [Head_stats_monthly(:,2)./Head_stats_monthly(:,1) Head_stats_monthly(:,3)./Head_stats_monthly(:,1)];
%%
prct_thres = 50;
thres = prctile(Coef_Var_monthly(:,1),prct_thres);
% find one node that is above the threshold to be used as starting node
i_start = find(Coef_Var_monthly(:,1) > thres,1,'first');
indices = [i_start:size(mesh_bnd_nd,1) 1:i_start-1];
%%
i_start_line = [];
i_end_line = [];
BC_LINES = {};
cnt_ln = 1;
for ii = 1:length(indices)
    if isempty(i_start_line)
        if Coef_Var_monthly(indices(ii),1) < thres
            i_start_line = ii;
        end
    else
        if Coef_Var_monthly(indices(ii),1) > thres
            i_end_line = ii - 1;
        end
    end
    if ~isempty(i_start_line) && ~isempty(i_end_line)
        if i_start_line ~= i_end_line
            BC_LINES{cnt_ln,1} = indices(i_start_line):indices(i_end_line);
            cnt_ln = cnt_ln + 1;
        end
        i_start_line = [];
        i_end_line = [];
    end
end
%% print boundary files
fid = fopen(['inputfiles' filesep 'c2vsimBC_02_12_p' num2str(prct_thres) '.npsat'], 'w');
fprintf(fid, '%d\n',length(BC_LINES));
for ii = 1:length(BC_LINES)
    temp_name = ['bnd_02_12_p' num2str(prct_thres) '_' num2str(ii) '.npsat'];
    fprintf(fid, 'EDGETOP 0 %s\n',['BC_files/' temp_name]);
    fid1 = fopen(['inputfiles' filesep 'BC_files' filesep temp_name],'w');
    fprintf(fid1, 'BOUNDARY_LINE\n');
    fprintf(fid1, '%d 1 1\n', length(BC_LINES{ii,1}));
    fprintf(fid1, '%0.2f %0.2f %0.2f\n', [mesh_bnd_nd(BC_LINES{ii,1},:) 0.3048*Head_stats_monthly(BC_LINES{ii,1},1)]');
    fclose(fid1);
end
fclose(fid);
%% Prepare Initial top head
% Read the expanded mesh of the original element mesh
c2vsimNDxpnd = shaperead('C2VsimMeshORExpanded_3310');
NDxpnd = [[c2vsimNDxpnd.X]' [c2vsimNDxpnd.Y]'];
% Read the C2Vsim nodes
c2vsim_nodes = shaperead('C2VsimNodes_3310.shp');
ND = [[c2vsim_nodes.X]' [c2vsim_nodes.Y]'];
%
HeadAv = zeros(size(ND,1),1);
for ii = start_tm_id:end_tm_id
    HeadAv = HeadAv + ALLHEADS{ii,2}(:,1);
end
HeadAv = HeadAv./length(start_tm_id:end_tm_id);
HeadAv = HeadAv*0.3048;
Fhead = scatteredInterpolant(ND(:,1), ND(:,2), HeadAv,'linear','nearest');
Hxpnd = Fhead(NDxpnd(:,1), NDxpnd(:,2));
%% write top
writeScatteredData(['inputfiles' filesep 'c2vsimTop_init.npsat'],...
    struct('PDIM', 2, 'TYPE', 'HOR', 'MODE', 'SIMPLE'),...
    [NDxpnd Hxpnd]);


