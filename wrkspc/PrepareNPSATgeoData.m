%% C2Vsim main Path
c2vsim_path = ['..' filesep 'c2vsimfg_beta2_publicrelease' filesep 'C2VSimFG_BETA2_PublicRelease' filesep];
%% Mesh file
% Read the obj file
[p, MSH] = readOBJmesh('../NPSAT/InputMesh_Modif.obj', 4);
%% print mesh to shapefile
for ii = 1:size(MSH,1)
   S(ii,1).Geometry = 'Polygon';
   x = p(MSH(ii,:),1);
   y = p(MSH(ii,:),2);
   S(ii,1).BoundingBox = [min(x) min(y); max(x) max(y)];
   S(ii,1).X = x;
   S(ii,1).Y = y;
   S(ii,1).id = ii;
end
shapewrite(S,'C2VsimMesh')
%% Convert to EPGS 3310 and read it back
SS = shaperead('C2VsimMesh_3310');
%% Write mesh file 
p = [];
MSH = nan(length(SS),4);
for ii = 1:length(SS)
    for jj = 1:4
        xp = SS(ii,1).X(jj);
        yp = SS(ii,1).Y(jj);
        if isempty(p)
            p = [xp yp];
            id = 1;
        else
            dst = sqrt((p(:,1) - xp).^2 + (p(:,2) - yp).^2);
            [cc, dd] = min(dst);
            if cc < 1
                id = dd;
            else
                p = [p;xp yp];
                id = size(p,1);
            end
        end
        MSH(ii,jj) = id;
    end
    
end
%%
writeMeshfile('inputfiles/c2vsimMesh_3310.npsat',p,MSH);
%% Convert buffered file to 3310
[px, MSHx] = readOBJmesh('../NPSAT/ExpandedMesh.obj', 4);
%%
for ii = 1:size(MSHx,1)
   S(ii,1).Geometry = 'Polygon';
   x = px(MSHx(ii,:),1);
   y = px(MSHx(ii,:),2);
   S(ii,1).BoundingBox = [min(x) min(y); max(x) max(y)];
   S(ii,1).X = x;
   S(ii,1).Y = y;
   S(ii,1).id = ii;
end
shapewrite(S,'C2VsimMeshExpanded')
%% Convert mesh and nodes to 3310
% First read the data 
%% Read Node coordinates
c2vsim_path = ['..' filesep 'c2vsimfg_beta2_publicrelease' filesep 'C2VSimFG_BETA2_PublicRelease' filesep];
fid = fopen([c2vsim_path 'Preprocessor' filesep 'C2VSimFG_Nodes.dat'],'r');
temp = textscan(fid, '%f %f %f', 30179, 'HeaderLines',90);
fclose(fid);
XY = [temp{1,2} temp{1,3}];
ND_ID = temp{1,1};
for ii = 1:length(ND_ID)
    C2Vsim_nodes(ii,1).Geometry = 'Point';
    C2Vsim_nodes(ii,1).X = XY(ii,1);
    C2Vsim_nodes(ii,1).Y = XY(ii,2);
    C2Vsim_nodes(ii,1).ID = ND_ID(ii);
end
%% Read stratigraphy
fid = fopen([c2vsim_path 'Preprocessor' filesep 'C2VSimFG_Stratigraphy.dat'],'r');
strat = textscan(fid, '%f %f %f %f %f %f %f %f %f %f', 30179, 'HeaderLines',105);
fclose(fid);
ELEV = strat{1,2};
ELEV(:,2) = ELEV(:,1) - strat{1,4};
ELEV(:,3) = ELEV(:,2) - max(strat{1,5},0.1);
ELEV(:,4) = ELEV(:,3) - strat{1,6};
ELEV(:,5) = ELEV(:,4) - strat{1,8};
ELEV(:,6) = ELEV(:,5) - strat{1,10};
ELEV = ELEV*0.3048;
for ii = 1:length(ND_ID)
    C2Vsim_nodes(ii,1).L1 = ELEV(ii,1);
    C2Vsim_nodes(ii,1).L2 = ELEV(ii,2);
    C2Vsim_nodes(ii,1).L3 = ELEV(ii,3);
    C2Vsim_nodes(ii,1).L4 = ELEV(ii,4);
    C2Vsim_nodes(ii,1).L5 = ELEV(ii,5);
    C2Vsim_nodes(ii,1).L6 = ELEV(ii,6);
end
shapewrite(C2Vsim_nodes, 'C2VsimNodes')
%% Read Elements of parametric grid for hydraulic conductivity
fid = fopen([c2vsim_path 'Simulation' filesep 'Groundwater' filesep 'C2VSimFG_Groundwater1974.dat'],'r');
temp = textscan(fid, '%f %f %f %f %f', 1419, 'HeaderLines',354);
fclose(fid);
%% Read Nodes of parametric grid for hydraulic conductivity
fid = fopen([c2vsim_path 'Simulation' filesep 'Groundwater' filesep 'C2VSimFG_Groundwater1974.dat'],'r');
temp = textscan(fid, '%f %f %f %f %f %f %f %f', 4*1419, 'HeaderLines',1793);
fclose(fid);
for ii = 1:1403
    paramND(ii,1).Geometry = 'Point';
    paramND(ii,1).X = temp{1,2}((ii-1)*4 + 1);
    paramND(ii,1).Y = temp{1,3}((ii-1)*4 + 1);
    paramND(ii,1).ID = temp{1,1}((ii-1)*4 + 1);
    paramND(ii,1).PKH1 = temp{1,4}((ii-1)*4 + 1);
    paramND(ii,1).PKH2 = temp{1,1}((ii-1)*4 + 2);
    paramND(ii,1).PKH3 = temp{1,1}((ii-1)*4 + 3);
    paramND(ii,1).PKH4 = temp{1,1}((ii-1)*4 + 4);
    paramND(ii,1).PL1 = temp{1,8}((ii-1)*4 + 1);
    paramND(ii,1).PL2 = temp{1,5}((ii-1)*4 + 2);
    paramND(ii,1).PL3 = temp{1,5}((ii-1)*4 + 3);
    paramND(ii,1).PL4 = temp{1,5}((ii-1)*4 + 4);
    paramND(ii,1).PV1 = temp{1,7}((ii-1)*4 + 1);
    paramND(ii,1).PV2 = temp{1,4}((ii-1)*4 + 2);
    paramND(ii,1).PV3 = temp{1,4}((ii-1)*4 + 3);
    paramND(ii,1).PV4 = temp{1,4}((ii-1)*4 + 4);
    paramND(ii,1).PS1 = temp{1,5}((ii-1)*4 + 1);
    paramND(ii,1).PS2 = temp{1,2}((ii-1)*4 + 2);
    paramND(ii,1).PS3 = temp{1,2}((ii-1)*4 + 3);
    paramND(ii,1).PS4 = temp{1,2}((ii-1)*4 + 4);
    paramND(ii,1).PN1 = temp{1,6}((ii-1)*4 + 1);
    paramND(ii,1).PN2 = temp{1,3}((ii-1)*4 + 2);
    paramND(ii,1).PN3 = temp{1,3}((ii-1)*4 + 3);
    paramND(ii,1).PN4 = temp{1,3}((ii-1)*4 + 4);
end
shapewrite(paramND, 'parametricGrid');
%% Read the expanded mesh of the original C2Vsim mesh and convert it to 3310
[p_or_xp, MSH_or_xp] = readOBJmesh(['..' filesep 'NPSAT' filesep 'C2VsimMeshORExpanded.obj'], 4);
%% write to shapefile to convert it to 3310
for ii = 1:size(p_or_xp,1)
   S(ii,1).Geometry = 'Point';
   S(ii,1).X = p_or_xp(ii,1);
   S(ii,1).Y = p_or_xp(ii,2);
   S(ii,1).id = ii;
end
shapewrite(S, 'C2VsimMeshORExpanded');
%% Read the expanded original mesh where the interpolated values should be assigned
c2vsimNDxpnd = shaperead('C2VsimMeshORExpanded_3310');
NDxpnd = [[c2vsimNDxpnd.X]' [c2vsimNDxpnd.Y]'];

%% Read the shapefile with the elevation information
c2vsimNodes_3310 = shaperead('C2VsimNodes_3310');
%% Interpolate the bottom elevation to the expanded node set
Fbot = scatteredInterpolant([c2vsimNodes_3310.X]', [c2vsimNodes_3310.Y]', [c2vsimNodes_3310.L6]', 'linear', 'nearest');
%% write data to file
bot_data = [NDxpnd Fbot(NDxpnd(:,1), NDxpnd(:,2))];
writeScatteredData(['inputfiles' filesep 'c2vsimBottom_init_3310.npsat'],...
    struct('PDIM', 2, 'TYPE', 'HOR', 'MODE', 'SIMPLE'), bot_data);
%% interpolate hydraulic properties for each layer
paramND_3310 = shaperead('parametricGrid_3310');
%% initialize data
DATA_HK = NDxpnd;
DATA_VK = NDxpnd;
DATA_SY = NDxpnd;
%% Layer 1
% interpolate the data from the parametric grid
% create interpolants
FHK = scatteredInterpolant([paramND_3310.X]', [paramND_3310.Y]', [paramND_3310.PKH1]', 'linear', 'nearest');
FVK = scatteredInterpolant([paramND_3310.X]', [paramND_3310.Y]', [paramND_3310.PL1]', 'linear', 'nearest');
FSY = scatteredInterpolant([paramND_3310.X]', [paramND_3310.Y]', [paramND_3310.PN1]', 'linear', 'nearest');
DATA_HK = [DATA_HK FHK(NDxpnd(:,1), NDxpnd(:,2))];
DATA_VK = [DATA_VK FVK(NDxpnd(:,1), NDxpnd(:,2))];
DATA_SY = [DATA_SY FSY(NDxpnd(:,1), NDxpnd(:,2))];
% interpolate the elevation between layer 1 - 2
FLay12 = scatteredInterpolant([c2vsimNodes_3310.X]', [c2vsimNodes_3310.Y]', [c2vsimNodes_3310.L2]', 'linear', 'nearest');
Lay12 = FLay12(NDxpnd(:,1), NDxpnd(:,2));
DATA_HK = [DATA_HK Lay12];
DATA_VK = [DATA_VK Lay12];
DATA_SY = [DATA_SY Lay12];
%% Layer 2 aquiclude is the aquiclude which is the only aquiclude that has non zero height
FVK = scatteredInterpolant([paramND_3310.X]', [paramND_3310.Y]', [paramND_3310.PL2]', 'linear', 'nearest'); % this is the vertical conductivity of the layer
FLK = scatteredInterpolant([paramND_3310.X]', [paramND_3310.Y]', [paramND_3310.PV2]', 'linear', 'nearest'); % this is the vertical in the aquiclude
% The vertical hydraulic conductivity gets the PL values everywhere with
% zero height (0.1) and PV where there is some height in the aquiclude
tmp_vk = FVK(NDxpnd(:,1), NDxpnd(:,2));
tmp_lk = FLK(NDxpnd(:,1), NDxpnd(:,2));
Flay2height = scatteredInterpolant([c2vsimNodes_3310.X]', [c2vsimNodes_3310.Y]', [c2vsimNodes_3310.L2]'-[c2vsimNodes_3310.L3]', 'linear', 'nearest');
lay2height = Flay2height(NDxpnd(:,1), NDxpnd(:,2));
id_aquiclude = find(lay2height > 1);
tmp_vk(id_aquiclude) = tmp_lk(id_aquiclude);
DATA_VK = [DATA_VK tmp_vk];
% interpolate the elevation between layer 2 - 3
FLay23 = scatteredInterpolant([c2vsimNodes_3310.X]', [c2vsimNodes_3310.Y]', [c2vsimNodes_3310.L3]', 'linear', 'nearest');
Lay23 = FLay23(NDxpnd(:,1), NDxpnd(:,2));
DATA_VK = [DATA_VK Lay23];
%% Layer 2 Aquifer
FHK = scatteredInterpolant([paramND_3310.X]', [paramND_3310.Y]', [paramND_3310.PKH2]', 'linear', 'nearest');
FSY = scatteredInterpolant([paramND_3310.X]', [paramND_3310.Y]', [paramND_3310.PN2]', 'linear', 'nearest');

DATA_HK = [DATA_HK FHK(NDxpnd(:,1), NDxpnd(:,2))];
DATA_VK = [DATA_VK FVK(NDxpnd(:,1), NDxpnd(:,2))];
DATA_SY = [DATA_SY FSY(NDxpnd(:,1), NDxpnd(:,2))];
% interpolate the elevation between layer 3 - 4
FLay34 = scatteredInterpolant([c2vsimNodes_3310.X]', [c2vsimNodes_3310.Y]', [c2vsimNodes_3310.L4]', 'linear', 'nearest');
Lay34 = FLay34(NDxpnd(:,1), NDxpnd(:,2));
DATA_HK = [DATA_HK Lay34];
DATA_VK = [DATA_VK Lay34];
DATA_SY = [DATA_SY Lay34];
%% Layer 3 Aquifer
FHK = scatteredInterpolant([paramND_3310.X]', [paramND_3310.Y]', [paramND_3310.PKH3]', 'linear', 'nearest');
FVK = scatteredInterpolant([paramND_3310.X]', [paramND_3310.Y]', [paramND_3310.PL3]', 'linear', 'nearest');
FSY = scatteredInterpolant([paramND_3310.X]', [paramND_3310.Y]', [paramND_3310.PN3]', 'linear', 'nearest');
DATA_HK = [DATA_HK FHK(NDxpnd(:,1), NDxpnd(:,2))];
DATA_VK = [DATA_VK FVK(NDxpnd(:,1), NDxpnd(:,2))];
DATA_SY = [DATA_SY FSY(NDxpnd(:,1), NDxpnd(:,2))];
% interpolate the elevation between layer 4 - 5
FLay45 = scatteredInterpolant([c2vsimNodes_3310.X]', [c2vsimNodes_3310.Y]', [c2vsimNodes_3310.L5]', 'linear', 'nearest');
Lay45 = FLay45(NDxpnd(:,1), NDxpnd(:,2));
DATA_HK = [DATA_HK Lay45];
DATA_VK = [DATA_VK Lay45];
DATA_SY = [DATA_SY Lay45];
%% Layer 4 Aquifer
FHK = scatteredInterpolant([paramND_3310.X]', [paramND_3310.Y]', [paramND_3310.PKH4]', 'linear', 'nearest');
FVK = scatteredInterpolant([paramND_3310.X]', [paramND_3310.Y]', [paramND_3310.PL4]', 'linear', 'nearest');
FSY = scatteredInterpolant([paramND_3310.X]', [paramND_3310.Y]', [paramND_3310.PN4]', 'linear', 'nearest');
DATA_HK = [DATA_HK FHK(NDxpnd(:,1), NDxpnd(:,2))];
DATA_VK = [DATA_VK FVK(NDxpnd(:,1), NDxpnd(:,2))];
DATA_SY = [DATA_SY FSY(NDxpnd(:,1), NDxpnd(:,2))];
%% write data 
writeScatteredData(['inputfiles' filesep 'c2vsimHK_3310.npsat'], ...
    struct('PDIM', 3, 'TYPE', 'FULL', 'MODE', 'STRATIFIED'), DATA_HK);
writeScatteredData(['inputfiles' filesep 'c2vsimVK_3310.npsat'], ...
    struct('PDIM', 3, 'TYPE', 'FULL', 'MODE', 'STRATIFIED'), DATA_VK);
writeScatteredData(['inputfiles' filesep 'c2vsimSY_3310.npsat'], ...
    struct('PDIM', 3, 'TYPE', 'FULL', 'MODE', 'STRATIFIED'), DATA_SY);
%% Read river file
fid = fopen([c2vsim_path 'Preprocessor' filesep 'C2VSimFG_StreamsSpec.dat'],'r');
% read NRH, NR, NRTB
for ii = 1:79; tline = fgetl(fid); end
tline = fgetl(fid);
temp = textscan(tline, '%f / %s');
NRH = temp{1,1}(1); %NUmber of stream reaches
tline = fgetl(fid);
temp = textscan(tline, '%f / %s');
NRTB = temp{1,1}(1); %Number of data points in tables per stream node
% read comments
for ii = 1:22; tline = fgetl(fid);end
% read rivers one by one
clear C2Vsim_rivers
for ii = 1:NRH
    ii
    % for each stream read
    % a comment line
    while 1
        tline = fgetl(fid);
        if strcmp('C-------------------------------------------------------------------------------',deblank(tline))
            break;
        end
    end
    for jj = 1:5; tline = fgetl(fid);end
    temp = textscan(fid, '%f %f %f %f', 1);
    name = fgetl(fid);
    C2Vsim_rivers(ii,1).ID = temp{1,1};
    C2Vsim_rivers(ii,1).NAME = deblank(name);
    C2Vsim_rivers(ii,1).IBUR = temp{1,2};
    %C2Vsim_rivers(ii,1).IBDR = temp{1,3};
    C2Vsim_rivers(ii,1).IDWN = temp{1,3};
    for jj = 1:5; tline = fgetl(fid);end
    temp = textscan(fid, '%f %f', C2Vsim_rivers(ii,1).IBUR);
    C2Vsim_rivers(ii,1).IRV = temp{1,1};
    C2Vsim_rivers(ii,1).IGW = temp{1,2};
end
fclose(fid);
%% Create a shapefile with the rivers
% Read XY coordinates from the section 'Read Node coordinates'
clear S
for ii = 1:length(C2Vsim_rivers)
    S(ii,1).Geometry = 'Line';
    S(ii,1).X = [XY(C2Vsim_rivers(ii,1).IGW,1)' nan];
    S(ii,1).Y = [XY(C2Vsim_rivers(ii,1).IGW,2)' nan];
    S(ii,1).BoundingBox = [min(XY(C2Vsim_rivers(ii,1).IGW,1)) min(XY(C2Vsim_rivers(ii,1).IGW,2)); ...
        max(XY(C2Vsim_rivers(ii,1).IGW,1)) max(XY(C2Vsim_rivers(ii,1).IGW,2))];
    S(ii,1).ID = C2Vsim_rivers(ii,1).ID;
    S(ii,1).Name = C2Vsim_rivers(ii,1).NAME;
    S(ii,1).IBUR = C2Vsim_rivers(ii,1).IBUR;
    S(ii,1).IDWN = C2Vsim_rivers(ii,1).IDWN; 
end
shapewrite(S,'C2Vsim_Rivers');
%% Read small watersheds
fid = fopen([c2vsim_path 'Simulation' filesep 'C2VSimFG_SWatersheds_v1.dat'],'r');
% read NSW
for ii = 1:101; tline = fgetl(fid); end
tline = fgetl(fid);
temp = textscan(tline, '%f / %s');
NSW = temp{1,1}; %NUmber of stream reaches
for ii = 1:19; tline = fgetl(fid); end

clear C2Vsim_SWatersheds
for ii = 1:NSW
    tline = [];
    while isempty(tline)
        tline = fgetl(fid);
    end
    temp = textscan(tline, '%f %f %f %f %f %f', 1);
    C2Vsim_SWatersheds(ii,1).ID = temp{1,1};
    C2Vsim_SWatersheds(ii,1).AREAS = temp{1,2};
    C2Vsim_SWatersheds(ii,1).IWBTS = temp{1,3};
    NWB = temp{1,4};
    iwb = temp{1,5};
    temp = textscan(fid, '%f %f', NWB-1);
    C2Vsim_SWatersheds(ii,1).IWB = [iwb; temp{1,1}];
end
fclose(fid);
%%
id = [];
for ii = 1:length(C2Vsim_SWatersheds)
    if length(C2Vsim_SWatersheds(ii,1).IWB) == 1
        id = [id;ii];
    end
end
%% Create a shapefile for the small watersheds that discharge at least in 2 nodes
clear S
cnt = 1;
for ii = 1:length(C2Vsim_SWatersheds)
    if length(C2Vsim_SWatersheds(ii,1).IWB) == 1
        continue;
    end
    S(cnt,1).Geometry = 'Line';
    S(cnt,1).X = [XY(C2Vsim_SWatersheds(ii,1).IWB,1)' nan];
    S(cnt,1).Y = [XY(C2Vsim_SWatersheds(ii,1).IWB,2)' nan];
    S(cnt,1).BoundingBox = [min(XY(C2Vsim_SWatersheds(ii,1).IWB,1)) min(XY(C2Vsim_SWatersheds(ii,1).IWB,2)); ...
        max(XY(C2Vsim_SWatersheds(ii,1).IWB,1)) max(XY(C2Vsim_SWatersheds(ii,1).IWB,2))];
    S(cnt,1).ID = C2Vsim_SWatersheds(ii,1).ID;
    S(cnt,1).AREAS = C2Vsim_SWatersheds(ii,1).AREAS;
    S(cnt,1).IWBTS = C2Vsim_SWatersheds(ii,1).IWBTS;
    cnt = cnt + 1;
end
shapewrite(S, 'C2VsimSWatersheds')
%% make a unique list of segments
% save also for each segment the ids of rivers that contribute to each
% segment recharge
% Run "Read river file"
riv_segm = nan(10000,4);
IDS(10000,1).riv = [];
IDS(10000,1).sw = [];
ind = 1;
for ii = 1:length(C2Vsim_rivers)
    for jj = 1:length(C2Vsim_rivers(ii,1).IGW)-1
        a = C2Vsim_rivers(ii,1).IGW(jj);
        b = C2Vsim_rivers(ii,1).IGW(jj+1);
        c = C2Vsim_rivers(ii,1).IRV(jj);
        d = C2Vsim_rivers(ii,1).IRV(jj+1);
        if isempty(riv_segm)
            riv_segm = [a b c d];
        else
            id = find(riv_segm(:,1) == a & riv_segm(:,2) == b);
            if isempty(id)
                id = find(riv_segm(:,1) == b & riv_segm(:,2) == a);
                if isempty(id)
                    riv_segm(ind,:) = [a b c d];
                    IDS(ind,1).riv = [IDS(ind,1).riv; ii];
                    ind = ind + 1;
                else
                    IDS(id,1).riv = [IDS(id,1).riv; ii];
                end
            else
                IDS(id,1).riv = [IDS(id,1).riv; ii];
            end
        end
    end
end
%% repeat for smallwatersheds
wshed_segm = nan(10000,2);
ind = 1;
for ii = 1:length(C2Vsim_SWatersheds)
    if length(C2Vsim_SWatersheds(ii,1).IWB) == 1
        continue
    end
    for jj = 1:length(C2Vsim_SWatersheds(ii,1).IWB)-1
        a = C2Vsim_SWatersheds(ii,1).IWB(jj);
        b = C2Vsim_SWatersheds(ii,1).IWB(jj+1);
        %if isempty(wshed_segm)
        %    wshed_segm = [a b];
        %else
            id = find(wshed_segm(:,1) == a & wshed_segm(:,2) == b);
            if isempty(id)
                id = find(wshed_segm(:,1) == b & wshed_segm(:,2) == a);
                if isempty(id)
                    wshed_segm(ind,:) = [a b];
                    %IDS(ind,1).sw = [IDS(ind,1).sw; ii];
                    ind = ind + 1;
                else
                    %IDS(id,1).sw = [IDS(id,1).sw; ii];
                end
            else
                %IDS(id,1).sw = [IDS(id,1).sw; ii];
            end
        %end
    end
end
%% Make a list of the SW that dischagre to one node only
SW_single_nodes(100,1).nd = [];
SW_single_nodes(100,1).id = [];
ind1 = 1;
for ii = 1:length(C2Vsim_SWatersheds)
    if length(C2Vsim_SWatersheds(ii,1).IWB) == 1
        id = find([SW_single_nodes.nd]' == C2Vsim_SWatersheds(ii,1).IWB);
        if isempty(id)
            SW_single_nodes(ind1,1).nd = C2Vsim_SWatersheds(ii,1).IWB;
            SW_single_nodes(ind1,1).id = [SW_single_nodes(ind1,1).id; ii];
            ind1 = ind1 + 1;
        else
            SW_single_nodes(id,1).id = [SW_single_nodes(id,1).id; ii];
        end
    end
end








