%{
%% Groundwater levels
% For the groundwater levels use the following line to read the HeadAll.out
% file. The name of this file is specified under the
% Simulation\Groundwater\C2VSimFG_Groundwater1974.dat file GWALLOUTFL
%
% h = readC2Vsim_headalloutput('..\..\C2VSimFG-BETA_PublicRelease\Results\C2VSimFG_GW_HeadAll.out');
% This take alot of time so read once and save it. Here we load this from
% the C2VSimHead_OR.mat file

load('..\wrkspc\C2VSimHead_OR.mat')
% read the meshnodes
ND = shaperead('..\gis_data\C2Vsim_Nodes');
p = [[ND.X]',[ND.Y]'];
% load the PLSS data
% First load the csv in QGIS to convert it to a shapefile with EPSG:26910
%read the shapefile
plss = shaperead('Liam_legal_to_latlong');

waterLevel = nan(3606,505);
% load the C2Vsim Mesh elements by running the Read Elements section of the
% ExtractData.m script

% Separate the quads from the triangles
id_qd = find(MSH(:,4) ~= 0 );
id_tr = find(MSH(:,4) == 0 );
msh_qd = MSH(id_qd,:);
msh_tr = MSH(id_tr,1:3);
% find in which element is point belong. This is using an mSim function
% find_elem_id_point <http://subsurface.gr/joomla/msim_doc/find_elem_id_point_help.html>
el_id_qd = find_elem_id_point(p, msh_qd,  [[plss.X]' [plss.Y]'], 5);
el_id_tr = find_elem_id_point(p, msh_tr,  [[plss.X]' [plss.Y]'], 5);
%
for ii = 1:length(plss)
    p_id = [];
    if ~isnan(el_id_tr(ii))
        BC = BaryCoord2D([plss(ii,1).X plss(ii,1).Y], p(msh_tr(el_id_tr(ii),:),:));
        p_id = msh_tr(el_id_tr(ii),:);
    end
    if ~isnan(el_id_qd(ii))
        BC = BaryCoord2D([plss(ii,1).X plss(ii,1).Y], p(msh_qd(el_id_qd(ii),:),:));
        p_id = msh_qd(el_id_qd(ii),:);
    end
    if ~isempty(p_id)
        for jj = 1:length(h)
            waterLevel(ii,jj) = sum(h{jj,2}(p_id,2).*BC);
        end
    end
end

%% Water levels Write to ascii format
frmt = '%s';
for ii = 1:size(waterLevel,2) + 4
    frmt = [frmt ' %.3f'];
end
frmt = [frmt '\n'];

fid = fopen('WaterLevelPLSSLiam.dat', 'w');
for ii = 1:size(waterLevel,1)
    fprintf(fid, frmt, plss(ii,1).CO_MTRS, [plss(ii,1).X1 plss(ii,1).Y1 plss(ii,1).X plss(ii,1).Y waterLevel(ii,:)]);
end
fclose(fid);
%}
%
%% Groundwater Pumping
GBinfo = h5info('..\C2VSimFG-BETA_PublicRelease\Results\C2VSimFG_GW_ZBudget.hdf');
IWFM_DatasetNames = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(1).Name GBinfo.Name GBinfo.Groups(1).Datasets(5).Name]);
for ii = 1:4
    IWFMelemIDS{ii,1} = double(h5read(GBinfo.Filename, ...
        [GBinfo.Groups(1).Name GBinfo.Name GBinfo.Groups(1).Datasets(5+ii).Name]));
end
PumpElemCol = 20;
PumpWellCol = 22;
for ii = 1:4
    PumpElem{ii,1} = h5read(GBinfo.Filename, ...
        [GBinfo.Groups(ii+1).Name GBinfo.Name GBinfo.Groups(ii+1).Datasets(13).Name])';
    PumpWell{ii,1} = h5read(GBinfo.Filename, ...
        [GBinfo.Groups(ii+1).Name GBinfo.Name GBinfo.Groups(ii+1).Datasets(14).Name])';
end

ND = shaperead('..\gis_data\C2Vsim_Nodes');
p = [[ND.X]',[ND.Y]'];
% load the C2Vsim Mesh elements by running the Read Elements section of the
% ExtractData.m script
% Calculate element barycenters
for ii = 1:size(MSH,1)
    n = 4;
    if MSH(ii,4) == 0
        n = 3;
    end
    elcc(ii,:) = [mean(p(MSH(ii,1:n),1)) mean(p(MSH(ii,1:n),2))];
end
load('..\..\mat_data\C2VsimWells');

plss = shaperead('Liam_legal_to_latlong');
plssPumping = nan(length(plss), size(PumpElem{1,1},1));
mile = 1609.344;
thereAreWells = [];
for ii = 1:length(plss)
    ii
    % Assume that the plss centers correspond to 1sq mi rectangles aligned
    % with the X, Y coordinates
    % create rectangle
    cc = [plss(ii,1).X plss(ii,1).Y];
    poly = [cc(1) - mile/2 cc(2) - mile/2; cc(1) + mile/2 cc(2) - mile/2; ...
            cc(1) + mile/2 cc(2) + mile/2; cc(1) - mile/2 cc(2) + mile/2];
    % sort elements by the distance
    dst = sqrt((cc(1) - elcc(:,1)).^2 + (cc(2) - elcc(:,2)).^2);
    [D, I] = sort(dst);
    cum_area = 0;
    Qpump = zeros(1,size(PumpElem{1,1},1));
    for jj = 1:length(I)
        n = 4;
        if MSH(I(jj),4) == 0; n = 3; end
        polyEl = [p(MSH(I(jj),1:n),1) p(MSH(I(jj),1:n),2)];
        [bx, by] = polybool('&',poly(:,1), poly(:,2), polyEl(:,1), polyEl(:,2));
        if ~isempty(bx)
            a = polyarea(bx, by);
            Elarea = polyarea(polyEl(:,1), polyEl(:,2));
            for kk = 1:4
                lay_el_id = find(IWFMelemIDS{kk,1}(:,PumpElemCol)~=0);
                rowid = find(lay_el_id == I(jj) );
                if ~isempty(rowid)
                    Qpump = Qpump + (a/Elarea).*PumpElem{kk,1}(:,rowid)';
                end
            end
            cum_area = cum_area + a;
        end
        if abs(cum_area - mile^2) < 0.1 || jj > 20
            plssPumping(ii,:) = Qpump;
            break
        end
    end
    
    % Check if any iwfm well is inside the poly
    inwell = inpolygon(iwfmwells(:,2), iwfmwells(:,3), poly(:,1), poly(:,2));
    well_id = find(inwell);
    if isempty(well_id)
        thereAreWells = [thereAreWells;well_id];
    end
end
%% Write Pumping to ascii format
frmt = '%s';
for ii = 1:size(plssPumping,2) + 4
    frmt = [frmt ' %.3f'];
end
frmt = [frmt '\n'];

writeasciiDataLiam('PumpingPLSSLiam.dat', frmt, plss, plssPumping);
%
%% Ag and Urban deliveries
GBinfo = h5info('..\C2VSimFG-BETA_PublicRelease\Results\C2VSimFG_L&WU_ZBudget.hdf');
% NonPond
NonPondArea = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(1).Name]);
NonPondAgDeliv = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(2).Name]);
NonPondAgPump = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(8).Name]);

% Refuge
RefugeArea = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(12).Name]);
RefugeDeliv = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(13).Name]);
RefugePump = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(20).Name]);

% Rice
RiceArea = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(23).Name]);
RiceDeliv = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(24).Name]);
RicePump = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(31).Name]);

% Urban
UrbanArea = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(34).Name]);
UrbanDeliv = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(35).Name]);
UrbanPump = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(37).Name]);
%%
ND = shaperead('..\gis_data\C2Vsim_Nodes');
p = [[ND.X]',[ND.Y]'];
% load the C2Vsim Mesh elements by running the Read Elements section of the
% ExtractData.m script
% Calculate element barycenters
for ii = 1:size(MSH,1)
    n = 4;
    if MSH(ii,4) == 0
        n = 3;
    end
    elcc(ii,:) = [mean(p(MSH(ii,1:n),1)) mean(p(MSH(ii,1:n),2))];
end

plss = shaperead('Liam_legal_to_latlong');
[Nel, Nt] = size(NonPondArea);
%%
%}
plssNonPondArea = nan(length(plss), Nt);
plssNonPondDeliv = nan(length(plss), Nt);
plssNonPondPump = nan(length(plss), Nt);
plssRefugeArea = nan(length(plss), Nt);
plssRefugeDeliv = nan(length(plss), Nt);
plssRefugePump = nan(length(plss), Nt);
plssRiceArea = nan(length(plss), Nt);
plssRiceDeliv = nan(length(plss), Nt);
plssRicePump = nan(length(plss), Nt);
plssUrbanArea = nan(length(plss), Nt);
plssUrbanDeliv = nan(length(plss), Nt);
plssUrbanPump = nan(length(plss), Nt);

mile = 1609.344;
for ii = 1:length(plss)
    ii
    % Assume that the plss centers correspond to 1sq mi rectangles aligned
    % with the X, Y coordinates
    % create rectangle
    cc = [plss(ii,1).X plss(ii,1).Y];
    poly = [cc(1) - mile/2 cc(2) - mile/2; cc(1) + mile/2 cc(2) - mile/2; ...
        cc(1) + mile/2 cc(2) + mile/2; cc(1) - mile/2 cc(2) + mile/2];
    % sort elements by the distance
    dst = sqrt((cc(1) - elcc(:,1)).^2 + (cc(2) - elcc(:,2)).^2);
    [D, I] = sort(dst);
    cum_area = 0;
    
    tempNPArea = zeros(1,Nt);
    tempNPDeliv = zeros(1,Nt);
    tempNPPump = zeros(1,Nt);
    tempRFArea = zeros(1,Nt);
    tempRFDeliv = zeros(1,Nt);
    tempRFPump = zeros(1,Nt);
    tempRCArea = zeros(1,Nt);
    tempRCDeliv = zeros(1,Nt);
    tempRCPump = zeros(1,Nt);
    tempURArea = zeros(1,Nt);
    tempURDeliv = zeros(1,Nt);
    tempURPump = zeros(1,Nt);
    
    for jj = 1:length(I)
        n = 4;
        if MSH(I(jj),4) == 0; n = 3; end
        polyEl = [p(MSH(I(jj),1:n),1) p(MSH(I(jj),1:n),2)];
        [bx, by] = polybool('&',poly(:,1), poly(:,2), polyEl(:,1), polyEl(:,2));
        if ~isempty(bx)
            a = polyarea(bx, by);
            Elarea = polyarea(polyEl(:,1), polyEl(:,2));
            r = a/Elarea;
            tempNPArea = tempNPArea + r.*NonPondArea(I(jj),:);
            tempNPDeliv = tempNPDeliv + r.*NonPondAgDeliv(I(jj),:);
            tempNPPump = tempNPPump + r.*NonPondAgPump(I(jj),:);
            tempRFArea = tempRFArea + r.*RefugeArea(I(jj),:);
            tempRFDeliv = tempRFDeliv + r.*RefugeDeliv(I(jj),:);
            tempRFPump = tempRFPump + r.*RefugePump(I(jj),:);
            tempRCArea = tempRCArea + r.*RiceArea(I(jj),:);
            tempRCDeliv = tempRCDeliv + r.*RiceDeliv(I(jj),:);
            tempRCPump = tempRCPump + r.*RicePump(I(jj),:);
            tempURArea = tempURArea + r.*UrbanArea(I(jj),:);
            tempURDeliv = tempURDeliv + r.*UrbanDeliv(I(jj),:);
            tempURPump = tempURPump + r.*UrbanPump(I(jj),:);
            
            cum_area = cum_area + a;
        end
        
        if abs(cum_area - mile^2) < 0.1 || jj > 20
            plssNonPondArea(ii,:) = tempNPArea;
            plssNonPondDeliv(ii,:) = tempNPDeliv;
            plssNonPondPump(ii,:) = tempNPPump;
            plssRefugeArea(ii,:) = tempRFArea;
            plssRefugeDeliv(ii,:) = tempRFDeliv;
            plssRefugePump(ii,:) = tempRFPump;
            plssRiceArea(ii,:) = tempRCArea;
            plssRiceDeliv(ii,:) = tempRCDeliv;
            plssRicePump(ii,:) = tempRCPump;
            plssUrbanArea(ii,:) = tempURArea;
            plssUrbanDeliv(ii,:) = tempURDeliv;
            plssUrbanPump(ii,:) = tempURPump;
            break
        end
    end
end
%% Write data
frmt = '%s';
for ii = 1:size(tempNPArea,2) + 4
    frmt = [frmt ' %.3f'];
end
frmt = [frmt '\n'];
writeasciiDataLiam('NonPondArea.dat', frmt, plss, plssNonPondArea);
writeasciiDataLiam('NonPondDeliv.dat', frmt, plss, plssNonPondDeliv);
writeasciiDataLiam('NonPondPump.dat', frmt, plss, plssNonPondPump);
writeasciiDataLiam('RefugeArea.dat', frmt, plss, plssRefugeArea);
writeasciiDataLiam('RefugeDeliv.dat', frmt, plss, plssRefugeDeliv);
writeasciiDataLiam('RefugePump.dat', frmt, plss, plssRefugePump);
writeasciiDataLiam('RiceArea.dat', frmt, plss, plssRiceArea);
writeasciiDataLiam('RiceDeliv.dat', frmt, plss, plssRiceDeliv);
writeasciiDataLiam('RicePump.dat', frmt, plss, plssRicePump);
writeasciiDataLiam('UrbanArea.dat', frmt, plss, plssUrbanArea);
writeasciiDataLiam('UrbanDeliv.dat', frmt, plss, plssUrbanDeliv);
writeasciiDataLiam('UrbanPump.dat', frmt, plss, plssUrbanPump);
