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
%%
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

fid = fopen('PLSSLiam.dat', 'w');
for ii = 1:size(waterLevel,1)
    fprintf(fid, frmt, plss(ii,1).CO_MTRS, [plss(ii,1).X1 plss(ii,1).Y1 plss(ii,1).X plss(ii,1).Y waterLevel(ii,:)]);
end
fclose(fid);


