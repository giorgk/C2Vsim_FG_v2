%% C2Vsim main Path
c2vsim_path = ['..' filesep 'c2vsimfg_beta2_publicrelease' filesep 'C2VSimFG_BETA2_PublicRelease' filesep];
%% Set time stamp
start_date = datetime(1973,10,1);
tm = start_date + calmonths(0:503);
%% Read groundwater Budget
GB_bud_info = h5info([c2vsim_path 'Results/C2VSimFG_GW_ZBudget.hdf']);
% This is how to obtain the IDS and the column names. However this is not
% used here
%Exampl1_GW_Budget.hdf
%h5read(GB_bud_info.Filename,[GB_bud_info.Groups(1).Name GB_bud_info.Name GB_bud_info.Groups(1).Datasets(6).Name])';
%% Read and calculate cumulative storage
CVstorage = zeros(504, 32537);
for ii = 1:4
    CVstorage = CVstorage + h5read(GB_bud_info.Filename,...
        [GB_bud_info.Groups(1+ii).Name GB_bud_info.Name GB_bud_info.Groups(1+ii).Datasets(22).Name])';
end
%% 
subregions = shaperead('F:\UCDAVIS\C2VSIM_FG_OR\C2Vsim_FG_v2\gis_data\C2Vsim_Elements.shp');
Basin1_ids = find([subregions.IRGE]' <= 7);
Basin2_ids = find([subregions.IRGE]' >= 8 & [subregions.IRGE]' <= 13);
Basin3_ids = find([subregions.IRGE]' >= 14);
Basin1_storage = sum(CVstorage(:,Basin1_ids),2);
Basin2_storage = sum(CVstorage(:,Basin2_ids),2);
Basin3_storage = sum(CVstorage(:,Basin3_ids),2);
Basin1_storage = Basin1_storage - Basin1_storage(1);
Basin2_storage = Basin2_storage - Basin2_storage(1);
Basin3_storage = Basin3_storage - Basin3_storage(1);

plot(tm,(Basin1_storage*2.29569e-5)/1000000,'DisplayName','SAC')
hold on
plot(tm,(Basin2_storage*2.29569e-5)/1000000,'DisplayName','SJV')
plot(tm,(Basin3_storage*2.29569e-5)/1000000,'DisplayName','TLB')
%% Groundwater Pumping
fid = fopen([c2vsim_path 'Simulation' filesep 'Groundwater' filesep 'C2VSimFG_ElemPump.dat']);
C = textscan(fid,'%d %d %d %d %f %f %f %f %d %d %d %d %d %d /%s',65074,'HeaderLines',119);
fclose(fid);
%% Pumping
GW_pump = zeros(504,1);
for ii = 1:4
    ii
    temp =  h5read(GB_bud_info.Filename,...
        [GB_bud_info.Groups(ii+1).Name GB_bud_info.Name GB_bud_info.Groups(ii+1).Datasets(13).Name])';
    GW_pump = GW_pump + sum(temp,2);
    
    temp =  h5read(GB_bud_info.Filename,...
        [GB_bud_info.Groups(ii+1).Name GB_bud_info.Name GB_bud_info.Groups(ii+1).Datasets(15).Name])';
    GW_pump = GW_pump + sum(temp,2);

end
%% Small Watersheds
ii = 1;
GW_swtrshds = zeros(504,1);
temp =  h5read(GB_bud_info.Filename,...
        [GB_bud_info.Groups(ii+1).Name GB_bud_info.Name GB_bud_info.Groups(ii+1).Datasets(18).Name])';
GW_swtrshds = GW_swtrshds + sum(temp,2);
temp =  h5read(GB_bud_info.Filename,...
        [GB_bud_info.Groups(ii+1).Name GB_bud_info.Name GB_bud_info.Groups(ii+1).Datasets(20).Name])';
GW_swtrshds = GW_swtrshds + sum(temp,2);
%%
GB_bud_info1 = h5info([c2vsim_path 'Results/C2VSimFG_GW_Budget.hdf']);
ColNames = GB_bud_info1.Groups.Attributes(12).Value(2:end);
%ttt = h5read(GB_bud_info1.Filename, [GB_bud_info1.Name GB_bud_info1.Datasets(1).Name]);
%% Average Period
start_tm_id = find(tm == datetime(2002,10,1));
end_tm_id = find(tm == datetime(2012,9,1));
ndays = sum(eomday(year(tm(start_tm_id:end_tm_id)),month(tm(start_tm_id:end_tm_id))));
%% Budget calculation
% Deep percolation + Gain from stream + boundary inflow - Pumping
Subbasin_ids{1,1} = 1:7;
Subbasin_ids{2,1} = 8:13;
Subbasin_ids{3,1} = 14:21;
for ii = 1:length(Subbasin_ids)
    in_flows = zeros(504,1);
    out_flows = zeros(504,1);
    for jj = 1:length(Subbasin_ids{ii,1})
        for k = 1:22
            if strcmp(GB_bud_info1.Datasets(k).Name, ['Subregion ' num2str(Subbasin_ids{ii,1}(jj)) ' (SR' num2str(Subbasin_ids{ii,1}(jj)) ')'])
                SR_bud = h5read(GB_bud_info1.Filename, [GB_bud_info1.Name GB_bud_info1.Datasets(k).Name]);
                break;
            end
        end
        in_flows = in_flows + SR_bud(4,:)' + SR_bud(5,:)' + SR_bud(8,:)';
        out_flows = out_flows + SR_bud(12,:)';
    end
    ratio(ii,1) = sum(in_flows(start_tm_id:end_tm_id))/sum(out_flows(start_tm_id:end_tm_id));
end

%% Read C2Vsim mesh and calculate the element area and barycenters
c2vsim_elem = shaperead('F:\UCDAVIS\C2VSIM_FG_OR\C2Vsim_FG_v2\wrkspc\C2Vsim_Elements_3310.shp');
elem_cc = zeros(length(c2vsim_elem),2);
elem_area = zeros(length(c2vsim_elem),1);
for ii = 1:length(c2vsim_elem)
    id = c2vsim_elem(ii,1).IE;
    elem_area(ii,1) = polyarea(c2vsim_elem(ii,1).X(1:end-1), c2vsim_elem(ii,1).Y(1:end-1));
    elem_cc(ii,:) = [mean(c2vsim_elem(ii,1).X(1:end-2)) mean(c2vsim_elem(ii,1).Y(1:end-2))];
end
%% Read nodes of the expanded mesh
c2vsim_nodes_expanded = shaperead('F:\UCDAVIS\C2VSIM_FG_OR\C2Vsim_FG_v2\wrkspc\C2VsimMeshORExpanded_3310.shp');

%% Groundwater recharge
%=======================================================

GW_recharge = h5read(GB_bud_info.Filename,...
        [GB_bud_info.Groups(2).Name GB_bud_info.Name GB_bud_info.Groups(2).Datasets(5).Name])';
GW_recharge = GW_recharge(start_tm_id:end_tm_id,:)*0.028316877004293; % ft^3 -> m^3 / Month
GW_recharge = [sum(GW_recharge,1)/ndays]'./elem_area; % m/day
% Createinterpolant
F = scatteredInterpolant(elem_cc(:,1), elem_cc(:,2), GW_recharge,'linear', 'nearest');
rch = F([c2vsim_nodes_expanded.X]', [c2vsim_nodes_expanded.Y]');

writeScatteredData(['inputfiles' filesep 'c2vsim_rch.npsat'],...
    struct('PDIM', 2, 'TYPE', 'HOR', 'MODE', 'SIMPLE'), ...
    [[c2vsim_nodes_expanded.X]', [c2vsim_nodes_expanded.Y]' rch]);
%__________________________________________________________



%% Read stream node budget
%=======================================================

Stream_node_bud_info = h5info([c2vsim_path 'Results/C2VSimFG_Stream_Node_Budget.hdf']);
Stream_col_names = Stream_node_bud_info.Groups.Attributes(11).Value(2:end);
%% Calculate the gain from inside the model on each node
Stream_node_Rch = zeros(length(Stream_node_bud_info.Datasets),1);
for ii = 1:length(Stream_node_bud_info.Datasets)
    ii
    id = textscan(Stream_node_bud_info.Datasets(ii).Name,'NODE%s');
    id = str2double(id{1,1}{1,1});
    tmp = h5read(Stream_node_bud_info.Filename,[Stream_node_bud_info.Name Stream_node_bud_info.Datasets(ii).Name]);
    
    tmp_rch = tmp(7, start_tm_id:end_tm_id)*0.028316877004293; % ft^3 -> m^3 / Month;
    tmp_rch = sum(tmp_rch)/ndays;
    if sum(tmp(8, start_tm_id:end_tm_id)) ~= 0
        tmp_rch1 = tmp(8, start_tm_id:end_tm_id)*0.028316877004293; % ft^3 -> m^3 / Month;
        tmp_rch1 = sum(tmp_rch1)/ndays;
    else
        tmp_rch1 = 0;
    end
    Stream_node_Rch(id,1) = tmp_rch + tmp_rch1;
end
%% Sum the budget for the nodes that coincide
gw_strmid = unique(riv_segm(:,1:2));
rr = riv_segm(:,1:2);
tt = riv_segm(:,3:4);
for ii = 1:length(gw_strmid)
    ind = find(rr == gw_strmid(ii));
    id_coincide = unique(tt(ind));
    if length(id_coincide) > 1
        sumq_tmp = sum(Stream_node_Rch(id_coincide));
        Stream_node_Rch(id_coincide) = sumq_tmp;
        tt(ind) = tt(ind(1));
    end
end
riv_segm(:,3:4) = tt;
%% Distribute the stream node recharge to segments
c2vsim_nodes_3310 = shaperead('F:\UCDAVIS\C2VSIM_FG_OR\C2Vsim_FG_v2\wrkspc\C2VsimNodes_3310.shp');
XY = [[c2vsim_nodes_3310.X]' [c2vsim_nodes_3310.Y]'];
%%
riv_segm(isnan(riv_segm(:,1)),:) =[];
riv_segm_Q_len_width = zeros(size(riv_segm,1),3);
for ii = 1:size(riv_segm,1)
    % find segment length
    segm_len = sqrt(sum((XY(riv_segm(ii,1),:) - XY(riv_segm(ii,2),:)).^2));
    Q1 = -Stream_node_Rch(riv_segm(ii,3));
    Q2 = -Stream_node_Rch(riv_segm(ii,4));
    %check with how many segments each node is connected
    % for Q1
    [II, JJ] = find(riv_segm(:,3:4) == riv_segm(ii,3));
    II(II == ii,:) = [];
    if ~isempty(II)
        dst = sqrt(sum((XY(riv_segm(II,1),:) - XY(riv_segm(II,2),:)).^2,2));
        Q1 = Q1*segm_len/(segm_len + sum(dst));
    end
    
    % for Q2
    [II, JJ] = find(riv_segm(:,3:4) == riv_segm(ii,4));
    II(II == ii,:) = [];
    if ~isempty(II)
        dst = sqrt(sum((XY(riv_segm(II,1),:) - XY(riv_segm(II,2),:)).^2,2));
        Q2 = Q2*segm_len/(segm_len + sum(dst));
    end
    riv_segm_Q_len_width(ii,1:2) = [Q1+Q2 segm_len];
end
%% assign width
for ii = 1:size(riv_segm_Q_len_width,1)
    w = 50;
    while true
        q_temp = riv_segm_Q_len_width(ii,1)/(riv_segm_Q_len_width(ii,2)*w);
        if q_temp < 0.3
            break
        else
            w = w + 50;
            if w > 400
                w = 400;
                break;
            end
        end
    end
    riv_segm_Q_len_width(ii,3) = w;
end
%% Gather data
clear STRM_DATA
STRM_DATA(1,1).p1 = [];
STRM_DATA(1,1).p2 = [];
STRM_DATA(1,1).Q = [];
STRM_DATA(1,1).W = [];
cnt = 1;
for ii = 1:length(riv_segm_Q_len_width)
    if abs(riv_segm_Q_len_width(ii,1)) > 0
        STRM_DATA(cnt,1).p1 = XY(riv_segm(ii,1),:);
        STRM_DATA(cnt,1).p2 = XY(riv_segm(ii,2),:);
        STRM_DATA(cnt,1).Q = riv_segm_Q_len_width(ii,1)/(riv_segm_Q_len_width(ii,2)*riv_segm_Q_len_width(ii,3));
        STRM_DATA(cnt,1).W = riv_segm_Q_len_width(ii,3);
        cnt = cnt + 1;
    end
end

%% Read Small watershed budget
%=======================================================
Wsheds_bud_info = h5info([c2vsim_path 'Results/C2VSimFG_SWatersheds_Budget.hdf']);
Wsheds_col_names = Wsheds_bud_info.Groups.Attributes(12).Value(2:end);
Wsheds_Perc2GW = zeros(length(Wsheds_bud_info.Datasets),1);
for ii = 1:length(Wsheds_bud_info.Datasets)
    ii
    id = textscan(Wsheds_bud_info.Datasets(ii).Name,'WATERSHED%s');
    id = str2double(id{1,1}{1,1});
    tmp = h5read(Wsheds_bud_info.Filename,[Wsheds_bud_info.Name Wsheds_bud_info.Datasets(ii).Name]);
    tmp_rch = tmp(16, start_tm_id:end_tm_id)*0.028316877004293; % ft^3 -> m^3 / Month;
    tmp_rch = sum(tmp_rch)/ndays;
    Wsheds_Perc2GW(id,1) = tmp_rch;
end
%%
Wsheds_Q = zeros(size(XY,1),1);
for ii = 1:length(C2Vsim_SWatersheds)
    Wsheds_Q(C2Vsim_SWatersheds(ii,1).IWB,1) = Wsheds_Q(C2Vsim_SWatersheds(ii,1).IWB,1) + ...
        Wsheds_Perc2GW(ii,1)/length(C2Vsim_SWatersheds(ii,1).IWB);
end
%%
wshed_segm(isnan(wshed_segm(:,1)),:) =[];
Wsheds_segm_Q_len_width = zeros(size(wshed_segm,1),3);
for ii = 1:size(wshed_segm,1)
    segm_len = sqrt(sum((XY(wshed_segm(ii,1),:) - XY(wshed_segm(ii,2),:)).^2));
    Q1 = Wsheds_Q(wshed_segm(ii,1));
    Q2 = Wsheds_Q(wshed_segm(ii,2));
    
    [II, JJ] = find(wshed_segm == wshed_segm(ii,1));
    II(II == ii,:) = [];
    if ~isempty(II)
        dst = sqrt(sum((XY(wshed_segm(II,1),:) - XY(wshed_segm(II,2),:)).^2,2));
        Q1 = Q1*segm_len/(segm_len + sum(dst));
    end
    
    % for Q2
    [II, JJ] = find(wshed_segm == wshed_segm(ii,2));
    II(II == ii,:) = [];
    if ~isempty(II)
        dst = sqrt(sum((XY(wshed_segm(II,1),:) - XY(wshed_segm(II,2),:)).^2,2));
        Q2 = Q2*segm_len/(segm_len + sum(dst));
    end
    Wsheds_segm_Q_len_width(ii,1:2) = [Q1+Q2 segm_len];
    Wsheds_segm_Q_len_width(ii,3) = 50;
end
Wsheds_segm_Q_len_width(:,1) = Wsheds_segm_Q_len_width(:,1)./(Wsheds_segm_Q_len_width(:,2)*Wsheds_segm_Q_len_width(ii,3));
%%
cnt = length(STRM_DATA) + 1;
for ii = 1:length(Wsheds_segm_Q_len_width)
    if abs(Wsheds_segm_Q_len_width(ii,1)) > 0
        STRM_DATA(cnt,1).p1 = XY(wshed_segm(ii,1),:);
        STRM_DATA(cnt,1).p2 = XY(wshed_segm(ii,2),:);
        STRM_DATA(cnt,1).Q = Wsheds_segm_Q_len_width(ii,1)/(Wsheds_segm_Q_len_width(ii,2)*Wsheds_segm_Q_len_width(ii,3));
        STRM_DATA(cnt,1).W = Wsheds_segm_Q_len_width(ii,3);
        cnt = cnt + 1;
    end
end
%% Print file
fid = fopen(['inputfiles' filesep 'c2vsim_Streams.npsat'], 'w');
fprintf(fid, '%d\n',length(STRM_DATA));
for ii = 1:length(STRM_DATA)
   fprintf(fid, '2 %.7f %.1f\n', [STRM_DATA(ii,1).Q STRM_DATA(ii,1).W]);
   fprintf(fid, '%.2f %.2f\n', STRM_DATA(ii,1).p1);
   fprintf(fid, '%.2f %.2f\n', STRM_DATA(ii,1).p2); 
end

fclose(fid);
% 
%__________________________________________________________






%% Read streams - OBSOLETE
%Stream_bud_info = h5info([c2vsim_path 'Results/C2VSimFG_Stream_Budget.hdf']);
%Stream_col_names = Stream_bud_info.Groups.Attributes(11).Value(2:end);
GB_bud_info = h5info([c2vsim_path 'Results/C2VSimFG_GW_ZBudget.hdf']);
Stream_rch_in = h5read(GB_bud_info.Filename,[GB_bud_info.Groups(2).Name GB_bud_info.Name GB_bud_info.Groups(2).Datasets(23).Name])';
Stream_rch_out = h5read(GB_bud_info.Filename,[GB_bud_info.Groups(2).Name GB_bud_info.Name GB_bud_info.Groups(2).Datasets(24).Name])';
Stream_elem_id = h5read(GB_bud_info.Filename,[GB_bud_info.Groups(1).Name GB_bud_info.Name GB_bud_info.Groups(1).Datasets(6).Name])';
Stream_id_names = h5read(GB_bud_info.Filename,[GB_bud_info.Groups(1).Name GB_bud_info.Name GB_bud_info.Groups(1).Datasets(5).Name]);
% Run "make a unique list of segments" sections of the Prepare script
Stream_elem_id = find(Stream_elem_id(3,:) ~= 0)';
msh = h5read(GB_bud_info.Filename,[GB_bud_info.Groups(1).Name GB_bud_info.Name GB_bud_info.Groups(1).Datasets(17).Name])';
%%
% Associate stream elements with stream nodes
riv_segm(isnan(riv_segm(:,1)),:) =[];
stream_nodes = unique(riv_segm);
clear stream_nodes_elem_assoc stream_elem_nodes_assoc
stream_nodes_elem_assoc{length(stream_nodes),1} = [];
stream_elem_nodes_assoc{length(Stream_elem_id),1} = [];
for ii = 1:length(Stream_elem_id)
    for j = 1:4
       if  msh(Stream_elem_id(ii),j) == 0
           continue;
       end
       igw = msh(Stream_elem_id(ii),j);
       id = find(stream_nodes == igw);
       if ~isempty(id)
           stream_nodes_elem_assoc{id,1} = [stream_nodes_elem_assoc{id,1};Stream_elem_id(ii)];
           stream_elem_nodes_assoc{ii,1} = [stream_elem_nodes_assoc{ii,1}; stream_nodes(id)];
       end
    end
end
%% Assign rates to groundwater nodes
for ii = 1:length(stream_nodes)
    % Find out from how many elements this node receives stream water
    n_el = length(stream_nodes_elem_assoc{ii,1});
    for j = 1:n_el
        id = find(Stream_elem_id == stream_nodes_elem_assoc{ii,1}(j));
        % find out in how many nodes the element contributes
        n_nd = length(stream_elem_nodes_assoc{id,1});
        
    end
    
end




%% Read Small watershed budget
Wsheds_bud_info = h5info([c2vsim_path 'Results/C2VSimFG_SWatersheds_Budget.hdf']);
Wsheds_col_names = Wsheds_bud_info.Groups.Attributes(12).Value(2:end);

stream_bud = h5read(Stream_bud_info.Filename, [Stream_bud_info.Name Stream_bud_info.Datasets(1).Name]);

ttt = h5read(GB_bud_info.Filename,[GB_bud_info.Groups(2).Name GB_bud_info.Name GB_bud_info.Groups(2).Datasets(23).Name])';
%%
tmp = h5readatt(GB_bud_info.Filename, ...
    [GB_bud_info.Groups(2).Name GB_bud_info.Name GB_bud_info.Groups(2).Datasets(23).Name],...
    GB_bud_info.Groups(2).Datasets(23).Attributes(1).Name)';
%%
ttt = h5read(GB_bud_info.Filename,[GB_bud_info.Groups(1).Name GB_bud_info.Name GB_bud_info.Groups(1).Datasets(6).Name])';
