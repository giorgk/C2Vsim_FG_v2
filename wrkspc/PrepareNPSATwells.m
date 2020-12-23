%% Ratios
% Scenario 0
scenName = 'bud0';
ratio = [0.8231; 0.7437; 0.6960];
%% C2Vsim main Path
c2vsim_path = ['..' filesep 'c2vsimfg_beta2_publicrelease' filesep 'C2VSimFG_BETA2_PublicRelease' filesep];
%% Set time stamp
start_date = datetime(1973,10,1);
tm = start_date + calmonths(0:503);
%% Average Period
start_tm_id = find(tm == datetime(2002,10,1));
end_tm_id = find(tm == datetime(2012,9,1));
ndays = sum(eomday(year(tm(start_tm_id:end_tm_id)),month(tm(start_tm_id:end_tm_id))));
%% Load required data
welldata = shaperead(['..' filesep '..' filesep '..' filesep 'CVHM_DATA' filesep 'WellData' filesep 'CVwells_30y_3310']);
Subregions = shaperead('C2Vsim_Subregions_3310');
Subregions_withoutRivers = shaperead('SubregionsWithoutRivers_3310');
C2Vsim_elements = shaperead('C2Vsim_Elements_3310');
%% Read GW pumping Info
PumpOut = readC2VsimPumpOut([c2vsim_path filesep 'Results' filesep 'C2VSimFG_Pumping.out'], []);
% Read the well coordinates
fid = fopen([c2vsim_path filesep 'Simulation' filesep 'Groundwater' filesep 'C2VSimFG_WellSpec.dat'],'r');
C = textscan(fid, '%f %f %f %f %f %f %s',610, 'HeaderLines', 95, 'delimiter', '\n');
fclose(fid);
c2vsimWells = [C{1,2} C{1,3} zeros(length(C{1,2}),1)];
% read elements in original coordinate system
elem_26910 = shaperead(['..' filesep 'gis_data' filesep 'C2Vsim_Elements']);
for ii = 1:length(elem_26910)
    c2vsimWells(inpolygon(c2vsimWells(:,1), c2vsimWells(:,2), elem_26910(ii,1).X, elem_26910(ii,1).Y),3) = ii;
end
clear elem_26910

ElemPumpAg = PumpOut.DATA(PumpOut.Type == 11,:);
ElemPumpUrb = PumpOut.DATA(PumpOut.Type == 12,:);
WellPumpAg = PumpOut.DATA(PumpOut.Type == 21,:);
WellPumpUrb = PumpOut.DATA(PumpOut.Type == 22,:);

%GB_bud_info = h5info([c2vsim_path 'Results/C2VSimFG_GW_ZBudget.hdf']);
%Nelem = double(GB_bud_info.Groups(1).Attributes(17).Value);
% This is how to obtain the IDS and the column names. However this is not
% used here
%ColNames = h5read(GB_bud_info.Filename,[GB_bud_info.Groups(1).Name GB_bud_info.Name GB_bud_info.Groups(1).Datasets(5).Name]);
%h5read(GB_bud_info.Filename,[GB_bud_info.Groups(1).Name GB_bud_info.Name GB_bud_info.Groups(1).Datasets(6).Name])';
%% Main Loop

WELLS = nan(25000,6); % X Y D SL Q type
cnt_well = 0;
do_plot = true;
for isub = 21:21%1:length(Subregions)
    cnt_surgeg_wells = 0;
    id_sub = Subregions(isub,1).IRGE;
    disp(['Subregion: ' num2str(id_sub) ' [' num2str(isub) ']']);
    id_elem_in_sub = find([C2Vsim_elements.IRGE]' == id_sub);
    
    [C, ia, ib] = intersect(id_elem_in_sub, c2vsimWells(:,3));
    
    % find the total Ag and Urb pumping for the subregion
    SubregAgPump = sum(sum(ElemPumpAg(id_elem_in_sub, start_tm_id:end_tm_id)))*1233.48/ndays + ...
        sum(sum(WellPumpAg(ib, start_tm_id:end_tm_id)))*1233.48/ndays;
    SubregUrbPump = sum(sum(ElemPumpUrb(id_elem_in_sub, start_tm_id:end_tm_id)))*1233.48/ndays + ...
        sum(sum(WellPumpUrb(ib, start_tm_id:end_tm_id)))*1233.48/ndays;
    % Scale pumping accordin to the ratio
    if id_sub <= 7
        SubregAgPump = SubregAgPump*ratio(1);
        SubregUrbPump = SubregUrbPump*ratio(1);
    elseif id_sub >= 8 && id_sub <=13
        SubregAgPump = SubregAgPump*ratio(2);
        SubregUrbPump = SubregUrbPump*ratio(2);
    else
        SubregAgPump = SubregAgPump*ratio(3);
        SubregUrbPump = SubregUrbPump*ratio(3);
    end
    
    % find the sample wells that are inside the current subregion
    [Xs, Ys] = polysplit(Subregions(isub,1).X, Subregions(isub,1).Y);
    in_sample_wells = false(length(welldata),1);
    for k = 1:length(Xs)
        if ispolycw(Xs{k,1}, Ys{k,1})
            in_sample_wells(inpolygon([welldata.X]', [welldata.Y]', Xs{k,1}, Ys{k,1})) = true;
        end
    end
    in_sample_wells = find(in_sample_wells);
    urbWelldata = welldata(in_sample_wells,1);
    urbWelldata = urbWelldata([urbWelldata.type]' == 1,1);
    AgWelldata = welldata(in_sample_wells,1);
    AgWelldata = AgWelldata([AgWelldata.type]' == 0,1);
    
    Ag_loc_pdf = Calc_2D_PDF([AgWelldata.X]',[AgWelldata.Y]', 15);
    Pb_loc_pdf = Calc_2D_PDF([urbWelldata.X]',[urbWelldata.Y]', 15);
    
    % plot the subregion boundaries and the distribution of wells
    % plot3(Subregions(isub,1).X, Subregions(isub,1).Y,2*ones(length(Subregions(isub,1).X),1),'k')
    if do_plot
        figure(1); clf
        figure(1);plot(Subregions(isub,1).X, Subregions(isub,1).Y)
        figure(1); hold on
        figure(1);title(['Subregion ' num2str(id_sub)]);
        figure(1);plot([urbWelldata.X]',[urbWelldata.Y]','or');
        figure(1);plot([AgWelldata.X]',[AgWelldata.Y]','.b');
        figure(1); axis equal

        % plot the PDF for the Ag well locations
        figure(2); clf
        figure(2);surf(Ag_loc_pdf.X, Ag_loc_pdf.Y, Ag_loc_pdf.V, 'edgecolor', 'none')
        figure(2);hold on
        figure(2); plot([AgWelldata.X]',[AgWelldata.Y]','+r')
        figure(2); view(0,90)
        figure(2); alpha(0.5)
        figure(2); axis equal
        drawnow

        % plot the PDF for the Urban wells locations
        figure(3); clf
        figure(3); surf(Pb_loc_pdf.X, Pb_loc_pdf.Y, Pb_loc_pdf.V, 'edgecolor', 'none')
        figure(3); hold on
        figure(3); plot([urbWelldata.X]',[urbWelldata.Y]','+r')
        figure(3); view(0,90)
        figure(3); alpha(0.5)
        figure(3); axis equal
        drawnow
    end
    
    Qag = [AgWelldata.Q]';
    Qpb = [urbWelldata.Q]';
    Qag(isnan(Qag) | Qag <= 0,:) = [];
    Qpb(isnan(Qpb) | Qpb <= 0,:) = [];
    % In some farms there are very few public supply wells
    if length(Qpb) < 5
        Qpb = [Qpb;Qag];
    end
    %we have to sample from a log10 distribution
    Qag = log10(Qag);
    Qpb = log10(Qpb);
    
    % Pumping - depth distribution
    Qtmp = [welldata(in_sample_wells,1).Q]';
    Dtmp = [welldata(in_sample_wells,1).depth]';
    test = find(~isnan(Qtmp) & ~isnan(Dtmp) & Qtmp > 0 & Dtmp > 0);
    Qtmp = Qtmp(test);
    Dtmp = Dtmp(test);
    QD_pdf = Calc_2D_PDF(log10(Qtmp), log10(Dtmp), 15);
    
    if do_plot
        figure(5); clf
        figure(5);surf(QD_pdf.X, QD_pdf.Y, QD_pdf.V, 'edgecolor','none')
        figure(5);hold on
        figure(5); plot3(log10(Qtmp), log10(Dtmp), 1.2*ones(length(Dtmp),1),'+r')
        figure(5); view(0,90)
        figure(5); alpha(0.5)
        figure(5); axis equal
        drawnow
    end
    
    % Depth - screen length
    Dtmp = [welldata(in_sample_wells,1).depth]';
    SLtmp = [welldata(in_sample_wells,1).bot]' - [welldata(in_sample_wells,1).top]';
    test = find(~isnan(Dtmp) & ~isnan(SLtmp) & Dtmp > 0 & SLtmp > 0);
    Dtmp = Dtmp(test);
    SLtmp = SLtmp(test);
    DS_pdf = Calc_2D_PDF(log10(Dtmp), log10(SLtmp), 15);
    
    if do_plot
        figure(6); clf
        figure(6);surf(DS_pdf.X, DS_pdf.Y, DS_pdf.V, 'edgecolor','none')
        figure(6);hold on
        figure(6); plot3(log10(Dtmp), log10(SLtmp),1.2*ones(length(SLtmp),1),'+r')
        figure(6); view(0,90)
        figure(6); alpha(0.5)
        figure(6); axis equal
        drawnow
    end
    
    Qag_subreg = 0;
    Qpb_subreg = 0;
    
    % The data in Rich well data set are in gpm while the units of C2Vsim
    % have been converted to m^3/day. 
    % To convert the pumping to m^3/day Q = Q*5.450992992.
    % However this correpond to the design pumping. The argument is that
    % each year only six months operate with this rate therefore we divide
    % the pumping by 6 months Q = Q*5.451/6.
    % As we dont want rates lower than 400 m^3/day this would be
    % log10(400*(6/5.451)) = 2.6437 gpm.
    % 10^2.6437*5.451/6
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AG WELLS
    BBox = Subregions(isub,1).BoundingBox;
    [Xs, Ys] = polysplit(Subregions_withoutRivers(isub,1).X, Subregions_withoutRivers(isub,1).Y); 
    exit_this_loop = false;
    [Qag_f, Qag_x] = ecdf(Qag);
    while Qag_subreg < abs(SubregAgPump) 
        % generate a valid well location
        % print for debug purpose
        %if true
        %    figure(10);clf
        %    figure(10);hold on
        %    for iii = 1:length(Xs)
        %        figure(10);plot(Xs{iii,1},Ys{iii,1})
        %    end
        %end
        while true
            xr = BBox(1,1) + (BBox(2,1) - BBox(1,1))*rand;
            yr = BBox(1,2) + (BBox(2,2) - BBox(1,2))*rand;
            %figure(10);plot(xr,yr,'xk')
            point_in = false;
            for jj = 1:length(Xs)
                is_in = inpolygon(xr, yr, Xs{jj,1}, Ys{jj,1});
                if is_in
                    if ispolycw(Xs{jj,1}, Ys{jj,1})
                        point_in = true;
                    else
                        point_in = false;
                        break;
                    end
                end
            end
            if ~point_in; continue; end
            %check if the generated point is closer to a previously generated well
            if cnt_well > 0
                if min(sqrt((WELLS(1:cnt_well,1) - xr).^2 + (WELLS(1:cnt_well,2) - yr).^2)) < 400; continue;end
            end
        
            % accept the point with certain probability
            r_accept = rand;
            if r_accept < Ag_loc_pdf.F(xr,yr)
                break
            end
        end
        
        % assign a random Pumping with minimum 400 m^3/day
        while true
            Qr = interp1(Qag_f,Qag_x,rand);
            if Qr > 2.6437
                break;
            end
        end
        %Qr = Qag(Qag > 2.6437);
        %Qr = Qr(randperm(length(Qr),1));
        if abs(SubregAgPump) - (Qag_subreg + (10^Qr)*5.450992992/6) < 400
            Qr = log10((abs(SubregAgPump) - Qag_subreg)*6/5.450992992);
            exit_this_loop = true;
        end
        
        % assign Depth and screen length
        [D, SL, out] = AssignDS_v2(QD_pdf,DS_pdf, Qr);
        if ~out; continue; end
        cnt_well = cnt_well + 1;
        cnt_surgeg_wells = cnt_surgeg_wells + 1;
        WELLS(cnt_well,:) = [xr yr 10^Qr*5.450992992/6 10^D, 10^SL, 1];
        
        Qag_subreg = Qag_subreg + 10^Qr*5.450992992/6;
        
        if do_plot
            figure(1);title(['Subregion ' num2str(id_sub) ', Allwells#: ' num2str(cnt_well) ' Subreg# ' num2str(cnt_surgeg_wells)]);
            figure(2);plot(xr,yr,'.b','MarkerSize',13)
            figure(2);title([num2str(100*Qag_subreg/abs(SubregAgPump)) '%  of ' num2str(abs(SubregAgPump))])
            figure(5);plot(Qr,D,'.b','MarkerSize',13)
            figure(6);plot(D,SL,'.b','MarkerSize',13)
            drawnow
        end
        
        if exit_this_loop
            break;
        end
    end
    disp(['Ag: ' num2str(cnt_surgeg_wells),]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % URBAN WELLS
    exit_this_loop = false;
    [Qpb_f, Qpb_x] = ecdf(Qpb);
    while Qpb_subreg < abs(SubregUrbPump) 
        % generate a valid well location
        while true
            xr = BBox(1,1) + (BBox(2,1) - BBox(1,1))*rand;
            yr = BBox(1,2) + (BBox(2,2) - BBox(1,2))*rand;
            point_in = false;
            for jj = 1:length(Xs)
                is_in = inpolygon(xr, yr, Xs{jj,1}, Ys{jj,1});
                if is_in
                    if ispolycw(Xs{jj,1}, Ys{jj,1})
                        point_in = true;
                    else
                        point_in = false;
                        break;
                    end
                end
            end
            if ~point_in; continue; end
            %check if the generated point is closer to a previously generated well
            if cnt_well > 0
                if min(sqrt((WELLS(1:cnt_well,1) - xr).^2 + (WELLS(1:cnt_well,2) - yr).^2)) < 400; continue;end
            end
        
            % accept the point with certain probability
            r_accept = rand;
            if r_accept < Pb_loc_pdf.F(xr,yr)
                break
            end
        end
        
        % assign a random Pumping with minimum 400 m^3/day
        
        while true
            Qr = interp1(Qpb_f,Qpb_x,rand);
            if Qr > 2.6437
                break;
            end
        end
        %Qr = Qpb(Qpb > 2.6437);
        %Qr = Qr(randperm(length(Qr),1));
        if abs(SubregUrbPump) - (Qpb_subreg + (10^Qr)*5.450992992/6) < 400
            Qr = log10((abs(SubregUrbPump) - Qpb_subreg)*6/5.450992992);
            exit_this_loop = true;
        end
        
        % assign Depth and screen length
        [D, SL, out] = AssignDS_v2(QD_pdf,DS_pdf, Qr);
        if ~out; continue; end
        cnt_well = cnt_well + 1;
        cnt_surgeg_wells = cnt_surgeg_wells + 1;
        WELLS(cnt_well,:) = [xr yr 10^Qr*5.450992992/6 10^D, 10^SL, 0];
        
        Qpb_subreg = Qpb_subreg + 10^Qr*5.450992992/6;
        
        if do_plot
            figure(1);title(['Subregion ' num2str(id_sub) ', Allwells#: ' num2str(cnt_well) ' Subreg# ' num2str(cnt_surgeg_wells)]);
            figure(3);plot(xr,yr,'.g','MarkerSize',13)
            figure(3);title([num2str(100*Qpb_subreg/abs(SubregUrbPump)) '%  of ' num2str(abs(SubregUrbPump))])
            figure(5);plot(Qr,D,'.g','MarkerSize',13)
            figure(6);plot(D,SL,'.g','MarkerSize',13)
            drawnow
        end
        
        if exit_this_loop
            break;
        end
    end
    disp(['Ag + Urb: ' num2str(cnt_surgeg_wells),]);
    disp(['Total: ' num2str(cnt_well),]);
end
%% Remove nans
WELLS(cnt_well+1:end,:) = [];
save(['WELLS_QDS_' scenName], 'WELLS'); 
%% STATS
sum(sum(ElemPumpAg(:, start_tm_id:end_tm_id)))*1233.48/ndays + ...
    sum(sum(WellPumpAg(:, start_tm_id:end_tm_id)))*1233.48/ndays + ...
    sum(sum(ElemPumpUrb(:, start_tm_id:end_tm_id)))*1233.48/ndays + ...
    sum(sum(WellPumpUrb(:, start_tm_id:end_tm_id)))*1233.48/ndays
%% Fit wells between top and bottom
% Get data
top = read_Scattered(['inputfiles' filesep 'c2vsimTop_init.npsat'],2);
bot = read_Scattered(['inputfiles' filesep 'c2vsimBottom_initModif_3310.npsat'],2);
Ftop = scatteredInterpolant(top.p(:,1), top.p(:,2), top.v, 'linear', 'nearest');
Fbot = scatteredInterpolant(bot.p(:,1), bot.p(:,2), bot.v, 'linear', 'nearest');
load(['WELLS_QDS_' scenName])
%% Do fitting
lower_thres = 10;
upper_thres = 10;
clear wellstbQ
wellstbQ = nan(size(WELLS,1),5);
for ii = 1:size(WELLS,1)
    xw = WELLS(ii,1);
    yw = WELLS(ii,2);
    wellstbQ(ii,1:2) = [xw yw];
    tw = Ftop(xw, yw);
    bw = Fbot(xw, yw);
    total_length = tw - bw;
    depth = WELLS(ii,4)*0.3048;
    sl = WELLS(ii,5)*0.3048;
    remain_dist = total_length - (depth + lower_thres + sl);
    if remain_dist < 0
        remain_dist = abs(remain_dist);
        udepth = depth/(depth+sl);
        usl = sl/(depth+sl);
        depth = depth - remain_dist*udepth;
        if depth < 0
            warning(["well " ii " has negative depth" ])
        end
        
        sl = sl - remain_dist*usl;
        if depth < 0
            warning(["well " ii " has negative screen length" ])
        end
        
    end
    wellstbQ(ii, 3) = tw - depth;
    wellstbQ(ii, 4) = wellstbQ(ii, 3) - sl;
    wellstbQ(ii, 5) = -WELLS(ii,3);
end
%% Write wells
fid = fopen(['inputfiles' filesep 'wellsModif' scenName 'init.npsat'],'w');
fprintf(fid, '%d\n', size(wellstbQ,1));
fprintf(fid, '%.2f %.2f %.2f %.2f %.2f\n', wellstbQ');
fclose(fid);