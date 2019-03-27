%
%% Paths
Bdata_path = '/media/giorgk/DATA/giorgk/Documents/C2Vsim_FG_v2/'; % UCD Ubuntu
Bdata_path = 'D:\giorgk\Documents\C2Vsim_FG_v2\'; % UCD windows
mat_data = [Bdata_path 'mat_data' filesep];
%% load shapefiles
load([mat_data 'ESR_shapefiles'])
load([mat_data 'C2Vsim_Elements'],'C2Vsim_elem');
%% compute the barycenter for the elements
for ii = 1:length(C2Vsim_elem)
    Elems(ii,:) = [ C2Vsim_elem(ii,1).ID ...
                    C2Vsim_elem(ii,1).IRGE ...
                    mean(C2Vsim_elem(ii,1).X) ...
                    mean(C2Vsim_elem(ii,1).Y)];
end
%% Unique list of study areas
for ii = 1:length(Kern_only)
    in = inpolygon(Elems(:,3), Elems(:,4), Kern_only(ii,1).X, Kern_only(ii,1).Y);
    elem_id = Elems(in,1);
    Kern_only(ii,1).Elem_ids = elem_id;
end
%% GROUP DATA using a structure similar to Kings
for ii = 1:length(Kern_only)
    GROUPS(ii,1).name = deblank(Kern_only(ii,1).KernDistri);
    GROUPS(ii,1).EL_ID = Kern_only(ii,1).Elem_ids;
end
%
%% Append precipitation
load([mat_data 'PrecipTimeSeries.mat'])
load([mat_data 'C2Vsim_Elements_RTZ_Precip.mat'])
for ii = 1:length(GROUPS)
    Precip = zeros(1131,1);
    weights = zeros(length(GROUPS(ii,1).EL_ID),1);
    for jj = 1:length(GROUPS(ii,1).EL_ID)
        weights(jj,1) = polyarea(C2Vsim_elem(GROUPS(ii,1).EL_ID(jj,1),1).X, C2Vsim_elem(GROUPS(ii,1).EL_ID(jj,1),1).Y);
        Precip = Precip + C2Vsim_elem(GROUPS(ii,1).EL_ID(jj,1),1).RTZ.Precip * weights(jj,1);
    end
    GROUPS(ii,1).Precip = Precip/sum(weights);
end
%% Write precipitation to excel spreadsheet
clear A
A{1,1} = 'Precipitation [INCH/MONTH]';
A{2,1} = 'Time';
for ii = 1:length(GROUPS)
    A{2,ii+1} = GROUPS(ii,1).name;
    for jj = 1:length(GROUPS(ii,1).Precip)
        if ii == 1
            A{2+jj,1} = PRC.time{jj,1};
        end
        A{2+jj,ii+1} = GROUPS(ii,1).Precip(jj,1);
    end
end
[status,message] = xlswrite('KernSubregionsData.xlsx',A,'Precipitation','A1');
%%
load([mat_data 'C2Vsim_elem_LUarea.mat']);
for ii = 1:length(GROUPS)
    GROUPS(ii,1).LU_area = zeros(94,28);
    for jj = 1:length(GROUPS(ii,1).EL_ID)
        GROUPS(ii,1).LU_area = GROUPS(ii,1).LU_area + C2Vsim_elem_LUarea(GROUPS(ii,1).EL_ID(jj,1),1).LUArea;
    end
    GROUPS(ii,1).LU_group = zeros(94,3);
    GROUPS(ii,1).LU_group(:,1) = sum(GROUPS(ii,1).LU_area(:,1:2),2);%Native vegetation
    GROUPS(ii,1).LU_group(:,2) = sum(GROUPS(ii,1).LU_area(:,3),2);%Urban
    GROUPS(ii,1).LU_group(:,3) = sum(GROUPS(ii,1).LU_area(:,4:end),2);%Urban
end
%% Write data to excel
load([mat_data 'C2VsimLU_Area'],'UrbanArea');% load to use just the dates

A = ['Time' 'Native & Riparian', 'Urban', 'Agricultural' NVfield_names Urbfield_names Pondfield_names NonPondfield_names];
for ii = 1:94
    A{ii+1,1} = UrbanArea.time{ii};
end
for ii = 1:length(GROUPS)
    for jj = 1:94
        for kk = 1:3
            A{jj+1,kk+1} = GROUPS(ii,1).LU_group(jj,kk);
        end
        
        for kk = 1:28
            A{jj+1,kk+4} = GROUPS(ii,1).LU_area(jj,kk);
        end
    end
    xlswrite('KernSubregionsData.xlsx',A,GROUPS(ii,1).name,'A1');
end
%% Aggregate Input ET
% The input ET is monthly step that are repeated each year
load([Bdata_path 'mat_data' filesep 'C2Vsim_elem_ET.mat'])
fieldNames = {'NonPonded', 'Ponded', 'NatVeg', 'Urban'};
for ii = 1:length(GROUPS)
    for jj = 1:length(fieldNames)
        cat = fields(C2Vsim_elem_ET(1,1).(fieldNames{jj}));
        for kk = 1:length(cat)
            weight = zeros(length(GROUPS(ii,1).EL_ID),1);
            tempET = zeros(12,1);
            for mm = 1:length(GROUPS(ii,1).EL_ID)
                weight(mm,1) = polyarea(C2Vsim_elem_ET(GROUPS(ii,1).EL_ID(mm,1),1).X, C2Vsim_elem_ET(GROUPS(ii,1).EL_ID(mm,1),1).Y);
                tempET = tempET + C2Vsim_elem_ET(GROUPS(ii,1).EL_ID(mm,1),1).(fieldNames{jj}).(cat{kk}) * weight(mm,1);
            end
            GROUPS(ii,1).(fieldNames{jj}).(cat{kk}) = tempET/sum(weight);
        end
    end
end
%%
clear A
ET_times = {...
    '10/31/4000'...
    '11/30/4000'...
    '12/31/4000'...
    '01/31/4000'...
    '02/29/4000'...
    '03/31/4000'...
    '04/30/4000'...
    '05/31/4000'...
    '06/30/4000'...
    '07/31/4000'...
    '08/31/4000'...
    '09/30/4000'};
A{1,1} = 'Time';
cnt = 2;
for ii = 1:length(fieldNames)
    cat = fields(GROUPS(1,1).(fieldNames{ii}));
    for jj = 1:length(cat)
        A{1,cnt} = [cat{jj} '_' fieldNames{ii}];
        cnt = cnt + 1;
    end
end
for ii = 1:length(ET_times)
    A{1+ii,1} = ET_times{ii};
end
for ii = 1:length(GROUPS)
    cnt_c = 2;
    for jj = 1:length(fieldNames)
        cat = fields(GROUPS(1,1).(fieldNames{jj}));
        for kk = 1:length(cat)
            for mm = 1:length(GROUPS(ii,1).(fieldNames{jj}).(cat{kk}))
                A{mm+1,cnt_c} = GROUPS(ii,1).(fieldNames{jj}).(cat{kk})(mm);
            end
            cnt_c = cnt_c + 1;
        end
    end
    xlswrite('KernInputET.xlsx', A, GROUPS(ii,1).name, 'A1');
end
%% Aggregate Output ET
c2vsim_path = 'D:\giorgk\Documents\C2Vsim_FG_v2\C2VSimFG-BETA_PublicRelease\'; % UCD Windows
c2vsim_path = '/media/giorgk/DATA/giorgk/Documents/C2Vsim_FG_v2/C2VSimFG-BETA_PublicRelease/'; % UCD Ubuntu
%%
GBinfo = h5info([c2vsim_path 'Results' filesep 'C2VSimFG_RZ_ZBudget.hdf']);
load([Bdata_path 'mat_data' filesep 'C2Vsim_Elements.mat'])

NatRipPET = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(12).Name])';
NonPondPET = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(28).Name])';
RefugePET = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(46).Name])';
RicePET = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(64).Name])';
UrbanPET = h5read(GBinfo.Filename, ...
    [GBinfo.Groups(2).Name GBinfo.Name GBinfo.Groups(2).Datasets(81).Name])';
%% 
for ii = 1:length(GROUPS)
    natrip = zeros(504,1);
    nonpond = zeros(504,1);
    refuge = zeros(504,1);
    rice = zeros(504,1);
    urban = zeros(504,1);
    weights = zeros(length(GROUPS(ii,1).EL_ID),1);
    for jj = 1:length(GROUPS(ii,1).EL_ID)
        weights(jj,1) = polyarea(C2Vsim_elem(GROUPS(ii,1).EL_ID(jj),1).X, C2Vsim_elem(GROUPS(ii,1).EL_ID(jj),1).Y);
        natrip = natrip + NatRipPET(:, GROUPS(ii,1).EL_ID(jj));
        nonpond = nonpond + NonPondPET(:, GROUPS(ii,1).EL_ID(jj));
        refuge = refuge + RefugePET(:, GROUPS(ii,1).EL_ID(jj));
        rice = rice + RicePET(:, GROUPS(ii,1).EL_ID(jj));
        urban = urban + UrbanPET(:, GROUPS(ii,1).EL_ID(jj));
    end
    GROUPS(ii,1).NativRipar = natrip./sum(weights);
    GROUPS(ii,1).NonPonded = nonpond./sum(weights);
    GROUPS(ii,1).Refuge = refuge./sum(weights);
    GROUPS(ii,1).Rice = rice./sum(weights);
    GROUPS(ii,1).Urban = urban./sum(weights);
end

A = {'Time', 'NativeRiparian', 'NonPonded', 'Refuge', 'Rice', 'Urban'} ;


yy = 1973;mm = 10;
for ii = 1:504
    A{ii+1,1} = datestr(datevec(datetime(yy,mm,eomday(yy,mm))),'mm/yyyy');
    mm = mm + 1;
    if mm > 12
        mm = 1;
        yy = yy+1;
    end
end


for ii = 1:length(GROUPS)
    for jj = 1:504
        A{jj+1,2} = GROUPS(ii,1).NativRipar(jj);
        A{jj+1,3} = GROUPS(ii,1).NonPonded(jj);
        A{jj+1,4} = GROUPS(ii,1).Refuge(jj);
        A{jj+1,5} = GROUPS(ii,1).Rice(jj);
        A{jj+1,6} = GROUPS(ii,1).Urban(jj);
    end
    xlswrite('KingsOutputPET.xlsx', A, GROUPS(ii,1).name, 'A1');
end

