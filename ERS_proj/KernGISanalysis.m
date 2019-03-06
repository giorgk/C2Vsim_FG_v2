%
%% Paths
mat_data = '/media/giorgk/DATA/giorgk/Documents/C2Vsim_FG_v2/mat_data/';
%% load shapefil
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
    for jj = 1:length(GROUPS(ii,1).EL_ID)
        Precip = Precip + C2Vsim_elem(GROUPS(ii,1).EL_ID(jj,1),1).RTZ.Precip;
    end
    GROUPS(ii,1).Precip = Precip;
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