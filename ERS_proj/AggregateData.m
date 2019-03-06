%% read the grouping from the csv file
%
load('KingsSubregions');
%% unique groups
temp = unique(Subregion);
%% find the element ids for each group
for ii = 1:length(temp)
   GROUPS(ii,1).name = deblank(temp{ii,1});
   GROUPS(ii,1).EL_ID = [];
   for jj = 1:length(Subregion)
       if strcmp(deblank(Subregion{jj,1}), GROUPS(ii,1).name)
           GROUPS(ii,1).EL_ID = [GROUPS(ii,1).EL_ID; ElementID(jj,1)];
       end
   end
end
%% Append precipitation
load('F:\UCDAVIS\C2VSIM_FG_OR\mat_data\PrecipTimeSeries.mat')
load('F:\UCDAVIS\C2VSIM_FG_OR\mat_data\C2Vsim_Elements_RTZ_Precip.mat')
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
        if ii == 1;
            A{2+jj,1} = PRC.time{jj,1};
        end
        A{2+jj,ii+1} = GROUPS(ii,1).Precip(jj,1);
    end
end
[status,message] = xlswrite('KingsSubregionsData.xlsx',A,'Precipitation','A1');
%%
load('F:\UCDAVIS\C2VSIM_FG_OR\mat_data\C2Vsim_elem_LUarea.mat');
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
%}
%% Write data to excel
%load('F:\UCDAVIS\C2VSIM_FG_OR\mat_data\C2VsimLU_Area','UrbanArea');% load to use just the dates

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
    xlswrite('KingsSubregionsData.xlsx',A,GROUPS(ii,1).name,'A1');
end
