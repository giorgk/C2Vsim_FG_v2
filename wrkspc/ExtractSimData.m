%% relative Paths
gis_data = ['..' filesep 'gis_data' filesep];
%% Paths at UCDAVIS LINUX
c2vsim_path = '/media/giorgk/DATA/giorgk/Documents/C2Vsim_FG_v2/C2VSimFG-BETA_PublicRelease/';
mat_data = '/media/giorgk/DATA/giorgk/Documents/C2Vsim_FG_v2/mat_data/';

%% Precipitation
PRC = readC2VsimPrecip([c2vsim_path 'Simulation' filesep 'C2VSimFG_Precip.dat']);
save([mat_data 'PrecipTimeSeries'],'PRC');
%% Read Root zone parameters. They are defined per element
load([mat_data 'C2Vsim_Elements.mat']);
load([mat_data 'PrecipTimeSeries.mat']);
fid = fopen([c2vsim_path 'Simulation' filesep 'RootZone' filesep 'C2VSimFG_RootZone.dat'],'r');
temp = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f', 32537, 'HeaderLines',191);
fclose(fid);
for ii = 1:length(temp{1,1})
    C2Vsim_elem(ii,1).RTZ.WP = temp{1,2}(ii,1);
    C2Vsim_elem(ii,1).RTZ.FC = temp{1,3}(ii,1);
    C2Vsim_elem(ii,1).RTZ.TN = temp{1,4}(ii,1);
    C2Vsim_elem(ii,1).RTZ.LAMBDA = temp{1,5}(ii,1);
    C2Vsim_elem(ii,1).RTZ.K = temp{1,6}(ii,1);
    C2Vsim_elem(ii,1).RTZ.RHC = temp{1,7}(ii,1);
    C2Vsim_elem(ii,1).RTZ.CPRISE = temp{1,8}(ii,1);
    C2Vsim_elem(ii,1).RTZ.IRNE = temp{1,9}(ii,1);
    C2Vsim_elem(ii,1).RTZ.FRNE = temp{1,10}(ii,1);
    C2Vsim_elem(ii,1).RTZ.IMSRC = temp{1,11}(ii,1);
    C2Vsim_elem(ii,1).RTZ.TYPDEST = temp{1,12}(ii,1);
    C2Vsim_elem(ii,1).RTZ.DEST = temp{1,13}(ii,1);
    C2Vsim_elem(ii,1).RTZ.Precip = PRC.ARAIN(:,C2Vsim_elem(ii,1).RTZ.IRNE);
    
end
RTZ_metadata = writeMetaData([],'RTZ.WP','Wilting point');
RTZ_metadata = writeMetaData(RTZ_metadata,'RTZ.FC','Field capacity');
RTZ_metadata = writeMetaData(RTZ_metadata,'RTZ.TN','Total porosity');
RTZ_metadata = writeMetaData(RTZ_metadata,'RTZ.LAMBDA','Pore size distribution index');
RTZ_metadata = writeMetaData(RTZ_metadata,'RTZ.K','Saturated hydraulic conductivity');
RTZ_metadata = writeMetaData(RTZ_metadata,'RTZ.RHC','1 = Campbell equation, 2 = van Genucten-Mualem equation');
RTZ_metadata = writeMetaData(RTZ_metadata,'RTZ.CPRISE', 'Capillary rise');
RTZ_metadata = writeMetaData(RTZ_metadata,'RTZ.IRNE','Precipitation data column in the Precipitation file');
RTZ_metadata = writeMetaData(RTZ_metadata,'RTZ.FRNE','Factor to convert rainfall');
RTZ_metadata = writeMetaData(RTZ_metadata,'RTZ.IMSRC','Generic source of moisture data column');
RTZ_metadata = writeMetaData(RTZ_metadata,'RTZ.TYPDEST','Destination type for the surface flow from element IE');
RTZ_metadata = writeMetaData(RTZ_metadata,'RTZ.DEST','Destination for the surface flow from element IE');
save([mat_data 'C2Vsim_Elements_RTZ_Precip'],'C2Vsim_elem','RTZ_metadata');
%% Read Land use Area Data
NVfield_names = {'ALANDNV','ALANDRV'};
Urbfield_names = {'ALANDU'};
Pondfield_names = {'ALANDRI_FL','ALANDRI_NFL','ALANDRI_NDC','ALANDRF_SL','ALANDRF_PR'};
NonPondfield_names = {'GR','CO','SB','CN','DB','SA','FL','AL','PA','TP','TF','CU','OG','PO','TR','AP','OR','CS','VI','ID'};
%%
NativeVegArea = LU_Area_per_element([c2vsim_path 'Simulation' filesep 'RootZone' filesep 'C2VSimFG_NativeVeg_Area.dat'], 32537, NVfield_names, 104);
PondedCropArea = LU_Area_per_element([c2vsim_path 'Simulation' filesep 'RootZone' filesep 'C2VSimFG_PondedCrop_Area.dat'], 32537, Pondfield_names, 111);
UrbanArea = LU_Area_per_element([c2vsim_path 'Simulation' filesep 'RootZone' filesep 'C2VSimFG_Urban_Area.dat'], 32537, Urbfield_names, 101);
NonPondedCropArea = LU_Area_per_element([c2vsim_path 'Simulation' filesep 'RootZone' filesep 'C2VSimFG_NonPondedCrop_Area.dat'], 32537, NonPondfield_names, 126);
save([mat_data 'C2VsimLU_Area'],'NativeVegArea','PondedCropArea','UrbanArea','NonPondedCropArea');
%
%% Aggregate Land use data per elements
load([mat_data 'C2Vsim_Elements.mat']);
load([mat_data 'C2VsimLU_Area']);
for ii = 1:length(C2Vsim_elem)
    ii
    C2Vsim_elem(ii,1).ElemArea = polyarea(C2Vsim_elem(ii,1).X,C2Vsim_elem(ii,1).Y);
    C2Vsim_elem(ii,1).LUArea = zeros(94,2+5+1+20);
    for jj = 1:length(NVfield_names)
        C2Vsim_elem(ii,1).LUArea(:,jj) = NativeVegArea.(NVfield_names{jj})(:,ii);
    end
    for jj = 1:length(Urbfield_names)
        C2Vsim_elem(ii,1).LUArea(:,jj+2) = UrbanArea.(Urbfield_names{jj})(:,ii);
    end
    for jj = 1:length(Pondfield_names)
        C2Vsim_elem(ii,1).LUArea(:,jj+3) = PondedCropArea.(Pondfield_names{jj})(:,ii);
    end
    for jj = 1:length(NonPondfield_names)
        C2Vsim_elem(ii,1).LUArea(:,jj+8) = NonPondedCropArea.(NonPondfield_names{jj})(:,ii);
    end
end
%%
C2Vsim_elem_LUarea = C2Vsim_elem;
save([mat_data 'C2Vsim_elem_LUarea'],'C2Vsim_elem_LUarea');
%% compare element area with the Land use area
% The ElemArea is square meter while the Land use areas are in ACRE.
% Divide m^3 by 4046.86 to obtain ACRE
for ii = 1:length(C2Vsim_elem)
    CalcArea = C2Vsim_elem(ii,1).ElemArea/4046.86;
    SumArea = sum(C2Vsim_elem(ii,1).LUArea,2);
    maxError(ii,1) = max(abs(CalcArea - SumArea));
end
%% Read ET
% load Element file
load([mat_data 'C2Vsim_Elements.mat']);
%% read the ET data
fid = fopen([c2vsim_path 'Simulation' filesep 'C2VSimFG_ET.dat'],'r');
for ii = 1:136
    temp = fgetl(fid);
end
clear ET
for ii = 1:12
    temp = fgetl(fid);
    C = strsplit(temp,'\t');
    c = textscan(C{1,1},'%f/%f/%f/_%s');
    ET.Time{ii,1} = [num2str(c{1,1}) '/' num2str(c{1,2}) '/' num2str(c{1,3})];
    for jj = 2:length(C)
        ET.data(ii,jj-1) = cell2mat(textscan(C{1,jj},'%f'));
    end
end

fclose(fid);
%% Non Ponded Crops
fid = fopen([c2vsim_path 'Simulation' filesep 'RootZone' filesep 'C2VSimFG_NonPondedCrop.dat'],'r');
% read all comments up to the line with the crop codes
for ii = 1:32767
    temp = fgetl(fid);
end
crops = strsplit(temp,' ');
crops(:,1) = [];
frmt = '%f';
for ii = 1:length(crops)-1
    frmt = [frmt ' %f'];
end
temp = fgetl(fid);
cols = textscan(fid, frmt, 32537);
fclose(fid);
%% Map non ponded crop ET to elements
for ii = 1:length(C2Vsim_elem)
    for jj = 2:length(crops)
        C2Vsim_elem(ii,1).NonPonded.(crops{1,jj}) = ET.data(:,cols{1,jj}(ii)); 
    end
end
%% Ponded Crops
fid = fopen([c2vsim_path 'Simulation' filesep 'RootZone' filesep 'C2VSimFG_PondedCrop.dat'],'r');
for ii = 1:32702
    temp = fgetl(fid);
end
crops = strsplit(temp,' ');
crops(:,1) = [];
frmt = '%f';
for ii = 1:length(crops)-1
    frmt = [frmt ' %f'];
end
temp = fgetl(fid);
cols = textscan(fid, frmt, 32537);
fclose(fid);
%% Map non ponded crop ET to elements
for ii = 1:length(C2Vsim_elem)
    for jj = 2:length(crops)
        C2Vsim_elem(ii,1).Ponded.(crops{1,jj}) = ET.data(:,cols{1,jj}(ii)); 
    end
end
%% Native Vegetation
fid = fopen([c2vsim_path 'Simulation' filesep 'RootZone' filesep 'C2VSimFG_NativeVeg.dat'],'r');
for ii = 1:102
    temp = fgetl(fid);
end
props = strsplit(temp,' ');
props(:,1) = [];
frmt = '%f';
for ii = 1:length(props)-1
    frmt = [frmt ' %f'];
end
temp = fgetl(fid);
data = textscan(fid, frmt, 32537);
fclose(fid);
%% Map native vegetation ET to elements
for ii = 1:length(C2Vsim_elem)
    for jj = [2 3 6]
        C2Vsim_elem(ii,1).NatVeg.(props{1,jj}) = data{1,jj}(ii);
    end
    for jj = [4 5]
        C2Vsim_elem(ii,1).NatVeg.(props{1,jj}) = ET.data(:,data{1,jj}(ii)); 
    end
end
%% Urban 
fid = fopen([c2vsim_path 'Simulation' filesep 'RootZone' filesep 'C2VSimFG_Urban.dat'],'r');
for ii = 1:126
    temp = fgetl(fid);
end
props = strsplit(temp,' ');
props(:,1) = [];
frmt = '%f';
for ii = 1:length(props)-1
    frmt = [frmt ' %f'];
end
temp = fgetl(fid);
data = textscan(fid, frmt, 32537);
fclose(fid);
%% Map Urban ET to elements
for ii = 1:length(C2Vsim_elem)
    for jj = [7]
        C2Vsim_elem(ii,1).Urban.(props{1,jj}) = ET.data(:,data{1,jj}(ii)); 
    end
end
%% save ET per element
C2Vsim_elem_ET = C2Vsim_elem;
save([mat_data 'C2Vsim_elem_ET'], 'C2Vsim_elem_ET');