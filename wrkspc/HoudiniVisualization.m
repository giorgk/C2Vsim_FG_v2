%% paths
mat_data = '..\mat_data\';
hou_data = '..\hou_data\';
mat_data = '/media/giorgk/DATA/giorgk/Documents/C2Vsim_FG_v2/mat_data/';
hou_data = '/media/giorgk/DATA/giorgk/Documents/C2Vsim_FG_v2/hou_data/';
%% print the nodes with elevations and the mesh
load([mat_data 'C2Vsim_Nodes.mat']);
load([mat_data 'C2Vsim_Elements.mat']);
AA = zeros(length(C2Vsim_nodes), 3+4*2);
% in houdini the z axis is the y
for ii = 1:length(C2Vsim_nodes)
    AA(ii,1) = C2Vsim_nodes(ii,1).X;
    AA(ii,2) = C2Vsim_nodes(ii,1).GSE*0.3048;
    AA(ii,3) = C2Vsim_nodes(ii,1).Y;
    cc = 4:2:10;
    EL = AA(ii,2);
    for jj = 1:4
        a = C2Vsim_nodes(ii,1).A(jj)*0.3048;
        l = C2Vsim_nodes(ii,1).L(jj)*0.3048;
        if a == 0
            a = 0.1;
        end
        ELa = EL - a;
        ELl = ELa - l;
        AA(ii,cc(jj):cc(jj)+1) = [ELa ELl];
        EL = ELl;
    end
end
%% 
MSH = zeros(length(C2Vsim_elem),4);
for ii = 1:length(C2Vsim_elem)
    MSH(ii,1:length(C2Vsim_elem(ii,1).ND_ID)) = C2Vsim_elem(ii,1).ND_ID;
end
   %% 
fid = fopen([hou_data 'C2VsimMeshElev.txt'],'w');
fprintf(fid, '%d %d\n', [length(C2Vsim_nodes) length(C2Vsim_elem)]);
fprintf(fid, '%f %f %f %f %f %f %f %f %f %f %f\n',AA');
fprintf(fid, '%d %d %d %d\n', MSH');
fclose(fid);
%% print land use history
load([mat_data 'C2Vsim_elem_LUarea.mat']);
for ii = 1:size(C2Vsim_elem_LUarea(1,1).LUArea,1)
    ii
    AA = zeros(length(C2Vsim_elem_LUarea),4);
    
    for jj = 1:length(C2Vsim_elem_LUarea)
        AA(jj,:) = [jj sum(C2Vsim_elem_LUarea(jj,1).LUArea(ii,1:2)) C2Vsim_elem_LUarea(jj,1).LUArea(ii,3) sum(C2Vsim_elem_LUarea(jj,1).LUArea(ii,4:end))];
        AA(jj,2:4) = AA(jj,2:4)/sum(AA(jj,2:4));
    end
    fid = fopen([hou_data 'LandUseHistory\LUdist_' num2str(ii) '.txt'],'w');
    fprintf(fid, '%d %.5f %.5f %.5f\n', AA');
    fclose(fid);
end
%% load hydraulic head data
load([mat_data 'C2VsimHead.mat']);
for ii = 1:size(C2VsimHead,1)
    ii
   fid = fopen([hou_data 'waterTable/lay1_' num2str(ii, '%04d') '.dat'],'w'); 
   fprintf(fid, '%f\n', C2VsimHead{ii,2}(:,1).*0.3048);
   fclose(fid);
end
%% Deep percolation
% Read Deep percolation first then write it to houdini
yy = 1973;
mm = 9;
for ii = 1:505
    ii
    if ii == 1
        DP = zeros(size(DeepPerc,2),1);
    else
        DP = DeepPerc(ii-1,:)*1233.48/eomday(yy, mm);  %ACFT-> *1233.48 m^3 / days -> m^3/day
    end
    
   fid = fopen([hou_data 'DeepPerc/DeepPerc_' num2str(ii, '%04d') '.dat'],'w'); 
   fprintf(fid, '%f\n', DP);
   fclose(fid);
   mm = mm +1;
   if mm > 12
       yy = yy + 1;
       mm = 1;
   end
end
%% Write time file
yy = 1973;
mm = 9;
fid = fopen([hou_data 'TimeInfo.dat'],'w');
for ii = 1:505
    fprintf(fid, '%d/%d\n',[mm yy]);
    mm = mm +1;
   if mm > 12
       yy = yy + 1;
       mm = 1;
   end
end
fclose(fid);