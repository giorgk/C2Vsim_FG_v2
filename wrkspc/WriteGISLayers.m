%% Under octave load the required packages
pkg load mapping
%% Paths
mat_data = '/media/giorgk/DATA/giorgk/Documents/C2Vsim_FG_v2/mat_data/';
%% load node file
load([mat_data 'C2Vsim_Nodes.mat']);
%% Write it as shapefile
for ii = 1:length(C2Vsim_nodes)
    S(ii,1).Geometry = 'Point';
    S(ii,1).X = C2Vsim_nodes(ii,1).X;
    S(ii,1).Y = C2Vsim_nodes(ii,1).Y;
    S(ii,1).GSE = C2Vsim_nodes(ii,1).GSE;
    for jj = 1:4
        S(ii,1).(['A' num2str(jj)]) = C2Vsim_nodes(ii,1).A(jj);
        S(ii,1).(['L' num2str(jj)]) = C2Vsim_nodes(ii,1).L(jj);
    end
end
%%
C2Vsim_Nodes_shp = S;
%%
save([mat_data 'C2Vsim_Nodes.mat'],'C2Vsim_Nodes_shp', '-append')
%% load mesh file
load([mat_data 'C2Vsim_Elements.mat']);
%% Write it as shapefile
for ii = 1:length(C2Vsim_elem)
    S(ii,1).Geometry = 'Polygon';
    S(ii,1).BoundingBox = [min(C2Vsim_elem(ii,1).X) min(C2Vsim_elem(ii,1).Y);max(C2Vsim_elem(ii,1).X) max(C2Vsim_elem(ii,1).Y)];
    S(ii,1).X = [C2Vsim_elem(ii,1).X C2Vsim_elem(ii,1).X(1) nan];
    S(ii,1).Y = [C2Vsim_elem(ii,1).Y C2Vsim_elem(ii,1).Y(1) nan];
    S(ii,1).IRGE = C2Vsim_elem(ii,1).IRGE;
end
C2Vsim_elem_shp = S;
%%
save([mat_data 'C2Vsim_Elements.mat'],'C2Vsim_elem_shp', '-append')
%% write nodes to shapefile
load([mat_data 'C2Vsim_Nodes.mat'], 'C2Vsim_Nodes_shp');
load([mat_data 'C2Vsim_Elements.mat'], 'C2Vsim_elem_shp');
shapewrite(C2Vsim_Nodes_shp, '../gis_data/C2Vsim_Nodes')
shapewrite(C2Vsim_elem_shp, '../gis_data/C2Vsim_Elements')