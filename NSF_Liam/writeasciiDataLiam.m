function writeasciiDataLiam(filename, frmt, plss, data)

fid = fopen(filename, 'w');
for ii = 1:size(data,1)
    fprintf(fid, frmt, plss(ii,1).CO_MTRS, [plss(ii,1).X1 plss(ii,1).Y1 plss(ii,1).X plss(ii,1).Y data(ii,:)]);
end
fclose(fid);