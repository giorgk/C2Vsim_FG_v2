function LU = LU_Area_per_element(filename, NE, Cols, Nheaderlines)
fid = fopen(filename,'r');
for ii = 1:Nheaderlines
    temp = fgetl(fid);
end
LU = [];
Ncol = length(Cols);
for ii = 1:Ncol
    LU.(Cols{ii}) = zeros(2000,NE);
end

cnt_per = 0;
while  ~feof(fid)
    temp = fgetl(fid);
    C = strsplit(temp, {' ','\t'});
    if ~isempty(C{1,1})
        c = textscan(C{1,1},'%f/%f/%f/_%s');
        cnt_per = cnt_per + 1;
        LU.time{cnt_per,1} = [num2str(c{1,1}) '/' num2str(c{1,2}) '/' num2str(c{1,3})];
        display(LU.time{cnt_per,1});
        IE = str2double(C{1,2});
        for jj = 1:Ncol
            LU.(Cols{jj})(cnt_per,IE) = str2double(C{1,2+jj});
        end
    else
        IE = str2double(C{1,2});
        for jj = 1:Ncol
            LU.(Cols{jj})(cnt_per,IE) = str2double(C{1,2+jj});
        end
    end
end
fclose(fid);

cnt_per = cnt_per + 1;
for jj = 1:Ncol
    LU.(Cols{jj})(cnt_per:end,:) = [];
end