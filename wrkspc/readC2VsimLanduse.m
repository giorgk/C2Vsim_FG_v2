function LU = readC2VsimLanduse(filename)

LU.time = [];
LU.AG = [];
LU.UR = [];
LU.NV = [];
LU.RV = [];

fid = fopen(filename, 'r');
cnt_per = 0;
while 1
    try
    temp = fgetl(fid);
        if isempty(temp)
            continue
        end
        if strcmp(temp(1),'C')
            continue
        end
        
        C = strsplit(temp,' ');
        if strcmp(C{1,3},'/') || strcmp(C{1,2},'/')
            continue
        end
        if ~isempty(C{1,1})
            c = textscan(C{1,1},'%f/%f/%f/_%s');
        else
            c{1,2} = [];
        end
        if ~isempty(c{1,2})
            cnt_per = cnt_per + 1;
            LU.time{cnt_per,1} = [num2str(c{1,1}) '/' num2str(c{1,2}) '/' num2str(c{1,3})];
            display(LU.time{cnt_per,1})
            LU.AG(str2double(C{1,2}),cnt_per) = str2double(C{1,3});
            LU.UR(str2double(C{1,2}),cnt_per) = str2double(C{1,4});
            LU.NV(str2double(C{1,2}),cnt_per) = str2double(C{1,5});
            LU.RV(str2double(C{1,2}),cnt_per) = str2double(C{1,6});
        else
            LU.AG(str2double(C{1,2}),cnt_per) = str2double(C{1,3});
            LU.UR(str2double(C{1,2}),cnt_per) = str2double(C{1,4});
            LU.NV(str2double(C{1,2}),cnt_per) = str2double(C{1,5});
            LU.RV(str2double(C{1,2}),cnt_per) = str2double(C{1,6});
        end
       
   catch
       break;
   end
end
fclose(fid);

