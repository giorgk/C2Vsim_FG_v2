function PRC = readC2VsimPrecip(filename)

NRAIN = 33570;


fid = fopen(filename, 'r');
cnt_per = 0;
PRC.ARAIN = zeros(1,NRAIN);
while 1
    try
        temp = fgetl(fid);
        if isempty(temp)
            continue;
        end
        if strcmp(temp(1), 'C') || strcmp(temp(1), '*')
            continue
        end
      
        C = strsplit(temp, {' ','\t'});
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
            PRC.time{cnt_per,1} = [num2str(c{1,1}) '/' num2str(c{1,2}) '/' num2str(c{1,3})];
            display(PRC.time{cnt_per,1})
            for j = 1:NRAIN
                PRC.ARAIN(cnt_per,j) = str2double(C{1,1+j});
                
            end
        else
            for j = 1:NRAIN
                PRC.ARAIN(cnt_per,NRAIN) = str2double(C{1,1+j});
            end
        end
   end
    if strcmp(PRC.time{cnt_per,1}, '12/31/2015')
        break
    end
end

fclose(fid);