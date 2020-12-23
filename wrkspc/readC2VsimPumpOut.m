function PumpOut = readC2VsimPumpOut(filename, rowInfo)
if isempty(rowInfo)
    rowInfo = [5; ...% The row where the 1st section with the variables starts
               11;... % The row where the Types are set
               16; ... % The row where the data start
               ];
end
fid = fopen(filename);
cnt_row = 0;

Time = [];
DATA = [];
while 1
   try
       temp = fgetl(fid); 
       if isempty(temp)
          continue
       end
       cnt_row = cnt_row + 1;
       if cnt_row == rowInfo(1)
           NCOLPUMP = textscan(temp, '%f /%s');
           NCOLPUMP = NCOLPUMP{1,1};
       end
       if cnt_row == rowInfo(2)
           ttt = textscan(temp,'%s');
           ttt{1,1}(1:2,:) = [];
           Type = zeros(length(ttt{1,1}),1);
           for ii = 1:length(ttt{1,1})
               tp = strsplit(ttt{1,1}{ii,1},'_');
               if strcmp(tp{1,1},'Elem')
                   Type(ii,1) = Type(ii,1) + 10;
               elseif strcmp(tp{1,1},'Well')
                   Type(ii,1) = Type(ii,1) + 20;
               end
               if strcmp(tp{1,2},'Ag')
                   Type(ii,1) = Type(ii,1) + 1;
               elseif strcmp(tp{1,2},'Urb')
                   Type(ii,1) = Type(ii,1) + 2;
               end
           end
       end
       if cnt_row == rowInfo(2)+1
           ttt = textscan(temp,'%s');
           ttt{1,1}(1:2,:) = [];
           ID = zeros(length(ttt{1,1}),1);
           for ii = 1:length(ttt{1,1})
               ID(ii,1) = str2double(ttt{1,1}{ii,1});
           end
       end
       if cnt_row == rowInfo(2)+2
           ttt = textscan(temp,'%s');
           ttt{1,1}(1:2,:) = [];
           COL = zeros(length(ttt{1,1}),1);
           for ii = 1:length(ttt{1,1})
               COL(ii,1) = str2double(ttt{1,1}{ii,1});
           end
       end
       
       if cnt_row >= rowInfo(3)
           C = strsplit(temp,' ');
           c = textscan(C{1,1},'%f/%f/%f/_%s');
           display([num2str(c{1,1}) '/' num2str(c{1,2}) '/' num2str(c{1,3})]);
           Time = [Time; c{1,1} c{1,2} c{1,3}];
           V = zeros(NCOLPUMP,1);
           for ii = 1:NCOLPUMP
               V(ii,1) = str2double(C{1,ii+1});
           end
           DATA = [DATA V];
       end
           
      
   catch
       break;
   end
    
end


fclose(fid);

PumpOut.Type = Type;
PumpOut.ID = ID;
PumpOut.Time = Time;
PumpOut.DATA = DATA;



