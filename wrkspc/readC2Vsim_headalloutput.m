function h = readC2Vsim_headalloutput(filename)

fid = fopen(filename, 'r');
cnt_per = 0;
cnt_lay = 1;
h{1,2} = [];

while 1 
   try 
      temp = fgetl(fid); 
      if isempty(temp)
          continue
      end
      if strcmp(temp(1),'*')
          continue
      end
      C = strsplit(temp,' ');
      for i = 1:length(C)
          if isempty(C{1,i})
              continue
          end
          c = textscan(C{1,i},'%f/%f/%f/_%s');
          if isempty(c{1,2})
              %Then its a head value
              h{cnt_per,2}(cnt_nodes, cnt_lay) = c{1,1};
              cnt_nodes = cnt_nodes + 1;
              

          else
              %its time stamp
              cnt_per = cnt_per + 1;
              h{cnt_per,1} = [num2str(c{1,1}) '/' num2str(c{1,2}) '/' num2str(c{1,3})];
              display(h{cnt_per,1})
              cnt_nodes = 1;
              cnt_lay = 1;
          end
      end
      cnt_lay = cnt_lay + 1;
      cnt_nodes = 1;
   catch
       break;
   end
end
fclose(fid);

