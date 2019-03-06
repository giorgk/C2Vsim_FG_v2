function V = writeMetaData(V,field1,field2)
if isempty(V)
    V{1,1} = field1;
    V{1,2} = field2;
else
    n = size(V,1);
    V{n+1,1} = field1;
    V{n+1,2} = field2;
end