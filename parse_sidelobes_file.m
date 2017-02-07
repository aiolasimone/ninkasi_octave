function[sl_red,sl_det]=parse_sidelobes_file(filename)

sl = myload(filename);
[s,i]=sortrows(sl,2);
ii=find(diff([s(:,2); 0])!=0);
sl_red = s(ii,2:end);

sl_det={};
for j=1:size(sl_red,1),
   ii=find(sl(:,2)==sl_red(j,1));
   sl_det(j)=sl(ii,1);
end

