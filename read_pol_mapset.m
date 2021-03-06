function[big_map]=read_pol_mapset(froot,pols,ftail)
if ~exist('pols')
  pols={'I','Q','U'};
end
if isempty(pols)
  pols={'I','Q','U'};
end
if ~exist('ftail')
  ftail='';
end


npol=length(pols);
for j=1:npol,
  %fname=[froot pols{j} '.fits'];
  fname=[froot pols{j} ftail '.fits'];
  try
    %map=fits_image_read([froot pols{j} '.fits']);
    map=fits_image_read(fname);
  catch
    try
      fname=[froot pols{j} ftail];
      map=fits_image_read(fname);
    catch
      disp(['unable to find ' fname ' in read_pol_mapset.']);
      assert(1==0);
    end
  end
  if j==1,
    big_map=zeros([npol size(map,1) size(map,2)]);
  end
  big_map(j,:,:)=map;
end


