function[value]=write_skymap_cea_fits(map,fname,varargin)
%write a cea fits file out, but set crval1 to be in the center of the map
%flips Q/U signs for HWP data.
%To be able to write multiple comments needs to be varargin

map_vec=get_keyval_default('rescale_map_vec',[],varargin{:}); %read the vector to rescale maps
comments=varargin;
isok=true(length(comments),1); %even if map_vec is empty, let's see if the keyword is there
for j=1:length(comments),
   foo=0;
   if(ischar(comments{j})),
      foo=strcmp(comments{j},'rescale_map_vec'); 
   end
   if(foo),
      isok(j)=false;
      isok(j+1)=false;
   end
end
comments=comments(isok);

if isstruct(map)
  assert(isfield(map,'map'));
    assert(isfield(map,'mapptr'));
    octave2skymap(map);
    map=map.mapptr;
end
assert(class(map)=='int64');

mm=skymap2octave(map);
if(~isempty(map_vec)),
  mdisp('Here')
  mm(1,:,:)*=map_vec(1);
  mm(2,:,:)*=map_vec(2);
  mm(3,:,:)*=map_vec(3);
end

[rapix,decpix,radelt,decdelt,pv]=get_skymap_cea_params_c(map);
keys={};
vals={};


%now force the reference pixel to be in the center of the map
ra_ref_pix=round(0.5*size(mm,1));
ra_ref=radelt*(ra_ref_pix-rapix);
rapix=ra_ref_pix;

%[keys,vals]=set_keyval_val('SIMPLE','T',keys,vals);
[keys,vals]=set_keyval_val('SIMPLE',true,keys,vals);
[keys,vals]=set_keyval_val('BITPIX',-64,keys,vals);
[keys,vals]=set_keyval_val('NAXIS',2,keys,vals);
[keys,vals]=set_keyval_val('NAXIS1',size(map,1),keys,vals);
[keys,vals]=set_keyval_val('NAXIS2',size(map,2),keys,vals);
%[keys,vals]=set_keyval_val('EXTEND','T',keys,vals);
[keys,vals]=set_keyval_val('EXTEND',true,keys,vals);

%[keys,vals]=set_keyval_val('CTYPE1','''RA---CEA''',keys,vals);
%[keys,vals]=set_keyval_val('CTYPE2','''DEC--CEA''',keys,vals);
%[keys,vals]=set_keyval_val('CUNIT1','''        ''',keys,vals);
%[keys,vals]=set_keyval_val('CUNIT2','''        ''',keys,vals);

[keys,vals]=set_keyval_val('CTYPE1','RA---CEA',keys,vals);
[keys,vals]=set_keyval_val('CTYPE2','DEC--CEA',keys,vals);
[keys,vals]=set_keyval_val('CUNIT1','        ',keys,vals);
[keys,vals]=set_keyval_val('CUNIT2','        ',keys,vals);

[keys,vals]=set_keyval_val('CRVAL1',ra_ref,keys,vals);
[keys,vals]=set_keyval_val('CRVAL2',0,keys,vals);
[keys,vals]=set_keyval_val('CDELT1',radelt,keys,vals);
[keys,vals]=set_keyval_val('CDELT2',decdelt,keys,vals);
[keys,vals]=set_keyval_val('CRPIX1',rapix,keys,vals);
[keys,vals]=set_keyval_val('CRPIX2',decpix,keys,vals);
[keys,vals]=set_keyval_val('PV2_1',pv,keys,vals);
[keys,vals]=set_keyval_val('EQUINOX',2000,keys,vals);
[keys,vals]=set_keyval_val('PC1_1',1.0,keys,vals);
[keys,vals]=set_keyval_val('PC1_2',0.0,keys,vals);
[keys,vals]=set_keyval_val('PC2_1',0.0,keys,vals);
[keys,vals]=set_keyval_val('PC2_2',1.0,keys,vals);
if exist('comments')
  if ischar(comments)
    comments={comments};
  end
  for j=1:length(comments),
    keys(end+1)='comment';
    vals(end+1)=comments{j};
  end
end

%this is a test.  this is only a test.  please enjoy this.  I really cannot say how ling this line is.  I do hope it exceeds 80 characters, though!';


if ~is_map_polarized(map)
  write_fits_cell(fname,mm,keys,vals);
else
  for j=1:size(mm,1)
    ii=max(strfind(fname,'.fits'));
    mytag=['_' get_map_poltag(map,j)];
    if ~isempty(ii)
      fname_use=[fname(1:ii-1) mytag fname(ii:end)];
    else
      fname_use=[fname mytag '.fits'];
    end
    write_fits_cell(fname_use,squeeze(mm(j,:,:)),keys,vals);
  end
end

      

