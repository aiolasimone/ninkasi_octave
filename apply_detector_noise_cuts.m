function[ll]=apply_detector_noise_cuts(tods,cuts_file,do_purge)
%cut detectors explicitly listed in the cuts_file.  format is from Sigurd:
%  1378863997.1378864029.ar1   3: 271 367 687
%sample file: /project2/r/rbond/sigurdkn/actpol/mapdata/noise_cut/noise_whiteness_cut_c7v5b_150712.txt

if ~exist('do_purge')
  do_purge=true;
end


fid=fopen(cuts_file,'r');
ll=fread(fid,inf,'char');
fclose(fid);
nl=find(ll==sprintf('\n'));
ntod=length(nl);
nl=[0; nl];
mylen=diff(nl)-1;
maxlen=max(mylen);
bigstr(1:ntod,1:maxlen)=' ';
for j=1:ntod,  
  bigstr(j,1:mylen(j))=ll(nl(j)+1:nl(j+1)-1);
end

for j=1:length(tods),
  flub=get_tod_name(tods(j));ii=max(find(flub=='/'));flub=flub(ii+1:end);
  nn=length(flub);
  for k=1:ntod
    if strcmp(flub,bigstr(k,1:nn))
      %disp(strtrim(bigstr(k,:)));
      ii=find(bigstr(k,:)==':');
      bad_dets=str2num(bigstr(k,ii+1:end));
      cc_tmp=rem(bad_dets,32);
      rr_tmp=(bad_dets-cc_tmp)/32;
      for kk=1:length(cc_tmp),
        cut_detector_c(tods(j),rr_tmp(kk),cc_tmp(kk));
      end
      if (do_purge)
        purge_cut_detectors(tods(j));
      end
      break;
    end
  end
end
