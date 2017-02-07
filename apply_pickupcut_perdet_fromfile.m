function[ll]=apply_pickupcut_perdet_fromfile(tod,cuts_file)
%azimuth cut defined for each detector. We assume that we know which file we need to send inside the function (i.e. check on alt done outside).

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
%bigstr has dimension (all_det)x(some)

% TOD info
[alt,az]=get_tod_altaz(tod);
az=az*180/pi;
az(az>180)=az(az>180)-360;
[rows,cols]=get_tod_rowcol(tod);

mdisp('Applying az cut from Sigurd files.');
%Reading and applying the cuts for each detector
for j=1:length(rows),
   uid=cols(j)+32*rows(j);
   isbad=false(size(az));
   if (str2num(bigstr(uid+1,1:4))==uid)
      cuts=bigstr(uid+1,floor(uid/1000)+6:max(find(bigstr(uid+1,:)!=' ')));
      foo=find(cuts==' ');

      if (length(cuts)>1)
         foo=[0, foo, length(cuts)+1];
         for k=1:length(foo)-1,
            cut=cuts(foo(k)+1:foo(k+1)-1);
            crap=find(cut==':');
            az_cut_min=str2num(cut(1:crap-1));
            az_cut_max=str2num(cut(crap+1:end));
            ii=(az>=az_cut_min)&(az<=az_cut_max);
            isbad(ii)=true;
         end
         [istart,istop]=cutvec2inds(isbad);
         for i=1:length(istart),
            cuts_extend_c(tod,istart(i)-1,istop(i),rows(j),cols(j));
         end
      end
   end
end

