function[mytod_names]=distribute_tods_from_file(tod_names,file,myid,nproc)
%This function distributes TODs based on central RA, to avoid large maps in partitioned mapping

%Read tods limits
lines_in=read_lines(file);
names_in=cell(size(lines_in));
lims_in=cell(size(lines_in));

nline=length(names_in);
for j=1:nline,
  [a,b]=strtok(lines_in{j});
  names_in(j)={a};
  lims_in(j)={b};
end



%Getting tags for match
names_in=get_tod_tags_from_names(names_in);
tags=get_tod_tags_from_names(tod_names);

tt(length(tags),length(tags{1}))=' ';
for j=1:length(tags),
  tt(j,:)=tags{j};
end

nn(length(names_in),length(names_in{1}))=' ';
for j=1:length(names_in),
  nn(j,:)=names_in{j};
end



%Checking were our tod_names are located in file, cutting missing TODs.
ntag=length(tags);
tod_ra_min=zeros(length(ntag));
tod_ra_max=zeros(length(ntag));
tod_ra_mean=zeros(length(ntag));
tod_dec_min=zeros(length(ntag));
tod_dec_max=zeros(length(ntag));

for j=1:ntag,
  ind=strmatch(tt(j,:),nn);
  if numel(ind)==1,
     [a,b]=strtok(lims_in{ind});
     [b,c]=strtok(b);
     tod_ra_min(j)=str2num(a);
     tod_ra_max(j)=str2num(b);
     tod_ra_mean(j)=(tod_ra_min(j)+tod_ra_max(j))/2.;
     [a,b]=strtok(c);
     tod_dec_min(j)=str2num(a);
     tod_dec_max(j)=str2num(b);
  end
end


if sum(tod_ra_mean==0)>0,
  mdisp(['warning - in load_balance_tods_from_file, have TODs with missing time info. Cutting!!']);
  isok=true(size(tod_ra_min));
  isok(tod_ra_mean==0)=false;
  tod_ra_min=tod_ra_min(isok);
  tod_ra_max=tod_ra_max(isok);
  tod_ra_mean=tod_ra_mean(isok);
  tod_dec_min=tod_dec_min(isok);
  tod_dec_max=tod_dec_max(isok);
  tod_names=tod_names(isok);
end


%Computing myid for distribution
[s,i]=sort(tod_ra_mean);
tod_ra_mean=tod_ra_mean(i);
tod_ra_min=tod_ra_min(i);
tod_ra_max=tod_ra_max(i);
tod_dec_min=tod_dec_min(i);
tod_dec_max=tod_dec_max(i);
tod_names=tod_names(i);


ntods_per_node=floor(length(tod_names)/nproc);
if (ntods_per_node<1),
   warning(['Less than a TOD per node!']);
   return;
end
if (myid!=nproc),
   mytod_names=tod_names(((myid-1)*ntods_per_node)+1:(myid)*ntods_per_node);
   ra_min=min(tod_ra_min(((myid-1)*ntods_per_node)+1:(myid)*ntods_per_node));
   ra_max=min(tod_ra_max(((myid-1)*ntods_per_node)+1:(myid)*ntods_per_node));
   dec_min=min(tod_dec_min(((myid-1)*ntods_per_node)+1:(myid)*ntods_per_node));
   dec_max=min(tod_dec_max(((myid-1)*ntods_per_node)+1:(myid)*ntods_per_node));
else
   mytod_names=tod_names(((myid-1)*ntods_per_node)+1:end);
   ra_min=min(tod_ra_min(((myid-1)*ntods_per_node)+1:end));
   ra_max=min(tod_ra_max(((myid-1)*ntods_per_node)+1:end));
   dec_min=min(tod_dec_min(((myid-1)*ntods_per_node)+1:end));
   dec_max=min(tod_dec_max(((myid-1)*ntods_per_node)+1:end));
end
disp(['For myid:' num2str(myid) ', the size of the C map in RA, dec is: [' num2str(ra_min) ':' num2str(ra_max) ',' num2str(dec_min) ':' num2str(dec_max) ']'])

mpi_barrier;
return 
