function[map]=make_weightmap_octave(tods,map,do_window,varargin)

do_new_pointing=get_keyval_default('do_new_pointing',false,varargin{:});
have_pointing=get_keyval_default('have_pointing',false,varargin{:}); %pointing is already saved in TODs, and you don't want to throw it away
do_reduce=get_keyval_default('mpi_reduce',false,varargin{:});
free_2gamma=get_keyval_default('free_2gamma',true,varargin{:});
do_noise=get_keyval_default('do_noise',false,varargin{:});  %added ability to do noise instead of hitcounts.
do_ivar2std=get_keyval_default('ivar2std',true,varargin{:});  %for noise: we want to mpi_reduce before taking the sqrt
noise_file=get_keyval_default('noise_file','',varargin{:});

if isstruct(map)
  if isfield(map,'partition')
    clear_map(map.mapptr);
    make_weightmap_octave(tods,map.mapptr,do_window,'mpi_reduce',false,'ivar2std',false,varargin{:}); %specified twice but seems to do what set right now
    if do_reduce
      map=mpi_reduce_map(map); 
      octave2skymap(map);
    end
    if (do_noise&do_ivar2std)
       map=f_ivar2std(map);
    end
    return;
  end
end

clear_map(map);

if isempty('do_window')
  clear do_window;
end


if ~exist('do_window')
  do_window=true;
end

for j=1:length(tods),
  mytod=tods(j);
  allocate_tod_storage(mytod);
  if (do_noise)
    tod_name=get_tod_name(mytod);
    if iscell(tod_name), tod_name=tod_name{1}; end
    tod_name=tod_name(max(find(tod_name=='/'))+1:end);
    [rows,cols]=get_tod_rowcol(mytod);
    uid=32*rows+cols;
    mynoise=get_detector_average_noise_banded_projvecs(mytod);
    nn=get_tod_ndata(mytod);
    wt_per_samp=1./mynoise.^2;
    dat=repmat(wt_per_samp',[nn 1]);
    push_tod_data(dat,mytod);
    clear dat;
    if (~isempty(noise_file)),
       array_noise=1/sqrt(sum(wt_per_samp))*sqrt(get_tod_dt(mytod));
       mynoise=mynoise*sqrt(get_tod_dt(mytod));
       mat=zeros(2*length(rows),1);mat([1:2:2*length(rows)])=uid;mat([2:2:2*length(rows)])=mynoise;
       to_write{j}=sprintf(['%s %5.3f :: ' repmat(' %4d:%5.3f ',1,length(rows))],tod_name,array_noise,mat(:));
    end
  else
    assign_tod_value(mytod,1.0);
  end
  if do_window,
    window_data(mytod);
    window_data(mytod);
  end

  
  if (do_new_pointing)
    mdisp('doing pointing')
    precalc_actpol_pointing_exact(mytod);
    convert_saved_pointing_to_pixellization(mytod,map);
    free_tod_pointing_saved(mytod,free_2gamma);
  end
  if (have_pointing)
    mdisp('converting_pointing')
    convert_saved_pointing_to_pixellization(mytod,map);
  end  

  mdisp('assigned value');
  if (is_map_polarized(map))
    tod2polmap(mytod,map);
  else
    tod2map(mytod,map);
  end
  mdisp('projected');
  
  if (do_new_pointing)
    %free_tod_pointing_saved(mytod);
    free_saved_pixellization(mytod);
    
  end


  free_tod_storage(mytod);    
  %disp('freed');
end

if (do_noise&(~isempty(noise_file))),
   myid=mpi_comm_rank+1;
   for j=1:mpi_comm_size,
      if (myid==j)
        if (myid==1)
          fid=fopen([noise_file '.noise'],'w');
        else
          fid=fopen([noise_file '.noise'],'a');
        end
        fprintf(fid,'%s\n',to_write{:});
        fclose(fid);
      end
      mpi_barrier;
   end
end

if (do_reduce)
  mpi_barrier;
  mdisp('reducing')
  weight=skymap2octave(map);
  %JLS: scinet mpi_allreduce started to act up for large matrices, so split up polarization ones
  if length(size(weight))==3,
    for j=1:size(weight,1),
      weight(j,:,:)=mpi_allreduce(weight(j,:,:));
    end
    myweight=sum(sum(sum(abs(weight))));
    disp(['node ' num2str(mpi_comm_rank+1) ' has weight ' num2str(myweight)]);
  else
    myweight=sum(sum(sum(abs(weight))));
    disp(['node ' num2str(mpi_comm_rank+1) ' has weight ' num2str(myweight)]);
    
    weight=mpi_allreduce(weight);
  end
  mdisp(['new total weight is ' num2str(sum(sum(sum(abs(weight)))))]);
  octave2skymap(weight,map);
  clear weight;
end
  
if (do_noise&do_ivar2std)
   f_ivar2std(map);
end


function[map]=f_ivar2std(map)
  if isfield(map,'partition');
    weight=map.map;
  else
    weight=skymap2octave(map);
  end

  if numel(size(weight))==3,
    asdf=weight(1,:,:);
    asdf(asdf>0)=1./sqrt(asdf(asdf>0));
    weight(1,:,:)=asdf;
  else
    weight(weight>0)=1./sqrt(weight(weight>0));
  end

  if isfield(map,'partition');
    map.map=weight;
    octave2skymap(map);
  else
    octave2skymap(weight,map);
  end
return;
