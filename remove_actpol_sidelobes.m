function[value]=remove_actpol_sidelobes(mytod,myopts)

do_actpol_pointing=get_struct_mem(myopts,'do_actpol_pointing',false);
sl=get_struct_mem(myopts,'sl','');
sl_det=get_struct_mem(myopts,'sl_det','');
T_sl=get_struct_mem(myopts,'T_sl',[]);

[dy,dx,ang] = get_detector_offsets_actpol(mytod);
[rows,cols] = get_tod_rowcol(mytod);
sl_corr=dat*0.;
comps = size(sl,1);

for k=1:comps, %loop over sidelobe components
   %destroying poiting
   assign_tod_value(mytod,0);
   free_tod_pointing_saved(mytod,free_2gamma);
   free_saved_pixellization(mytod);
   destroy_actpol_pointing(mytod);
   destroy_actpol_pointing_fit(mytod);
   
   %creating shifted pointing, with angle correction from JLS
   initialize_actpol_pointing(mytod,dy+sl(k,2),dx+sl(k,1),ang,148.0,1);
   get_actpol_pointing_fit(mytod,'downsamp',20);
   set_tod_twogamma_fit(mytod,'downsamp',20,'rescale_az_2gamma',true);
   precalc_actpol_pointing_exact(mytod);
   convert_saved_pointing_to_pixellization(mytod,T_sl);
   
   %estimating buddy signal
   map2tod(T_sl,mytod);
   det_list=sl_det{k};
   buddy=(sl(k,3)+sl(k,4)*cos(2.*ang)+sl(k,5)*sin(2.*ang));
   buddy=repmat(buddy,[1 size(dat,1)])';
   T_sig=get_tod_data(mytod);
 
   %Applying the correction to detectors
   for h=1:size(rows,1),
      uid=32*rows(h)+cols(h);
      ii=find(det_list==uid);
      if ~(isempty(ii))
         sl_corr(:,h)+=T_sig(:,h).*buddy(:,h);
      end
   end
   free_tod_pointing_saved(mytod,free_2gamma);
   free_saved_pixellization(mytod);
   clear buddy; clear T_sig;
end

assign_tod_value(mytod,0);
push_tod_data(sl_corr,mytod);
clear sl_corr; clear rows; clear cols; clear comps;

mdisp('Restoring original pointing');
destroy_actpol_pointing(mytod);
destroy_actpol_pointing_fit(mytod);
initialize_actpol_pointing(mytod,dy,dx,ang,148.0,1);
clear dy; clear dx; clear ang;
get_actpol_pointing_fit(mytod,'downsamp',20);
set_tod_twogamma_fit(mytod,'downsamp',20,'rescale_az_2gamma',true);
if (do_actpol_pointing)
   precalc_actpol_pointing_exact(mytod,1);
   if isfield(mapset.skymap,'mapptr')
      convert_saved_pointing_to_pixellization(mytod,mapset.skymap.mapptr)
   end
   free_tod_pointing_saved(mytod,free_2gamma);
   if (~does_tod_have_twogamma_fit(mytod))
      precalc_actpol_pointing_exact(mytod,2);
      set_tod_twogamma_fit(mytod,'npoly_2gamma',3,'rescale_az_2gamma',true);
      free_tod_pointing_saved(mytod,free_2gamma);
   end
end
