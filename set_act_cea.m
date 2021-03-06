function[value]=set_act_cea(map,ar)
if ~exist('ar')
  ar='ar1';
end
mdisp(ar)
if strcmp(ar,'ar1')
  set_skymap_cea_predef_c(map,-0.0140108961842663, 0.0140108961842663,8554.0,10034.0,0.353814121044493,14055,1714);
  return
end
if strcmp(ar,'ar1_strip_3year')
    set_skymap_cea_predef_c(map,-0.0140108961842663, 0.0140108961842663,10493,9804,0.353814121044493,18381,1284);
    return
end
if strcmp(ar,'ar2')
  set_skymap_cea_predef_c(map,-0.0140108961842663, 0.0140108961842663,8495.0,9909.0,0.353814121044493,13063,1596);
  return
end
if strcmp(ar,'ar2_new')
  set_skymap_cea_predef_c(map,-0.0140108961842663, 0.0140108961842663,9922.0,9750.0,0.353814121044493,14205,1196);
  return
end
if strcmp(ar,'ar1equ')
  set_skymap_cea_predef_c(map,-0.0008334, 0.0008334,7584.0,323.0,1.0,28830,687);
  return
end

if strcmp(ar,'ar1equ_2010')
  %set_skymap_cea_predef_c(map,-0.0083336141, 0.0083336141,8558.0,265.0,0.99993261,28672,644);
  set_skymap_cea_predef_c(map,-0.008334, 0.008334,8558.0,265.0,1.0,28672,644);
  return
end

if strcmp(ar,'allequ')  %Hasselfield's preferred header for all equatorial maps
  %set_skymap_cea_predef_c(map,-0.0083336141, 0.0083336141,8558.0,265.0,0.99993261,28672,644);
  set_skymap_cea_predef_c(map,-0.00825, 0.00825,7972.0,292.0,1.0,38275,611);
  return
end

if strcmp(ar,'allequ_ar3')  %Hasselfield's preferred header for all equatorial maps
  %set_skymap_cea_predef_c(map,-0.0083336141, 0.0083336141,8558.0,265.0,0.99993261,28672,644);
  set_skymap_cea_predef_c(map,-0.00825, 0.00825,8150.0,292.0,1.0,38500,611);
  return
end

if strcmp(ar,'allequ_round')  %Hasselfield's preferred header for all equatorial maps
  %set_skymap_cea_predef_c(map,-0.0083336141, 0.0083336141,8558.0,265.0,0.99993261,28672,644);
  set_skymap_cea_predef_c(map,-0.00825, 0.00825,7968.0,288.0,1.0,38304,640);
  return
end

if strcmp(ar,'allequ_round_ar3')  %Hasselfield's preferred header for all equatorial maps
  %set_skymap_cea_predef_c(map,-0.0083336141, 0.0083336141,8558.0,265.0,0.99993261,28672,644);
  set_skymap_cea_predef_c(map,-0.00825, 0.00825,8150.0,288.0,1.0,38500,640);
  return
end


if strcmp(ar,'allsouth')  %Hasselfield's preferred header for all equatorial maps
  set_skymap_cea_predef_c(map,-0.0138696776687, 0.0138696776687,10582,9833,0.35381412,18787,1177);
  return
end

if strcmp(ar,'allsouth_round')  %Hasselfield's preferred header for all equatorial maps
  set_skymap_cea_predef_c(map,-0.0138696776687, 0.0138696776687,10592,9824,0.35381412,18784+32,1184);
  return
end

if strcmp(ar,'newsouth')  %patching up the header since some data maps out of bounds
  set_skymap_cea_predef_c(map,-0.0138696776687, 0.0138696776687,10582,9893+40,0.35381412,18787,1177+60+80+120);
  return
end

if strcmp(ar,'newsouth_round')  %let the above header be multiple of 32 pixels
  set_skymap_cea_predef_c(map,-0.0138696776687, 0.0138696776687,10582,9893+40,0.35381412,18816,1440);
  return
end

if strcmp(ar,'ar1_dunner_as_das')
  set_skymap_cea_predef_c(map,-0.0140108961842663, 0.0140108961842663,10478.0,9734.0,0.353814121044493,18600,1170);
  return
end


if strcmp(ar,'abs_field_A')
  set_skymap_cea_predef_c(map,-0.10962491,0.10962491,920.0,816.0,0.57785513,1256,561);
  return
end

if strcmp(ar,'abs_field_A_fine')
  fac=1.5/5;
  set_skymap_cea_predef_c(map,-0.10962491*fac,0.10962491*fac,round(920.0/fac),round(816.0/fac),0.57785513,round(1256/fac),round(561/fac));
  return
end



if strcmp(ar,'newsouth_fine')  %patching up the header since some data maps out of bounds
  fac=1.5;
  set_skymap_cea_predef_c(map,-0.0138696776687/fac, 0.0138696776687/fac,10582*fac,(9893+40)*fac,0.35381412,18787*fac,(1177+60+80+120)*fac);
  return
end


error(['error in set_act_cea - unrecognized array ' ar]);
