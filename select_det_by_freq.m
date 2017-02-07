function[det_offsets_good]=select_det_by_freq(det_offsets_good,det_specs_file,freq)
%This function selects detectors of given freq from spreadsheet. Hopefully will be merged in array_info.

   det_freq = read_text_file_comments(det_specs_file);
   for j=1:length(det_freq),
      tmp = strsplit(det_freq{j});
      tmp_uid(j) = str2num(tmp{2});
      tmp_freq(j) = str2num(tmp{5});
   end

   tmp_uid=tmp_uid(tmp_freq==freq);

   isok = false(length(det_offsets_good));

   for j=1:length(tmp_uid),
      myind = find(tmp_uid(j)==det_offsets_good(:,1));
     isok(myind) = true;
   end

   det_offsets_good = det_offsets_good(isok,:);

   return;
