
%% modify annotation for FreeSurfer anatomical parcellation
% removes corpuscallosum 

H = { 'l', 'r'};

for ih = 1:2
    [vertices, label, colortable] = read_annotation('Freesurfer/fsaverage/label/',H{ih},'h.aparc.annot');

    out_label = 'corpuscallosum';
    out_inx =strcmp(colortable.struct_names, 'corpuscallosum');

    label(find( label == colortable.table(out_inx,end))) = colortable.table(1,end);

    colortable.numEntries = nnz(~out_inx);
    colortable.struct_names(out_inx) = [];
    colortable.table(out_inx,:) = [];

    write_annotation('Freesurfer/',H{ih},'h.FsAnat.annot',vertices, label, colortable);
end
