function str_cmd = surf_plot_FS_H(vals_2plot, fname, colors )
%% fressurfer plot
curr_path = pwd;
fs_path = [curr_path, '/Freesurfer/'];
fs_path_lb = [fs_path, 'fsaverage/label/'];
fs_path_sf = [fs_path, 'fsaverage/surf/'];
Results_folder = [pwd,'/Results/'];

load('Freesurfer/labels_names.mat')


labels = names_ANAT;

Z_vals_thrd = ones(size(labels));
[s, i] = sort(Z_vals_thrd, 'descend');

labs_plot = labels;

inx_l = 0;
inx_r = 0;
fs_area_r = {};
fs_area_l = {};
val_r = [];
val_l = [];
for l = 1:numel(labs_plot)
    if labs_plot{l}(1) =='r'
        inx_r =inx_r+1;
        fs_area_r{inx_r} = ['--l ',fs_path_lb,'/rh.', labs_plot{l}(3:end), '.label '];
        val_r(inx_r) = l;
    else
        inx_l =inx_l+1;
        fs_area_l{inx_l} = ['--l ',fs_path_lb,'/lh.', labs_plot{l}(3:end), '.label '];
        val_l(inx_l) = l;
    end
end

val_sort = Z_vals_thrd(i);


subj_var = ['SUBJECTS_DIR=',fs_path];

hems = {'r', 'l'};
for h = 1:2
    if exist( [fs_path_lb,hems{h},'h.',fname,'.annot'], 'file')
        delete([fs_path_lb,hems{h},'h.',fname,'.annot'])
    end
    
    val = repmat(colors(h,:),numel(val_l),1);
    n = 1: numel(val_l);
    txt_n = [n',n',round(val*255) ];
    
    lutname = sprintf('/Applications/freesurfer/LUT/LUT_H_%s.txt', fname );
    fileID = fopen(lutname,'w');
    formatSpec = '    %d seg%d    %d    %d    %d    0\n';
    fprintf(fileID,formatSpec,txt_n');    
    
    eval(['fs_area = fs_area_',hems{h},';'])
    
    system([subj_var,'; mris_label2annot --ctab ', lutname, ' --s fsaverage --h ',hems{h},'h  --a ', fname,' ', [fs_area{:}]])
    
    eval(['val_srt = val_sort(val_',hems{h},');'])
    [~, label, colortable] = read_annotation([fs_path_lb,hems{h},'h.' fname,'.annot']);
    
    str_cmd{h} = [subj_var,'; tksurfer fsaverage ',hems{h},'h inflated -annot ', [fs_path_lb,hems{h},'h.' fname,'.annot'],  ' &' ];
    
    system(  str_cmd{h} )
    
end
