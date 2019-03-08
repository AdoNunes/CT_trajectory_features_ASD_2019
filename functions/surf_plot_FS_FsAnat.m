
% Creates a Freesurfer annotation 
function str_cmd = surf_plot_FS_FsAnat(vals_2plot, thr, fname, plot_borders, randcolor )

% vals_2plot: values for each FS anatomical parcellation area ( n areas x 1)
% thr: threshold for the values, if negative only smaller values will be
%     plot, can be  0 or empty
% fname : filename for the output annotation
% plot_borders: 1 or 0, if 1 the outlay

%% fressurfer plot
curr_path = pwd;
fs_path = [curr_path, '/Freesurfer/'];
fs_path_lb = [fs_path, 'fsaverage/label/'];
fs_path_sf = [fs_path, 'fsaverage/surf/'];
Results_folder = [pwd,'/Results/'];

load('Freesurfer/labels_names.mat')

if ~isempty(thr) && thr>0
    thr_bin = vals_2plot> thr;
elseif ~isempty(thr) && thr<0
    thr_bin = vals_2plot< thr;
else
    thr_bin = true(size(vals_2plot));
end

labels = names_ANAT;

Z_vals_thrd = vals_2plot(thr_bin);
[s, i] = sort(Z_vals_thrd, 'descend');

labs_plot = labels(thr_bin)';
labs_plot = labs_plot(i);

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
    
    eval(['lut = make_LUT(val_sort(val_',hems{h},'), ''',fname,hems{h},''', ',num2str(randcolor), ');'])
    eval(['lut = make_LUT(val_sort(val_',hems{h},'), ''',fname,hems{h},''', ',num2str(randcolor), ');'])
    eval(['fs_area = fs_area_',hems{h},';'])
    
    system([subj_var,'; mris_label2annot --ctab ', lut, ' --s fsaverage --h ',hems{h},'h  --a ', fname,' ', [fs_area{:}]])
    
    eval(['val_srt = val_sort(val_',hems{h},');'])
    [~, label, colortable] = read_annotation([fs_path_lb,hems{h},'h.' fname,'.annot']);
    colortable.table(:,6) = [0; val_srt ];
    for n = 1:colortable.numEntries
        label(label==colortable.table(n,5)) = colortable.table(n,6);
    end
    
    if ~exist([fs_path_sf,hems{h},'h.white.mgh'], 'file')
        system(['mris_convert ',fs_path_sf,hems{h},'h.white ' ,fs_path_sf, hems{h},'h.white.mgh'])
    end
    
    fs = MRIread([fs_path_sf,hems{h},'h.white.mgh']);
    fs.vol(1,:,1,1) = label;
    fname_out= [fs_path_sf,hems{h},'h.' fname,'.mgz'];
    MRIwrite(fs,fname_out)
    
    if randcolor
        str_cmd{h} = [subj_var,'; tksurfer fsaverage ',hems{h},'h inflated -annot ', [fs_path_lb,hems{h},'h.' fname,'.annot'],  ' &' ];
    elseif plot_borders
        str_cmd{h} = [subj_var,'; tksurfer fsaverage ',hems{h},'h inflated  -o ', fname_out,' -annot ', [fs_path,hems{h},'h.FsAnat.annot'],' -fminmax ',num2str(min(abs(val_sort))*.98),' ',num2str(max(abs(val_sort))*.98),  '&' ];
        
    else
        str_cmd{h} = [subj_var,'; tksurfer fsaverage ',hems{h},'h inflated  -o ', fname_out,' -fminmax ',num2str(min(abs(val_sort))*.98),' ',num2str(max(abs(val_sort))*.98), '&' ];
    end
    system(  str_cmd{h} )
    
end
