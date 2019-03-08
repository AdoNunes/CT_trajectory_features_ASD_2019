
function surf_plot_matlab(vals_2plot, atlas_2plot, hem, thr)
%
%Plot cortcial surface with values in matlab
%
%       vals_2plot = n labels x 1, the labels to plot
%       atlas_2plot = annotatio tha will be used: 'HCP-MMAll' or 'FsAnat'
%      hem = hemishpere 'r' or 'l'
%

if nargin ==4
    if thr>0
        thr_bin = vals_2plot> thr;
    elseif  thr<0
        
        thr_bin = vals_2plot< thr;
    end
    vals_2plot(~thr_bin) = 0;
    figure;
    myColorMap = colormap(jet); close
    myColorMap(end,:,:) = [ .9 .9 .9 ];
else
    myColorMap = colormap(jet); %
end

surface_type = 'inflated';

% Read the surface file
g = gifti( ['Freesurfer/', hem, 'h.',surface_type,'.gii']);

[~, label, colortable] = read_annotation(['Freesurfer/', hem, 'h.', atlas_2plot,'.annot']);

parcels_vals = zeros(size(label));

for p = 2:size(colortable.table,1) % first element colortable is medial wall
    parcels_vals(label == colortable.table(p,end)) = vals_2plot(p-1);
end  


gg = g;
gg.cdata = parcels_vals;

figure,
plot(g,gg);
colormap(myColorMap);
material dull
if hem == 'l'
    view(-90, 10); 
    camlight
camlight(-80,-10);

else
    view(90, 10);
    camlight
 camlight(80,-10);
end



