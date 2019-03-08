% Creates a color table for creating annotation files
% it creates a gradiant of colors based on max and min values

function  fname = make_LUT(values, filename, randcl)
% values: list of values
% filename: to be saved
% randcl: 1 or 0, if 0 random colors are assigned
 
if randcl
    grad = colormap(lines(numel(unique(values))));
else
    grad = colormap(hot(100));
end
close

grad = grad* 255;
% unit variance
list =  abs(values);

list = (list - min(list)) / ( max(list) - min(list) );
if isnan(list)
    list(isnan(list)) = 1;
end
list = round(list*size(grad,1));
list(list==0) = 1;
val = grad(list, :, :);

n = 1: numel(list);

txt_n = [n',n',round(val) ];

fname = sprintf('/Applications/freesurfer/LUT/LUT_%d_values_%s.txt',numel(list), filename );
fileID = fopen(fname,'w');
formatSpec = '    %d seg%d    %d    %d    %d    0\n';

 fprintf(fileID,formatSpec,txt_n');

