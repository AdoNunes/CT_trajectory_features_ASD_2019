%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script removes:
%   - Center variability for a given model and group (the same model used in 
%       the subsequent analyses)
%   - Areas from an atlas where a given model is not a good fit
%
% Adonay Nunes, SFU, Vancouver, Feb 2019
% adonay.s.nunes@gmail.com
% from github: AdoNunes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear 
addpath('functions')
load('S01_data.mat')

do_plots_FS= 0;
do_plots_MAT = 0;

str_md = {'linear', 'quadratic', 'cubic'};
str_at = {'MSALL', 'FsAnat'};

%% Remove center variance 
T.siteID = removecats(T.siteID);
centers = unique(T.siteID);
centers(centers == 'NYU') =[];

dummy_site = zeros(size(CT.MSALL.raw,1), numel(centers) );
for c = 1:size(dummy_site,2)
    dummy_site(:,c) = +(T.siteID == centers(c));
end

a0 = ones(size(T.age));
a1=  T.age;
a2 = T.age.^2;
a3 = T.age.^3;
g1 = T.group =='asd';

mdl_site.linear    = [                    a1 a1.*g1 a0 g1 dummy_site];
mdl_site.quadratic = [          a2 a2.*g1 a1 a1.*g1 a0 g1 dummy_site]; 
mdl_site.cubic     = [a3 a3.*g1 a2 a2.*g1 a1 a1.*g1 a0 g1 dummy_site];

for at = 1:numel(str_at) % atlas  
    for md = 1:numel(str_md) % model  
        ct = CT.(str_at{at}).raw;
        
        for s = 1:size(ct,2)  % atlas areas          
            coeff  = regress(ct(:,s), mdl_site.linear);
            site_var = dummy_site* coeff(end-size(dummy_site,2)+1:end);
            
            CT.(str_at{at}).(str_md{md})(:,s) = ct(:,s)- site_var;
        end        
    end
end


str_grp = {'asd', 'ctr'};
 
group{1} = find(T.group== str_grp{1});
group{2} = find(T.group== str_grp{2});

ages{1} = T.age(group{1});
ages{2} = T.age(group{2});

mds_form = {'y ~  x1 ', 'y ~  x1^2 ', 'y ~  x1^3 '};

mdls_fit = struct;
for at = 1:numel(str_at) % atlas
    ct = CT.(str_at{at}).(str_md{md});
    
    for s = 1:size(ct,2)  % atlas areas
        for g = 1:numel(str_grp)% groups
            for md = 1:numel(str_md) % model
                
                GenMdl  = fitglm( ages{g} , ct(group{g},s), mds_form{md});
                
                mdls_fit.(str_at{at}).(str_grp{g}).pVal(s,md) = GenMdl.devianceTest.pValue(2);
                mdls_fit.(str_at{at}).(str_grp{g}).AIC(s,md) = GenMdl.ModelCriterion.AIC;
            end
        end
    end
end


for at = 1:numel(str_at) % atlas
    for md = 1:numel(str_md) % model
        for g = 1:numel(str_grp)% groups
            
            pvls = mdls_fit.(str_at{at}).(str_grp{g}).pVal(:,md);
            
            mdls_fit.(str_at{at}).(str_grp{g}).pVal_FDR.(str_md{md}) = fdr_bh(pvls ,0.05,'dep', 'yes');
        end
        mdls_fit.(str_at{at}).pVal_FDR_all.(str_md{md}) =  mdls_fit.(str_at{at}).asd.pVal_FDR.(str_md{md}) & mdls_fit.(str_at{at}).ctr.pVal_FDR.(str_md{md}); 
    end
end



pt_md= {'linear', 'quadratic', 'cubic'};
pt_atlas = {'MSALL', 'FsAnat'};
for at = 1:numel(pt_atlas)
    for  md = 1:numel(pt_md)
        
        vals_2plot = 1+mdls_fit.(pt_atlas{at}).pVal_FDR_all.(pt_md{md}) ;
        
        if do_plots_FS
            if strcmp(pt_atlas{at},'MSALL')
                str_cmd.(pt_atlas{at}).(pt_md{md}) = surf_plot_FS_MSALL(vals_2plot, [], ['FDR_sig_mdl',num2str(md),'.',pt_atlas{at}], 1,0 );                
            elseif strcmp(pt_atlas{at},'FsAnat')                
                str_cmd.(pt_atlas{at}).(pt_md{md}) = surf_plot_FS_FsAnat(vals_2plot, [], ['FDR_sig_mdl',num2str(md),'.',pt_atlas{at}], 1,0 );                
            end
        end
    
        if do_plots_MAT
            surf_plot_matlab(vals_2plot(1:end/2)+1, pt_atlas{at}, 'r'); set(gcf,'color','w'); title([pt_atlas{at}, 'FDR sign. mdl ', pt_md{md}])
            surf_plot_matlab(vals_2plot(1+end/2:end)+1, pt_atlas{at}, 'l'); set(gcf,'color','w'); title([pt_atlas{at}, ' FDR sign. mdl ', pt_md{md}])
        end    
        
    end
end

save('S02_data.mat', 'mdls_fit', 'CT')

%%%%%%%%%
%% Report
%%%%%%%%%

disp('MSALL singificant areas for: ')
for  md = 1:numel(pt_md)
    disp(['     order',num2str(md),': ', num2str(nnz(mdls_fit.MSALL.pVal_FDR_all.(pt_md{md}))), ' out of 360 areas']);
end
disp('FsAnat singificant areas for: ')
for  md = 1:numel(pt_md)
    disp(['     order',num2str(md),': ', num2str(nnz(mdls_fit.FsAnat.pVal_FDR_all.(pt_md{md}))), ' out of 64 areas']);
end
