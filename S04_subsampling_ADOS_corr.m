%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script correlates the rates of change extracted from the ASD group 
%  with the ADOS scores
%   - The CT for this script is the averaged CT from atlases
%   
% To run statistical testing it is required the PLS toolbox:
%   https://www.rotman-baycrest.on.ca/index.php?section=345
%
% To visualize Freesurfer annotations it is necessary to have Freesurfer in
%   the environment path
%
% To visualize the surface plots it is necessary the gifti toolbox:
%   https://www.artefact.tk/software/matlab/gifti/
%
% Adonay Nunes, SFU, Vancouver, Feb 2019
% adonay.s.nunes@gmail.com
% from github: AdoNunes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
addpath('functions')
load('S01_data.mat')
load('S02_data.mat')

do_plots_FS= 0;
do_plots_MAT = 0;

str_md = {'linear', 'quadratic', 'cubic'};
str_at = {'MSALL', 'FsAnat'};


%% Get ADOS scores

str_ados = {'ADOS_COMM', 'ADOS_SOCIAL','ADOS_STERO_BEH'};

ADOS_nans  = isnan(T.ADOS_COMM);
ADOS_zeros = T.ADOS_COMM;
ADOS_zeros(ADOS_nans) = 0;
ADOS = logical(ADOS_zeros);

grp_asd = find(T.group== 'asd' & ADOS);

%% get subsamples
Idx_max = [];


age{1}= T.age(grp_asd);
group{1} = grp_asd;
g = 1;

binranges = 5.9:2.5:30;
K = 20; % num subj in each subsample

m = 100000; % num subsamples to generate
num_subs = 400;


E = zeros(m,1); % Entropy
Idx = zeros(m,K); % Indices
for kk = 1:m
    n = length(age{g});
    idx = randsample(1:n,K,false);
    age_sub = age{g}(idx);
    p = histc(age_sub,binranges);
    p = p./(sum(p));
    e = -nansum(p.*log2(p));
    Idx(kk,:) = idx;
    E(kk) = e;
end

%take highest E
[Emax, ii] = sort(E, 'descend');
Idx_max(:,:,g) = Idx(ii(1:num_subs),:);
 
age_subs = age{1}(Idx_max);
  
mdls_coeff = struct;
ados_subs = [];

for s = 1:num_subs % group subsamples
    tmp = Idx_max(s,:,g);
    age_sub = age{g}(tmp);
    
    mdl_age.linear    = [                      age_sub ones(K,1) ];
    mdl_age.quadratic = [           age_sub.^2 age_sub ones(K,1) ];
    mdl_age.cubic     = [age_sub.^3 age_sub.^2 age_sub ones(K,1) ];
    
    ados_subs(s,:) = [mean(T.ADOS_total(group{g}(tmp))), mean(T.ADOS_COMM(group{g}(tmp))), mean(T.ADOS_SOCIAL(group{g}(tmp))), nanmean(T.ADOS_STERO_BEH(group{g}(tmp)))];
    
    for at = 1:numel(str_at) % atlas
        for md = 1:numel(str_md)% model
            for k = 1:size(CT.(str_at{at}).(str_md{md}),2) % atlas areas
                mdls_coeff.(str_at{at}).(str_md{md})(s,k,g,:) = regress(CT.(str_at{at}).(str_md{md})(group{g}(tmp),k), mdl_age.(str_md{md}));
            end
        end
    end
end


%% exclude areas without good model fit

for md = 1:numel(str_md)
    for at = 1:numel(str_at)
        mdls_coeff.(str_at{at}).(str_md{md})(:,~mdls_fit.(str_at{at}).pVal_FDR_all.(str_md{md}),:,:)= [];
    end
end


%% statistical testing

str_md ={'linear' 'quadratic' 'cubic'};%   {'quadratic'}; % linear' 'quadratic' 'cubic'

for ados = 1:numel(str_ados)
    for md = 1:numel(str_md)
        
        C_MSALL  =  mdls_coeff.MSALL.(str_md{md})(:,:,:, 1);
        C_FsAnat =  mdls_coeff.FsAnat.(str_md{md})(:,:,:, 1);
        
        option.num_boot = 500;
        option.num_perm = 1000;
        option.method = 3;
        option.meancentering_type = 1;
        option.stacked_behavdata = ados_subs(:,ados);%(:,ados);%(randperm(100),2);
        
        
        dmat{1} = C_MSALL;
        num_subj(1) = size(dmat{1},1);
        out.ADOS_MSALL.(str_md{md}).(str_ados{ados}) = pls_analysis(dmat,num_subj,1,option);
        p_ADOS.MSALL.(str_md{md}).(str_ados{ados}) = out.ADOS_MSALL.(str_md{md}).(str_ados{ados}).perm_result.sprob;
        
        p_msall(md,ados) =p_ADOS.MSALL.(str_md{md}).(str_ados{ados})(1);
        
        dmat{1} = C_FsAnat;
        num_subj(1) = size(dmat{1},1);
        out.ADOS_FsAnat.(str_md{md}).(str_ados{ados}) = pls_analysis(dmat,num_subj,1,option);
        p_FsAnat.ANAT.(str_md{md}).(str_ados{ados}) = out.ADOS_FsAnat.(str_md{md}).(str_ados{ados}).perm_result.sprob;
        
        p_fsanat(md,ados) = p_FsAnat.ANAT.(str_md{md}).(str_ados{ados})(1);
        
    end
end

% p-vals
p_msall < 0.05/numel(p_msall)
p_fsanat < 0.05/numel(p_fsanat)

%% brain plots





pt_md= {'quadratic'};
pt_atlas = {'MSALL', 'FsAnat'};


ados = 1;

Z_vals_all = struct;
for at = 1:numel(pt_atlas)
    inx = 1;
    figure,
    for  md = 1:numel(pt_md)
 
        at_size =  size(mdls_coeff.(pt_atlas{at}).(pt_md{md})(:,:,:, 1),2); 
        Z_vals_all.(pt_atlas{at}).(pt_md{md}) = out.(['ADOS_',pt_atlas{at}]).(str_md{md}).(str_ados{ados}).boot_result.compare_u(inx:inx+at_size-1);
        inx = inx + at_size;
        
        subplot(1,3,md), hist( Z_vals_all.(pt_atlas{at}).(pt_md{md})), title([pt_atlas{at},' ',pt_md{md}])
        
    end
end


for at = 1:numel(pt_atlas)
    for  md = 1:numel(pt_md)
        
        Z_vals =  Z_vals_all.(pt_atlas{at}).(pt_md{md});
        if max(Z_vals)< abs(min(Z_vals))
            thr= - percentile(Z_vals,99);
        else
            thr= - percentile(Z_vals,.1);
        end
        
        vals_2plot = zeros(size(mdls_fit.(pt_atlas{at}).pVal_FDR_all.(pt_md{md})));
        vals_2plot(mdls_fit.(pt_atlas{at}).pVal_FDR_all.(pt_md{md}))=Z_vals;
 
        
        if do_plots_FS
            if strcmp(pt_atlas{at},'MSALL')
                str_cmd.(pt_atlas{at}).(pt_md{md}) = surf_plot_FS_MSALL(vals_2plot, thr, ['Zscores_mdl',num2str(md),'.',pt_atlas{at}], 1,0 );
            elseif strcmp(pt_atlas{at},'FsAnat')
                str_cmd.(pt_atlas{at}).(pt_md{md}) = surf_plot_FS_FsAnat(vals_2plot, thr, ['Zscores_mdl',num2str(md),'.',pt_atlas{at}], 1,0 );
            end
        end
        if do_plots_MAT
            surf_plot_matlab(vals_2plot(1:end/2)+1, pt_atlas{at}, 'r',thr); set(gcf,'color','w'); title([pt_atlas{at}, ' Zscores mdl ', pt_md{md}])
            surf_plot_matlab(vals_2plot(1+end/2:end)+1, pt_atlas{at}, 'l',thr); set(gcf,'color','w'); title([pt_atlas{at}, 'Zscores mdl ', pt_md{md}])
        end
    end
end


