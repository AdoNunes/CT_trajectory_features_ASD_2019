%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script estimates group differences between rates of change extracted
%   from the higest order age coefficient of a given model.
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

%% for one site only
one_site=1;
if one_site
    cin = T.siteID=='NYU';

    for md = 1:numel(str_md)
        for at = 1:numel(str_at)
            CT.(str_at{at}).(str_md{md})(~cin,:)= [];
        end
    end

    T(~cin,:) = [];
end


%% get subsamples
age= {};
age{1}= T.age(T.group== 'asd');
age{2}= T.age(T.group== 'ctr');

group{1} = find(T.group== 'asd');
group{2} = find(T.group== 'ctr');

binranges = 5.9:2.5:30;
K = 50; % num subj in each subsample

m = 100000; % num subsamples to generate
num_subs = 40;


for g = 1:2  % group
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
    Idx_max(:,:,g) = Idx(ii,:);
    
end
 
        
A = unique(Idx_max(1:num_subs,1:K,1));
B = unique(Idx_max(1:num_subs,1:K,2));
out1 = [A, histc(Idx_max(1:num_subs,1:K,1), A)];
out2 = [B, histc(Idx_max(1:num_subs,1:K,2), B)];
rep1 = sum(out1(:, 2:end), 2)/num_subs;
rep2 = sum(out2(:, 2:end), 2)/num_subs;
[h p ] = ttest2(rep1, rep2);


%% Fit models
mdls_coeff = struct;

for g = 1:2 % group
    for s = 1:num_subs % group subsamples
        tmp = Idx_max(s,:,g);
        age_sub = age{g}(tmp);
        
        mdl_age.linear    = [                      age_sub ones(K,1) ];
        mdl_age.quadratic = [           age_sub.^2 age_sub ones(K,1) ];
        mdl_age.cubic     = [age_sub.^3 age_sub.^2 age_sub ones(K,1) ];
        
        for at = 1:numel(str_at) % atlas  
            for md = 1:numel(str_md)% model
                for k = 1:size(CT.(str_at{at}).(str_md{md}),2) % atlas areas
                 
                    mdls_coeff.(str_at{at}).(str_md{md})(s,k,g,:) = regress(CT.(str_at{at}).(str_md{md})(group{g}(tmp),k), mdl_age.(str_md{md}) );
                end
            end
        end
    end
end

%% exclude areas without a good model fit

for md = 1:numel(str_md)
    for at = 1:numel(str_at)
        mdls_coeff.(str_at{at}).(str_md{md})(:,~mdls_fit.(str_at{at}).pVal_FDR_all.(str_md{md}),:,:)= [];
    end
end


%% statistical testing

C_MSALL   = [mdls_coeff.MSALL.linear(:,:,:, 1), mdls_coeff.MSALL.quadratic(:,:,:, 1), mdls_coeff.MSALL.cubic(:,:,:, 1)];
C_FsAnat  = [mdls_coeff.FsAnat.linear(:,:,:, 1), mdls_coeff.FsAnat.quadratic(:,:,:, 1), mdls_coeff.FsAnat.cubic(:,:,:, 1)];

     
option.num_boot             = 5000;
option.num_perm             = 10000; 
option.method               = 1; 
option.meancentering_type   = 1;
option.stacked_designdata   = [1 -1]';

dmat{1}     = C_MSALL(:,:,1);
dmat{2}     = C_MSALL(:,:,2);
num_subj(1) = size(dmat{1},1);
num_subj(2) = size(dmat{2},1);
out.MSALL   = pls_analysis(dmat,num_subj,1,option);
p_MSALL     = out.MSALL.perm_result.sprob
figure, hist(out.MSALL.boot_result.compare_u,20);

dmat{1}     = C_FsAnat(:,:,1);
dmat{2}     = C_FsAnat(:,:,2);
num_subj(1) = size(dmat{1},1);
num_subj(2) = size(dmat{2},1);
out.FsAnat  = pls_analysis(dmat,num_subj,1,option);
p_ANAT      = out.FsAnat.perm_result.sprob
figure, hist(out.FsAnat.boot_result.compare_u);

%% brain plots


pt_md= {'linear', 'quadratic', 'cubic'};
pt_atlas = {'MSALL', 'FsAnat'};


Z_vals_all = struct;
for at = 1:numel(pt_atlas)
    inx = 1;
    figure,
    for  md = 1:numel(pt_md)
 
        at_size =  size(mdls_coeff.(pt_atlas{at}).(pt_md{md})(:,:,:, 1),2); 
        Z_vals_all.(pt_atlas{at}).(pt_md{md}) = out.(pt_atlas{at}).boot_result.compare_u(inx:inx+at_size-1);
        inx = inx + at_size;
        
        subplot(1,3,md), hist( Z_vals_all.(pt_atlas{at}).(pt_md{md})), title([pt_atlas{at},' ',pt_md{md}])
        
    end
end


for at = 1:numel(pt_atlas)
    for  md = 1:numel(pt_md)
        
        Z_vals =  Z_vals_all.(pt_atlas{at}).(pt_md{md});
        if max(Z_vals)< abs(min(Z_vals))
            thr= - percentile(Z_vals,99);
            % thr =  -5.7; % 5.7 - 11 thr color
        else
            thr= - percentile(Z_vals,.1);
           % thr =  5.7; % 5.7 - 11 thr color
        end
        
        vals_2plot = zeros(size(mdls_fit.(pt_atlas{at}).pVal_FDR_all.(str_md{md})));
        vals_2plot(mdls_fit.(pt_atlas{at}).pVal_FDR_all.(str_md{md}))=Z_vals;
 
        
        if do_plots_FS
            if strcmp(pt_atlas{at},'MSALL')
                str_cmd.(pt_atlas{at}).(pt_md{md}) = surf_plot_FS_MSALL(vals_2plot,  thr, ['Zscores_mdl',num2str(md),'.',pt_atlas{at}], 1,0 );
            elseif strcmp(pt_atlas{at},'FsAnat')
                str_cmd.(pt_atlas{at}).(pt_md{md}) = surf_plot_FS_FsAnat(vals_2plot, thr, ['Zscores_mdl',num2str(md),'.',pt_atlas{at}], 1,0 );
            end
        end
        if do_plots_MAT
            surf_plot_matlab(vals_2plot(    1:end/2)+1, pt_atlas{at}, 'r',thr); set(gcf,'color','w'); title([pt_atlas{at}, ' Zscores mdl ', pt_md{md}])
            surf_plot_matlab(vals_2plot(1+end/2:end)+1, pt_atlas{at}, 'l',thr); set(gcf,'color','w'); title([pt_atlas{at}, ' Zscores mdl ', pt_md{md}])
        end
    end
end


%% plot curvatures


pt_md    = {'linear', 'quadratic', 'cubic'};
pt_atlas = {'MSALL', 'FsAnat'};
pt_age   = [6:.1:30]';
at       = 1;
md       = 1;

Z_vals = Z_vals_all.(pt_atlas{at}).(pt_md{md});

if max(Z_vals)< abs(min(Z_vals))
    [v,sort_inx] = sort(Z_vals);
else
    [v,sort_inx] = sort(-Z_vals);
end

vv = 9;
pt_area  = sort_inx(vv); % MSMALL 1: 84(1), 174(6); 2: 78(1), 247(5); 3: 3(4), 181(9)

ensamble = squeeze(mdls_coeff.(pt_atlas{at}).(pt_md{md})(:,pt_area,:, :));

ptage.linear    = [                      pt_age, ones(size(pt_age,1),1)];
ptage.quadratic = [           pt_age.^2, pt_age, ones(size(pt_age,1),1)];
ptage.cubic     = [pt_age.^3, pt_age.^2, pt_age, ones(size(pt_age,1),1)];

ensam_traj        = [];
ensam_traj(:,:,1) =  squeeze(ensamble(:,1,:)) * ptage.(pt_md{md})';
ensam_traj(:,:,2) =  squeeze(ensamble(:,2,:)) * ptage.(pt_md{md})';     
    

ensam_traj_m    = squeeze(mean(ensam_traj,1));
ensam_traj_std  = squeeze(std(ensam_traj,[],1));

CI_u =   ensam_traj_std+ ensam_traj_m;
CI_l =  -ensam_traj_std+ ensam_traj_m;

figure, hold on,
 
[ph,msg]=jbfill(pt_age',CI_u(:,1)',CI_l(:,1)',[155 175 228]/256,'b',1,1);
alpha(0.1)
[ph,msg]=jbfill(pt_age',CI_u(:,2)',CI_l(:,2)',[237 139 135]/256,'r',1,1);
hold on
alpha(0.2)

plot(pt_age,ensam_traj_m(:,2), 'r', 'LineWidth', 5)
plot(pt_age,ensam_traj_m(:,1),'b', 'LineWidth', 5)
xlim([6 30])
 
set(gcf,'Position', [1348 884 353 173])
set(gca,'LineWidth',2)
set(gcf,'color','w')

set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])

plot_allsbj = 1;
if plot_allsbj    
    g1 = T.group == 'asd';    
     
    pt_CT = CT.(pt_atlas{at}).(pt_md{md})(:,mdls_fit.(pt_atlas{at}).pVal_FDR_all.(pt_md{md}));
    pt_CT = pt_CT(:,pt_area);
    
    plot(a1(g1),pt_CT(g1), 'b.'), plot(a1(~g1),pt_CT(~g1), 'r.'),
end



  