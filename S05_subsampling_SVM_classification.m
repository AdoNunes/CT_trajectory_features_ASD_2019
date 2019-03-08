%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates a traning and a testing set of ASD and CTR
% participants. A SVM is trained on the rate of change of a given model, and
% the accuracy is tested on the testing set. The SVM is trained for each of
% the 3 models, and accuracy is calculated on 500 cross-validation runs. 
%   - The CT for this script is the averaged CT from atlases
%
% Adonay Nunes, SFU, Vancouver, Feb 2019
% adonay.s.nunes@gmail.com
% from github: AdoNunes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
clear 
load('S01_data.mat')
load('S02_data.mat')

do_plots_FS= 0;
do_plots_MAT = 0;

str_md = {'linear', 'quadratic', 'cubic'};
str_at = {'MSALL', 'FsAnat'};


%% get subsamples
accu =[];
good_samples = {};
inx = 0;
for rep = 1:500
    binranges = 5.9:2.5:30;
    n_asd = nnz(T.group== 'asd');
    n_ctr = nnz(T.group== 'ctr');
    
    
   age_a = T.age(T.group== 'asd');
    
   [ n, edges, bins] = histcounts(age_a,binranges);
   inbin_A1 =[];
   inbin_A2 =[];
    for bn = 1:numel(n)
        aux = find(bins == bn);
        inx_rnd = randperm(size(aux,1));
        inbin_A1 = cat(1,inbin_A1,  aux(inx_rnd(1:round(end/2))));
        inbin_A2 = cat(1,inbin_A2,  aux(inx_rnd(1+round(end/2):end)));
    end
    
    age_c = T.age(T.group== 'ctr');
    [ n, edges, bins] = histcounts(age_c,binranges);
      inbin_C1 =[];
   inbin_C2 =[];
    for bn = 1:numel(n)
        aux = find(bins == bn);
        inx_rnd = randperm(size(aux,1));
        inbin_C1 = cat(1,inbin_C1,aux(inx_rnd(1:round(end/2))));
        inbin_C2 = cat(1,inbin_C2,aux(inx_rnd(1+round(end/2):end)));
    end
    
    age= {};
    age{1}= age_a(inbin_A1);
    age{2}= age_a(inbin_A2);
    age{3}= age_c(inbin_C1);
    age{4}= age_c(inbin_C2);
    
    
    
    group_a = find(T.group== 'asd');
    group_c = find(T.group== 'ctr');
    
    group= {};
    group{1}= group_a(inbin_A1);
    group{2}= group_a(inbin_A2);
    group{3}= group_c(inbin_C1);
    group{4}= group_c(inbin_C2);
    
    binranges = 5.9:2.5:30;
    K = 80; % num subj in each subsample
    
    m = 100000; % num subsamples to generate
    num_subs = 70;
    
    
    for g = 1:numel(group)  % group
        E = zeros(m,1); % Entropy
        Idx = zeros(m,K); % Indeces
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
    

    mdls_coeff = struct;
    
    for g = 1:numel(group)  % group
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
    
    
    %% SVM training 
    
    
    test_atl= 'MSALL';%'FsAnat';%
    
    for md = 1:numel(str_md)
        
        test_mdl = str_md{md};
       
        train_vals = [mdls_coeff.(test_atl).(test_mdl)(:,:,1,1); mdls_coeff.(test_atl).(test_mdl)(:,:,3,1)];
        nsbs = size(train_vals,1)/2;
        grp =[];
        grp(1:nsbs) = -1;
        grp(nsbs+1:nsbs*2) = 1;
        
        SVMModel = fitcsvm(train_vals,grp, 'Standardize',0);
        
        pred_vals= [mdls_coeff.(test_atl).(test_mdl)(:,:,2,1); mdls_coeff.(test_atl).(test_mdl)(:,:,4,1)];
        
        label = predict(SVMModel,pred_vals);
        sensi(rep,md) = nnz(label(1:end/2) ==-1)/nsbs;
        spesi(rep,md) = nnz(label(1+end/2:end) ==1)/nsbs;
        accu(rep,md) = (sensi(rep,md) + spesi(rep,md))/2
        
    end
     
end

v_str= 'spesi';%'sensi';% 'accu';
eval(['vplot = ',v_str,';'])
figure, hold on, 

histogram(vplot(:,2),.1:.05:1,'FaceColor','b')
histogram(vplot(:,1),.1:.05:1,'FaceColor','r')
histogram(vplot(:,3),.1:.05:1,'FaceColor','g')

set(gca,'LineWidth',1, 'FontSize', 14)
set(gca,'Position', [0.13 0.11 0.77 0.815])
set(gcf,'color','w')
legend( 'quadratic', 'linear','cubic')

title(v_str)
mean(vplot)
max(vplot)

save('s05_accu_sens_spes_fsanat_550rep.mat', 'accu','sensi','spesi' )

  