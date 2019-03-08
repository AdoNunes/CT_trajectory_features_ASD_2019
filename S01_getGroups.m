%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates an ASD and CTR groups with:
%   more than 10 subjects per center
%   without CT outliers per center
%   non-significant age differences
%
% Data is obtained from a table with sections:
%   ID (subj ID), age( numeric), group(categorical: asd-ctr, 
%   sex(categorical:male-female), siteID(categorical: center_names)
%
% Adonay Nunes, SFU, Vancouver, Feb 2019
% adonay.s.nunes@gmail.com
% from github: AdoNunes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
load('Table_CT_ABIDE.mat')
 
T(T.age>30,:) = [];
T(T.age< 5.9,:) = [];  
 
centers = unique(T.siteID);

sbj=[];
for c= 1 : numel(centers) 
    sbj(c,1) = numel(T.age(T.siteID ==centers(c) & T.group=='asd' ));
    sbj(c,2) = numel(T.age(T.siteID ==centers(c) & T.group=='ctr' ) );
end

min_subj = 10;
too_few = find(sum(sbj < min_subj,2));

for tf = 1:numel(too_few)
    T(T.siteID == centers(too_few(tf)),:) = [];
end

%% get CT

tmp = cat(3, T.cortstats{:}); 
CT.MSALL.raw = squeeze(tmp(:,4,:))';

tmp = cat(3, T.anatstats{:}); 
CT.FsAnat.raw = squeeze(tmp(:,4,:))';

tmp = cat(3,T.meanthick{:});
CT.H.raw = squeeze(tmp)';

%% Remove outliers
out_all = nan(size(T,1),2);
centers = unique(T.siteID);
Z_thres = 3.5;
for c= 1 : numel(centers) 
    
    cent_vec =  T.siteID == centers(c);
    
    ct_1 = zscore(CT.MSALL.raw(cent_vec,:));
    out_1 = sum(abs(ct_1)>Z_thres,2);
    
    ct_2 = zscore(CT.FsAnat.raw(cent_vec,:));
    out_2 = sum(abs(ct_2)>Z_thres,2);
    
    out_all(cent_vec,:) = [out_1,out_2];
end

out_all_ind = find(out_all(:,1)>1 |out_all(:,2)>1);

T(out_all_ind,:) = [];
CT.MSALL.raw(out_all_ind,:) = [];
CT.FsAnat.raw(out_all_ind,:)   = [];
CT.H.raw(out_all_ind,:)  = [];
 
%% check minimum
centers = unique(T.siteID);

sbj = [];
for c= 1 : numel(centers) 
    sbj(c,1) = numel(T.age(T.siteID ==centers(c) & T.group=='asd' ));%& ~isnan(T.ADOS_total) & ~isnan(T.ADOS_SOCIAL)) )
    sbj(c,2) = numel(T.age(T.siteID ==centers(c) & T.group=='ctr' ) );
end

min_subj = 10;
too_few = find(sum(sbj < min_subj,2));

for tf = 1:numel(too_few)
    center_out = T.siteID == centers(too_few(tf));
    T(center_out,:) = [];
    CT.MSALL.raw(center_out,:) = [];
    CT.FsAnat.raw(center_out,:)   = [];
    CT.H.raw(center_out,:)  = [];
    
end
 
%% equalize groups
centers = unique(T.siteID);

sbj = [];
for c= 1 : numel(centers) 
    sbj(c,1) = numel(T.age(T.siteID ==centers(c) & T.group=='asd' ));%& ~isnan(T.ADOS_total) & ~isnan(T.ADOS_SOCIAL)) )
    sbj(c,2) = numel(T.age(T.siteID ==centers(c) & T.group=='ctr' ) );
    
    sbj_fm(c,1) = numel(T.age(T.siteID ==centers(c) & T.group=='asd' & T.sex=='female' ));%& ~isnan(T.ADOS_total) & ~isnan(T.ADOS_SOCIAL)) )
    sbj_fm(c,2) = numel(T.age(T.siteID ==centers(c) & T.group=='ctr' & T.sex=='female' ) );
end

grp_ratio = sbj(:,1)./sbj(:,2);
min_tol = .13;
sbj_rm = [];
for c= 1 : numel(centers) 
    if grp_ratio(c)>1+min_tol % rm asd
        sbj_rm(c,1) = round(sbj(c,1) - sbj(c,2)* (1+ min_tol));        
    elseif grp_ratio(c)<1-min_tol % rm ctr
        sbj_rm(c,2) = ceil(sbj(c,2) - sbj(c,1)* (1 +min_tol/2));
    end    
end

 
% chi-squared test for ASD and CTR ratios
c1 = sbj(:,1) -sbj_rm(:,1); 
c2 = sbj(:,2) -sbj_rm(:,2);

n1 = c1'; n2 = c2';
N1 = sum(n1); N2 = sum(n2);
p0 = (n1+n2)./(N1+N2);% Pooled estimate of proportion
n10 = N1 * p0; n20 = N2 * p0;% Expected counts under H0 
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
chi2stat = sum((observed-expected).^2 ./ expected,2);
p = 1 - chi2cdf(chi2stat,1) % check if different proportions



%% take out extra subjects - females if possible

fem_4_rej = sbj_fm(:,2) - sbj_fm(:,1);
fem_4_rej(fem_4_rej<0) = 0;
rej_center = {};
for c= 1 : numel(centers)
    if any (sbj_rm(c,:))
        if sbj_rm(c,2)>0 % ctr reject
            fm_rej = find(T.siteID ==centers(c) & T.group=='ctr' & T.sex=='female');
            m_rej = find(T.siteID ==centers(c) & T.group=='ctr' & T.sex=='male' );
            
            nrej = sbj_rm(c,2);
            
            if fem_4_rej(c)>0 &&  fem_4_rej(c)>=nrej % rej females only
                rng(9,'twister')
                aux = fm_rej(randperm(size(fm_rej,1)));
                rej_center{c} =  aux(1:nrej);
                
            elseif fem_4_rej(c)>0 &&  fem_4_rej(c)<nrej % rej females then males
                rng(9,'twister')
                aux_1 = fm_rej(randperm(size(fm_rej,1)));
                rng(9)
                aux_2 = m_rej(randperm(size(m_rej,1)));
                
                rej_center{c} =  [aux_1(1:fem_4_rej(c));aux_2(1:nrej-fem_4_rej(c)) ];
                
            elseif fem_4_rej(c)<=0 % rej males only
                rng(9,'twister')
                aux = m_rej(randperm(size(m_rej,1)));
                rej_center{c} =  aux(1:nrej);
            end
            
        else % asd reject
            nrej = sbj_rm(c,1);
            m_rej = find(T.siteID ==centers(c) & T.group=='asd' & T.sex=='male' ); % rej males only
            rng(9,'twister')
            aux = m_rej(randperm(size(m_rej,1)));
            rej_center{c} =  aux(1:nrej);
        end
    end
end
 
rej_center_all = cell2mat(rej_center');
 
 
CT.MSALL.raw(rej_center_all,:) = []; % subjects  x ROIS
CT.FsAnat.raw(rej_center_all,:) = [];
CT.H.raw(rej_center_all,:) = [];
T(rej_center_all,:) = [];

sbj = [];
for c= 1 : numel(centers) 
   
    sbj(c,1) = numel(T.age(T.siteID ==centers(c) & T.group=='asd' ));%& ~isnan(T.ADOS_total) & ~isnan(T.ADOS_SOCIAL)) )
    sbj(c,2) = numel(T.age(T.siteID ==centers(c) & T.group=='ctr' ) );
    
    ageT(c) = ttest2(T.age(T.siteID ==centers(c) & T.group=='asd'), T.age(T.siteID ==centers(c) & T.group=='ctr'));
    
    sbj_fm(c,1) = numel(T.age(T.siteID ==centers(c) & T.group=='asd' & T.sex=='female' ));%& ~isnan(T.ADOS_total) & ~isnan(T.ADOS_SOCIAL)) )
    sbj_fm(c,2) = numel(T.age(T.siteID ==centers(c) & T.group=='ctr' & T.sex=='female' ) );
end


%% remove any age diff

abins = 5:.5:30;

ageT_inx = find(ageT);

rm_age_out= [];
   
for c= 1 : nnz(ageT)
    
    site = centers(ageT_inx(c));
    
    g1c = T.siteID ==site & T.group=='asd';
    g2c = T.siteID ==site & T.group=='ctr';
    
    aa = histc( T.age(g1c),abins);
    ac = histc( T.age(g2c),abins);
    
    a_diff = aa/sum(aa) - ac/sum(ac);
    a_diff_bin = find( abs(a_diff)> 0.049);
    
 
    for aout = 1:numel(a_diff_bin)
        abin_OI = abins(a_diff_bin(aout));
        
        if a_diff(a_diff_bin(aout)) > 0
            nout = round(a_diff(a_diff_bin(aout))*sum(aa)*(1-min_tol));
            inx_OI = find( T.age>abin_OI-.5 & T.age<abin_OI+.5 & T.siteID ==site & T.group=='asd' );
            
        else
            nout = round(abs(a_diff(a_diff_bin(aout))*sum(ac))*(1-min_tol));
            inx_OI = find( T.age>abin_OI-.5 & T.age<abin_OI+.5  & T.siteID ==site & T.group=='ctr');
            
        end
        rm_age_out =  cat(1, rm_age_out,  inx_OI(1:nout));
        
    end
    T2 = T;
    T2(rm_age_out,:) = [];
    
    [h p ] = ttest2(T2.age(T2.siteID ==site & T2.group=='asd'), T2.age(T2.siteID ==site & T2.group=='ctr'));
end

[h p ] = ttest2(T2.age( T2.group=='asd'), T2.age( T2.group=='ctr')) % age diferent between groups?


CT.MSALL.raw(rm_age_out,:) = []; % subjects  x ROIS
CT.FsAnat.raw(rm_age_out,:) = [];
CT.H.raw(rm_age_out,:) = [];
T(rm_age_out,:) = [];


sbj = [];
for c= 1 : numel(centers) 
   
    sbj(c,1) = numel(T.age(T.siteID ==centers(c) & T.group=='asd' ));%& ~isnan(T.ADOS_total) & ~isnan(T.ADOS_SOCIAL)) )
    sbj(c,2) = numel(T.age(T.siteID ==centers(c) & T.group=='ctr' ) );
    
    sbj_fm(c,1) = numel(T.age(T.siteID ==centers(c) & T.group=='asd' & T.sex=='female' ));%& ~isnan(T.ADOS_total) & ~isnan(T.ADOS_SOCIAL)) )
    sbj_fm(c,2) = numel(T.age(T.siteID ==centers(c) & T.group=='ctr' & T.sex=='female' ) );
end


c1 = sbj(:,1) ; 
c2 = sbj(:,2) ;

% chi-squared test for ASD and CTR ratios
n1 = c1'; n2 = c2';
N1 = sum(n1); N2 = sum(n2);
p0 = (n1+n2)./(N1+N2);% Pooled estimate of proportion
n10 = N1 * p0; n20 = N2 * p0;% Expected counts under H0 
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
chi2stat = sum((observed-expected).^2 ./ expected,2);
p = 1 - chi2cdf(chi2stat,1) % age

% chi-squared test for female proportions per center
N1 = c1; N2 = c2;
n1 = sbj_fm(:,1); n2 = sbj_fm(:,2);
p0 = (n1+n2)./(N1+N2);% Pooled estimate of proportion
n10 = N1 .* p0; n20 = N2 .* p0;% Expected counts under H0 
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
s_0 = ~logical(sum([observed==0;expected==0]));
chi2stat = sum((observed-expected).^2 ./ expected,2);
p = 1 - chi2cdf(chi2stat,1); % female proportions

save('S01_data.mat','T', 'CT'  )

%%%%%%%%%
%% Report
%%%%%%%%%

disp(['Total subjects: ', num2str(size(T,1))])
disp(['ASD subjects: ', num2str(nnz(T.group =='asd'))])
disp(['CTR subjects: ', num2str(nnz(T.group =='ctr'))])
disp(['Number centers: ', num2str(numel(unique(T.siteID)))])


