# Cort_Thick_age-changes_Autism_2018
Age-related changes in Cortical Thickness in Autism


This set of scprips takes cortical thickness (CT) values from the ABIDE repository and :

[ in S01_getGroups.m ] 
- Excludes outliers per center.
- Creates a group of ASD and TD, balanced by number of subjects, age and the ratio male/female per center.


[ in S02_ModelFitting.m ]
- Center variability is taken out while accounting for group, age and age-group interactions in a linear, quadratic and cubic models. In the following steps, CT values corrected for center variability are used. 
- Significance of a model fit using a deviance test is measured for each area of the CT parcellation, areas that do not fit on a group level are masked out. ASD and TD are tested separately.


[ in S03_subsampling_analysis_Hemi.m & S03_subsampling_analysis_Atlases.m ]
- For each group, subsamples of subjects are created. This subsamples have the flattest age distribution possible. 
- For each subgroup, a linear, quadratic and cubic models are fitted and the highest order coefficient (linear, quadratic and cubic trends) are extracted for statistical testing between groups. These coefficients are invariant across age.
- Using Partial Least Squares (PLS), the coefficients at each area of teh parcellation are tested for significance. 


[ in S04_subsampling_ADOS_corr.m ]
- The coefficients derived from ASD subgroups and the mean ADOS scores from the subgrup of subjects are tested for correlations. Each model is correlated with mean ADOS with PLS, the p-vals are bonferroni corrected. 
- The three models were analyzed separately, for simplicity, as only one highly correlated.


[ in S05_subsampling_SVM_classification.m ]
- A support vector machine(SVM) was trained and tested with 500 cross-validations. 
- For each run, the ASD and TD groups where splitted in half. In the training set, coefficients were extracted using one parcellation, and for each model coefficient a SVM was trained to classify age-related changes velonging to ASD or TD samples. Then, with the other dataset, the accuracy, sensistivity and specificity was tested. 
