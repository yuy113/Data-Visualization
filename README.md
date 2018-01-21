# Data-Visualization

## Introduction and description of Data set 
Explore and try identifying important metabolites potentially associated with cardiovascular disease in one matched case control study(CATHGEN) through machine learning algorithms and data visualization techniques.


CATHGEN medical study was a matched case-control study to discover prognostic biomarkersin blood plasma for near-term cardiovascular events. Subjects were selected from the CATHGEN project which collected peripheral blood samples from consenting research subjects undergoing cardiac catheterization at Duke University Medical Center from 2001 through 2011.  68 cases were selected from among individuals who had a major adverse cardiac event (MACE) within two years following the time of their sample collection. In a 1:1 matched study design, 68 controls were selected from individuals who were MACE-free for the two years following sample collection and were matched to cases on age, gender, race/ethnicity and severity of coronary artery disease. High-content mass spectrometry and multiplexed immunoassay-based techniques were employed to quantify 625 proteins and metabolites from each subject serum specimen. Comprehensive metabolite profiling of the individual samples was based on a combination of four platforms employing mass spectrometry (MS) based techniques to profile lipids, fatty acids, amino acids, sugars and other metabolites. 


This dataset contains 136 subjects, 472 features which are the concentrations of metabolites where the number of observations(N) is smaller than the number of features(D), and also the class labeling for patients(1= cardiovascular disease; 2=Healthy subjects). In order to reduce selection bias in clinical study, the exact names of those 472 metabolite features are not given in the dataset, instead named by order number-Metabolite1-472.


## Data preprocessing and Missing Data Imputation


Assume that the values of each 33 metabolite feature with percentage of 0 values less than or equal to 5% are following log normal distribution, and that mean and variance of log Normal distribution equal to sample mean, sample variance from nonzero values for each feature respectively. Thus for each of these 33 features, the values of 0 are imputed by the random values with less than lower 2.5% quantile of log normal distributions.


## Computation and Exploratory Data Analysis

The 406 imputed metabolite features will be firstly log transformed and then normalize by subtracting their means divided by corresponding standard deviations. The dimensional reduction techinques will be applied to reduce the size of features for the data with more features(406) than number of observations(136) such as sparse principal component analysis (sPCA), t-SNE technique, Partial Least Square(PLS) etc. t-SNE analysis is performed on Python Scikit-learn manifold TSNE function(\url{http://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html}). The outputs from t-SNE are either tightly clustered together or very sparse distributed without separation across diseased and healthy groups by tuning the values of two parameters in TSNE function-perplexity and learning rate. Thus the results from t-SNE will not shown in the webpages.

The first 20 principal components from sparse PCA analysis on all 406 imputed normalized metabolite features are saved and run data visualzation in RadViz and interactive visualization between scatterplot and Histograms on any two of these 20 principal components across diseased and healthy groups. Similary the first 20 components from PLS analysis will also saved and are presented in RadViz plot and interactive visualization between Scatterplot and Histograms on any two of these 20 principal components across diseased and healthy groups.

Also paired t tests between paired diseased and healthy groups are performed on each of those 406 imputed normalized metabolite features. The 20 features with smallest P values in paired t tests are selected. 

The R code file for both data preprocessing and imputation and the analysis results for sPCA, PLS and paired t tests is data_processing_project_yyb.R.


## Interactive Data Visualization


RadViz visualization will be used to explore the association patterns among all normalized imputed 406 metabolite features in all the observations, and in diseased and healthy groups. Identify several strong association groups of metabolite features from 20 metabolite features with lowest paired t test P values(the data combined with diseased group labels saved in \url{feature_20lowestpvalue_diseasegroup.csv}) or associated components from first 20 components in PLS method (the data combined with diseased group labels saved in \url{PLS_20components_diseasegroup.csv}), and associated principal components from first 20 principal components in sPCA method (the data combined with diseased group labels saved in \url{sPCA_20components_diseasegroup.csv}), by RadVis visualization.

The RadViz visualization on 20 metabolite features with lowest paired t test P values across diseased and healthy groups are shown in the webpage-http://www-edlab.cs.umass.edu/~yyao/Radviz_feature_20lowestpvalue.html


The RadViz visualization on first 20 components in Partial Least Square(PLS) method across diseased and healthy groups are shown in the webpage-http://www-edlab.cs.umass.edu/~yyao/Radviz_20components_pls.html


The RadViz visualization on first 20 principal components in sparse principal component analysis (sPCA) method across diseased and healthy groups are shown in the webpage-\url{http://www-edlab.cs.umass.edu/~yyao/Radviz_20components_sPCA.html}


Association and clustering patterns among 20 metabolite features with lowest paired t test p values will be presented by parallel coordinates. This visualization is presented in the webpage--\url{http://www-edlab.cs.umass.edu/~yyao/parallel_20features_ttest.html}.


The interactive techniques in this parallel coordinate visualization include selection of lines and moving the 20 axises of metabolite features in order to explore the details of any clustering and association patterns among selected 20 metabolite features.


The interactive scatterplot and histograms of two variables in the scatterplot and also piechart of numbers of subjects in both diseased and healthy groups will be presented to explore the distributions and associations among either metabolite features or principal components from sPCA method or components from PLS method, and their association with diseased labels.


The interactive scatterplot and two histograms and piechart for 20 metabolite features with lowest paired t test P values across diseased and healthy groups is shown in the webpage-\url{http://www-edlab.cs.umass.edu/~yyao/hist_scatter_20features_ttest.html}.

The interactive scatterplot and two histograms and piechart for first 20 components in Partial Least Square(PLS) method across diseased and healthy groups is shown in the webpage-\url{http://www-edlab.cs.umass.edu/~yyao/hist_scatter_20components_pls.html}.


The interactive scatterplot and two histograms and piechart for first 20 principal components in sparse principal component analysis (sPCA) method across diseased and healthy groups is shown in the webpage-\url{http://www-edlab.cs.umass.edu/~yyao/hist_scatter_20components_sPCA.html}.



### Summary Results from Data Visualization


Combined the visualization results from RadViz, parallel coordinate, and interactive scatterplot, two histograms and piechart for 20 metabolite features with lowest paired t test P values across diseased and healthy groups, each of the 20 metabolite features contributes very weakly to prediction/classification of diseased groups from not sepearating results from both RadViz and scatterplot between diseased and healthy groups. It appears to exist weakly positively or negatively association within some subgroups of these metabolite features based on parallel coordinate visualization. The values of most transformed metabolite features are not regularly distributed, either bimodal or with skewed distribution.


Combined the visualization results from RadViz and interactive scatterplot, two histograms and piechart for first 20 components in Partial Least Square(PLS) method across diseased and healthy groups, each of the 20 PLS components contributes very weakly to prediction/classification of diseased groups from not sepearating results from both RadViz and scatterplot between diseased and healthy groups. The first few PLS components have relatively well separation in scatterplot between diseased and healthy groups suggesting they have relatively more effect on classification of diseased groups than other PLS components. This observation can be explained by partial least square method which focuses on the degree of correlation beween the diseased labels and the metabolite features. The values of most PLS components are not regularly distributed, either bimodal or with skewed distribution.


Combined the visualization results from RadViz and interactive scatterplot, two histograms and piechart for first 20 principal components in sparse principal component analysis (sPCA) method across diseased and healthy groups, each of the 20 sPCA components contributes very weakly to prediction/classification of diseased groups from not sepearating results from both RadViz and scatterplot between diseased and healthy groups.The values of most sPCA components are not regularly distributed, either bimodal or with skewed distribution.


Overally with such a large number of metabolite features(406), due to the very small number of observation of subjects participated in the clinical study and the very weak or no association between each metabolite feature and diseased labels, the metabolite features themselves may be not sufficient to predict/classify the diseased groups, more observation or more features such as prognostic features of the patients such as age, specific disease clinical measures need to be added and combine with those biological features to have enough power to predict the diseased labels for the subjects.

