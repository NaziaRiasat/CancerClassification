#Molecular Classification of Cancer Class Discovery and Class Prediction by Gene Expression Monitoring

The project examined a small gene expression dataset and attempted to classify leukemia patients into one of two classes. Golub et. al (1999) measured the expression levels of 7,129 genes for 72 patients diagnosed with acute lymphoblastic leukemia (ALL) or acute myeloid leukemia (AML).  

Firstly, Principal Component Analysis was applied to the data to reduce its dimensionality. There are 7129 numerical features, The results suggested to use of 36 principal components that define around 90 percent of the variance among the data. 

Secondly, the Naïve Bayes Algorithm trained on data from 38 of these patients to predict the class of acute leukemia for unseen patients. This model correctly classified 33 of 34 patients in an independent test. As a research tool, the model selected genes with accession numbers M19507_at, M27891_at, M96326_rna1_at, and Y00787_s_at as being positively associated with AML and negatively associated with ALL. These associations were also the most prominent in PCA analysis.
