# Infrared-fingerprint-Health-diagnosis

## ABSTRACT

**Background** NAFLD are characterized by lipid accumulation in hepatocytes that may evolve toward NASH with the appearance of inflammation and hepatocyte suffering that contribute to the development of hepatic fibrosis, cirrhosis and/or hepatocellular carcinoma (HCC). Histological evaluation of liver biopsy is the gold standard for NASH diagnosis but cannot be performed in all patients, and is not repeatable for a follow-up. Our objective, was to evaluate whether a non-invasive and easy to perform test, based on mid-infrared vibrational spectroscopy (MIR) and machine learning, may give significant information on liver status in a NAFLD mice model.

**Methods** Serum samples and blood smears were obtained from two independent sets of mice receiving control and/or high fat and/or high carbohydrate diet and/or parenteral iron. They have been respectively used as calibration set and validation set for machine learning algorithms.

**Results** Algorithm obtained by training a part of the set 1 (75% randomly chosen) led to predict significant hepatic triglyceride accumulation with an AUROC of 0.88 in the other part of set 1 (25%) and 0.90 in the second set, used as external validation set. Chemical functional groups identified in the spectra and contributing to the signature were associated to lipid ester, carbohydrate content and phospholipids. Regarding histologically quantified hepatic steatosis the reproducibility was limited in the external set. Hepatic inflammation prediction by MIR, evaluated on blood smear, led to an AUROC of 0.81 in the first set and 0.67 in the external set. Chemical groups identified during evaluation of inflammation were associated to lipid peroxidation, protein secondary structure variations and protein aggregation.

**Conclusions** MIR spectroscopy performed in serum, combined with machine learning, is an efficient tool to evaluate lipid accumulation in the liver. Using MIR for the evaluation of hepatic inflammation on blood smear is promising.

## CODE

1.1 PREPROCESSING OF SERUM: describe the preprocess used (normalisation and data compression methods) \n
1.2 Present the classification of the hepatic triglyceride using random forest and associated figures. \n
1.3 Present the classification of the hepatic steatosis using random forest and associated figures. 

2.1 PREPROCESSING OF BLOOD SMEAT: describe the preprocess used (normalisation and data compression methods)
2.2 Present the determination of the hepatic liver statu using clustering method and associated figures
2.3 Present the classification of the hepatic inflammation using random forest and associated figures.
