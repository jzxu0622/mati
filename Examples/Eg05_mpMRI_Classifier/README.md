This example shows how to perform a two-step approach (Lasso and Ridge regressions) on mpMRI (multi-parametric MRI) with arbitrary number of parameters. The results return (1) the dominant (most important) MRI parameters and (2) the best combination of the dominant MRI parameters as a binary classifier. This is suitable to find out the most important parameters from an array of parameters to perform a binary classification, such as the differentiation of brain tumors from radiation-induced necrosis. 

Usage: Run Eg05_mpMRI_Classifier.m

Note: the key function mpMRI_Classifier.m takes input data in a condensed table. 


----------------------------------------------------------------------------
If you use this function in your work, please cite the following paper
1.	Devan SP, Jiang X, Kang H, Luo G, Xie J, Zu Z, Stokes AM, Gore JC, McKnight CD, Kirschner AN, Xu J. Towards differentiation of brain tumor from radiation necrosis using multi-parametric MRI: Preliminary results at 4.7 T using rodent models. Magn Reson Imaging. 2022;94:144-50. Epub 20221006. doi: 10.1016/j.mri.2022.10.002. PubMed PMID: 36209946.
