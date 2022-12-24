This Examples folder contains examples about how to use the *mati* package for specific applications. Each example contains the shared code that was used in one or more published research papers. 

# Eg01_ADCSpectrum_1compt
This example shows how to calculate theoretical ADC spectra of restricted diffusion inside some simple geometries such as parallel planes, cylinders, spheres, and spherical shells. 

Reference: 
1.	Jiang X, Li H, Xie J, Zhao P, Gore JC, Xu J. Quantification of cell size using temporal diffusion spectroscopy. Magn Reson Med. 2016;75(3):1076-85. Epub 20150404. doi: 10.1002/mrm.25684. PubMed PMID: 25845851; PMCID: PMC4592783.
2.	Xu J, Jiang X, Li H, Arlinghaus LR, McKinley ET, Devan SP, Hardy BM, Xie J, Kang H, Chakravarthy AB, Gore JC. Magnetic resonance imaging of mean cell size in human breast tumors. Magn Reson Med. 2020;83(6):2002-14. Epub 20191125. doi: 10.1002/mrm.28056. PubMed PMID: 31765494; PMCID: PMC7047520.

# Eg02_IMPULSED
This example shows how to (1) synthesize dMRI signal data using the IMPULSED model and (2) how to fit the IMPULSED model to dMRI signals to extract microstructural parameters such as d (mean cell size), vin (apparent intracellular volume fraction), Din (intracellular diffusion coefficient), and Dex (extracellular diffusion coefficient). 

Reference: 
1.	Jiang X, Li H, Xie J, Zhao P, Gore JC, Xu J. Quantification of cell size using temporal diffusion spectroscopy. Magn Reson Med. 2016;75(3):1076-85. Epub 20150404. doi: 10.1002/mrm.25684. PubMed PMID: 25845851; PMCID: PMC4592783.
2.	Xu J, Jiang X, Li H, Arlinghaus LR, McKinley ET, Devan SP, Hardy BM, Xie J, Kang H, Chakravarthy AB, Gore JC. Magnetic resonance imaging of mean cell size in human breast tumors. Magn Reson Med. 2020;83(6):2002-14. Epub 20191125. doi: 10.1002/mrm.28056. PubMed PMID: 31765494; PMCID: PMC7047520.

# Eg03_MRICytometry
This example shows how to (1) synthesize dMRI signal data using the MRI-cytometry model and (2) how to fit the MRI-cytometry model to dMRI signals to extract microstructural parameters such as distributions of d (cell size), Din (intracellular diffusion coefficient), and Dex (extracellular diffusion coefficient); and other parameters such as vin (apparent intracellular volume fraction). 

Reference:
1.	Xu J, Jiang X, Devan SP, Arlinghaus LR, McKinley ET, Xie J, Zu Z, Wang Q, Chakravarthy AB, Wang Y, Gore JC. MRI-cytometry: Mapping nonparametric cell size distributions using diffusion MRI. Magn Reson Med. 2021;85(2):748-61. Epub 20200916. doi: 10.1002/mrm.28454. PubMed PMID: 32936478; PMCID: PMC7722100.

# Eg04_JOINT
This example shows how to use the JOINT model to simultaneously fit d (mean cell size), Din (intracellular diffusion coefficient), vin (apparent intracellular volume fraction), and transcytolemmal water exchange described with parameters such as tau_in (mean intracellular water lifetime) and/or Pm (membrane permeability). Sample data from computer simulations and cultured cell experiments in vitro are included in this example. 

Reference:
1.	Jiang X, Devan SP, Xie J, Gore JC, Xu J. Improving MR cell size imaging by inclusion of transcytolemmal water exchange. NMR Biomed. 2022;35(12):e4799. Epub 20220728. doi: 10.1002/nbm.4799. PubMed PMID: 35794795. 

# Eg05_mpMRI_Classifier
This example shows how to perform a two-step approach (Lasso and Ridge regressions) on mpMRI (multi-parametric MRI) with any arbitrary number of parameters. The results return (1) the dominant (most important) MRI parameters and (2) the best combination of the dominant MRI parameters as a binary classifier. This is suitable to find out the most important parameters from an array of parameters to perform a binary classification, such as the differentiation of brain tumors from radiation-induced necrosis. 

Reference:
1.	Devan SP, Jiang X, Kang H, Luo G, Xie J, Zu Z, Stokes AM, Gore JC, McKnight CD, Kirschner AN, Xu J. Towards differentiation of brain tumor from radiation necrosis using multi-parametric MRI: Preliminary results at 4.7 T using rodent models. Magn Reson Imaging. 2022;94:144-50. Epub 20221006. doi: 10.1016/j.mri.2022.10.002. PubMed PMID: 36209946.

# Eg06_random_packing
This example creates a 3D structure of randomly packed spheres with a distribution of sphere diameters. The intra-spherical volume fraction can exceed 70% or 80% by adjusting sphere swelling ratios. 

Reference:
1.	Xu J, Xie J, Semmineh NB, Devan SP, Jiang X, Gore JC. Diffusion time dependency of extracellular diffusion. (under review)


# Comments or questions?
 
Please send your comments or questions to Junzhong (JZ) Xu (junzhong.xu@vanderbilt.edu) or Xiaoyu Jiang (xiaoyu.jiang@vumc.org)

