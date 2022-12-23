# MATI (Microstructural Analysis of Tissues by Imaging)
---------------------------------------------------------------------------------------------------------------
This is a MATLAB package that includes microstructural analysis code developed by Vanderbilt University Institute of Imaging Science. 

Authors:  Junzhong Xu, Xiaoyu Jiang, Sean P. Devan 
---------------------------------------------------------------------------------------------------------------
# Installation 
1. Create a local sandbox folder named "+mati". The "+" indicates a MATLAB namespace. 
2. Create a new MATLAB project using the github repository https://github.com/jzxu0622/mati.git with the local sandbox folder "+mati". If an error occurs, try "git://github.com/jzxu0622/mati" instead. 
3. Make sure the parent folder of +mati is addded to the MATLAB path. 
4. Run the examples in the +mati/Examples to show how mati works in different specific applications

# Contents
The MATI package was written using Objective Oriented Programming and contains a few easy-to-use classes:
1.	mati.PulseSequence is the parent class for all pulse sequence classes. It contains basic pulse sequence parameters such as Nacq (number of data points acquired), TR (repetition time) and TE (echo time). 
2.	mati.SignalModel is the parent class for all model classes. It contains basic signal model parameters such as parmsLabels (name tags of fitting parameters), FcnSignal (function to calculate signals), and defaultFitopts (the default options and parameters for data fitting). 
3.	mati.ImageData is the class for image data. It contains basic properties of the image data such as img (data matrix), sigma (optional noise standard deviation), mask (mask for data fitting). 
4.	mati.FitPars is the class for general data fitting. It contains model (signal model) and fitopts (fitting options and parameters). mati.FitPars.Fit() is the main function to fit data. 

There are subclasses are defined for specific applications. 
1.	mati.DiffusionPulseSequence is a subclass of mati.PulseSequence specifically for diffusion pulse sequences. It includes more diffusion MRI parameters such as b value, δ (each diffusion gradient duration), Δ (separation of two diffusion gradients), and shape (the diffusion gradient waveform shape such as trapezoid or trapezoidal-cosine). 
2.	There are three subclasses of mati.SignalModel, corresponding to three different methods that either published or under review: 
    1.	mati.IMPULSED for estimation of mean cell size using the IMPULSED method
    2.	mati.MRICytometry for estimation of cell size distribution using the MRI-cytometry method
    3.	mati.JOINT for simultaneous estimation of mean cell size and transcytolemmal water exchange using a joint method (IMPULSED + water exchange)

Examples about how to use the MATI package are provided in the mati/Examples folder. All necessary functions, data, and help files are included with the main script for each example. Users are encouraged to run examples first to get familiar with the MATI package. 

More examples have been added including Eg04 - Eg06. Each of these examples contains the shared code of a recently-published paper. Check them out! 

# Reference 
To quote the relevant references for the MATI code you can look in the individual example readme.txt files, and also please reference the review papers about OGSE and MR cell size imaging:
1.	Gore JC, Xu J, Colvin DC, Yankeelov TE, Parsons EC, Does MD. Characterization of tissue structure at varying length scales using temporal diffusion spectroscopy. NMR Biomed, 2010; 23: 745-56.
2.	Jiang X, Li H, Devan SP, Gore JC, Xu J. MR cell size imaging with temporal diffusion spectroscopy. Magn Reson Imaging, 2021; 77: 109-23.

# Comments or questions? 
Please send your comments or questions to Junzhong (JZ) Xu (junzhong.xu@vanderbilt.edu) or Xiaoyu Jiang (xiaoyu.jiang@vumc.org) 


