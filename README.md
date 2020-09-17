#Objective
Here we describe our algorithm develop on The PhysioNet/Computing in Cardiology Challenge 2020 [1]. The goal of the 2020 Challenge is to identify the clinical diagnosis from 12-lead ECG recordings. The poster presents features that we used in classification. We also test features classification  power base on simple filter approach [2].


##Algorithm 
The algorithm has 3 steps: ECG pre-processing features extraction and classification. A sketchof the algorithm is shown in figure on the right. Our solution is based on bootstrap-aggregated (bagged) decision trees.
Signal preprocessing 
First, we load data from a file, we used a similar approach as an example. From the header file, we utilize age and sex data as classification features and gain and sample frequency for signal calibration. After calibration, we perform signal filtration. We used a median filter to remove some noise and the Butterworth high pass filter at cut-off frequency = 1 Hz to remove isoline’s floating. 

##Features 
In our algorithm we base on PhysioNet-Cardiovascular-Signal-Toolbox [3]. We used the following set of features: 
Global Electrical Heterogeneity  (GEH features) – such as in example code, this group contains 22 parameters based on spatial ventricular gradient vector (SVG) such as SVG magnitude, SVG elevation, SVG azimuth, etc. 
AF features – function AF_features.m in ECG_Analysis_Tools – this subset of features provide analysis of variability in RR interval and some sophisticated features such as coefficient if fuzzy measure entropy etc. [4]
Ratio of PVC beat – we use modified PVC_detect.m function to detect premature ventricular contraction beats[3]. As feature we use ratio of PVC beats to all recorded ECG beats. 
ECG periods –wavedet_3D_ECGKit function return time of characteristic points in ECG – based on it we write function that calculate some ECG periods, selected from literature [5, 6, 7]: PR, QS, QR, PT, TP. Mentioned periods are shown on image (on the left). We also calculate RAPR which is a ratio of PR and RR.
ECG morphology parameters – we also calculate QRS area [8] for each lead (as integrated ECG signal between Q and S points)  ST elevation as difference between value of ECG in J-point and isoline. J-point was defined as local extremum after QRSend points.R elevation as difference of value of ECG in R points and isoline. Methods shown on the left.

##Classification algorithm
We use build-in MATLAB, The Classification Learner app to choose the best classifier. Forests classification supported vector machine; k-nearest neighbor classifier was tested. As training data, we use all described features obtained for all samples in the challenge data set. Only valid recording, with one scored diagnose, is evaluated.Bootstrap-aggregated (bagged) decision trees shown the best accuracy and was chosen to perform the classification task.
Our algorithm was scored several times in the official phase, the best challenge score obtained be our team is 0,308

### Literature
1. Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220
2. John G. Kohavi R. (1997) "Wrappers for feature subset selection", Artificial Intelligence, Vol.97, No.1-2, pp.272-324.
3. Vest A, Da Poian G, Li Q, Liu C, Nemati S, Shah A, Clifford GD, "An Open Source Benchmarked Toolbox for Cardiovascular Waveform and Interval Analysis", Physiological measurement 39, no. 10 (2018): 105004. DOI:10.5281/zenodo.1243111; 2018. 
4. Q. Li, C. Y. Liu, J. Oster and G. D. Clifford. Chapter: Signal processing and feature selection preprocessing for classificationin noisy healthcare data. In book: Machine Learning for Healthcare Technologies, Edition: 1st, Publisher: IET, Editors: David A. Clifton, 2016.
5. Mao, L., Chen, H., Bai, J., Wei, J., Li, Q., & Zhang, R. (2019, October). Automated Detection of First-Degree Atrioventricular Block Using ECGs. In International Conference on Health Information Science (pp. 156-167). Springer, Cham.
6. Elgendi, M., Jonkman, M., & De Boer, F. (2008, August). Premature atrial complexes detection using the Fisher Linear Discriminant. In 2008 7th IEEE International Conference on Cognitive Informatics (pp. 83-88). IEEE.
7. Rafa³ Baranowski 25 EKG na XXV-lecie SENIT (ECG atlas in Polish, from Kasprowisko 2019 conference)
7. Krasteva, V. T., Jekova, I. I., & Christov, I. I. (2006). Automatic detection of premature atrial contractions in the electrocardiogram. Electrotechniques Electronics E & E, 9, 10

## Links
https://physionetchallenges.github.io/2020/ - PhysioNet/CinC challenges main pagecin
https://www.cinc2020.org/ - Computing in Cardiology 2020 - link to full paper

The code uses three main toolboxes:
- HRV toolbox to compute the RR intervals. https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox.git. 
  "An Open Source Benchmarked Toolbox for Cardiovascular Waveform and Interval Analysis", 
   Physiological measurement 39, no. 10 (2018): 105004. DOI:10.5281/zenodo.1243111; 2018. 
 - ECGkit to find the ECG fiducial points: https://github.com/marianux/ecg-kit.git
  Demski AJ, Llamedo Soria M. "ecg-kit a Matlab Toolbox for Cardiovascular Signal Processing".  
  Journal of Open Research Software. 2016;4(1):e8. DOI: http://doi.org/10.5334/jors.86
- GEH parameter extraction and origin point: https://github.com/Tereshchenkolab/Global-Electrical-Heterogeneity.git and https://github.com/Tereshchenkolab/Origin.git. 
  Perez-Alday,et al; "Importance of the Heart Vector Origin Point Definition for an ECG analysis: 
  The Atherosclerosis Risk in Communities (ARIC) study". Comp Biol Med, Volume 104, January 2019, 
  pages 127-138. https://doi.org/10.1016/j.compbiomed.2018.11.013
  Waks JW, et al. "Global Electric Heterogeneity Risk Score for Prediction of Sudden Cardiac Death in the General Population: 
  The Atherosclerosis Risk in Communities (ARIC) and Cardiovascular Health (CHS) Studies". Circulation. 2016;133:2222-2234.
