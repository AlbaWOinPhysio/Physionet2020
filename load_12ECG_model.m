function model = load_12ECG_model()

    GEH=load('GEH_dAC.mat');
    AF_model = load ('AF_model_balanced.mat');
    model.GEH= GEH.discriminantAnalysisClassifier;
    model.AF = AF_model.SVMModel;
        
end


