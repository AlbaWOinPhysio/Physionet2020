function model = load_12ECG_model()

        filename='finalized_model.mat';
        A=load(filename);
        AF_model = load ('AF_model_balanced.mat');
        model.A= A.model;
        model.AF = AF_model.SVMModel;
        
end


