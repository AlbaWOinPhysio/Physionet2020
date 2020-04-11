function [score, label] = run_12ECG_classifier(data,header_data,classes, ~)

    
    num_classes = length(classes);

    label = zeros([1,num_classes]);
    
    
    % Use your classifier here to obtain a label and score for each class.
    [Result]  = get_12ECG_features(data,header_data);
    
    prepared_data = [Result.GEH, Result.AF_param([2,3,4,5,6,7,12,14]), Result.PVC ,Result.resultVector];
    score = myNeuralNetworkBayesian (prepared_data);	
    [~,idx] = max (score);

    label(idx)=1;
end



