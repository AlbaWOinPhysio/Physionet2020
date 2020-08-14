function [score, label,classes] = run_12ECG_classifier(data,header_data, loaded_model)


	model=loaded_model.model;
	classes=loaded_model.classes;

    num_classes = length(classes);

    label = zeros([1,num_classes]);
    score = zeros([1,num_classes]);
    
    % Use your classifier here to obtain a label and score for each class.
    features = get_12ECG_features(data,header_data);

    
    result =predict(model,features);	
    idx=strcmp(classes,result);
    
    
    score(idx) = 1;
    label(idx) = 1;
end



