function [score, label] = run_12ECG_classifier(data,header_data,classes, model)

    num_classes = length(classes);

    label = zeros([1,num_classes]);
    %score = ones([1,num_classes]);
    
    % Use your classifier here to obtain a label and score for each class.
    [GEH, AF_param,PVC]  = get_12ECG_features(data,header_data);
    
    dAC_predicts = predict(model.GEH,GEH);
    list = strsplit (dAC_predicts{1}, ',');
    
    for i=1:length(list)
        label(strcmp (classes, list{i})) = 1;
    end
    
    AF = predict(model.AF, AF_param);
    
    if (AF == 1)
        label(strcmp (classes, 'AF')) = 1;
    end
    
    if (PVC > 0.1)
        label(strcmp (classes, 'PVC')) = 1;
    end
    		
    score = label;
end



