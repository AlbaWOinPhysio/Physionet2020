function [score, label] = run_12ECG_classifier(data,header_data,classes, model)

    num_classes = length(classes);

    label = zeros([1,num_classes]);
    %score = ones([1,num_classes]);
    
    % Use your classifier here to obtain a label and score for each class.
    [features, AF_param]  = get_12ECG_features(data,header_data);
    
    score = mnrval(model.A,features);
    
    AF = predict(model.AF, AF_param);
    
    if (AF == 1)
        label(strcmp (classes, 'AF')) = 1;
    end
    
    		
    [~,idx] = max (score);

    label(idx)=1;
    
end



