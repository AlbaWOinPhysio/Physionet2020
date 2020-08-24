function [score, label,classes] = run_12ECG_classifier(data,header_data, loaded_model)
	model=loaded_model.model;

	classes=loaded_model.classes;

    num_classes = length(classes);

    label = zeros([1,num_classes]);
    score = zeros([1,num_classes]);
    
    % Use your classifier here to obtain a label and score for each class.
    features = get_12ECG_features(data,header_data);
    
    [~,values] =predict(model.TreeBagger,features);	
    
    score = addResultToScore(score,classes,values,model.TreeBagger.ClassNames);
    
    result = predictForest (features, model.forest);
    score = addResultToScore(score,classes,result,model.forest.classNames);
    score = score./2;
    [~,I]  = max(score);
  
    score(I) = 1;
    label(I) = 1;
end

function newScore = addResultToScore(scores,classes,results,modelClassesNames)
    newScore =  scores;
    for i=1:length(modelClassesNames)
        idx=strcmp(classes,modelClassesNames{i});
        
        newScore(idx)=scores(idx)+results(i);
    end
end

function result = predictForest (features, Forest)
    Models=Forest.Models;
    selectedFeatures = Forest.selectedFeatures;
    result = zeros (length(Models),1);

    for m=1:length(Models)
        if (~isnumeric(Models{m}))
            [LABEL,POSTERIOR,~,~] = predict (Models{m},features(selectedFeatures(m,:)));
            if ((LABEL~= 0)&&(length(POSTERIOR)==2))
                result(m) = POSTERIOR(2);
            end
        end
    end
end

