function [score, label] = run_12ECG_classifier(data,header_data,classes, model)

    num_classes = length(classes);

    label = zeros([1,num_classes]);
    %score = ones([1,num_classes]);
    
    % Use your classifier here to obtain a label and score for each class.
    [Result]  = get_12ECG_features(data,header_data);
    
    dAC_predicts = predict(model.GEH,Result.GEH);
    list = strsplit (dAC_predicts{1}, ',');
        
    for i=1:length(list)
        label(strcmp (classes, list{i})) = 1;
    end
    
    AF = predict(model.AF, Result.AF_param);
    
    if (AF == 1)
        label(strcmp (classes, 'AF')) = 1;
    end
    
    if (Result.PVC > 0.1)
        label(strcmp (classes, 'PVC')) = 1;
    end
    
    if (~isnan(Result.ST_elevation))
        sex = Result.GEH(2);
        age = Result.GEH(1);
        if (length(Result.ST_elevation)== 12)
            t1 = (sum(Result.ST_elevation(1:3) <= (-50000)));
            t2 = (sum(Result.ST_elevation(4:6) <= (-50000)));
            t3 = (sum(Result.ST_elevation(7:9) <= (-50000)));
            t4 = (sum(Result.ST_elevation(10:12) <= (-50000)));

            if (sex==1)
                V2V3_treshold  = 150000;
            else
                if (age >= 40)
                    V2V3_treshold  = 200000;
                else
                    V2V3_treshold  = 250000;
                end
            end

            if any(Result.ST_elevation(8:9)>V2V3_treshold)||any(Result.ST_elevation(1:7)>100000)||any(Result.ST_elevation(10:12)>100000)
                label(strcmp (classes, 'STE')) = 1;
            else
                if (t1>=2||t2>=2||t3>=2||t4>=2)
                    label(strcmp (classes, 'STD')) = 1;
                else
                    label(strcmp (classes, 'STD')) = 0;
                    label(strcmp (classes, 'STE')) = 0;
                end
            end
        end

        if (isempty (rmmissing(Result.ECG_periods.PR)))
            Result.ECG_periods.PR = 0;
        end
        if (isempty (rmmissing(Result.ECG_periods.QR)))
            Result.ECG_periods.QR = inf;
        end
        if (isempty (rmmissing(Result.ECG_periods.QS)))
            Result.ECG_periods.QS = 0;
        end
        
        if (median(rmmissing(Result.ECG_periods.PR))>0.2)
            label(strcmp (classes, 'I-AVB')) = 1;
        else
            label(strcmp (classes, 'I-AVB')) = 0;
        end
        
        if (label(strcmp (classes, 'RBBB')))
            if ((max(rmmissing(Result.ECG_periods.QS))<0.120) && (max(rmmissing(Result.ECG_periods.QR))<0.050))
                label(strcmp (classes, 'RBBB')) = 0;
            end
        end
        if (label(strcmp (classes, 'LBBB')))
            if ((max(rmmissing(Result.ECG_periods.QS))<0.120) && (max(rmmissing(Result.ECG_periods.QR))<0.060))
                label(strcmp (classes, 'LBBB')) = 0;
            end
        end
    end
    
    
    if (sum(label) == 0)
        label(strcmp (classes, 'Normal')) = 1;
    end
    
    		
    score = label;
end



