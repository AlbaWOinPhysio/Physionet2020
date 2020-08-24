function [forrest] = BuildForestfromSelectedFeatures (features, labels, classes)
sl = size (labels);
nr = 1; 
for i=1:sl(2)
    classObservation = logical(labels(:,i));
    if (any(classObservation))
        fromClassFeatures= features(classObservation,:);
        notFromClassFeatures= features(~classObservation,:);
try%For each class chose significant features
        [~,p,~,~] = ttest2(rmmissing(fromClassFeatures),rmmissing(notFromClassFeatures),'Vartype','unequal');
        selectedFeatures(nr,:) = (p<0.05);

        classNumber(nr) = i;
        nr = nr + 1;
catch
end
%         figure()
%         ecdf(p);
%         xlabel('P value');
%         ylabel('CDF value')
    end
end

sf = size (selectedFeatures);

disp ('Trainning models');
classNames(sf(1),1) = "";
Models= cell(sf(1),1);
for i=1:sf(1)
     disp(['    ', num2str(i), '/', num2str(sf(1)), '...']);
    try 
    classNames(i) =  string(classes{classNumber(i)});
    classObservation = logical(labels(:,classNumber(i)));

    if (any(selectedFeatures(i,:))&&any(classObservation))
        
        fromClassFeatures= features(classObservation,:);
        sfCF = size(fromClassFeatures);

       notFromClassFeatures= features(~classObservation,:);
       sNfCF = size(notFromClassFeatures);
       if ((sfCF(1)>0)&&(sNfCF(1)>0))
            %select observation for trainning, balance the proportion
            n = round(sNfCF(1)/sfCF(1));
            if (n>2)
                n=2;
            end
            K=sfCF(1)*n;
            if (K>sNfCF(1))
                K= sNfCF(1);
            end
            c=randperm(sNfCF(1),K);
            %search for best clasifier
            X = [fromClassFeatures(:,selectedFeatures(i,:)); notFromClassFeatures(c,selectedFeatures(i,:))];
            y = [];
            y(1:sfCF(1),1) = string(classes(classNumber(i)));
            y((sfCF(1)+1):K+sfCF(1),1) = string(0);
            %divide data in trainning and test dataset
            sX = size(X);
            trainningPart = round (sX(1)*0.8);
            trainX = X(1:trainningPart,:);
            testX = X(trainningPart:sX(1),:);
            
            trainY = y(1:trainningPart);
            testY = y(trainningPart:sX(1),:);
            
            %train and validiate classifiers
            testModels{1} = fitctree(trainX,trainY);
            yM1 = predict (testModels{1},testX);
            score(1) = sum(yM1==testY);
            
            testModels{2} = fitcnb(trainX,trainY);
            yM2 = predict (testModels{2},testX);
            score(2) = sum(yM2==testY);
            
            testModels{3} = fitcsvm(trainX,trainY);
            yM3 = predict (testModels{3},testX);
            score(3) = sum(yM3==testY);
            
            [~,Idx] = max(score)
            %get the best 
            Models{i} = testModels{Idx};
        end
%         CVSVMModel = crossval(Models{i});        
%         classLoss = kfoldLoss(CVSVMModel);
%         disp (['class Loss: ',  num2str(classLoss)])
    end
    catch
        
    end
end
%Return result
forrest.Models = Models;
forrest.selectedFeatures = selectedFeatures;
forrest.ClassNames = classNames;

