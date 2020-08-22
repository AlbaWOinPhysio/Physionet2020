function [forrest] = BuildForestfromSelectedFeatures (features, labels, classes)
sl = size (labels);
nr = 1; 
for i=1:sl(2)
    classObservation = logical(labels(:,i));
    if (any(classObservation))
        fromClassFeatures= features(classObservation,:);
        notFromClassFeatures= features(~classObservation,:);
try
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
    classNames(i) =  string(classes{classNumber(i)});
    classObservation = logical(labels(:,classNumber(i)));

    if (any(selectedFeatures(i,:))&&any(classObservation))
        
        fromClassFeatures= features(classObservation,:);
        sfCF = size(fromClassFeatures);
        notFromClassFeatures= features(~classObservation,:);
        sNfCF = size(notFromClassFeatures);
        
        n = round(sNfCF(1)/sfCF(1));
        if (n>2)
            n=2;
        end
        c=randperm(sNfCF(1),sfCF(1)*n);

        X = [fromClassFeatures(:,selectedFeatures(i,:)); notFromClassFeatures(c,selectedFeatures(i,:))];
        y = [];
        y(1:sfCF(1),1) = string(classes(classNumber(i)));
        y((sfCF(1)+1):(sfCF(1)*(n+1)),1) = string(0);
     
        Models{i} = fitctree(X,y);
%         CVSVMModel = crossval(Models{i});        
%         classLoss = kfoldLoss(CVSVMModel);
%         disp (['class Loss: ',  num2str(classLoss)])
    end
end
forrest.Models = Models;
forrest.selectedFeatures = selectedFeatures;
forrest.classNames = classNames;

