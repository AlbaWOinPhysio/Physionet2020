function [selectedFeatures] = SelectFeatures (features, labels,classes)
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
        selectedFeaturesp(nr,:) = p;
        classNumber(nr) = i;
        nr = nr + 1;
catch
end
%         figure()
%         ecdf(p);
%         xlabel('P value');
%         ylabel('CDF value')
    end
     %build table of p-values (presented in paper) 
    build_grouped_table (selectedFeaturesp,classes,classNumber)
end

function build_grouped_table (sf_values,classes,classNumber)

    n= 1;
    result = {"Name", "age", "sex", "GEH", "AF feature", "% PVC beat", "PR", "QS", "QR", "PT", "TP", "RAPR", "ST elevation", "QRS area", "R elevation"};
    s = size(sf_values);
    for n=1:s(1)
        result{n+1,1} = classes(classNumber(n));
        result{n+1,2} = sf_values (n,1);%age
        result{n+1,3} = sf_values (n,2);%sex
        result{n+1,4} = (min(sf_values(n,3:24)));%GEH
        result{n+1,5} = min(sf_values(n,25:32));%AF
        result{n+1,6} = min(sf_values(n,33)) ;%PVC
        result{n+1,7} = min(sf_values(n,33:39)) ;%PR
        result{n+1,8} = min(sf_values(n,40:46)) ;%QS
        result{n+1,9} = min(sf_values(n,47:53)) ;%QR
        result{n+1,10} = min(sf_values(n,54:60)) ;%PT
        result{n+1,11} = min(sf_values(n,61:67)) ;%TP
        result{n+1,12} = min(sf_values(n,68:74)) ;%RAPR
        result{n+1,13} = min(sf_values(n,75:86)) ;%ST
        result{n+1,14} = min(sf_values(n,87:170));%QRS
        result{n+1,15} = min(sf_values(n,171:183));%R
    end
    xlswrite('result.xls',string(result))
