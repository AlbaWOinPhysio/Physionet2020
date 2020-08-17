function  model = train_12ECG_classifier(input_directory,output_directory)

disp('Parsing data...')

% Find files.
input_files = {};
for f = dir(input_directory)'
    if exist(fullfile(input_directory, f.name), 'file') == 2 && f.name(1) ~= '.' && all(f.name(end - 2 : end) == 'mat')
        input_files{end + 1} = f.name;
    end
end

% read number of unique classes
classes = get_classes(input_directory,input_files);

supportedClasses = {'270492004';'164889003';'164890007';'426627000';'713427006';'713426002';'445118002';'39732003';'164909002';'251146004';'698252002';'10370003';'284470004';'427172004';'164947007';'111975006';'164917005';'47665007';'59118001';'427393009';'426177001';'426783006';'427084000';'63593006';'164934002';'59931005';'17338001'};


num_files = length(input_files);
%index - nr - number of valid sample (from supported classes)
nr = 1;
labels = zeros (nr, length(classes));

% Iterate over files.
for i = 1:num_files
    disp(['    ', num2str(i), '/', num2str(num_files), '...'])
    
    % Load data.
    file_tmp=strsplit(input_files{i},'.');
    tmp_input_file = fullfile(input_directory, file_tmp{1});
    
    [data,hea_data] = load_challenge_data(tmp_input_file);
    
%     Total_data{i}=data;
%     Total_header{i}=hea_data;
    validObservation = false;
    for j = 1 : length(hea_data)
        if startsWith(hea_data{j},'#Dx')
            tmp = strsplit(hea_data{j},': ');
            tmp_c = strsplit(tmp{2},',');
            if (length(tmp_c)==1)
                if any(strcmp(supportedClasses,tmp_c{1}))
                    idx=strcmp(classes,tmp_c{1});
                    if sum(idx)==1
                        labels(nr,idx)=1;
                        validObservation = true;
                    end
                end
            end
            break
        end
    end
    if (validObservation)
        tmp_features = get_12ECG_features(data,hea_data);
        features(nr,:)=tmp_features;
        nr = nr + 1;
    else
        disp("Sample are not valid, not supported or multiple class")
    end
end
%remove classes with no samples
%classes = supportedClasses (any(labels,1));
%labels = labels(:,any(labels,1));

save ('data.mat', 'features', 'labels','classes');

disp('Training model..')
names = cell (length (labels),1); 
for j = 1:length (labels)
    names(j,1) = classes (logical(labels (j,:)));
end

model = TreeBagger(100,features,names);
save_12_ECG_model(model,output_directory,classes);
end

function save_12_ECG_model(model,output_directory,classes)
% Save results.
tmp_file = 'finalized_model.mat';
filename=fullfile(output_directory,tmp_file);
save(filename,'model','classes','-v7.3');


disp('Done.')
end


% find unique number of classes
function classes = get_classes(input_directory,files)

classes={};
num_files = length(files);
k=1;
for i = 1:num_files
    g = strrep(files{i},'.mat','.hea');
    input_file = fullfile(input_directory, g);
    fid=fopen(input_file);
    tline = fgetl(fid);
    tlines = cell(0,1);
    
    while ischar(tline)
        tlines{end+1,1} = tline;
        tline = fgetl(fid);
        if startsWith(tline,'#Dx')
            tmp = strsplit(tline,': ');
            tmp_c = strsplit(tmp{2},',');
            for j=1:length(tmp_c)
                idx2 = find(strcmp(classes,tmp_c{j}));
                if isempty(idx2)
                    classes{k}=tmp_c{j};
                    k=k+1;
                end
            end
            break
        end
    end
    
    fclose(fid);
    
end
classes=sort(classes);
end

function [data,tlines] = load_challenge_data(filename)

% Opening header file
fid=fopen([filename '.hea']);

if (fid<=0)
    disp(['error in opening file ' filename]);
end

tline = fgetl(fid);
tlines = cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline = fgetl(fid);
end
fclose(fid);

f=load([filename '.mat']);

try
    data = f.val;
catch ex
    rethrow(ex);
end

end
