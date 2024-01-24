clear all;
clc;

%% Creating noise covariance matrix
% A known noise covariance file (square matrix) for the WSVD can be input into the
% variable "Covariance" and saved into a format that the McMRSGUI accepts. 
Covariance = ones(8,8);

filter = {'*.mat';'*.txt'};
[file, path] = uiputfile(filter);
complete_file = fullfile(path,file);
[filepath,name,ext] = fileparts(file);

if strcmp(ext,'.mat')
    save(complete_file,"Covariance")

elseif strcmp(ext,'.txt')
    T = table(Covariance);
    writetable(T,complete_file)

else % User didn't select a file type
    errordlg('A file type was not loaded and so no noise covariance matrix was exported.')
    return
end