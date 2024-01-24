clear all;
clc;
close all;

%% Inputs: This script reads in the Results.txt produced by jMRUI in order to grab the amplitude values for the reference peak being used for the S/N2 method and outputs them in a .txt or .mat file to be loaded in to McMRSGUI
% The 8 channel data was exported from the McMRSGUI into .txt jMRUI Textfiles. Then it was loaded in to jMRUI and processed. 
% filename = the name of AMARES results output from jMRUI in .txt format
% Num_channels = number of channels of data (for the original case this was 8 channels)
% Peak_selection = column corresponding to desired reference peak in jMRUI results .txt file. The second peak was the CH3 peak in the AMARES program.  
filename = 'Power_45_results.txt';
Num_channels = 8;
Peak_selection = 2;


% Calculations 
Y = readcell(filename,'Range','15:80000'); 
[A,B] = size(Y);

for i = 1:A
    switch Y{i,1}
%         case 'Frequencies (ppm)' % Frequencies position
        case 'Frequencies (ppm)' % Frequencies position
            Freq_pos = i;
        case 'Standard deviation of Frequencies (ppm)' % Standard deviation of frequencies
            Frequency_pos_std = i;
        case 'Amplitudes (-)' % Amplitudes of integration for peaks
            Amplitude_pos = i;
        case 'Standard deviation of Amplitudes (-)'
            Amplitude_pos_std = i;
        case 'Linewidths (Hz)'
            Linewidth_pos = i;
        case 'Standard deviation of Linewidths (Hz)'
            Linewidth_pos_std = i;
        case 'Phases (degrees)'
            Phases_pos = i;
        case 'Standard deviation of Phases (degrees)'
            Phases_pos_std = i;
        case 'Noise :'
            Noise_pos = i;
        case 'pH'%	Standard deviation of pH	[Mg2+]	Standard deviation of [Mg2+]'
            Ph_pos = i;            
        otherwise
    end
end

Num_acquisitions = (Frequency_pos_std-Freq_pos-1)/Num_channels;
Amplitudes = cell2mat(Y(Amplitude_pos+1:Amplitude_pos_std-1,1:2));


%% Create sections based on Monte Carlo simulation

vec = [1:Num_acquisitions:Num_acquisitions*Num_channels];
Weighting_amplitude = Amplitudes(:,Peak_selection);
[row2,col2] = find(isinf(Weighting_amplitude));

if isempty(row2) == 0
    for i = 1:length(row2)
        Weighting_amplitude(row2(i)) = NaN;
    end
else
end

for i = 1:Num_channels
    if i ~= Num_channels
        Weighting_amplitude_matrix(i,:) = Weighting_amplitude(vec(i):vec(i+1)-1);
    else
        Weighting_amplitude_matrix(i,:) = Weighting_amplitude(vec(i):end);
    end
end

filter = {'*.mat';'*.txt'};
[file, path] = uiputfile(filter);

complete_file = fullfile(path,file);
if file ~= 0
    [filepath,name,ext] = fileparts(file);

    if strcmp(ext,'.mat')
        save(complete_file,"Weighting_amplitude_matrix")

    elseif strcmp(ext,'.txt')
        T = table(Weighting_amplitude_matrix);
        writetable(T,complete_file)

    else % User didn't select a file type
    end
else
    errordlg('No weighting factors were exported for the S_N^2 method.')
    return
end
