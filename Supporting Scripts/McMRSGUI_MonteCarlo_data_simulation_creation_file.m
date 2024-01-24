clear all;
clc;
close all;

%% Creating the model signal
% Inputs
a = [15 10];  % amplitude of the respective peaks (3:2)
theta = [0 0]; % angle in rad/s for respective peaks
d = [10 10]; % dampings for respective peaks
f = [132 56]; % frequency in Hz for respective peaks
fs = 2500;   % sampling frequency
num_pts = 2048; % number of data points 
F0 = 108.57; % Set center frequency of signal (MHz)

% Calculations
dwell = 1/fs; % dwell time
BW = 1/dwell; % Bandwidth of signal (Hz)

tacq = num_pts/fs;
t_axis = [0:1/fs:tacq-1/fs]; % Set the time axis based on inputs

% Create FID data
Fid_signal = zeros(1,num_pts);
for i = 1:length(a)
        Fid_signal(i,:) = a(i)*exp(1j*theta(i))*exp((-d(i)+1j*2*pi*f(i))*t_axis);
end
time_domain = sum(Fid_signal,1);

%% Creating multiple channel data
% Inputs
Num_coils = 8;

% Set scaling values for Num_coil channels
scale_factor = [1, 0.8, 0.7, 0.6, 0.5, 0.65, 0.9, 0.93]';

% Initial Scaling factors in degrees
phase_shift_deg = [0,0,0,0,0,0,0,0]';

% Calculations
phase_shift_rad = phase_shift_deg*pi/180;

% Generate Num_coils channel data
Ch_TD = zeros(Num_coils, num_pts);
for i = 1:Num_coils
    Ch_TD(i,:) = scale_factor(i)*time_domain*exp(1j*phase_shift_rad(i));
end

%% Noise calculations between channels
% Inputs
% A = Num_coil x Num_coil matrix, off-diagonal elements associated with noise correlation

A = zeros(Num_coils,Num_coils);
% Guess and Test Found values That worked for normalized correlation values
Vec1 = [0, -0.031578947368421, 0.194736842105263, -0.036842105263158, -0.036842105263158, 0.163157894736842, 0.110526315789474, 0.163157894736842];
Vec2 = [0, -0.021052631578947, 0.131052631578947, 0.205263157894737, 0.020526315789474, 0.128947368421053, -0.042105263157895];
Vec3 = [0,0.131578947368421, 0.015789473684211, 0.106842105263158, 0.011578947368421, -0.039473684210526];
Vec4 = [0, 0.010526315789474, 0.116842105263158, -0.042105263157895, -0.039473684210526];
Vec5 = [0, 0.019473684210526, 0.163157894736842, 0.007894736842105];
Vec6 = [0, 0.016315789473684, 0.136842105263158];
Vec7 = [0, -0.030526315789474];
Vec8 = [0];
Vec_name = {Vec1;Vec2;Vec3;Vec4;Vec5;Vec6;Vec7;Vec8};
for i = 1:Num_coils
    A(i,i:end) = Vec_name{i};
    A(i:end,i) = Vec_name{i};
end

% B = num_coil x num_coil matrix, diagonal elements denote noise variance on each coil element

% Guess and Test Found values That worked for normalized Covariance values
vec = [1, 0.818947368421053, 0.763157894736842, 0.814210526315790, 0.836842105263158, 0.797894736842105, 0.816315789473684, 0.751578947368421];
B = diag(vec);

% Num_coil x Num_coil matrix with the diagonal components representing the
% covariance of each channel and the off diagonal components representing
% the noise correlation between channels. Different full noise correlations
% can be input for analysis.
Full_noise_correlation_matrix = A+B;

%% Generating Gaussian independent noise
% Inputs:
Noise_level_dBm = [15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45];  % Define the noise level in dBM to be created
Num_Monte_Carlo_runs = 300; % Define the number of Monte Carlo simulations to run

% Calculations
for k = 1:length(Noise_level_dBm)
    for l = 1:Num_Monte_Carlo_runs
        % Generate the complex noise matrix with white gaussian noise for a 50 Ohm impedance at the specified Noise_level_dbM
        N_bar = wgn(num_pts, Num_coils,Noise_level_dBm(k),50,'dBm','complex').';

        % Noise Correlation Calculation
        R = corrcoef((N_bar.'));

        % Combine Total Noise Matrix with white gaussian noise
        N = (Full_noise_correlation_matrix)*(N_bar);
        R1 = corrcoef((N.'));

        % Adding the noise to the fids
        Simul_data_w_noise_TD(k,:,l,:) = (Ch_TD + N);
    end
end
 
%% Exporting the data to structure "Starting" for use with McMRSGUI
[m,n,o,p] = size(Simul_data_w_noise_TD);

for i = 1:m
    Filename = sprintf('Generated_FIDS_Power_%0.1f',Noise_level_dBm(i));
    Starting = struct("Filename",Filename,"Data",squeeze(Simul_data_w_noise_TD(i,:,:,:)),"F0",F0*10^6,"BW",BW);
    save(strcat(Filename,'.mat'), 'Starting')
end
