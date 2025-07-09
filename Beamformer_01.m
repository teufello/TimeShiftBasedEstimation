%% Init
field_init(-1);
clear all; close all;

%% Run parallel
if isempty(gcp('nocreate'))
    p = parpool('local',3); % cluster/local max. 20/4
end

%% Define the parameters of the transducer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_data = false; % choose true: simulated data; false: acquired data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sim_data
    fs = 100e6; % 100e6 for simulation
    f0 = 7e6; % Center freq in Hz
    c = 1540;
    multipl_factor = 1e23;
else
    fs = 20.833e6; 
    f0 = 5.2e6; % Center freq in Hz
    c = 1487.4544;
    multipl_factor = 1;
end
N_elements = 192;
Active_elements = 64;
lambda = c/f0; % wavelength
kerf = 0.03/1000;
width = 0.2/1000;
pitch = width + kerf; % 0.23 mm
F_transm = -1;
pixel_size = 0.1/1000; %  pixel size laterally/axially 

%%%%%%%%%%
%% Load %% 
%%%%%%%%%%
%% Truncate and filter acquired data
if ~sim_data
start_time = 2.689e-5; % in seconds
timeDelay = round(start_time * fs);

for seq = 1:2
    for i = 1:17
        % Construct the filename with zero-padding
        filename = sprintf('elem_data_em%04d_seq_%04d.mat', i, seq);
        
        % Load the data
        data = load(filename);
        
        % Extract the variable 'samples' from the loaded structure
        samples = double(data.samples);
        
        % Correct padding loop
        idx_start = timeDelay + 1;
        idx_end = timeDelay + size(samples,1);
        image_data(idx_start:idx_end, :, i+(seq-1)*17) = samples;
    end
end

% Filter
f_low = f0 - 1e6;          % Lower cutoff frequency in Hz
f_high = f0 + 1e6;         % Upper cutoff frequency in Hz
filter_order = 98;   % Filter order 
% Normalize cutoff frequencies
wn = [f_low, f_high] / (fs / 2);
% Design FIR bandpass filter
fir_coeff = fir1(filter_order, wn, 'bandpass');
image_data = filtfilt(fir_coeff, 1, image_data);
% freqz(fir_coeff, 1, 1024, fs); % validate filter 
end

%% Load phantom
if sim_data
% Plug flow
load img_LRI_sim_plug.mat % times x elements x transm*frames
load img_LRI_sim_plug_2scat.mat
% Parabolic flow
% load img_LRI_sim_parabolic.mat

image_data = image_LRI;
end

%% Parameters
% Define how many transmissions, frames
Nshoots = size(image_data,3);
Ntransm = 17;
Nframes = Nshoots/Ntransm;
% Compute the image width
image_width = (Ntransm-1)*1.84/1000; % vs spacing 1.84mm
% Compute the position of the left most line
x = -image_width/2;
% Define transmit focus
z_focus = (Active_elements*pitch)/F_transm;

% Define array
arrayPos = (-N_elements/2+0.5:N_elements/2-0.5)*pitch; % Make a pitched array
arrayPos = [arrayPos.', zeros(N_elements,1)]; % the z-coordinate of the array pos
% Define virtual sources
VS = [(-64:8:64)' * pitch, ones(Ntransm,1)*z_focus]; 
% Apodization
apod = hanning(N_elements);

% Create the pixel map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depth_mm = 60; % mm
lateral_mm = 30; % mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_depth = depth_mm/1000/pixel_size;
num_lateral = lateral_mm/1000/pixel_size;
depth_range = linspace(0.1, depth_mm, num_depth) / 1000; % m
lateral_range = linspace(-lateral_mm/2, lateral_mm/2, num_lateral) / 1000; % m 
[X,Z] = meshgrid(lateral_range, depth_range);
pixelMap = cat(3, X, Z);

%%%%%%%%%%%%%%%%%
%% Beamforming %%
%%%%%%%%%%%%%%%%%
% Compute LRIs
LRI = zeros(size(pixelMap,1),size(pixelMap,2),Nshoots);
rf_data_hilbert = hilbert(image_data*multipl_factor);
parfor i = 1:Nshoots
    line_idx = mod(i-1, Ntransm) + 1;
    LRI(:,:,i) = dyn_imageFormation(rf_data_hilbert(:,:,i),pixelMap, VS(line_idx,:), arrayPos, c, fs);
    fprintf("Loop no.: %d of %d\n", i, Nshoots);
end

% Create HRIs
HRI_all = zeros(size(LRI,1),size(LRI,2), Nframes); % Split up in the frames
HRI_display = zeros(size(LRI,1),size(LRI,2), Nframes); % Split up in the frames
for frame = 1:Nframes
    start_idx = frame*Ntransm-Ntransm+1;
    end_idx = frame*Ntransm;
    HRI_all(:,:,frame) = sum(LRI(:,:,start_idx:end_idx),3);
    HRI_display(:,:,frame) = abs(HRI_all(:,:,frame)) / max(abs(HRI_all(:,:,frame)), [], 'all');
end

%% Visu
set(groot, 'DefaultAxesFontSize', 14);       % Sets default font size for axes
set(groot, 'DefaultColorbarFontSize', 14);   % Sets default font size for colorbar
set(groot, 'DefaultTextFontSize', 14);       % Sets default font size for text elements


% Plot
for frame = 1:Nframes
    figure('Position',[100,100,600,800]);
    imagesc([-lateral_mm/1000/2*1000, lateral_mm/1000/2*1000], ...
        [0, depth_range*1000], 20*log10(HRI_display(:,:,frame)));
    xlabel('Lateral (mm)', 'FontSize',14);
    ylabel('Depth (mm)', 'FontSize',14);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = title(sprintf('Parabolic flow, Frame %d of %d', frame, Nframes), ...
        'Interpreter', 'tex');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t.Units = 'normalized';         % Use normalized coordinates (0 to 1)
    t.Position(2) = 1.05;           % Move title higher (default is around 1.01)
    axis ij equal tight;
    % ylim([20, 40]); % used for cyst phantom
    clim([-60, 0]);
    cb = colorbar; cb.Label.String = 'Received pressure (dB)'; cb.FontSize = 14; colormap(gray);
end

%% LRIs
LRI_display = zeros(size(LRI));
for k = 1:size(LRI, 3)
    frame = abs(LRI(:,:,k));
    LRI_display(:,:,k) = frame / max(frame(:));
end

for trans = 1:Ntransm
    figure('Position',[100,100,600,800]);
    imagesc([-lateral_mm/1000/2*1000, lateral_mm/1000/2*1000], ...
        [0, depth_range*1000], 20*log10(LRI_display(:,:,trans)));
    xlabel('Lateral (mm)', 'FontSize',14);
    ylabel('Depth (mm)', 'FontSize',14);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = title(sprintf('LRI, Transm %d', trans), 'Interpreter', 'tex');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t.Units = 'normalized';         % Use normalized coordinates (0 to 1)
    t.Position(2) = 1.05;           % Move title higher (default is around 1.01)
    axis ij equal tight;
    % ylim([20, 40]); % used for cyst phantom
    clim([-60, 0]);
    cb = colorbar; cb.Label.String = 'Received pressure (dB)'; cb.FontSize = 14; colormap(gray);
end

%% Save
% save the variable HRI_all to retain the complex data for velocity
% estimate, HRI_display is just subject to display the signal nice