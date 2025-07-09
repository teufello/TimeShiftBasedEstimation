clear all; close all;

%% Load
% load img_HRI_beamformed_sim_parabolic.mat
% load img_HRI_beamformed_sim_plug.mat
load img_HRI_beamformed_acq.mat

% Normalize
Nframes = size(HRI_all,3);
HRI_display = zeros(size(HRI_all));
for frame = 1:Nframes
    HRI_display(:,:,frame) = (HRI_all(:,:,frame)) / max(abs(HRI_all(:,:,frame)),[],'all');
end

%% Parameters
kerf = 0.03/1000;
width = 0.2/1000;
pitch = width + kerf; % Spacing between virtual sources = 0.23/1000 (m)
% Pixel map
pixel_size = 0.1/1000; %  pixel size laterally/axially 
axial = size(HRI_all,1)*pixel_size;
lateral = size(HRI_all,2)*pixel_size;
num_depth = axial/pixel_size;
depth_range = linspace(pixel_size, axial, num_depth); % Convert

%% Visu
fontSize=16;
set(groot, 'DefaultAxesFontSize', fontSize);       % Sets default font size for axes
set(groot, 'DefaultColorbarFontSize', fontSize);   % Sets default font size for colorbar
set(groot, 'DefaultTextFontSize', fontSize);       % Sets default font size for text elements
figure('Position',[100,100,600,800]);

% Plot
for frame = 1:Nframes
    figure('Position',[100*frame,100,600,800]);
    imagesc([-lateral/2*1000, lateral/2*1000], [0, depth_range*1000], 20*log10(abs(HRI_display(:,:,frame))));
    xlabel('Lateral (mm)');
    ylabel('Depth (mm)');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % t = title(sprintf('Parabolic flow, Frame %d of %d', frame, Nframes), 'Interpreter', 'tex');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t.Units = 'normalized';         % Use normalized coordinates (0 to 1)
    t.Position(2) = 1.05;           % Move title higher (default is around 1.01)
    axis ij equal tight;
    % ylim([20, 40]); % only vessel
    % ylim([25, 45]); % acq data% 
    clim([-60, 0]);
    cb = colorbar; cb.Label.String = 'Envelope amplitude (dB)'; cb.FontSize = fontSize+1; colormap(gray);
end

% In subplots
% HRI_display(:,:,1) = HRI_all(:,:,1) / max(abs(HRI_all(:,:,1)), [], 'all');
% HRI_display(:,:,2) = HRI_all(:,:,2) / max(abs(HRI_all(:,:,2)), [], 'all');
% figure('Position',[100,100,600,800]);
% for frame = 1:Nframes
%     subplot(Nframes,1,frame);
%     imagesc([-lateral/2*1000, lateral/2*1000], [0, depth_range*1000], 20*log10(abs(HRI_display(:,:,frame))));
%     xlabel('Lateral (mm)', 'FontSize',14);
%     ylabel('Depth (mm)', 'FontSize',14);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     t = title(sprintf('Parabolic flow, Frame %d of %d', frame, Nframes), 'Interpreter', 'tex');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     t.Units = 'normalized';         % Use normalized coordinates (0 to 1)
%     t.Position(2) = 1.05;           % Move title higher (default is around 1.01)
%     axis ij equal tight;
%     ylim([25, 35]); % only vessel
%     clim([-60, 0]);
%     cb = colorbar; cb.Label.String = 'Received pressure (dB)'; cb.FontSize = 14; colormap(gray);
% end