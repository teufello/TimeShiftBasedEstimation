% Initialize the field system
field_init(-1);
clear all;
close all;

% Generate the transducer apertures for send and receive
f0 = 7e6;               % Transducer center frequency [Hz]
M = 2;                  % Number of cycles in emitted pulse
fs = 100e6;             % Sampling frequency [Hz]
c = 1540;               % Speed of sound [m/s]
lambda = c/f0;          % Wavelength [m]
pitch = 0.23/1000;       % Pitch of transducer
element_height = 5/1000;% Height of element [m]
kerf = 0.03/1000;        % Kerf [m]
width = pitch - kerf;   % Width of element
focus = [0, 0, 30]/1000;% Fixed focal point [m]
N_elements = 192;         % N_elements = 192;
Active_elements = 64;

% Set the sampling frequency
set_sampling(fs);

% Generate aperture for emission
emit_aperture = xdc_linear_array(N_elements, width, element_height, kerf, 1, 1, focus);

% Set the impulse response and excitation of the emit aperture
impulse_response = sin(2*pi*f0*(0:1/fs:2/f0));
impulseresponse = impulse_response.*hann(numel(impulse_response))';
xdc_impulse(emit_aperture, impulse_response);

excitation = sin(2*pi*f0*(0:1/fs:M/f0));
xdc_excitation(emit_aperture, excitation);

% Generate the receive aperture
receive_aperture = xdc_linear_array(N_elements, width, element_height, kerf, 1, 1, focus);
% Set the impulse response for the receive aperture
xdc_impulse(receive_aperture, impulse_response);

% Set a hanning apodization on the apertures
apo = hann(Active_elements)';

%% Define the transmit sequence
vs_x = (-64:8:64)*pitch;
vs_z = -Active_elements*pitch;  % F = -1;
no_lines = numel(vs_x);

%% Determine how many frames to simulate
Nframes = 2;
Nshoots = no_lines*Nframes;         % Number of shots to be processed

%% Create the digital flow phantom
x_range = 15/1000;
z_range = 8/1000;
peak_velocity = 0.3; % [m/s]
PRF = 5000;          % [Hz]
Tprf = 1/PRF;       % [sec]
num_scatterers = 15000;
z_offset = 30/1000;
theta = 0/180*pi; % in deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,y,z,amp,velocity] = make_flow_phantom(x_range, z_range, peak_velocity, ...
    num_scatterers, 'plug'); % parabolic, plug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make the flow simulation
for i = 1:Nshoots

    currentVSPos = mod(i-1, no_lines) + 1;
    disp(vs_x(currentVSPos));
    % Set the transmit center by setting the focus center
    xdc_center_focus(emit_aperture, [vs_x(currentVSPos) 0 0]);
    % Set the position of the transmit focus
    xdc_focus(emit_aperture, 0, [vs_x(currentVSPos), 0, vs_z]);
    % Set the receive center by setting the focus center
    xdc_center_focus(receive_aperture, [0 0 0]);
    % Set the position of the receive focus and the time in which the Rx focus
    % is valid
    xdc_focus(receive_aperture, 0, [0, 0, 5]); % set the focus to 5 meters to ensure no receive delay is applied
    
    % Set the correct apodization (active the correct elements)
    actual_apo_vector = zeros(N_elements,1);
    actual_apo_vector( (1:Active_elements) + (currentVSPos-1)*8 ) = apo;
    % figure(1); plot(actual_apo_vector);
    % Apply the apodization and specify the time in which the apodization
    % is valid
    xdc_apodization(emit_aperture, 0, actual_apo_vector.');    % Only set the tx apodization. no rx apodization

    % Generate the rotated and offset block of sample
    xnew = x*cos(theta) + z*sin(theta);
    znew = z*cos(theta) - x*sin(theta) + z_offset;
    scatterers = [xnew; y; znew]';

    % visualize the scatterer movement
    % scatter(xnew, znew),drawnow;pause(0.02);  % You can uncomment this line to vsiualize the flow tube

    % Calculate the receive response
    [v, t1] = calc_scat_multi(emit_aperture, receive_aperture, scatterers, amp');

    % Store result
    image_data_raw(1:max(size(v)),:, i) = v;
    times(i) = t1;

    % Propagate the scatterers and aliaze them to lie within the correct
    % range
    x1 = x;
    x = x + velocity*Tprf;
    outside_range = (x>x_range/2);
    x = x-x_range*outside_range;
end

%% Zero padding
% Compute how much zeros to pad
timeDelay = times*fs;
						  
% Get raw data size
[samples, channels, Nshoots] = size(image_data_raw);
maxPad = max(timeDelay);

% Initialize padded image
image_data = zeros(samples + maxPad, channels, Nshoots);

% Correct padding loop
for i = 1:Nshoots
    idx_start = round(timeDelay(i)) + 1;
    idx_end = timeDelay(i) + samples;
    image_data(idx_start:idx_end, :, i) = image_data_raw(:, :, i);
end

%%
image_HRI = sum(image_data(:,:,1:no_lines), 3);
image_HRI = image_HRI / max(image_HRI, [], 'all');
image_LRI = image_data / max(image_data, [], 'all');

%% Visu
% HRI
z = [0, size(image_data,1)/fs * c/2];
image_width = (N_elements-1)*pitch;
figure('Position',[100,100,600,600]);
imagesc([-image_width/2*1000, image_width/2*1000], [0, z*1000], 20*log10(abs(hilbert(image_HRI))));
xlabel('Lateral (mm)');
ylabel('Depth (mm)');
title('High Resolution Image');
axis ij equal tight;
clim([-60, 0]);
cb = colorbar; cb.Label.String = 'Received pressure (dB)'; colormap(gray)

% Scatterers
figure;
scatter(scatterers(:,1)*1000, scatterers(:,3)*1000, 1, 'filled');
xlabel('Lateral (mm)');
ylabel('Depth (mm)');
title('Scatterer Distribution');
axis equal;
axis ij;

% Raw image data
hri_raw = sum(image_data_raw(:,:,1:no_lines), 3);
hri_raw = hri_raw / max(hri_raw, [], 'all');
figure;
imagesc([-image_width/2*1000, image_width/2*1000], [0, z*1000], 20*log10(abs(hilbert(hri_raw))));
xlabel('Lateral (mm)');
ylabel('Depth (mm)');
title('Raw Image Data');
axis ij equal tight;
clim([-60, 0]);
cb = colorbar; cb.Label.String = 'Received pressure (dB)'; colormap(gray);

%% Save the created LRI
% save("img_LRI_sim_parabolic_andrea.mat", "image_LRI")