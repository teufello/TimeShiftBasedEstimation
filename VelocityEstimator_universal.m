%% Init
close all; clear all;

%% Run in parallel
% if isempty(gcp('nocreate'))
%     p = parpool('local',20); % cluster/local max. 20/4
% end

%%%%%%%%%%
%% Load %% 
%%%%%%%%%%
% Load phantom
% Plug flow
% load img_HRI_beamformed_sim_plug.mat % times x elements x frames
% Parabolic flow
% load img_HRI_beamformed_sim_parabolic.mat
% load img_HRI_beamformed_sim_parabolic33.mat


% Acquired data 
% load img_HRI_beamformed_acq.mat
load img_HRI_beamformed_acq34.mat

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%% !!! SELECT !!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_data = false; % choose true: simulated data; false: acquired data
parabolic = true; % choose true: parabolic flow, false: plug flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleSize = 81;
searchSize = 121;
thres_R = 1e-3; % default: 1e-3
PRF = 5000;         % [Hz]
Tprf = 1/PRF;       % [sec]
no_lines = 17;
% Get size of beamformed image(s)
num_lateral_pixels= size(HRI_all,2);
% Define lateral pixel size [m]
pixel_size = 0.1e-3;
% Compute total image width
image_width = num_lateral_pixels * pixel_size;
% Create centered lateral axis [m]
x_axis = linspace(-image_width/2, image_width/2, num_lateral_pixels);

% Find index of x = 0
[~, x0_idx] = min(abs(x_axis));
centerSample = x0_idx;

%%%%%%%%%%%%%%%%%%%%%%%
%% Cross-Correlation %%
%%%%%%%%%%%%%%%%%%%%%%%
Nframes = size(HRI_all,3);
if sim_data
    depth_range = 200:400;  % define range to receive the velocity along depth
else
    depth_range = 250:450;  % define range to receive the velocity along depth
end
padded = round((sampleSize + searchSize)/2); % (searchSize - sampleSize)/2;
sample = -(sampleSize-1)/2:(sampleSize-1)/2;
search = -(searchSize-1)/2:(searchSize-1)/2;
sample_apo = tukeywin(length(sample), 0.25)'; % horizontal vec
search_apo = tukeywin(length(search), 0.25)'; % horizontal vec
% sample_apo = rectwin(length(sample))'; % horizontal vec
% search_apo = rectwin(length(search))'; % horizontal vec

profiles = zeros(size(depth_range,2),num_lateral_pixels,2);
dataVecs1 = zeros(size(depth_range,2),sampleSize);
dataVecs2 = zeros(size(depth_range,2),searchSize);
Rs = zeros(size(depth_range,2),2*searchSize-1);
R_max = zeros(size(depth_range,2),1);

% Perform global normalization
max_val = max(abs(HRI_all), [], 'all');
HRI_norm = HRI_all / max_val;

% Definition positions to analyze in x-direction 
if parabolic
    x_positions_mm = -5:0.5:5;
else
    x_positions_mm = 0;
end

x_positions_m = x_positions_mm / 1000;

% Find closest pixel indices for each desired x position
x_indices = arrayfun(@(x) find(abs(x_axis - x)...
    == min(abs(x_axis - x)), 1), x_positions_m);

%% Echo-canceling section
if ~sim_data
    % THIS APPROACH IS BASED ON FRAME DIFFERENCE
    % Perform echo cancellation by subtracting the previous frame
    HRI_echo_canceled = zeros(size(HRI_norm));

    % First frame has no previous frame to subtract from - optionally leave it unchanged or zero
    HRI_echo_canceled(:,:,1) = HRI_norm(:,:,1); 

    % Loop over time frames and subtract previous frame
    for k = 2:size(HRI_norm, 3)
        % HRI_echo_canceled(:,:,k) = 0.5* (HRI_norm(:,:,k-1) - HRI_norm(:,:,k)); % lecture
        HRI_echo_canceled(:,:,k) = HRI_norm(:,:,k) - HRI_norm(:,:,k-1);
    end

    HRI_norm = HRI_echo_canceled;
    
    
    % THIS APPROACH IS BASED ON MEAN SUBTRACTION
    % Preallocate output array
    % HRI_echo_mean = zeros(size(HRI_norm));
    % % Compute mean across all frames
    % HRI_mean = mean(HRI_norm, 3);
    % 
    % % Subtract the same mean from each frame
    % for k = 1:size(HRI_norm, 3)
    %     HRI_echo_mean(:,:,k) = HRI_norm(:,:,k) - HRI_mean;
    % end
    % HRI_norm = HRI_echo_mean;
    
end
%% Velocity estimation section
velocities_matrix =...
    zeros(length(depth_range), length(x_indices), size(HRI_norm,3)-1);
velocities_matrix_filtered =...
    zeros(length(depth_range), length(x_indices), size(HRI_norm,3)-1);
max_interp = 0;
for frame = 1:Nframes-1 % iterate over frames
    for xi = 1:length(x_indices)
        centerSample = x_indices(xi);
    
        for i = 1:length(depth_range)
            % Extract lateral profiles
            profile1 = HRI_norm(depth_range(i),:,frame);
            profile2 = HRI_norm(depth_range(i),:,frame+1);
            profiles(i,:,1) = profile1;
            profiles(i,:,2) = profile2;
        
            % Check boundaries
            if (centerSample + max(sample) > size(HRI_norm,2)) || ...
               (centerSample + min(sample) < 1) || ...
               (centerSample + max(search) > size(HRI_norm,2)) || ...
               (centerSample + min(search) < 1)
                warning('Check your boundaries in loop i=%d',i);
                continue;
            end
        
            dataVec1 = profile1(centerSample + (sample)) .* sample_apo; % horizontal vec
            dataVec1 = dataVec1 - mean(dataVec1);
            dataVec2 = profile2(centerSample + (search)) .* search_apo; % horizontal vec
            dataVec2 = dataVec2 - mean(dataVec2);
            dataVecs1(i,:) = dataVec1;
            dataVecs2(i,:) = dataVec2;
            % dimension of R: 2*max(length(dataVec1),length(dataVec2))-1
            % lag zero at max(length(dataVec1),length(dataVec2))
            [R,lag_xcorr] = xcorr(dataVec1, dataVec2, 'none'); 
            % norm_factor = norm(dataVec1) * norm(dataVec2); 
            % R = R / norm_factor;  % Normalized Cross Correlation
            Rs(i,:) = R;
            R_padded = R;
            % Peak detection only in the direction of flow
            if sim_data
                R_padded(padded+1:end) = 0; % positive flow assumed
                % R_padded(1:padded-11) = 0; % upper velocity constraint (0.2941 m/s)
            else
                R_padded(1:padded-1) = 0; % negative flow assumed
            end

            % Find maximum correlation peak
            [peak, pos] = max(R_padded);
            % Apply threshold
            if real(peak) < thres_R 
                lag = 0;
            else
                % Overcome quantification constraints - interpolate peak
                % Ensure the peak is not on the edge
                if (pos > 1) && (pos < length(R_padded)) 
                    % Exploit adjacent values
                    Rm1 = real(R_padded(pos-1));
                    Rp1 = real(R_padded(pos+1));
                    R0 = real(R_padded(pos));                
                    % Quadratic interpolation formula
                    interp_offset = (Rp1 - Rm1) / (2*(Rp1 - 2*R0 + Rm1)); 
                    % ensure expected value is in range [-1,1]
                    if abs(interp_offset) <= 1  
                    % Total lag with interpolation
                        lag = -(pos - padded) - interp_offset;
                    else
                        lag = -(pos - padded); % dismiss interpolation
                    end
                else
                    % No interpolation if the peak is at boundary
                    lag = -(pos - padded);
                end
            end

            displacement = lag * pixel_size; % use lag to obtain the relative distance in pixel
            velocities_matrix(i,xi,frame) = displacement / (Tprf * no_lines);
        end
    end
    % Filter
    velocities_matrix_filtered(:,:,frame) = medfilt2(velocities_matrix(:,:,frame), [5 1]);
end

%% Averaging the velocity profiles
% Compute the mean over selected frames
% 1st frame is faulty due to characteristics
% of echo cancelation, therefore the avg is build with frames 2-34 to
% obtain 33 estimates
if sim_data
    % since no echo cancelation is applied, all frame can be used for avg
    velocities_avg = mean(velocities_matrix_filtered, 3); 
else
    velocities_avg = mean(velocities_matrix_filtered(:, :, 2:end), 3);
end

%% Theoretical flow profile
% Spatial grid
x = linspace(-5e-3, 5e-3, 21);  % Lateral: -5 mm to 5 mm
z = depth_range(1)*0.1:0.1:depth_range(end)*0.1; % Axial

[X, Z] = meshgrid(x_positions_mm, depth_range * 0.1);

% Flow parameters
if sim_data
    z_min = 26;
    z_max = 34;
    v_max = 0.3;  % Maximum velocity m/s
else
    z_min = 30;
    z_max = 40;
    v_max = -0.3;  % Maximum velocity m/s
end
z_center = (z_min + z_max) / 2;
H = (z_max - z_min) / 2;  % Half height for parabola

% Initialize flow field
V = zeros(size(Z));

if parabolic
    % For each lateral x position, assign a vertical parabolic profile
    for x0 = x
        % Find indices in the mesh closest to current x0
        [~, x_idx] = min(abs(x - x0));
        
        % Get the axial slice
        for z_idx = 1:length(z)
            if z(z_idx) >= z_min && z(z_idx) <= z_max
                z_rel = z(z_idx) - z_center;
                V(z_idx, x_idx) = v_max * (1 - (z_rel / H)^2);
            end
        end
    end
else
    % For each lateral x position, assign a vertical plug profile
    for x0 = x
        % Find indices in the mesh closest to current x0
        [~, x_idx] = min(abs(x - x0));
        
        % Get the axial slice
        for z_idx = 1:length(z)
            if z(z_idx) >= z_min && z(z_idx) <= z_max
                V(z_idx, x_idx) = v_max;
            end
        end
    end

end

%% Visu
% Sets default font size
fontSize=14;
set(groot, 'DefaultAxesFontSize', fontSize);       
set(groot, 'DefaultColorbarFontSize', fontSize);   
set(groot, 'DefaultTextFontSize', fontSize);       

% Video of averaging
figure;
for t = 1:size(velocities_matrix_filtered, 3)
    imagesc(x_positions_mm, depth_range * 0.1,velocities_matrix_filtered(:,:,t))
    title(['Velocity Frame ', num2str(t)]);
    xlabel('Lateral (mm)'), ylabel('Depth (mm)')
    colorbar;
    if sim_data
        colormap(jet);
    else 
        colormap(flipud(jet));
    end
    clim([min(v_max,0) max(v_max,0)]);
    pause(0.05)
end

showDepth = 325; % px
% Profiles
lateral_axis = (-(num_lateral_pixels-1)/2:(num_lateral_pixels-1)/2) * 0.1;  % Convert to mm
figure;
plot(lateral_axis,real(profiles(showDepth-depth_range(1),:,1)),'LineWidth',1,'DisplayName','Profile 1'); hold on;
plot(lateral_axis,real(profiles(showDepth-depth_range(1),:,2)),'LineWidth',1,'DisplayName','Profile 2');
ylabel('Amplitude');
xlabel('Lateral (mm)');
title(['Lateral profile at depth: ',num2str(showDepth/10),' mm']);
legend; grid on;

% Data vector
lateral_sample = (-(sampleSize-1)/2:(sampleSize-1)/2) * 0.1;  % Convert to mm
lateral_search = (-(searchSize-1)/2:(searchSize-1)/2) * 0.1;  % Convert to mm
figure; 
plot(lateral_sample,real(dataVecs1(showDepth-depth_range(1),:)),'LineWidth',1,'DisplayName','Data vec 1'); hold on;
% plot(lateral_search,real(dataVecs2(showDepth-depth_range(1),:)),'--','LineWidth',1,'DisplayName','Data vec 2');
ylabel('Amplitude');
xlabel('Lateral (mm)');
title(['Windowed data vector at depth: ',num2str(showDepth/10),' mm']);
legend; grid on;

% Theoretical
depth_axis = (depth_range) * 0.1;  % Convert to mm
figure('Position',[100,100,400,550]);
V_idx = find(x==0);
plot(V(:,V_idx),depth_axis, 'LineWidth',1.5); hold on;
set(gca, 'YDir', 'reverse');
ylabel('Depth (mm)');
xlabel('Velocity (m/s)');
yticks(20:2:40);
t=title('Theoretical velocity');
t.Units = 'normalized';         % Use normalized coordinates (0 to 1)
t.Position(2) = 1.05;           % Move title higher (default is around 1.01)
grid on;

% Velocity
figure('Position',[100,100,400,550]);
% plot(velocities_matrix(:,1),depth_axis, 'LineWidth',1); hold on;
plot(velocities_matrix_filtered(:,1),depth_axis,'LineWidth',1.5); hold on;
set(gca, 'YDir', 'reverse');
ylabel('Depth (mm)');
xlabel('Velocity (m/s)');
t=title('Velocity estimation');
t.Units = 'normalized';         % Use normalized coordinates (0 to 1)
t.Position(2) = 1.05;           % Move title higher (default is around 1.01)
grid on;

%% Flow profiles
fontSize = 12;
% Theoretical
figure('Position',[100,100,400,550]);
imagesc(x*1000, depth_range * 0.1, V);
cb = colorbar;
if sim_data
    colormap(jet);
else 
    colormap(flipud(jet));
end
xlabel('Lateral Position (mm)');
ylabel('Depth (mm)');
ylabel(cb, 'Velocity (m/s)','FontSize', fontSize+1)
t=title('Theoretical Flow Profile');
t.Units = 'normalized';         % Use normalized coordinates (0 to 1)
t.Position(2) = 1.05;           % Move title higher (default is around 1.01)
axis ij equal tight;
set(gca, 'YDir', 'reverse','FontSize', fontSize);

% Estimate
figure('Position',[100,100,400,550]);
imagesc(x_positions_mm, depth_range * 0.1, velocities_avg);  % x in mm, depth in mm
set(gca, 'YDir', 'reverse');
xlabel('Lateral Position (mm)');
ylabel('Depth (mm)');
cb=colorbar;
if sim_data
    colormap(jet);
else 
    colormap(flipud(jet));
end
ylabel(cb, 'Velocity (m/s)','FontSize', fontSize+1)
clim([min(v_max,0) max(v_max,0)]);
t=title('Acquired Data - Parabolic Flow Profile');
t.Units = 'normalized';         % Use normalized coordinates (0 to 1)
t.Position(2) = 1.05;           % Move title higher (default is around 1.01)
axis ij equal tight;
set(gca, 'YDir', 'reverse','FontSize', fontSize);