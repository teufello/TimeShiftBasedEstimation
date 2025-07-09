function [delay, apod_vec] = dyn_delayCal(pixel, idx, vs, arrayPos, c, num_lateral)
    % dyn_delayCal: Computes full-length delay vector with NaNs for inactive elements
    F_recv = 1;
    pitch = 2.3e-4;
    N_elements = 192;
    
    % Determine current central element based on lateral index
    cur_element = ceil(idx(2)/num_lateral*N_elements);
    
    % Determine number of active elements based on depth
    active_elements = ceil(pixel(2)*F_recv/pitch);

    % Compute aperture limits
    start_idx = cur_element - floor((active_elements/2));
    if start_idx < 1
        start_idx = 1;
    end

    end_idx = cur_element + floor((active_elements/2));
    if end_idx > N_elements
        end_idx = N_elements;
    end

    idx_active = start_idx:end_idx;
    active_array = arrayPos(idx_active, :);

    % Compute transmit path
    tx_path = sqrt(sum((vs - pixel).^2)) + vs(2);

    % Compute receive path
    rx_path = sqrt(sum((active_array - pixel).^2, 2));

    % Fill full-length delay vector with NaNs for inactive elements
    delay = nan(N_elements,1);

    % Compute delay for active elements only
    delay(idx_active) = (tx_path + rx_path) / c;

    % Apodization vector: Hanning for active elements, 0 elsewhere
    apod_vec = zeros(N_elements,1);
    apod_vec(idx_active) = hann(length(idx_active));
end