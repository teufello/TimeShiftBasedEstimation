function LRI = dyn_imageFormation(rf_data, pixelMap, VS, arrayPos, c, fs)
    % Extract parameters
    si = size(pixelMap);
    LRI = zeros(si(1:2));
    t = (0:size(rf_data,1)-1) / fs;
    % Grid
    numElements = length(arrayPos);
    x = 1:numElements;
    [X, T] = meshgrid(x, t);
    X2 = repmat(x, [si(2), 1]);  % same x-axis for each row

    for i = 1:si(1)
        delay_line = nan(si(2), numElements);  % [lateral points x elements]
        apod_line = zeros(si(2), numElements); % corresponding apodization

        for j = 1:si(2)
            pixel = [pixelMap(i,j,1), pixelMap(i,j,2)];
            [full_delay, apod_vec] = dyn_delayCal(pixel, [i,j], VS, arrayPos, c, si(2));

            delay_line(j,:) = full_delay';
            apod_line(j,:) = apod_vec';
        end

        % Interpolation
        value = interp2(X, T, rf_data, X2, delay_line, 'cubic');

        % Apodization masks the relevant values
        value_weighted = value .* apod_line;

        % Row-wise summation -> LRI
        LRI(i,:) = sum(value_weighted, 2, 'omitnan');
    end
end
