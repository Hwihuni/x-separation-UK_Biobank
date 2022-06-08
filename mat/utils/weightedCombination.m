function combined = weightedCombination(image, sensitivity)

% This script was adapted from ASPIRE (https://github.com/korbinian90/ASPIRE)

% Please cite: Eckstein, K., et al. (2018), Computationally Efficient Combination of Multi-channel Phase Data From Multi-echo Acquisitions (ASPIRE). 
% Magn. Reson. Med, 79: 2996-3006. https://doi.org/10.1002/mrm.26963

    image = double(image);
    sensitivity = double(sensitivity);

    dimension = size(image);
    combined = zeros(dimension(1:(end - 1)), 'double');

    for iChannel = 1:dimension(4)
        combined = combined + image(:, :, :, iChannel) .* sensitivity(:, :, :, iChannel);
    end
    combined = combined ./ sqrt(abs(combined));

end
