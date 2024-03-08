function bin_map = phasevariance_nonlin2(mask, phase, radius)
    
% Description: Script to calculate phase reliability maps (to remove voxels in the vicinity of sinuses)
% 
% Authors: Chaoyue Wang, Benjamin C. Tendler & Karla L. Miller
% 
% Copyright 2021 University of Oxford
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
% http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% radius in mm e.g. 3mm
    

    dim = size(phase);

    dimX = dim(1)/13;
    dimY = dim(2)/10;

    phase2 = zeros(dim(1), dim(2), dim(3) + radius * 2);
    phase2(:, :, radius + 1:radius + dim(3)) = phase;
    phase = phase2;

    mask2 = zeros(dim(1), dim(2), dim(3) + radius * 2);
    mask2(:, :, radius + 1:radius + dim(3)) = mask;
    mask = mask2;
    clear mask2;

    dim=size(phase);
    [cy, cx, cz] = meshgrid(double(-dim(2) / 2:(dim(2) / 2 - 1)),...
                double(-dim(1) / 2:(dim(1) / 2 - 1)),...
                double(-dim(3) * 30 / 2:(dim(3) * 30 / 2 - 1)));
        
    index = (cy) .^ 2 + (cx) .^ 2 + (cz) .^ 2 <= (radius * 10) ^ 2;
    rho_temp = zeros(size(cx), 'double');
    rho_temp(index) = 1;
        
    X(1:dimX) = 13;
    Y(1:dimY) = 10;
    Z(1:dim(3)) = 30;

    A = mat2cell(rho_temp, [X], [Y], [Z]);
    rho_temp = cellfun(@mean2, A);

    rho = double(zeros(dim));
    rho(((dim(1)-dimX)/2+1):((dim(1)-dimX)/2+dimX),((dim(2)-dimY)/2+1):((dim(2)-dimY)/2+dimY),:)=rho_temp;

    cdata = exp(1i * phase) .* mask;
    cdatalp = fftshift(ifftn(fftn(cdata) .* fftn(rho)));
    fz = abs(cdatalp);
    
    fm = mask .* fftshift(ifftn(fftn(mask) .* fftn(rho)));
    bin_map = fz ./ (fm + eps * ones(dim));
    bin_map(mask == 0) = 1;
    bin_map = bin_map(:, :, (radius + 1):(dim(3) - radius));

end

