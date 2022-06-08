function [output] = resample_image(input,voxel_size,voxelsize_new)
    matrix_size=size(input);
    if size(matrix_size,2)>3
        matrix_size_new = [round(matrix_size(1:3).*voxel_size(1:3)./voxelsize_new(1:3)) matrix_size(4:end)];
    else
        matrix_size_new = round(matrix_size(1:3).*voxel_size(1:3)./voxelsize_new(1:3));
    end
    measc_k = fft3c(input,11);
    measc_k_resample = zeros(matrix_size_new);
    if voxel_size(1)<= voxelsize_new(1)
        x_ori = round(matrix_size(1)/2-matrix_size_new(1)/2+1:matrix_size(1)/2+matrix_size_new(1)/2);
        x_resample = 1:matrix_size_new(1);
    else
        x_ori = 1:matrix_size(1);
        x_resample =  round(matrix_size_new(1)/2-matrix_size(1)/2+1:matrix_size(1)/2+matrix_size_new(1)/2);
    end
    if voxel_size(2) <= voxelsize_new(2)
        y_ori =round( matrix_size(2)/2-matrix_size_new(2)/2+1:matrix_size(2)/2+matrix_size_new(2)/2);
        y_resample = 1:matrix_size_new(2);
    else
        y_ori = 1:matrix_size(2);
        y_resample = round(matrix_size_new(2)/2-matrix_size(2)/2+1:matrix_size(2)/2+matrix_size_new(2)/2);
    end
    if voxel_size(3) <= voxelsize_new(3)
        z_ori =round( matrix_size(3)/2-matrix_size_new(3)/2+1:matrix_size(3)/2+matrix_size_new(3)/2);
        z_resample = 1:matrix_size_new(3);
    else
        z_ori = 1:matrix_size(3);
        z_resample =  round(matrix_size_new(3)/2-matrix_size(3)/2+1:matrix_size(3)/2+matrix_size_new(3)/2);
    end
    measc_k_resample(x_resample,y_resample,z_resample,:) = measc_k(x_ori,y_ori,z_ori,:)...
        .*(tukeywin(size(x_ori,2),0.5)*tukeywin(size(y_ori,2),0.5)');%.*reshape(tukeywin(size(z_ori,2),1),[1 1 size(z_ori,2)]);
    if isreal(input) && min(input(:))>=0
        output = prod(voxel_size./voxelsize_new)*abs(ifft3c(measc_k_resample,11));
    elseif isreal(input) && min(input(:))<0
        output = prod(voxel_size./voxelsize_new)*real(ifft3c(measc_k_resample,11));
    else
        output = prod(voxel_size./voxelsize_new)*(ifft3c(measc_k_resample,11));
    end
end

% .*(tukeywin(matrix_size(1),0.5)*tukeywin(matrix_size(2),0.5)')