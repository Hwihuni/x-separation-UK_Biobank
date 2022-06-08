function im = fft2c(dat)
% FFT2C performs a centered fft2

if size(size(dat),2) <=3
    for i = 1:size(dat,3);
        im(:,:,i) = fftshift(ifft2(fftshift(dat(:,:,i))));
    end
elseif size(size(dat),2) ==4
    for j = 1:size(dat,4)
        for i = 1:size(dat,3);
            im(:,:,i,j) = fftshift(ifft2(fftshift(dat(:,:,i,j))));
        end
    end
    
elseif size(size(dat),2) ==5
    
    for k = 1:size(dat,5)
        for j = 1:size(dat,4)
            for i = 1:size(dat,3);
                im(:,:,i,j,k) = fftshift(ifft2(fftshift(dat(:,:,i,j,k))));
            end
        end
    end
end
im = (im);



