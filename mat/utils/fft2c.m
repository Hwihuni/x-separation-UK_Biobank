function [im] = fft2c(dat,window_on, alpha)
% FFT2C performs a centered fft2
%function [im] = fft2c(dat,window_on, alpha)

if nargin < 2;
    window_on = 0;
end
if nargin < 3;
    alpha = 0.5;
end



if window_on == 0
    if size(size(dat),2) <= 3
        for i = 1:size(dat,3);
            im(:,:,i) = fftshift(fft2(ifftshift(dat(:,:,i))));
        end
    elseif size(size(dat),2) ==4
        for j = 1:size(dat,4)
            for i = 1:size(dat,3);
                im(:,:,i,j) = fftshift(fft2(ifftshift(dat(:,:,i,j))));
            end
        end
        
    elseif size(size(dat),2) ==5
        
        for k = 1:size(dat,5)
            for j = 1:size(dat,4)
                for i = 1:size(dat,3);
                    im(:,:,i,j,k) = fftshift(fft2(ifftshift(dat(:,:,i,j,k))));
                end
            end
        end
    end
    im = sqz(im);
else
    window_matrix = tukeywin(size(dat,1),alpha)*(tukeywin(size(dat,2),alpha))';
    for i = 1:size(dat,3);
        im(:,:,i) = fftshift(fft2(ifftshift(dat(:,:,i).*window_matrix)));
    end
end

