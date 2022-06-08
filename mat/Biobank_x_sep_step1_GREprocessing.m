clear
rmpath(genpath('./STISuite_V3.0/'));
addpath(genpath('./STISuite_V3.0/'));
path = 'D:\1001613_3\QSM';
path_fsl = '/mnt/d/1001613_3/QSM';

TE1 = '0.00942';
TE2 = '0.01970';

phs1 = double(niftiread([path '/QSM_mcpc3Ds_PHA_TE1.nii.gz']));
phs2 = double(niftiread([path '/QSM_mcpc3Ds_PHA_TE2.nii.gz']));
M1 = double(niftiread([path '/QSM_mcpc3Ds_MAG_TE1.nii.gz']));
M2 = double(niftiread([path '/QSM_mcpc3Ds_MAG_TE2.nii.gz']));
S1 = M1.*exp(1i*phs1);
S2 = M2.*exp(1i*phs2);
voxel_size = [0.8, 0.8, 3];
matrix_size=size(S2(:,:,:,1));
voxelsize_new = [1.05 1 3];

matrix_size_new = round(matrix_size.*voxel_size./voxelsize_new);
start = matrix_size/2-matrix_size_new/2+1;
en = matrix_size/2+matrix_size_new/2;

measc_k = fft3c(S1,11);
measc_k_down = measc_k(start(1):en(1),start(2):en(2),start(3):en(3),:);
S1 = ifft3c(measc_k_down,11);

measc_k = fft3c(S2,11);
measc_k_down = measc_k(start(1):en(1),start(2):en(2),start(3):en(3),:);
S2 = ifft3c(measc_k_down,11);
Mag_brain = mean(sqrt(abs(S1).^2+abs(S2).^2),4);

nii_info = niftiinfo([path '/QSM_mcpc3Ds_MAG_TE1.nii.gz']);
nii_info.Datatype = 'double';
nii_info.ImageSize = size(Mag_brain);
nii_info.PixelDimensions = [1.05 1 3];

niftiwrite(Mag_brain, [path '/hip_abs.nii'], nii_info, ...
           'Compressed', true);

system(['bash -c "source ~/.profile && bet2 ', path_fsl, ...
        '/hip_abs.nii.gz ', path_fsl, ...
        '/QSM -m -f 0.35"']);
%% 
TEs(1)=str2double(TE1);
TEs(2)=str2double(TE2);

R2star_img = zeros([size(Mag_brain),2]);
for i = 1:size(S1,1)
    for j = 1:size(S1,2)
        for k = 1:size(S1,3)
            R2star_img(i,j,k,1)=abs(S1(i,j,k));
            R2star_img(i,j,k,2)=abs(S2(i,j,k));
        end
    end
end
r2star = -log(R2star_img(:,:,:,2)./R2star_img(:,:,:,1))/(TEs(2) - TEs(1));
r2star(isnan(r2star))=0;
mkdir([path(1:end-3) 'mat'])
save([path(1:end-3) 'mat\r2starimg_svd_2e.mat'],'R2star_img')
niftiwrite(r2star, [path '/r2star.nii'], nii_info, ...
           'Compressed', true);
mask = double(niftiread([path '/QSM_mask.nii.gz']));
%% 
phase1 = angle(S1);
phase2 = angle(S2);
map = phasevariance_nonlin2(mask, phase1, 2);
map2 = phasevariance_nonlin2(mask, phase2, 2);
dim = size(phase1);

mask(map < 0.6) = 0;
mask(map2 < 0.5) = 0;

for ii = 1:dim(3)
    mask(:, :, ii) =  bwareaopen(mask(:, :, ii), 300);
    mask(:, :, ii) =~ bwareaopen(~mask(:, :, ii), 50);
end
mask = imfill(mask,26,'holes');
%% 
% phase unwrap
voxsz = [1.05 1 3];

tmp = zeros([196 230 48]);
tmp(1:end-1,:,:) = phase1;
[tmp2, ~] = MRPhaseUnwrap(tmp, 'voxelsize', voxsz, 'padsize', [64 64 64]);
uwphase1 = tmp2(1:end-1,:,:);
tmp = zeros([196 230 48]);
tmp(1:end-1,:,:) = phase2;
[tmp2, ~] = MRPhaseUnwrap(tmp, 'voxelsize', voxsz, 'padsize', [64 64 64]);
uwphase2 = tmp2(1:end-1,:,:);
T2s = 40;

W1 = TEs(1) * exp(-TEs(1) / T2s);

W2 = TEs(2) * exp(-TEs(2) / T2s);

phs_comb = double((W1 * uwphase1 + W2 * uwphase2) / (W1 + W2));

TE = (W1 * TEs(1) + W2 * TEs(2)) / (W1 + W2);

clear phase1 uwphase1 uwphase2 W1 W2 TEs T2s combined phase2;
%% 
% v-SHARP and Dipole inversion

[dB_vsf, mask_vsf]=V_SHARP(phs_comb, mask, 'voxelsize', voxsz, 'smvsize', 12);


mask_vsf(map < 0.7) = 0;
mask_vsf(map2 < 0.6) = 0;

for ii = 1:dim(3)
    mask_vsf(:, :, ii) = bwareaopen(mask_vsf(:, :, ii), 200);
    mask_vsf(:, :, ii) =~ bwareaopen(~mask_vsf(:, :, ii), 30);
end
mask_vsf = imfill(mask_vsf,26,'holes');
r2s_k = fft2c(r2star.*mask_vsf);
L = tukeywin(size(r2s_k,1),0.5)*tukeywin(size(r2s_k,2),0.5)';
r2star = abs(ifft2c(r2s_k.*L));
voxel_size = [1.05 1 2];
voxelsize_new = [1.05 1 3];
di = [path '\15_20190928\'];
flist=dir([di '*.DCM']); 


dcm_info = dicominfo([di flist(1).name]);
ori = dcm_info.ImageOrientationPatient;
Xz = ori(3);
Yz = ori(6);
Zxyz = cross(ori(1:3), ori(4:6));
Zz = Zxyz(3);
H = [-Xz, -Yz, Zz];
B0 = dcm_info.MagneticFieldStrength;

mkdir([path '\result'])
save([path(1:end-3) 'mat\for_xsep.mat'],'dB_vsf','r2star','Mag_brain','mask_vsf','TE','B0','H','voxsz','path','voxelsize_new','voxel_size')
% pause
%% 
qsm_iLSQR_vsf = QSM_iLSQR(dB_vsf, mask_vsf, 'TE', TE, 'B0', B0, ...
            'H', H, 'padsize', [64 64 64], ...
            'voxelsize', voxsz);
niftiwrite(double(qsm_iLSQR_vsf * 1000), [path '/QSM.nii'], ...
           nii_info, 'Compressed', true);
figure(1);imshow_3df(fliplr(qsm_iLSQR_vsf)/2/pi/123,[-.1 .1])
x_sa = qsm_iLSQR_vsf/2/pi/123;
save([path '\result\QSM.mat'],'qsm_iLSQR_vsf','x_sa')
