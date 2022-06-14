clear
rmpath(genpath('./STISuite_V3.0/'));
addpath(genpath('./STISuite_V3.0/'));
path = 'D:\1001613_3\QSM\NII';
path_fsl = '/mnt/d/1001613_3/QSM/NII';

TE1 = '0.00942';
TE2 = '0.01970';
fn_pha1 = dir([path '/PHA_TE1/*.nii.gz']);

num_channels = length(fn_pha1);

for j = 1:num_channels
    phs1(:, :, :, j) = double(niftiread([path '/PHA_TE1/' fn_pha1(j).name]));
end

dim = size(phs1); 
M1 = zeros(dim);
fn_mag1 = dir([path '/MAG_TE1/*.nii.gz']);

for j = 1:num_channels
    M1(:, :, :, j) = double(niftiread([path '/MAG_TE1/' fn_mag1(j).name]));
end

phase1 = ((phs1 - 2048) / 2048) * pi;
S1 = M1 .* exp(1i * phase1);

clear M1 phase1 phs1 fn_pha1 fn_mag1;

phs2 = zeros(dim);
fn_pha2 = dir([path '/PHA_TE2/*.nii.gz']);

for j = 1:num_channels
    phs2(:, :, :, j)=double(niftiread([path '/PHA_TE2/' fn_pha2(j).name]));
end

M2 = zeros(dim);
fn_mag2 = dir([path '/MAG_TE2/*.nii.gz']);

for j = 1:num_channels
    M2(:, :, :, j) = double(niftiread([path '/MAG_TE2/' fn_mag2(j).name]));
end

phase2 = ((phs2 - 2048) / 2048) * pi;
clear phs2 j fn_mag2 fn_pha2;

S2 = M2 .* exp(1i * phase2);

clear M2 phase2;

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
clear voxel_size voxelsize_new start en matrix_size matrix_size_new measc_k measc_k_down

dim = size(S2); 
%%  Coil comb:MCPC-3D-S
hip = zeros(dim(1:3));
for iCha = 1:size(S1, 4)
    complexDifference = S2(:, :, :, iCha) .* conj(S1(:, :, :, iCha));
    hip = hip + complexDifference;
end

clear complexDifference iCha;

hip_abs = abs(hip);
hip_angle = angle(hip);

Mag_brain = mean(sqrt(abs(S1).^2+abs(S2).^2),4);
nii_info = niftiinfo([path(1:end-3) 'QSM_mcpc3Ds_MAG_TE1.nii.gz']);
nii_info.Datatype = 'double';
nii_info.ImageSize = dim(1:3);
nii_info.PixelDimensions = [1.05 1 3];
mkdir([path '/QSM'])
niftiwrite(hip_abs, [path '/QSM/hip_abs.nii'], nii_info, ...
           'Compressed', true);
niftiwrite(hip_angle, [path '/QSM/hip_angle.nii'], nii_info,...
           'Compressed', true);
%% prelude unwrapping
system(['bash -c "source ~/.profile && bet2 ', path_fsl, ...
        '/QSM/hip_abs.nii.gz ', path_fsl, ...
        '/QSM/QSM -m -f 0.35"']);
system(['bash -c "source ~/.profile && prelude -a ', path_fsl, '/QSM/hip_abs.nii.gz -p ', ...
        path_fsl, '/QSM/hip_angle.nii.gz -u ', path_fsl, ...
        '/QSM/hip_uw.nii.gz -m ', path_fsl, '/QSM/QSM_mask.nii.gz"']);
%% R2star mapping
unwrappedHip=double(niftiread([path '/QSM/hip_uw.nii.gz']));

TEs(1)=str2double(TE1);
TEs(2)=str2double(TE2);

R2star_img = zeros([dim(1:3) 2]);
pd = zeros(dim(1:3));
for i = 1:size(S1,1)
    for j = 1:size(S1,2)
        for k = 1:size(S1,3)
            [U,S,~] = svd(squeeze([S1(i,j,k,:) S2(i,j,k,:)]),'econ');
            R2star_img(i,j,k,:)=abs(S(1,1))*abs(U(:,1));
            pd(i,j,k) = abs(S(1,1));
        end
    end
end
r2star = -log(R2star_img(:,:,:,2)./R2star_img(:,:,:,1))/(TEs(2) - TEs(1));
r2star(isnan(r2star))=0;
save([path(1:end-7) 'mat\r2starimg_svd_2e.mat'],'R2star_img')
niftiwrite(r2star, [path '/QSM/r2star.nii'], nii_info, 'Compressed', true);
%% 
       
scale = TEs(1) / (TEs(2) - TEs(1));
unwrappedHip = unwrappedHip * scale;

hipComplex = exp(1i * unwrappedHip);

size_compl = size(hipComplex);

po = complex(zeros([size_compl(1:3) num_channels], 'double'));
for iCha = 1:num_channels
    po_ang=angle(exp(1i * (angle(S1(:, :, :, iCha)) - unwrappedHip)));
    po_double=double(abs(S1(:, :, :, iCha)) .* exp(1i * po_ang));
    po(:, :, :, iCha) = po_double;
end

clear hip hip_abs hip_angle hipComplex po_ang po_double unwrappedHip scale;

real_smmothed=double(zeros(size(po)));
imag_smmothed=double(zeros(size(po)));

for j=1:num_channels
    real_smmothed(:, :, :, j) = imgaussfilt3(real(po(:, :, :, j)), 4, ...
                                             'FilterDomain', 'spatial');
    imag_smmothed(:, :, :, j) = imgaussfilt3(imag(po(:, :, :, j)), 4, ...
                                             'FilterDomain', 'spatial');
end

clear po;


po_c = (real_smmothed + 1i * imag_smmothed) ./ abs(real_smmothed + 1i * imag_smmothed);

clear real_smmothed imag_smmothed;

S1_c = S1 .* squeeze(conj(po_c));

%% 
combined1 = weightedCombination(S1_c, abs(S1_c));
combined1(~isfinite(combined1)) = 0;
mask = double(niftiread([path '/QSM/QSM_mask.nii.gz']));
nii = angle(combined1) .* mask;

nii_info.Datatype = 'single';
niftiwrite(single(nii), [path '/QSM/PHASE_TE1.nii'], nii_info, ...
           'Compressed',true);

clear nii S1_c;

S2_c = S2 .* squeeze(conj(po_c));

clear po_c;

combined2 = weightedCombination(S2_c, abs(S2_c));

clear S1_c S2_c;

combined2(~isfinite(combined2)) = 0;
nii = angle(combined2) .* mask;
niftiwrite(single(nii), [path '/QSM/PHASE_TE2.nii'], nii_info, ...
           'Compressed',true);

clear hip hip_abs hip_angle nii po_c S1_c S2_c smoothed_weight ;
clear unwrappedHip ans iCha j num_channels size_compl dim;


% extract DICOM info and update mask
phase1 = mask .* angle(combined1);
phase2 = mask .* angle(combined2);

clear combined1 combined2;

di = [path(1:end-3) '15_20190928\'];
flist=dir([di '*.dcm']); 


dcm_info = dicominfo([di flist(1).name]);
ori = dcm_info.ImageOrientationPatient;
Xz = ori(3);
Yz = ori(6);
Zxyz = cross(ori(1:3), ori(4:6));
Zz = Zxyz(3);
H = [-Xz, -Yz, Zz];

voxsz = [1.05 1 3];
B0 = dcm_info.MagneticFieldStrength;

clear ori Xz Yz Zxyz Zz flist R2star_img siz U in ind formatSpec forname
%% 
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
%% phase unwrap
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
%%  v-SHARP and Dipole inversion

[dB_vsf, mask_vsf]=V_SHARP(phs_comb, mask, 'voxelsize', voxsz, 'smvsize', 12);

clear mask;

mask_vsf(map < 0.7) = 0;
mask_vsf(map2 < 0.6) = 0;

for ii = 1:dim(3)
    mask_vsf(:, :, ii) = bwareaopen(mask_vsf(:, :, ii), 200);
    mask_vsf(:, :, ii) =~ bwareaopen(~mask_vsf(:, :, ii), 30);
end
mask_vsf = imfill(mask_vsf,26,'holes');
clear map phs_comb map2 tmp tmp2 i ii k dcm_info S1 S2  TE1 TE2
r2s_k = fft2c(r2star.*mask_vsf);
L = tukeywin(size(r2s_k,1),0.5)*tukeywin(size(r2s_k,2),0.5)';
r2star = abs(ifft2c(r2s_k.*L));
voxel_size = [1.05 1 2];
voxelsize_new = [1.05 1 3];
mkdir([path(1:end-7) 'result'])
save([path(1:end-7) 'mat\for_xsep.mat'],'dB_vsf','r2star','Mag_brain','mask_vsf','TE','B0','H','voxsz','path','voxelsize_new','voxel_size')

%% 
qsm_iLSQR_vsf = QSM_iLSQR(dB_vsf, mask_vsf, 'TE', TE, 'B0', B0, ...
            'H', H, 'padsize', [64 64 64], ...
            'voxelsize', voxsz);
niftiwrite(single(qsm_iLSQR_vsf * 1000), [path '/QSM/QSM.nii'], ...
           nii_info, 'Compressed', true);
figure(1);imshow_3df(qsm_iLSQR_vsf/2/pi/123,[-.1 .1])
x_sa = qsm_iLSQR_vsf/2/pi/123;
save([path(1:end-7) 'result\QSM.mat'],'qsm_iLSQR_vsf','x_sa')
