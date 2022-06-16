clear
path = 'D:\UKB\Real_data\1004671_3\';
voxelsize = [1.05 1 2];

voxelsize_T1 = niftiinfo([path 'T1/T1.nii.gz']).PixelDimensions;
voxelsize_Flair= niftiinfo([path 'T2_FLAIR/T2_FLAIR.nii.gz']).PixelDimensions;

T1_pre =fliplr( double(niftiread([path 'T1/T1.nii.gz'])));
Flair_pre = fliplr(double(niftiread([path 'T2_FLAIR/T2_FLAIR.nii.gz'])));

T1_pre_vox = resample_image(T1_pre,voxelsize_T1,voxelsize);
Flair_pre_vox = resample_image(Flair_pre,voxelsize_Flair,voxelsize);

[T1w_reg_bet, ~] = fsl_bet(T1_pre_vox(:,:,20:end),0.35,voxelsize);
[Flair_reg_bet, ~] = fsl_bet(Flair_pre_vox(:,:,20:end),0.35,voxelsize);
Flair_bet_reg = Flair_reg_bet;
[ T1w_bet_reg, ~] = fsl_flirt2(Flair_reg_bet,T1w_reg_bet,6,voxelsize);
matrix_size = [198,256,size(Flair_bet_reg,3)];
padsize = [198-size(T1w_bet_reg,1),256-size(T1w_bet_reg,2),0];

T1w_bet_reg = padarray(T1w_bet_reg,floor(padsize/2),0,'post');
Flair_bet_reg = padarray(Flair_bet_reg,floor(padsize/2),0,'post');
T1w_bet_reg = padarray(T1w_bet_reg,padsize-floor(padsize/2),0,'pre');
Flair_bet_reg = padarray(Flair_bet_reg,padsize-floor(padsize/2),0,'pre');

mask = (T1w_bet_reg~=0).*(Flair_bet_reg~=0);


figure;imshow_3df(T1w_bet_reg.*mask,Flair_bet_reg.*mask)

mkdir([path 'mat'])
save([path 'mat\registered_images.mat'],'T1w_bet_reg', 'Flair_bet_reg','mask');


Mask = mask;
t1_norm = T1w_bet_reg/prctile(T1w_bet_reg(T1w_bet_reg~=0),99).*Mask;
flair_norm = Flair_bet_reg/prctile(Flair_bet_reg(Flair_bet_reg~=0),99).*Mask;

r2_norm = Mask;
figure;imshow_3df(t1_norm,flair_norm,r2_norm,'range',[0 1]);

save([path 'mat\Data_to_gpu_t2map_' path(end-9:end-1) '.mat'],'r2_norm','flair_norm','t1_norm','Mask')
%% 
voxel_size = [0.8 0.8 3];
voxelsize_new = [1.05 1 3];
load([path 'mat\r2starimg_svd_2e.mat'])

r2star_img = pd.*R2star_img;
[~,mask] = fsl_bet(mean(R2star_img,4),0.3,voxelsize_new);
r2star_norm = mask;
img = r2star_img.*mask/prctile(r2star_img(mask~=0),99);


Mask = mask;
save([path 'mat\Data_to_gpu_t2starmap_' path(end-9:end-1) '.mat'],'r2star_norm','img','Mask')
figure;imshow_3df(img)
mkdir([path 'mat\inf_from_gpu'])

