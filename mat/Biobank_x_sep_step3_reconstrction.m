clear
path = 'D:\UKB\Real_data\1004671_3\';
load([path 'mat\for_xsep.mat'],'dB_vsf','r2star','Mag_brain','mask_vsf','TE','B0','H','voxsz','voxelsize_new','voxel_size')
load([path 'result\QSM.mat'],'x_sa')
params.voxel_size = voxsz; % [mm unit]
params.CF = 123e6; 
params.b0_dir = H; 
params.Dr_pos = 137; 
params.Dr_neg = params.Dr_pos;
params.lambda = 1; 
params.lambda_CSF = 1; 

%% T2: MPRAGE+FLAIR Ds, T2*: DL Ds
dire = dir([path 'mat\inf_from_gpu\inf*']);
load([path 'mat\inf_from_gpu\' dire(1).name])
load([path 'mat\inf_from_gpu\' dire(2).name])
r2star_deep = permute(-10*(log(R2star+eps)),[2 3 1]).*mask_vsf;
load([path 'mat\registered_images.mat'])
R2_MPR_FLA = -10*log(permute(T2,[2 3 1])).*(T2 ~= 0);
R2_MPR_FLA(isinf(R2_MPR_FLA) | isnan(R2_MPR_FLA))=0;
R2_MPR_FLA_resample = resample_image(rot90(R2_MPR_FLA,-2),voxel_size,voxelsize_new);

[ r2_MPR_FLA_reg_r2s, ~] = fsl_flirt2( r2star_deep.*mask_vsf,R2_MPR_FLA_resample,6,voxelsize_new);
% [ r2_MPR_FLA_reg_r2s, ~] = fsl_flirt2( r2star_deep.*mask_vsf,rot90(R2_MPR_FLA,-2),6,voxelsize_new);
r2prime_MPR_FLA_r2sdeep= (r2star_deep-r2_MPR_FLA_reg_r2s).*mask_vsf;
r2prime_MPR_FLA_r2sdeep(r2prime_MPR_FLA_r2sdeep<0) = 0;
[x_pos_MPR_FLA_r2sdeep, x_neg_MPR_FLA_r2sdeep, x_tot_MPR_FLA_r2sdeep]= x_sep_SA(dB_vsf/2/pi/TE,r2prime_MPR_FLA_r2sdeep, mask_vsf,params,x_sa);
figure(1);imshow_3df(fliplr(x_pos_MPR_FLA_r2sdeep),[0 .1], fliplr(x_neg_MPR_FLA_r2sdeep),[0 .1],fliplr(x_tot_MPR_FLA_r2sdeep),[-.1 .1])
save([path 'result\x_sep_MPRAGE_FLAIR_r2stardeep.mat'],'r2star_deep','r2_MPR_FLA_reg_r2s','r2prime_MPR_FLA_r2sdeep','x_pos_MPR_FLA_r2sdeep','x_neg_MPR_FLA_r2sdeep','x_tot_MPR_FLA_r2sdeep')
clear R2 R2_MPR_FLA r2_MPR_FLA_reg_r2s R2_MPR_FLA_resample r2_zpd r2prime_MPR_FLA_r2sdeep R_mat2 T1w_bet_reg T2 t2w_reg t2w_resample t2w_zpd  x_neg_ds_r2s x_pos_ds_r2s x_tot_ds_r2s mask Flair_bet_reg R2star r2star_deep

%% T2: Cons, T2*: DL

load([path 'mat\inf_from_gpu\' dire(2).name])
r2star_deep = permute(-10*(log(R2star+eps)),[2 3 1]).*mask_vsf;

r2_con = 13.98;
r2prime_con = r2star_deep.*mask_vsf-r2_con;
r2prime_con(r2prime_con<0) = 0;
[x_pos_sa_con, x_neg_sa_con, x_tot_sa_con]= x_sep_SA(dB_vsf/2/pi/TE,r2prime_con, mask_vsf,params,x_sa);
figure(2);imshow_3df(fliplr(x_pos_sa_con),[0 .1], fliplr(x_neg_sa_con),[0 .1],fliplr(x_tot_sa_con),[-.1 .1])

save([path 'result\x_sep_constant_r2stardeep_cal.mat'],'r2_con','r2prime_con','x_pos_sa_con','x_neg_sa_con','x_tot_sa_con')

