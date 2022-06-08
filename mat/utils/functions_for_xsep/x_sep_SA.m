

function [x_pos, x_neg, x_tot]= x_sep_SA(local_f,r2prime, mask_qsm,params,...
    x_qsm, mask_Unrel_Mag, mask_Unrel_Phs)
%
% susceptibility source separation (chi-separation) by removing streaking artifacts
% [x_pos, x_neg, x_tot]= x_sep_SA(local_f,r2prime, mask_qsm,params,...
%    x_qsm, mask_Unrel_Mag, mask_Unrel_Phs)
%
% %%%%%%%%%%% Output %%%%%%%%%%%
%   -x_pos: positive susceptibility [ppm]
%   -x_neg: negative susceptibility [ppm]
%   -x_tot: bulk susceptibility of voxel (i.e., x_pos-x_neg) [ppm]
%
% %%%%%%%%%%% Input %%%%%%%%%%%
%   -local_f       : background-removed local frequency map (y, x, z) [Hz unit]
%   -r2prime       : R2' (= R2*-R2) (y, x, z) [Hz unit]
%   -mask_qsm      : binary mask for x-separation (e.g., brain mask)
%   -params        : parameter set for algorithm (see below)
%                   % essenstial params
%                   -> voxel_size: (y, x, z) [mm unit]
%                   -> CF                : central frequency [Hz unit]
%                   -> B0_dir            : (y, x, z)
%                   % optional params
%                   -> Dr_pos            : default = 137 for in-vivo human brain at 3T [Hz/ppm]
%                   -> Dr_neg            : default = Dr_pos for in-vivo human brain at 3T [Hz/ppm]
%                   -> flag_strong_x_correction: default = 1; % 1: turn on compensation
%                   -> cg_max_iter_LSQR  : default = 50;
%                   -> cg_tol_LSQR       : default = 0.02;% 0.1 
%                   -> cg_max_iter_SA    : default = default = 200;
%                   -> cg_tol_SA         : default = 0.02;
%                   -> D2_thr            : default = 0.2; 
%                   -> prc_W_I_min       : default = 60;
%                   -> prc_W_I_max       : default = 99.8;
%                   -> WG_prc_min        : default = 50;
%                   -> WG_prc_max        : default = 70;
%                   -> pad_size          : (y, x, z) integer number for padding size
%                   -> details : see "fill_params" function in this .m file
%   -x_qsm         : (optional) QSM map which offers initial guess to help fast convergence and 
%                               remove artifacts from large susceptibility
%                               sources (if x_qsm is empty matrix, automatically
%                               generate x_qsm)
%   -mask_Unrel_Mag: (optional) binary mask for unreliable values in magnitude or phase map
%                   e.g., low-SNR voxel, voxels with too fast or slow relaxation
%   -mask_Unrel_Phs: (optional) binary mask for unreliable values in magnitude or phase map
%                   e.g., voxels with large B0 field inhomogeneity
%
%
% %%%%%%%%%%% Info %%%%%%%%%%%
% Reference: Hyeong-Geol Shin et al., Neuroimage, 2021
% - https://doi.org/10.1016/j.neuroimage.2021.118371
%
% copyright @ Hyeong-Geol Shin (contact: sin4109@gmail.com)
% Laboratory for Imaging Science and Technology
% Department of Electrical and Computer Engineering
% Seoul National University
%
% version 1.0.0 : 2021.07.23
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%


params = fill_params(params); % check and fill params using default parameters

if nargin <5
    x_qsm = []; % empty matrix to automatically generate weighting matrix
    mask_Unrel_Mag  = zeros(size(r2prime)); 
    mask_Unrel_Phs = zeros(size(r2prime));  
elseif nargin < 6
    mask_Unrel_Mag  = zeros(size(r2prime)); 
    mask_Unrel_Phs = zeros(size(r2prime));  
end


disp('+++++++++ Starting SA x-sep  +++++++++')
tic
% LSQR solution 
% W_I: weights from LSQR reconstruction
disp('Step 1: Calculating LSQR solution')
[x_pos_LSQR, x_neg_LSQR, W_I, x_qsm1] = x_sep_SA_step1_func(local_f,r2prime, mask_qsm, x_qsm, mask_Unrel_Mag, mask_Unrel_Phs, params);


% estimating and removing streaking artifacts 
% W_G_p: weights determined from gradients of fast positive susceptibility (y,x,z, 3)
% W_G_n: weights determined from gradients of fast negative susceptibility (y,x,z, 3)
disp('Step 2: Calculating streaking artifacts')
[x_pos, x_neg, W_G_p, W_G_n]= x_sep_SA_step2_func(x_pos_LSQR, x_neg_LSQR,local_f,r2prime, mask_qsm, x_qsm1, mask_Unrel_Mag, mask_Unrel_Phs, params);
toc
disp('+++++++++++++++ Done ++++++++++++++++++')

% figure(12); 
% imshow(cat(2, W_I(:,:,end/2).*mask_qsm(:,:,end/2), sqrt(sum(abs(W_G_p(:,:,end/2,:)).^2,4)),sqrt(sum(abs(W_G_n(:,:,end/2,:)).^2,4))),[0 1])
% title('W_I      W_G_p_o_s       W_G_n_e_g');

x_tot = x_pos - x_neg; 
end



%%

function params_filled = fill_params(params)
% fill empty fields in params
gyro_ = 42.58e6;
params_filled = params;

% field which must be filled by manually
%Unit:                 mm        Hz     a.u
field_essential = {'voxel_size','CF','b0_dir'};

for ii = 1:length(field_essential)
    if ~isfield(params_filled, field_essential{ii})
        error(['params.',field_essential{ii},' is missed']);
    end
end


% optional parameters (please, refer to MEDI toolbox)
field_scale = params.CF/(gyro_*2.9); %relative to 2.9T
if ~isfield(params_filled, 'Dr_pos')
    params_filled.Dr_pos = 137*field_scale; % 137 Hz/ppm @ 3T for in-vivo human brain
end
if ~isfield(params_filled, 'Dr_neg')
    params_filled.Dr_neg = params_filled.Dr_pos; % default: copy Dr_pos
end
if ~isfield(params_filled, 'pad_size')% zero padding size
    params_filled.pad_size = [16 16 16];
end
if ~isfield(params_filled, 'flag_strong_x_correction')
    params_filled.flag_strong_x_correction = 0; % 0 [no corre.] / 1 [compensate for strong suscep.]
end
if ~isfield(params_filled, 'strong_x_thresh')
    params_filled.strong_x_thresh = 0.2; % Threshold for strong suscep. [ppm]
end

% parameters related to iLSQR algorithm (Li, Neuroimage, 2015). 
if ~isfield(params_filled, 'cg_max_iter_LSQR')
    params_filled.cg_max_iter_LSQR = 50;
end
if ~isfield(params_filled, 'cg_tol_LSQR')
    params_filled.cg_tol_LSQR = 0.02;
end
if ~isfield(params_filled, 'cg_max_iter_SA')
    params_filled.cg_max_iter_SA = 50;
end
if ~isfield(params_filled, 'cg_tol_SA')
    params_filled.cg_tol_SA = 0.02;
end
if ~isfield(params_filled, 'D2_thr')
    params_filled.D2_thr = 0.2;
end
if ~isfield(params_filled, 'prc_W_I_min')
    params_filled.prc_W_I_min = 60;
end
if ~isfield(params_filled, 'prc_W_I_max')
    params_filled.prc_W_I_max = 99.8;
end
if ~isfield(params_filled, 'WG_prc_min')
    params_filled.WG_prc_min = 50;
end
if ~isfield(params_filled, 'WG_prc_max')
    params_filled.WG_prc_max = 70;
end



end