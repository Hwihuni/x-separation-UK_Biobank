function [x_pos, x_neg, x_tot] = x_sep_l1norm(local_f, r2prime, mag, N_std, mask_qsm, params,...
    mask_CSF, wG, wG_r2p, x_qsm_init,...
    mask_FastRelax, mask_SlowRelax)
% susceptibility source separation (chi-separation) using L1 norm and CSF regularization
%  [x_pos, x_neg, x_tot] = x_sep_l1norm(local_f, r2prime, mag, N_std, mask_qsm, params,...
%    mask_CSF, wG, wG_r2p, x_qsm_init,...
%    mask_FastRelax, mask_SlowRelax)
%
% %%%%%%%%%%% Output %%%%%%%%%%%
%   -x_pos: positive susceptibility [ppm]
%   -x_neg: negative susceptibility [ppm]
%   -x_tot: bulk susceptibility of voxel (i.e., x_pos-x_neg) [ppm]
%
% %%%%%%%%%%% Input %%%%%%%%%%%
%   -local_f       : background-removed local frequency map (y, x, z) [Hz unit]
%   -r2prime       : R2' (= R2*-R2) (y, x, z) [Hz unit]
%   -Mag           : Magnitude imaged will be used to estimate wG
%   -N_std         : estimated noise standard deviation
%   -mask_qsm      : binary mask for x-separation (e.g., brain mask)
%   -mask_CSF      : binary mask of CSF regions for CSF regularization
%   -params        : set of parameters
%                   % essenstial params
%                  -> voxel_size : (y,x,z) [mm unit]
%                  -> lambda     : (regularization coefficient for l1 norm)
%                  -> lambda_CSF : (regularization coefficient for CSF)
%                  -> CF         : Central frequency [Hz unit]
%                  -> b0_dir     : (y,x,z)
%                   % optional params
%                  -> Dr_pos            : default = 137 for in-vivo human brain at 3T [Hz/ppm]
%                  -> Dr_neg            : default = Dr_pos for in-vivo human brain at 3T [Hz/ppm]
%                  -> r2p_highThres     : default = 30 
%                  -> r2p_lowThres      : default = 1 
%                  -> time_limit        : default = 1800 
%                  -> init_option       : default = 'TKD'  
%                  -> pad_size          : default = 5 
%                  -> tol_norm_ratio    : default = 0.03
%                  -> max_iter          : default = 30 
%                  -> cg_max_iter       : default = 100 
%                  -> cg_tol            : default = 0.05
%                  -> details : see "fill_params" function in this .m file
%   -wG            : (optional) [y x z 3] weighting matrix estimated from
%                  magnitude images
%   -wG_r2p        : (optional) [y x z 3] weighting matrix estimated from
%                  r2p maps
%   -x_qsm_init    : (optional) initial QSM map for fast convergence
%   -mask_FastRelax: (optional) mask for unreliable voxel with fast relaxation 
%                   e.g)r2p>30
%   -mask_SlowRelax: (optional) mask for unreliable voxel with slow relaxation 
%                   e.g)r2p<1
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

if nargin <7
    wG = []; % empty matrix to automatically generate weighting matrix
    wG_r2p = []; % empty matrix to automatically generate weighting matrix
    x_qsm_init = []; % empty matrix for not using initial point
    mask_CSF = []; % empty matrix if CSF regularization is not used
    
    % please check masks below are reliable
    mask_FastRelax = zeros(size(r2prime)); %or r2prime >params.r2p_highThres; % mask for unreliable voxel with fast relaxation
    mask_SlowRelax = zeros(size(r2prime)); %or r2prime<params.r2p_lowThres; % mask for unreliable voxel with slow relaxation
elseif nargin < 11
    % please check masks below are reliable
    mask_FastRelax = zeros(size(r2prime)); %or r2prime>params.r2p_highThres; % mask for unreliable voxel with fast relaxation
    mask_SlowRelax = zeros(size(r2prime)); %or r2prime<params.r2p_lowThres; % mask for unreliable voxel with slow relaxation
end


disp('+++++++ Starting l1 norm x-sep  +++++++')
tic
tic
[x_pos, x_neg, x_tot]  = x_sep_l1norm_func(local_f, r2prime,...
    mag, N_std, mask_qsm, mask_CSF,...
    wG, wG_r2p, mask_FastRelax,mask_SlowRelax, x_qsm_init, ...
    params);
toc
disp('+++++++++++++++ Done ++++++++++++++++++')

end

%%

function params_filled = fill_params(params)
% fill empty fields in params
gyro_ = 42.58e6;
params_filled = params;

% field which must be filled by manually
%Unit:                 mm          a.u.      a.u.      Hz     a.u
field_essential = {'voxel_size','lambda','lambda_CSF','CF','b0_dir'};

for ii = 1:length(field_essential)
    if ~isfield(params_filled, field_essential{ii})
        error(['params.',field_essential{ii},' is missed']);
    end
end


% optional parameters 
field_scale = params.CF/(gyro_*2.9); %relative to 2.9T
if ~isfield(params_filled, 'Dr_pos')
    params_filled.Dr_pos = 137*field_scale; % 137 Hz/ppm @ 3T for in-vivo human brain
end
if ~isfield(params_filled, 'Dr_neg')
    params_filled.Dr_neg = params_filled.Dr_pos; % default: copy Dr_pos
end
if ~isfield(params_filled, 'r2p_highThres') % less weights on unreliable relaxation estimation
    params_filled.r2p_highThres = 30*field_scale; % 30 Hz at 3T
end
if ~isfield(params_filled, 'r2p_lowThres') % less weights on unreliable relaxation estimation
    params_filled.r2p_lowThres = 1*field_scale; % 1 Hz at 3T
end 
if ~isfield(params_filled, 'time_limit')
    params_filled.time_limit = 1800; % iteration will stop when processing time exceed 1800s
end 
if ~isfield(params_filled, 'init_option')
    params_filled.init_option = 'TKD'; % estimate initial qsm map using TKD QSM algorithm
end 
if ~isfield(params_filled, 'pad_size')% zero padding size
    params_filled.pad_size = 5;
end 

% parameters related to MEDI algorhtim (Liu, MRM, 2018)
if ~isfield(params_filled, 'tol_norm_ratio')
    params_filled.tol_norm_ratio = 0.03; % iteration stops when l2 norm of gradient is smaller than 3% of original map
end
if ~isfield(params_filled, 'max_iter')
    params_filled.max_iter = 30; % maximum number of iteration in x-separation algorithm
end
if ~isfield(params_filled, 'cg_max_iter') % maximum number of iteration for cgsolver
    params_filled.cg_max_iter = 100;
end
if ~isfield(params_filled, 'cg_tol')% tolerance value for cgsolver
    params_filled.cg_tol = 0.05;
end



end