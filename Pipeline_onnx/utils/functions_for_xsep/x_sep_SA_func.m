

function [xp_iLSQR, xn_iLSQR]= x_sep_SA_func(local_f_in_HzUnit,r2p_in_HzUnit, mask_qsm, x_qsm_init, mask_Unrel_Mag, mask_Unrel_Phs,...
    params)
% x_sep_iLSQR_func v1.0.0 : modified from
% source_separation_iLSQR_Ver1_3_params with negative part swap

%% set params

voxel_size = params.voxel_size;
CF = params.CF;
B0_dir = params.b0_dir;
cg_max_iter_LSQR = params.cg_max_iter_LSQR;
cg_tol_LSQR =params.cg_tol_LSQR;
cg_max_iter_SA =params.cg_max_iter_SA;
cg_tol_SA =params.cg_tol_SA;
D2_thr =params.D2_thr;
verbose_option = 0 ;%params.verbose_option;
prc_W_I_min =params.prc_W_I_min;
prc_W_I_max =params.prc_W_I_max;
WG_prc_min =params.WG_prc_min;
WG_prc_max =params.WG_prc_max;
x_thr_strong = params.strong_x_thresh; % ppm unit
pad_size = params.pad_size;

if (isempty(x_qsm_init))&&(params.flag_strong_x_correction ==1)
    disp('Calculating conventional QSM to compensate for strong x...')
    [x_qsm_init] = iLSQR_myVersion_v1_4_(local_f_in_HzUnit, mask_qsm, params);
    disp('%% Done %%%')
    
end
local_f_in_HzUnit = padarray(local_f_in_HzUnit,pad_size);
r2p_in_HzUnit = padarray(r2p_in_HzUnit,pad_size);
mask_qsm = padarray(mask_qsm,pad_size);

if (~isempty(x_qsm_init))
    flag_x_init = 1;
    x_qsm_init = padarray(x_qsm_init,pad_size);
else
    flag_x_init = 0;
    x_qsm_init = zeros(size(mask_qsm));
end
mask_Unrel_Mag = padarray(mask_Unrel_Mag,pad_size);
mask_Unrel_Phs = padarray(mask_Unrel_Phs,pad_size);


wei_for_unrel_Mag = 1/5;
W_unrel_Mag = (wei_for_unrel_Mag*mask_Unrel_Mag+1*~mask_Unrel_Mag).*mask_qsm;

wei_for_unrel_Phs = 0;
W_unrel_Phs = (wei_for_unrel_Phs*mask_Unrel_Phs+1*~mask_Unrel_Phs).*mask_qsm;


wei_for_freq = 2*pi;

loc_f_norm = local_f_in_HzUnit./(CF); % normalized frequency [no unit ]
loc_f_ppm = loc_f_norm*1e6;
matrix_size = size(local_f_in_HzUnit);

r2p_in_HzUnit(r2p_in_HzUnit<0) = 0;
r2p_no_unit = r2p_in_HzUnit.*mask_qsm/CF; %for Dr in Hz/ppm unit
r2p_ppm = r2p_no_unit*1e6;


tol_norm_ratio = 0.02;
max_iter = 1; % one-iteration for LSQR solution




% let assume the Dm kernel is Dirac delta
% Dm_map_4d_no_unit = Dm_map_4d/CF*1e6; %
% [Dr_pos, Dr_neg,use_Dm_offset] = convert_Dm_map(Dm_map_4d_no_unit,matrix_size);
Dr_pos_no_unit = params.Dr_pos/CF*1e6;
Dr_neg_no_unit = params.Dr_neg/CF*1e6;
use_Dm_offset = 0;
if use_Dm_offset == 0
    Dr_offset_pos = 0;
    Dr_offset_neg = 0;
elseif use_Dm_offset ==1
    disp('Dm offset is used!!!');
    disp('Dm offset is correctly calculated..??');
    error('see source_separation_v_NIpaper.m function');
    Dr_offset_pos = Dm_map_offset_4d(1,1,1,1)*ones(matrix_size);
    Dr_offset_neg = Dm_map_offset_4d(1,1,1,2)*ones(matrix_size);
    
end


% inital guess to help fast convergence
if flag_x_init == 1
    x_pos_init = (r2p_in_HzUnit./params.Dr_neg+x_qsm_init)./(1+params.Dr_pos./params.Dr_neg);
    x_neg_init = (r2p_in_HzUnit./params.Dr_pos-x_qsm_init)./(1+params.Dr_neg./params.Dr_pos);
    x_init_cat = cat(2, x_pos_init, x_neg_init);
    x_init_cat(isnan(x_init_cat)) = 0;
    x_init_cat(isinf(x_init_cat)) = 0;
    x_init_cat(x_init_cat<=0) = 0;
elseif flag_x_init ==0
    x_init_cat = cat(2, zeros(size(r2p_in_HzUnit)), zeros(size(r2p_in_HzUnit)));
end

%% fast QSM process
% disp('++++++++++++++++++++++++++++++++++')
% disp('Calculating fast QSM data and gradient weighting...')
disp('Calculating x-separation .')

% see NI, Wei Li, 2015
percentile_1st = 1;
percentile_2nd = 30;
smv_size_for_LPF = 2.5; % mm unit
th_tkd = 1/8;



D_p=dipole_kernel_p(matrix_size, voxel_size, B0_dir);
Dconv_p = @(dx) (ifftn(D_p.*fftn(dx))); %  susceptibility [no unit not ppm] to normalized freq. domain [no_unit]

D_001= abs(D_p).^0.001;
a1 = prctile(D_001(:),percentile_1st);
b1 = prctile(D_001(:),percentile_2nd);
W_0 = (D_001-a1)/(b1-a1);
W_FS = W_0;
W_FS(W_0<0) = 0;
W_FS(W_0>1) = 1;

D_p_tkd = D_p;
mask_t = (D_p>=0).* (D_p<th_tkd);
D_p_tkd(mask_t>0) = th_tkd;
mask_t = (D_p<0).* (D_p>-th_tkd);
D_p_tkd(mask_t>0) = -th_tkd;
x_tkd = mask_qsm.*real(ifftn(D_p_tkd.^-1.*fftn(loc_f_ppm)));



x_f1_k = sign(D_p).*fftn(loc_f_ppm);
x_f1_r = real(ifftn(x_f1_k));

x_f1_kc = ifftshift(x_f1_k);
x_f1_kc_lpf = SMV(x_f1_kc, matrix_size,voxel_size, smv_size_for_LPF);

x_f2_kc = x_f1_kc.*W_FS + x_f1_kc_lpf.*(1-W_FS);
x_f2_k = fftshift(x_f2_kc);
x_f2_r = real(ifftn(x_f2_k));

x_f2_k_mask = fftn(mask_qsm.*x_f2_r);
x_f2_k_mask_lpf = SMV(x_f2_k_mask, matrix_size,voxel_size, smv_size_for_LPF);
x_f3_r = real(ifftn(x_f2_k_mask.*W_FS + x_f2_k_mask_lpf.*(1-W_FS))).*mask_qsm;

[a1, b1, rsq, yfit] = calc_linear_fit(x_f3_r(:),x_tkd(:), 1);

x_fs = a1.*x_f3_r+b1;

%x_f_ = a1*cat(2, x_f1_r,x_f2_r,x_f3_r)+b1;
% figure;imshow3Dfull(cat(2,x_f_,x_tkd),plot_range)
% title('x_f1, x_f2, x_f3, x_tkd')


r2p_ppm_lpf = real(calc_tukeywin(r2p_ppm,0.5,'3D'));
derv_x   = (fgrad_(            x_fs.*(mask_qsm>0), voxel_size));
derv_r2p = (fgrad_( r2p_ppm_lpf.*(mask_qsm>0), voxel_size)); % it can be modified  even in the case Dconv_m is defined by convolution (using deconv.)

%figure;imshow3d(sqz(derv_x(:,:,end/2,:)),[-10 10]*1e-8)
%figure;imshow3d(sqz(derv_r2p(:,:,end/2,:)),[-10 10]*1e-8)


grad_x_p = (derv_r2p + Dr_neg_no_unit.*derv_x)./(Dr_pos_no_unit+Dr_neg_no_unit);
grad_x_n = (derv_r2p - Dr_pos_no_unit.*derv_x)./(Dr_pos_no_unit+Dr_neg_no_unit);

[W_G_p,gmax, gmin] = gradient_mask_from_grad_(grad_x_p, mask_qsm, WG_prc_min, WG_prc_max);
[W_G_n,gmax, gmin] = gradient_mask_from_grad_(grad_x_n, mask_qsm, WG_prc_min, WG_prc_max);
%
% [W_G_x,gmax, gmin] = gradient_mask_from_grad(derv_x, mask_qsm, WG_prc_min, WG_prc_max);
% [W_G_r2p,gmax, gmin] = gradient_mask_from_grad(derv_r2p, mask_qsm, WG_prc_min, WG_prc_max);

figure;imshow3Dfull_toggle(sqrt(sum(abs(W_G_p).^2,4)),sqrt(sum(abs(W_G_n).^2,4)),[0 2]); title('grad weight pos and neg')
%
% figure;imshow3Dfull_toggle(sqrt(sum(abs(W_G_x).^2,4)),sqrt(sum(abs(W_G_r2p).^2,4)),[0 2]);



%% x-sep LSQR solution with matrix concat.

disp('Calculating x-separation ..')


if params.flag_strong_x_correction == 1
    % STAR-QSM using x_iLSQR
    x_strong  = x_qsm_init.*(abs(x_qsm_init)>x_thr_strong);
    f_strong = Dconv_p(x_strong).*mask_qsm;
    loc_f_ppm_weak = loc_f_ppm-f_strong;
    loc_f_ppm_new = loc_f_ppm_weak;
    
elseif params.flag_strong_x_correction == 0
    loc_f_ppm_new = loc_f_ppm;
    x_strong = zeros(size(loc_f_ppm));
end


% y = Ax problem

[W_I,lap_max, lap_min] = laplacian_mask_iLSQR_exact_(loc_f_ppm_new, mask_qsm, voxel_size, prc_W_I_min, prc_W_I_max);
figure;imshow3Dfull(W_I.*mask_qsm)

w_phs = wei_for_freq.*W_unrel_Phs.*W_I.*mask_qsm;
w_mag =               W_unrel_Mag.*W_I.*mask_qsm;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Kernel for Phase %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D_p=dipole_kernel_p(matrix_size, voxel_size, B0_dir);
Dconv_p =   @(dx) real(ifftn(D_p.*fftn(dx))); %  susceptibility [no unit not ppm] to normalized freq. domain [no_unit]
Af   =  @(xx) (Dconv_p(w_phs.*(Dconv_p(xx))));      %@(xx) (w_phs.*(Dconv_p(xx)));     %
Af_H =  @(xx) (Dconv_p(conj(w_phs).*(Dconv_p(xx))));%@(xx) (Dconv_p(conj(w_phs).*xx)); %
yf   =  Dconv_p(w_phs.*loc_f_ppm_new);  %real(w_phs.*loc_f_norm);%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Kernel for Magnitude %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dconv_m = @(x, Dr, offset) real(Dr.*x + offset);
A_r =   @(xx, Dr, offset) Dconv_p(w_mag.* real(Dconv_m(xx, Dr, offset)));          % w_mag.* real(Dconv_m(xx, Dr, offset));
A_r_H = @(xx, Dr, offset) real(Dconv_m( (conj(w_mag).* Dconv_p(xx)), Dr, offset)); % real(Dconv_m( (conj(w_mag).* xx), Dr, offset));
yr = Dconv_p(real(w_mag.*r2p_ppm)); % real(w_mag.*r2p_no_unit);

Arp = @ (xx) A_r(xx, Dr_pos_no_unit, Dr_offset_pos);
Arp_H = @ (xx) A_r_H(xx, Dr_pos_no_unit, Dr_offset_pos);
Arn = @ (xx) A_r(xx, Dr_neg_no_unit, Dr_offset_neg);
Arn_H = @ (xx) A_r_H(xx, Dr_neg_no_unit, Dr_offset_neg);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% derivative functions %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fidelity_p_deri_PandN =        @(xp, xn)        (Arp_H(Arp(xp)) + Af_H(Af(xp)) + Arp_H(Arn(xn)) - Af_H(Af(xn)));
fidelity_n_deri_PandN =        @(xp, xn)        (Arn_H(Arp(xp)) - Af_H(Af(xp)) + Arn_H(Arn(xn)) + Af_H(Af(xn)));
A_pn = @(x_cat) cat(2, fidelity_p_deri_PandN(x_cat(:,1:end/2,:),x_cat(:,end/2+1:end,:)),...
    fidelity_n_deri_PandN(x_cat(:,1:end/2,:),x_cat(:,end/2+1:end,:)));

derivative_remain_p_deri_PandN = @(xp, xn)  ( - Arp_H(yr) - Af_H(yf));
derivative_remain_n_deri_PandN = @(xp, xn)  ( - Arn_H(yr) + Af_H(yf));
derivative_remain_pn = @(x_cat) cat(2,derivative_remain_p_deri_PandN(x_cat(:,1:end/2,:),x_cat(:,end/2+1:end,:)),...
    derivative_remain_n_deri_PandN(x_cat(:,1:end/2,:),x_cat(:,end/2+1:end,:)) );


res_norm_ratio=inf;
res_norm_ratio_min = res_norm_ratio;
iter = 0;

x_p_iter = zeros(size(loc_f_ppm_new));%x_tkd.*(x_tkd>0);%x_fs;%zeros(size(x_tkd));
x_n_iter = zeros(size(loc_f_ppm_new));%x_tkd.*(x_tkd<0);%x_fs;%zeros(size(x_tkd));
x_n_iter = -x_n_iter;
x_cat_iter= cat(2,x_p_iter, x_n_iter);
matrix_size_cat = size(x_cat_iter);
mask_qsm_erode2 = imerode( mask_qsm, strel('sphere',2));

while (res_norm_ratio>tol_norm_ratio)&&((iter<max_iter)) % only one iterationa
    iter=iter+1;
    
    
    
    afunc_pn = @(dx_cat) reshape(A_pn(reshape(dx_cat,matrix_size_cat)),[prod(matrix_size_cat),1]);
    b_pn_cal = @(x_cat_k) derivative_remain_pn(x_cat_k);
    b_pn = b_pn_cal(x_cat_iter); % actually b_pn is independent on x_cat_iter
    [dx_cat_, ~] = pcg(afunc_pn, -b_pn(:), cg_tol_LSQR, cg_max_iter_LSQR,[],[], x_init_cat(:));
    dx_cat_ = reshape(real(dx_cat_),matrix_size_cat);
    
    x_cat_iter_prev = x_cat_iter;
    
    x_cat_iter = x_cat_iter + dx_cat_.*cat(2,mask_qsm, mask_qsm);
    
    
    x_p_iter = x_cat_iter(:,1:end/2,:);
    x_n_iter = x_cat_iter(:,end/2+1:end,:);
    
    %         if ~isempty(find(x_p_iter(:)<0))||~isempty(find(x_n_iter(:)<0))
    %             error('deal with negative value') % e.g., d_x_neg = -x_p.*(x_p<0);
    %         end
    
    if verbose_option > 0
        res_phs=w_phs.*(real(Dconv_p(x_p_iter-x_n_iter)) - loc_f_ppm_new);
        cost_pha =norm(res_phs(:),2);
        
        res_mag=w_mag.*(real(Dconv_m(x_p_iter,Dr_pos_no_unit,Dr_offset_pos)+Dconv_m(x_n_iter,Dr_neg_no_unit,Dr_offset_neg)) - r2p_ppm);
        cost_mag =norm(res_mag(:),2);
        
        %cost_phase_fidelity(iter) = norm(res_phs(:),2);
        fprintf('Data fidelity: Mag = %4.4f  /  Phase = %4.4f \n', cost_mag, cost_pha);
        
    end
    
    
end
%
% figure;imshow3Dfull(x_cat_iter,[-0.1 0.1]);colormap jet;
% figure;imshow3Dfull(x_cat_iter,[0 0.1]);


x_p_iter(isnan(x_p_iter)) = 0;
x_p_iter(isinf(x_p_iter)) = 0;
x_n_iter(isnan(x_n_iter)) = 0;
x_n_iter(isinf(x_n_iter)) = 0;


xp_LSQR = x_p_iter; % ppm
xn_LSQR = x_n_iter; % ppm
% figure;imshow3Dfull_toggle(xp_LSQR,xn_LSQR,[0 0.1]);

%fsleyes_mat(xp_LSQR,xn_LSQR);





%% estimation of streaking artifacts
disp('Calculating x-separation ...')


grad = @fgrad_; % discrete gradient
div = @bdiv_; % discrete divergence
D_p=dipole_kernel_p(matrix_size, voxel_size, B0_dir);
M_IC = abs(D_p)<D2_thr;
M_IC = double(M_IC);




res_norm_ratio=inf;
res_norm_ratio_min = res_norm_ratio;
W_G_cat = cat(5,W_G_p,W_G_n);
x_LSQR_cat = cat(4,xp_LSQR,xn_LSQR);
for iter_pn = 1:2
    if(iter_pn) == 1
        disp('Calculating x-separation ....')
    elseif (iter_pn) == 2
        disp('Calculating x-separation .....')
    end
    W_G_t = W_G_cat(:,:,:,:,iter_pn);
    x_0 = x_LSQR_cat(:,:,:,iter_pn);
    % init
    x_SA_iter_k = zeros(size(x_0));
    
    A_SA_           = @(xx_k)    (  conj(M_IC).*fftn(div(conj(W_G_t).*W_G_t.*grad(real(ifftn(M_IC.*xx_k))))));
    b_cal_SA = @(x_k, x_init)    (- conj(M_IC).*fftn(div(conj(W_G_t).*W_G_t.*grad(real(x_init)))));% +fidelity(x_k) )
    
    afunc_SA_ = @(dx) reshape(A_SA_(reshape(dx, matrix_size)),[prod(matrix_size),1]);
    b_SA = b_cal_SA(x_SA_iter_k, x_0);
    
    [x_SA_iter_k,~ ]= pcg(afunc_SA_, -b_SA(:), cg_tol_SA, cg_max_iter_SA);
    x_SA_iter_k = reshape(x_SA_iter_k,matrix_size);
    
    if(iter_pn) == 1
        xp_SA_indiv = real(ifftn(x_SA_iter_k.*M_IC)).*mask_qsm; % ppm
    elseif (iter_pn) == 2
        xn_SA_indiv = real(ifftn(x_SA_iter_k.*M_IC)).*mask_qsm; % ppm
    end
end


xp_iLSQR_indiv = (xp_LSQR -xp_SA_indiv).*mask_qsm;
xn_iLSQR_indiv = (xn_LSQR -xn_SA_indiv).*mask_qsm;

%%


xp_iLSQR = xp_iLSQR_indiv.*mask_qsm;
xn_iLSQR = xn_iLSQR_indiv.*mask_qsm;

xp_iLSQR(xp_iLSQR<0) = 0;
xn_iLSQR(xn_iLSQR<0) = 0;

xp_iLSQR = xp_iLSQR + x_strong.*(x_strong>0);
xn_iLSQR = xn_iLSQR - x_strong.*(x_strong<0);

% swap negative part
xp_iLSQR_n = xp_iLSQR.*(xp_iLSQR<0);
xn_iLSQR_n = xn_iLSQR.*(xn_iLSQR<0);
xp_iLSQR_swap = xp_iLSQR -xn_iLSQR_n;
xn_iLSQR_swap = xn_iLSQR -xp_iLSQR_n;


xp_iLSQR = xp_iLSQR_swap(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
xn_iLSQR = xn_iLSQR_swap(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
% xp_LSQR = xp_LSQR(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
% xn_LSQR = xn_LSQR(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));


disp('++++++++++++Done++++++++++++++++++')


end
%%


function img_filtered = calc_tukeywin(img, t_size, dim)



if strcmp(dim, '2D')
    rawd = ifft3c_(img, 3);
else
    rawd = ifft3c_(img, 11);
end

[ynt xnt znt cnt ent] = size(rawd);
hfilt2=tukeywin(ynt,t_size)*tukeywin(xnt,t_size)';
hfilt3 = zeros(ynt, xnt, znt);
if strcmp(dim, '2D')
    %disp('2D GRE DATA')
    hfilt_rep = repmat(hfilt2,[1 1 znt cnt]);
elseif strcmp(dim,'3D')
    %disp('3D GRE DATA')
    for yn=1:ynt
        for xn=1:xnt
            hfilt3(yn,xn,:)=hfilt2(yn,xn)*tukeywin(znt,t_size);
        end
    end
    hfilt_rep = repmat(hfilt3,[1 1 1 cnt]);
end

for en = 1:size(rawd,5)
    rawd(:,:,:,:,en) = rawd(:,:,:,:,en).*hfilt_rep;
end

if strcmp(dim, '2D')
    img_filtered=real(fft3c_(rawd,3)); % = fft2c(rawd)
else
    img_filtered = real(fft3c_(rawd,11));
end



end


%%

function [wI, G_max, G_min]=laplacian_mask_iLSQR_exact_(phs_, Mask, voxel_size, percentage_min, percentage_max)


wI = abs(del2(phs_.*(Mask>0),voxel_size(1), voxel_size(2), voxel_size(3)));

G_max = prctile(wI(Mask(:)>0),percentage_max);
G_min = prctile(wI(Mask(:)>0),percentage_min);
wG_t = (G_max-wI)./(G_max-G_min);
wG_t(wI<G_min) = 1;
wG_t(wI>G_max) = 0;% = wG;
wI = wG_t;

end


function [wG, G_max, G_min]=gradient_mask_from_grad_(derivated_value, Mask, percentage_min, percentage_max)


wG = abs(derivated_value);

wG(isnan(wG))= 0;
for iter_dim = 1:size(wG,4)
    temp =wG(:,:,:,iter_dim);
    G_max = prctile(temp(Mask(:)>0),percentage_max);
    G_min =  prctile(temp(Mask(:)>0),percentage_min);
    wG_t = (G_max-temp)./(G_max-G_min);
    wG_t(temp<G_min) = 1;
    wG_t(temp>G_max) = 0;% = wG;
    wG(:,:,:,iter_dim) = wG_t.*Mask;
end

end

function D=dipole_kernel_p(matrix_size, voxel_size, B0_dir)
domain = 'kspace';
if (B0_dir == 1)
    B0_dir = [1 0 0 ]';
elseif (B0_dir == 2)
    B0_dir = [0 1 0 ]';
elseif (B0_dir==3)
    B0_dir = [0 0 1]';
end

if strcmp(domain,'kspace')
    [Y,X,Z]=meshgrid(-matrix_size(2)/2:(matrix_size(2)/2-1),...
        -matrix_size(1)/2:(matrix_size(1)/2-1),...
        -matrix_size(3)/2:(matrix_size(3)/2-1));
    
    X = X/(matrix_size(1)*voxel_size(1));
    Y = Y/(matrix_size(2)*voxel_size(2));
    Z = Z/(matrix_size(3)*voxel_size(3));
    
    D = 1/3-  ( X*B0_dir(1) + Y*B0_dir(2) + Z*B0_dir(3) ).^2./(X.^2+Y.^2+Z.^2);
    D = fftshift(D);
    %     e=0.01;
    %     D(abs(D)<e)=e.*D(abs(D)<e)./abs(D(abs(D)<e));
    D(isnan(D)) = 0;
    
elseif strcmp(domain,'imagespace')
    [Y,X,Z]=meshgrid(-matrix_size(2)/2:(matrix_size(2)/2-1),...
        -matrix_size(1)/2:(matrix_size(1)/2-1),...
        -matrix_size(3)/2:(matrix_size(3)/2-1));
    
    X = X*voxel_size(1);
    Y = Y*voxel_size(2);
    Z = Z*voxel_size(3);
    
    D = (3*( X*B0_dir(1) + Y*B0_dir(2) + Z*B0_dir(3)).^2 - X.^2-Y.^2-Z.^2)./(4*pi*(X.^2+Y.^2+Z.^2).^2.5);
    
    D(isnan(D)) = 0;
    
    %     D = fftn(fftshift(d));
end

end

function Gx = fgrad_(chi, voxel_size)

if (nargin < 2)
    voxel_size = [1 1 1];
end

% chi = double(chi);

Dx = [chi(2:end,:,:); chi(end,:,:)] - chi;
Dy = [chi(:,2:end,:), chi(:,end,:)] - chi;
Dz = cat(3, chi(:,:,2:end), chi(:,:,end)) - chi;

Dx = Dx/voxel_size(1);
Dy = Dy/voxel_size(2);
Dz = Dz/voxel_size(3);

Gx = cat(4, Dx, Dy, Dz);

end





function div = bdiv_(Gx, voxel_size)

if (nargin < 2)
    voxel_size = [1 1 1];
end

% Gx = double(Gx);

Gx_x = Gx(:,:,:,1);
Gx_y = Gx(:,:,:,2);
Gx_z = Gx(:,:,:,3);

[Mx, My, Mz] = size(Gx_x);

Dx = [Gx_x(1:end-1,:,:); zeros(1,My,Mz)]...
    - [zeros(1,My,Mz); Gx_x(1:end-1,:,:)];

Dy = [Gx_y(:,1:end-1,:), zeros(Mx,1,Mz)]...
    - [zeros(Mx,1,Mz), Gx_y(:,1:end-1,:)];

Dz = cat(3, Gx_z(:,:,1:end-1), zeros(Mx,My,1))...
    - cat(3, zeros(Mx,My,1), Gx_z(:,:,1:end-1));

Dx = Dx/voxel_size(1);
Dy = Dy/voxel_size(2);
Dz = Dz/voxel_size(3);

div = -( Dx + Dy + Dz );

end




function [x_iLSQR, x_tkd, x_fs, x_LSQR] = iLSQR_myVersion_v1_4_(local_f_in_HzUnit, mask_qsm, params)
%
% function [x_iLSQR, x_tkd, x_fs, x_LSQR] = iLSQR_myVersion_v1_4(local_f_in_HzUnit, mask_qsm, params)
%
% input
%   loc_f_in_HzUnit: background-removed frequency map [Hz unit]
%   mask_qsm       : binary mask
%   voxel_size     : [y x z] in mm unit
%   CF             : central freq. [Hz]
%   B0_dir         : [y x z]
%   cg_max_iter_LSQR, cg_tol_LSQR, cg_max_iter_SA, cg_tol_SA, D2_thr_SA:
%   params for iLSQR (see literature)
%
% output
%   x_iLSQR        : iLSQR solution [ppm unit]
%   x_tkd          : TKD solution with threshould 1/8 [ppm unit]
%   x_fs           : fast QSM solution for edge mask [ppm unit]
%   x_LSQR         : LSQR solution with streaking artifact [ppm unit]
%
% Caution: Li's parameters are optimzed for COSMOS, which means the
% parameter can induce underestimaton of x due to anisotropy effects
%
% reference: Wei Li, NI, 2015
% generated by Hyeong-Geol Shin 2021.03.29 (sin4109@gmail.com)
%
% v1.0: use pcg instead of custom cgsolver (better than v1_1 (!) and v1_2 (?) )
% v1.4: Use params for easy parameter inputs (developed from v1.0)


%%
% if nargin <=5
%     % for step 1
%     cg_max_iter_LSQR = 5000;
%     cg_tol_LSQR = 0.01;% smaller means more artifacts but accurate, 0.02 is literature value (but can induce underestimation) (recommadation for quant. = 0.01)
%
%     % for step 2
%     cg_max_iter_SA = 1000;
%     cg_tol_SA = 0.02;
%     % large value make large SA but it can contain structure.
%     D2_thr_SA = 0.07; % smaller reduce smaller artifacts but less remove structure (literature value = 0.1, recommadation for quant. = 0.07)
%
%     % percentile for LSQR weighting
%     prc_W_I_min = 60;% see Wei LI, NI, 2015
%     prc_W_I_max = 99.8;
%     % for step 2
%     % percentil for edge weighting
%     WG_prc_min =50;
%     WG_prc_max =70;
% elseif nargin <=10
%     prc_W_I_min = 60;% see Wei LI, NI, 2015
%     prc_W_I_max = 99.8;
%     WG_prc_min =50;
%     WG_prc_max =70;
%
% end



voxel_size = params.voxel_size;
CF = params.CF;
B0_dir = params.b0_dir;
cg_max_iter_LSQR = params.cg_max_iter_LSQR;
cg_tol_LSQR =params.cg_tol_LSQR;
cg_max_iter_SA =params.cg_max_iter_SA;
cg_tol_SA =params.cg_tol_SA;
D2_thr_SA =params.D2_thr;
prc_W_I_min =params.prc_W_I_min;
prc_W_I_max =params.prc_W_I_max;
WG_prc_min =params.WG_prc_min;
WG_prc_max =params.WG_prc_max;





% for step 1
max_iter_LSQR = 1; % 1 for direct calculation using cgsolver
tol_norm_ratio = 0.02;
D2_thr_SA_high = 0.4; %for strong susce. distribution for STAR-QSM
x_thr_strong = 0.2; % ppm unit

div = @bdiv;

pad_size = [16 16 16];
local_f_in_HzUnit = padarray(local_f_in_HzUnit,pad_size);
mask_qsm = padarray(mask_qsm,pad_size);


loc_f_norm = local_f_in_HzUnit/(CF); % normalized frequency [no unit ]
loc_f_ppm = loc_f_norm*1e6;
matrix_size = size(local_f_in_HzUnit);



if max_iter_LSQR > 1
    error('use a single step solution')
    cg_max_iter_LSQR = 200;
    cg_tol_LSQR = 0.02;% 0.1
    tol_norm_ratio = 0.04;
end
%% fast QSM process
% disp('++++++++++++++++++++++++++++++++++')
% disp('Calculating fast QSM data...')

% see NI, Wei Li, 2015
percentile_1st = 1;
percentile_2nd = 30;
smv_size_for_LPF = 2.5; % mm unit
th_tkd = 1/8;



D_p=dipole_kernel_p(matrix_size, voxel_size, B0_dir);
Dconv_p = @(dx) (ifftn(D_p.*fftn(dx))); %  susceptibility [no unit not ppm] to normalized freq. domain [no_unit]

D_001= abs(D_p).^0.001;
a = prctile(D_001(:),percentile_1st);
b = prctile(D_001(:),percentile_2nd);
W_0 = (D_001-a)/(b-a);
W_FS = W_0;
W_FS(W_0<0) = 0;
W_FS(W_0>1) = 1;

D_p_tkd = D_p;
mask_t = (D_p>=0).* (D_p<th_tkd);
D_p_tkd(mask_t>0) = th_tkd;
mask_t = (D_p<0).* (D_p>-th_tkd);
D_p_tkd(mask_t>0) = -th_tkd;
x_tkd = mask_qsm.*real(ifftn(D_p_tkd.^-1.*fftn(loc_f_ppm)));



x_f1_k = sign(D_p).*fftn(loc_f_ppm);
x_f1_r = real(ifftn(x_f1_k));

x_f1_kc = ifftshift(x_f1_k);
x_f1_kc_lpf = SMV(x_f1_kc, matrix_size,voxel_size, smv_size_for_LPF);

x_f2_kc = x_f1_kc.*W_FS + x_f1_kc_lpf.*(1-W_FS);
x_f2_k = fftshift(x_f2_kc);
x_f2_r = real(ifftn(x_f2_k));

x_f2_k_mask = fftn(mask_qsm.*x_f2_r);
x_f2_k_mask_lpf = SMV(x_f2_k_mask, matrix_size,voxel_size, smv_size_for_LPF);
x_f3_r = real(ifftn(x_f2_k_mask.*W_FS + x_f2_k_mask_lpf.*(1-W_FS))).*mask_qsm;

[a1, b1, rsq, yfit] = calc_linear_fit(x_f3_r(:),x_tkd(:), 1);

x_fs = a1.*x_f3_r+b1;

%x_f_ = a1*cat(2, x_f1_r,x_f2_r,x_f3_r)+b1;
% figure;imshow3Dfull(cat(2,x_f_,x_tkd),plot_range)
% title('x_f1, x_f2, x_f3, x_tkd')


%% LSQR solution

% disp('Calculation...')
%

% y = Ax problem
% test prc_W_I_min = 70;prc_W_I_max = 95;
[W_I,lap_max, lap_min] = laplacian_mask_iLSQR_exact(loc_f_ppm, mask_qsm, voxel_size, prc_W_I_min, prc_W_I_max);
%figure;imshow3Dfull(W_I.*mask_qsm,[0 1]); title('W_I')

w_phs = W_I.*mask_qsm;
w_phs(isnan(w_phs))=0;
w_phs(isinf(w_phs)) = 0;
D_p=dipole_kernel_p(matrix_size, voxel_size, B0_dir);
Dconv_p = @(dx) real(ifftn(D_p.*fftn(dx))); %  susceptibility [no unit not ppm] to normalized freq. domain [no_unit]


res_norm_ratio=inf;
res_norm_ratio_min = res_norm_ratio;
iter = 0;

% init
x_iter = zeros(size(loc_f_ppm));%x_fs;%zeros(size(x_tkd));
mask_qsm_erode2 = imerode( mask_qsm, strel('sphere',2));


while (res_norm_ratio>tol_norm_ratio)&&((iter<max_iter_LSQR))
    iter = iter+1;
    precond = @(xx)   Dconv_p(xx);
    precond_H = @(xx) Dconv_p(xx);
    
    A = @(xx)  Dconv_p(real(conj(w_phs).*(precond_H(precond(w_phs.*(Dconv_p(xx)))))));
    afunc_LSQR = @(dx) reshape(A(reshape(dx,size(w_phs))),[numel(w_phs),1 ]);
    % direct inversion problem
    b_cal = @(x_k, target) ( - Dconv_p(real(conj(w_phs).*(precond_H(precond(w_phs.*target))))));%+fidelity(x_k));
    b = b_cal(x_iter, loc_f_ppm);
    [x_iter,~] = pcg(afunc_LSQR, -b(:), cg_tol_LSQR, cg_max_iter_LSQR);
    x_iter = reshape(real(x_iter),size(w_phs)).*mask_qsm;
    %figure;imshow3Dfull(x_iter,[-0.1 0.1])
    
    %res_norm_ratio = norm(dx(mask_qsm_erode2>0))/norm(x_iter(mask_qsm_erode2>0)); % gauss norm at stric mask due to the suscepbilityt errror around boundary
    
    %disp('+++++++++++++++++++++++');
    %fprintf('iter: %d; res_norm_ratio:%8.4f\n', iter,  res_norm_ratio);
    
    res_phs=w_phs.*(real(Dconv_p(x_iter)) - loc_f_ppm);
    cost_pha =norm(res_phs(:),2);
    
    %cost_phase_fidelity(iter) = norm(res_phs(:),2);
    fprintf('Data fidelity:   Phase = %4.4f \n', cost_pha);
    if 0%res_norm_ratio_min>res_norm_ratio
        res_norm_ratio_min=res_norm_ratio;
        x_min_res_norm = x_iter;
        
    end
    
end
%x_iter = x_min_res_norm;
x_iter(isnan(x_iter)) = 0;
x_iter(isinf(x_iter)) = 0;
x_LSQR = x_iter; % ppm

%figure;imshow3Dfull(x_LSQR,[-0.1 0.1])
%figure;imshow3Dfull(cat(2,x_medi+0.01,x_fs,x_LSQR),[-0.1 0.1])

%% Estimation of streaking artifacts
% disp('Step 2 (Estimation of streaking artifacts)...')


%gradient weight calculation
grad = @fgrad;% discrete gradient
% test => WG_prc_min= 60; WG_prc_max = 90;
[W_G,gmax, gmin] = gradient_mask_iLSQR_exact(x_fs, mask_qsm, grad, voxel_size, WG_prc_min, WG_prc_max);; title('W_G')
%figure;imshow3Dfull(cat(1,W_G(:,:,:,1),W_G(:,:,:,2),W_G(:,:,:,3)));


D_p=dipole_kernel_p(matrix_size, voxel_size, B0_dir);
M_IC = abs(D_p)<D2_thr_SA;
M_IC = double(M_IC);
M_IC_hth = double(abs(D_p)<D2_thr_SA_high);




res_norm_ratio=inf;
res_norm_ratio_min = res_norm_ratio;
iter = 0;

% init
x_SA_iter_k = zeros(size(x_LSQR));
mask_qsm_erode2 = imerode( mask_qsm, strel('sphere',2));

iter=iter+1;
%error('check location of real function');

sc_fft = numel(w_phs);
A_SA =     @(xx_k)         (  conj(M_IC).*fftn(div(conj(W_G).*W_G.*grad(real(ifftn(M_IC.*xx_k))))));
b_SA_cal = @(x_k, x_init)  (- conj(M_IC).*fftn(div(conj(W_G).*W_G.*grad(real(x_init)))));% +fidelity(x_k) )
b_SA_init = b_SA_cal(x_SA_iter_k*0, x_LSQR);
afunc_SA = @(xx_k)reshape(A_SA( reshape(xx_k,matrix_size)),[numel(w_phs),1]);

[dx_SA_smth, ~] = pcg(afunc_SA, -b_SA_init(:), cg_tol_SA, cg_max_iter_SA);
dx_SA_smth = reshape(dx_SA_smth,size(w_phs));
x_SA_sum_smth = real(ifftn(dx_SA_smth.*M_IC_hth)).*mask_qsm; % ppm


x_iLSQR_smth = x_LSQR-real(x_SA_sum_smth);
%figure;imshow3Dfull_toggle(x_iLSQR_smth,abs(x_iLSQR_smth)>x_thr_strong,[-0.1 0.1])
x_strong  = x_iLSQR_smth.*(abs(x_iLSQR_smth)>x_thr_strong);
f_strong = Dconv_p(x_strong).*mask_qsm;
loc_f_weak_ppm = loc_f_ppm-f_strong;


b_weak = double(b_cal(x_iter, loc_f_weak_ppm));
[x_LSQR_weak,~] = pcg(afunc_LSQR, -b_weak(:), cg_tol_LSQR, cg_max_iter_LSQR);
x_LSQR_weak = reshape(real(x_LSQR_weak),matrix_size).*mask_qsm;
x_LSQR_weak(isnan(x_LSQR_weak)) = 0;
x_LSQR_weak(isinf(x_LSQR_weak)) = 0;
%figure;imshow3Dfull(x_LSQR_weak,[-0.1 0.1])


%% Estimation of streaking artifacts

A_SA =     @(xx_k)         (  conj(M_IC).*fftn(div(conj(W_G).*W_G.*grad(real(ifftn(M_IC.*xx_k))))));
b_SA_cal = @(x_k, x_init)  (- conj(M_IC).*fftn(div(conj(W_G).*W_G.*grad(real(x_init)))));% +fidelity(x_k) )

afunc_SA = @(xx_k)reshape(A_SA( reshape(xx_k,matrix_size)),[numel(w_phs),1]);
b_SA = b_SA_cal(x_SA_iter_k, x_LSQR_weak);
[dx_SA, ~] = pcg(afunc_SA, -b_SA(:), cg_tol_SA, cg_max_iter_SA);
dx_SA = reshape(dx_SA,size(w_phs));
x_SA = real(ifftn(dx_SA.*M_IC)).*mask_qsm; % ppm
%figure;imshow3Dfull_toggle(x_LSQR_weak-real(x_SA),x_LSQR_weak,[-0.1 0.1])
%figure;imshow3Dfull(x_SA,[-0.1 0.1])




%% iLSQR solution

x_iLSQR = (x_LSQR -x_SA).*mask_qsm;

x_iLSQR = x_iLSQR(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
x_tkd = x_tkd(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
x_fs = x_fs(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
x_LSQR = x_LSQR(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
% disp('+++++++++++++++Done+++++++++++++++')


end



function im = ifft3c_(d,option)
% USAGE : im = ifft3c_(d,option)
%
% ifft3c_ performs a centered ifft3
%
% option :
%     0 -> all dir
%     1 -> y dir
%     2 -> x dir
%     8 -> z dir
%     3 -> y,x dir
%     9 -> y,z dir
%     10 -> x,z dir
%     11 -> x,y,z dir
%     
% coded by Sang-Young Zho
% last modified at 2009.05.27

if nargin==1
    option = 11;
end

switch option
    case 0 %     0 -> all dir
        im = ifftshift(ifftn(fftshift(d)));
        
    case 1 %     1 -> y dir
        im = ifftshift(ifft(fftshift(d),[],1));
    case 2 %     2 -> x dir
        im = ifftshift(ifft(fftshift(d),[],2));
    case 8 %     8 -> z dir
        im = ifftshift(ifft(fftshift(d),[],3));
        
    case 3 %     3 -> y,x dir
        im = fftshift(d);
        clear d;
        im = ifft(im,[],1);
        im = ifft(im,[],2);
        im = ifftshift(im);
    case 9 %     9 -> y,z dir
        im = fftshift(d);
        clear d;
        im = ifft(im,[],1);
        im = ifft(im,[],3);
        im = ifftshift(im);
    case 10 %     10 -> x,z dir
        im = fftshift(d);
        clear d;
        im = ifft(im,[],2);
        im = ifft(im,[],3);
        im = ifftshift(im);
        
    case 11 %     11 -> x,y,z dir
        im = fftshift(d);
        clear d;
        im = ifft(im,[],1);
        im = ifft(im,[],2);
        im = ifft(im,[],3);
        im = ifftshift(im);

    otherwise
        error();
        %disp('Error using "ifft3c_"...')
        %disp('Invalid option.')
        im = [];
        return;
end
end

function im = fft3c_(d,option)
% USAGE : im = fft3c(d,option)
%
% fft3c performs a centered fft3
%
% option :
%     0 -> all dir
%     1 -> y dir
%     2 -> x dir
%     8 -> z dir
%     3 -> y,x dir
%     9 -> y,z dir
%     10 -> x,z dir
%     11 -> x,y,z dir
%     
% coded by Sang-Young Zho
% last modified at 2009.05.27

if nargin==1
    option = 11;
end

switch option
    case 0 %     0 -> all dir
        im = ifftshift(fftn(fftshift(d)));
        
    case 1 %     1 -> y dir
        im = ifftshift(fft(fftshift(d),[],1));
    case 2 %     2 -> x dir
        im = ifftshift(fft(fftshift(d),[],2));
    case 8 %     8 -> z dir
        im = ifftshift(fft(fftshift(d),[],3));
        
    case 3 %     3 -> y,x dir
        im = fftshift(d);
        clear d;
        im = fft(im,[],1);
        im = fft(im,[],2);
        im = ifftshift(im);        
    case 9 %     9 -> y,z dir
        im = fftshift(d);
        clear d;
        im = fft(im,[],1);
        im = fft(im,[],3);
        im = ifftshift(im);
    case 10 %     10 -> x,z dir
        im = fftshift(d);
        clear d;
        im = fft(im,[],2);
        im = fft(im,[],3);
        im = ifftshift(im);
        
    case 11 %     11 -> x,y,z dir
        im = fftshift(d);
        clear d;
        im = fft(im,[],1);
        im = fft(im,[],2);
        im = fft(im,[],3);
        im = ifftshift(im);
        
    otherwise
        error();
%         disp('Error using "fft3c"...')
%         disp('Invalid option.')
        im = [];
        return;
end
end