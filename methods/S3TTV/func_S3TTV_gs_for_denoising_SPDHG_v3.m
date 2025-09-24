%% Spatio-Spectral Structure Tensor Total Variation for Hyperspectral Image Denoising and Destriping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f(U,S,T) = |P(D(Dl(U)))|_* + L1ball(S) + L1ball(T) + 
%               L2ball(U+S+T) + box constraint(U) + Dv(T)=0
%
% f1(U,S,T) = 0
% f2(U,S,T) = L2ball(U+S+T)
% f3(U,S,T) = |P(D(Dl(U)))|_* + box constraint(U) + L1ball(S) + L1ball(T) + Dv(T)=0
%
% A = (PDDl O O; I O O; O I O; O O I; O O Dv)
%
% Algorithm is based on SPDHG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HSI_restored, removed_noise, other_result] ...
     = func_S3TTV_gs_for_denoising_SPDHG_v3(HSI_clean, HSI_noisy, params)
fprintf('** Running func_S3TTV_gs_for_denoising_SPDHG_v3 **\n');
HSI_clean = single(HSI_clean);
HSI_noisy  = single(HSI_noisy);
HSI_noisy_gpu = gpuArray(single(HSI_noisy));
[n1, n2, n3] = size(HSI_noisy);

alpha       = gpuArray(single(params.alpha));
% beta        = gpuArray(single(params.beta));
epsilon     = gpuArray(single(params.epsilon));
blocksize   = gpuArray(single(params.blocksize));
maxiter     = gpuArray(single(params.maxiter));
stopcri     = gpuArray(single(params.stopcri));

prob_patch  = gpuArray(single(params.prob_patch));
eta         = gpuArray(single(params.eta));
clip_c      = log(2);
seed        = 'default';

%% Setting params
dispiter    = unique([1:10, 100:100:maxiter]);
dispband    = round(n3/2);


b1 = blocksize(1);
b2 = blocksize(2);

assert(mod(n1,b1)==0 && mod(n2,b2)==0, 'n1,n2 は blocksize で割り切れる必要があります。');

nb1 = n1/b1; nb2 = n2/b2; 
L    = b1*b2;           % シフト総数（= S3TTV の dual ブロック個数）
p    = prob_patch;   % 更新確率
m    = max(1, round(p*L));  % select shifts per iteration

if ~isempty(seed); rng(seed); end


%% Initializing primal and dual variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% primal variables
% U: clean HSI
% S: sparse noise(salt-and-pepper noise)
% T: stripe noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = zeros([n1, n2, n3], 'single', 'gpuArray');
S = zeros([n1, n2, n3], 'single', 'gpuArray');
% T = zeros([n1, n2, n3], 'single', 'gpuArray');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dual variables
% Y1: term of S3TTV
% Y2: term of box constraint
% Y3: term of sparse noise 
% Y4: term of stripe noise 
% Y5: stripe noise
%
% Y1 = Pk(D(Dl(U)))
% Y2 = U
% Y3 = S
% Y4 = T
% Y5 = Dv(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y1 = zeros([n1, n2, n3, 2, blocksize(1), blocksize(2)], 'single', 'gpuArray');
Y2 = zeros([n1, n2, n3], 'single', 'gpuArray');
Y3 = zeros([n1, n2, n3], 'single', 'gpuArray');
% Y4 = zeros([n1, n2, n3], 'single', 'gpuArray');
% Y5 = zeros([n1, n2, n3], 'single', 'gpuArray');

Y1_bar = Y1;
Y2_bar = Y2;
Y3_bar = Y3;


%% Setting operators
% Difference operators
D       = @(z) cat(4, z([2:end, 1],:,:) - z, z(:,[2:end, 1],:) - z);
Dt      = @(z) z([end,1:end-1],:,:,1) - z(:,:,:,1) + z(:,[end,1:end-1],:,2) - z(:,:,:,2);
% Dv      = @(z) z([2:end, 1],:,:) - z;
% Dvt     = @(z) z([end,1:end-1],:,:) - z(:,:,:);
Dl      = @(z) z(:, :, [2:end, 1], :) - z;
Dlt     = @(z) z(:,:,[end,1:end-1],:) - z(:,:,:,:);

% Expansion operators
P = @(z) func_PeriodicExpansion(z, blocksize);
Pt = @(z) func_PeriodicExpansionTrans(z);


%% Setting stepsize parameters for P-PDS
sq_opnorm_D = 8;
sq_opnorm_Dl = 4;
% sq_opnorm_Dv = 4;
sq_opnorm_P = (b1*b2)^2; 

% sigma_YL = gpuArray(single(1/(sq_opnorm_P * sq_opnorm_Dl * sq_opnorm_D + 1)));
sigma_Y1 = gpuArray(single(1/(sq_opnorm_P * sq_opnorm_Dl * sq_opnorm_D)));
sigma_Y2 = gpuArray(single(1));
sigma_Y3 = gpuArray(single(1));
% sigma_Y4 = gpuArray(single(1));
% sigma_Y5 = gpuArray(single(1/sq_opnorm_Dv));

tau = 0.9 * prob_patch;


%% main loop (P-PDS)
converge_rate_U = zeros([1, maxiter], "single");
converge_rate_S = zeros([1, maxiter], "single");
% converge_rate_T = zeros([1, maxiter], "single");
converge_rate_N = zeros([1, maxiter], "single");
move_mpsnr = zeros([1, maxiter], "single");
move_mssim = zeros([1, maxiter], "single");
running_time = zeros([1, maxiter], "single");
l2ball = zeros([1, maxiter], "single");
gamma = ones([1, maxiter], "single", "gpuArray");

r_primal = zeros(1,maxiter,'single');
r_dual   = zeros(1,maxiter,'single');

fprintf('~~~ P-PDS STARTS ~~~\n');

for i = 1:maxiter
    tic;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating U, S, T
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U_tmp   = U - tau*(Dlt(Dt(Pt(Y1_bar))) + Y2_bar);
    S_tmp   = S - tau*Y3_bar;
    % T_tmp   = T - tau*(Y4 + Dvt(Y5));

    Primal_sum = U_tmp + S_tmp;
    Primal_sum = ProjL2ball(Primal_sum, HSI_noisy, epsilon) - Primal_sum;

    U_next = U_tmp + 0.5 * Primal_sum;
    S_next = S_tmp + 0.5 * Primal_sum;
    % T_next = T_tmp + Primal_sum/3;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---（新）Y1更新：確率的ミニバッチ + SPDHG extrapolation ---
    % 1) ミニバッチ選択（b1*b2 の中から m 個）
    idx = randperm(L, m);
    sel_mask = false(b1,b2); sel_mask(idx)=true;
    mask = reshape(gpuArray(single(sel_mask)), [1,1,1,1,b1,b2]);

    % 選んだシフト部分だけ更新、それ以外はY1_new = Y1_tmp = Y1;
    Y1_tmp  = Y1 + sigma_Y1 * (P(D(Dl(U_next))).* mask);
    Y1_next  = Y1_tmp - sigma_Y1 * Prox_S3TTV_patch(Y1_tmp/sigma_Y1, 1/sigma_Y1, sel_mask);
    Y1_bar_new = Y1_next + (1/p)*(Y1_next - Y1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y2_tmp  = Y2 + sigma_Y2*(U_next);
    Y2_next = Y2_tmp - sigma_Y2*ProjBox(Y2_tmp/sigma_Y2, 0, 1);
    Y2_bar_new  = 2*Y2_next - Y2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y3_tmp  = Y3 + sigma_Y3*(S_next);
    Y3_next = Y3_tmp - sigma_Y3*ProjFastL1Ball(Y3_tmp/sigma_Y3, alpha);
    Y3_bar_new  = 2*Y3_next - Y3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Y4_tmp  = Y4 + sigma_Y4*(T_next);
    % Y4_new  = Y4_tmp - sigma_Y4*ProjFastL1Ball(Y4_tmp/sigma_Y4, beta);
    % Y4_next = 2*Y4_new - Y4;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Y5_next = Y5 + 2*sigma_Y5*Dv(T_next);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculating error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % N = HSI_noisy - U - S - T;
    % N_next = HSI_noisy - U_next - S_next - T_next;
    N = HSI_noisy - U - S;
    N_next = HSI_noisy - U_next - S_next;

    converge_rate_U(i) = norm(U_next(:) - U(:),2)/norm(U(:),2);
    converge_rate_S(i) = norm(S_next(:) - S(:),2)/norm(S(:),2);
    % converge_rate_T(i) = norm(T_next(:) - T(:),2)/norm(T(:),2);
    converge_rate_N(i) = norm(N_next(:) - N(:),2)/norm(N(:),2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating stepsizes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r_primal(i) = (converge_rate_U(i) + converge_rate_S(i)) / 2;


    mask_Y1_bar_new = Y1_bar_new .* mask;
    mask_Y1_bar = Y1_bar .* mask;
    dual1   = norm(mask_Y1_bar_new(:) - mask_Y1_bar(:), 2)/norm(mask_Y1_bar(:), 2);

    dual2   = norm(Y2_bar_new(:) - Y2_bar(:), 2)/norm(Y2_bar(:), 2);
    dual3   = norm(Y3_bar_new(:) - Y3_bar(:), 2)/norm(Y3_bar(:), 2);

    r_dual(i)  = (dual1 + dual2 + dual3) / 3;

    % gamma(i) = exp( max(-clip_c, min(clip_c, eta*(r_primal(i) - r_dual(i)))) );

    % gamma = 1;

    tau      = tau / gamma(i);
    sigma_Y1 = sigma_Y1 * gamma(i);
    sigma_Y2 = sigma_Y2 * gamma(i);
    sigma_Y3 = sigma_Y3 * gamma(i);
    % sigma_Y4 = sigma_Y4 * gamma(i);
    % sigma_Y5 = sigma_Y5 * gamma(i);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating all variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U   = U_next;
    S   = S_next;
    % T   = T_next;
    
    Y1  = Y1_next;
    Y2  = Y2_next;
    Y3  = Y3_next;
    % Y4  = Y4_next;
    % Y5  = Y5_next;

    Y1_bar = Y1_bar_new;
    Y2_bar = Y2_bar_new;
    Y3_bar = Y3_bar_new;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convergence checking
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Saving results per iter
    running_time(i) = toc;

    move_mpsnr(i) = calc_MPSNR(gather(U), HSI_clean);
    move_mssim(i) = calc_MSSIM(gather(U), HSI_clean);

    l2ball(i) = norm(gather(N(:)), 2);
    
    if i>=2 && converge_rate_U(i) < stopcri
        break
    end
    
    % Displaying progress 
    if ismember(i, dispiter)
        fprintf("Iter: %d, Error: %0.6f, MPSNR: %#.4g, MSSIM: %#.4g, Time: %0.2f.\n", ...
            i, converge_rate_U(i), move_mpsnr(i), move_mssim(i), sum(running_time));
        % fprintf("Iter: %d, Error: %0.6f, MPSNR: %#.4g, MSSIM: %#.4g, Gamma: %0.2f, Time: %0.2f.\n", ...
        %     i, converge_rate_U(i), move_mpsnr(i), move_mssim(i), gamma(i), sum(running_time));

        figure(1)
        subplot(2,3,1)
        imshow(HSI_clean(:,:,dispband));
        title("GT")
        
        subplot(2,3,2)
        imshow(HSI_noisy(:,:,dispband));
        title("Observed")
        
        subplot(2,3,3)
        imshow(gather(U(:,:,dispband)));
        title("Restored")

        subplot(2,3,4)
        semilogy(converge_rate_U(1:i));
        title("Converge rate U TV")
        
        subplot(2,3,5)
        semilogy(l2ball(1:i));
        title("L2ball")
        drawnow
    end
end

fprintf("Iter: %d, Error: %0.6f, MPSNR: %#.4g, MSSIM: %#.4g, Time: %0.2f.\n", ...
    i, converge_rate_U(i), move_mpsnr(i), move_mssim(i), sum(running_time));

fprintf('~~~ P-PDS ENDS ~~~\n');

%% Organizing results for output
HSI_restored                        = gather(U);
other_result.iteration              = gather(i);
removed_noise.sparse_noise          = gather(S);
% removed_noise.stripe_noise          = gather(T);
% removed_noise.gaussian_noise        = gather(HSI_noisy_gpu - U - S - T);
removed_noise.gaussian_noise        = gather(HSI_noisy_gpu - U - S);
removed_noise.all_noise             = gather(HSI_noisy_gpu - U);

other_result.converge_rate_U        = gather(converge_rate_U(1:other_result.iteration));
other_result.converge_rate_S        = gather(converge_rate_S(1:other_result.iteration));
% other_result.converge_rate_T        = gather(converge_rate_T(1:other_result.iteration));
other_result.converge_rate_N        = gather(converge_rate_N(1:other_result.iteration));

other_result.move_mpsnr             = gather(move_mpsnr(1:other_result.iteration));
other_result.move_mssim             = gather(move_mssim(1:other_result.iteration));
other_result.running_time           = gather(running_time(1:other_result.iteration));

other_result.l2ball                 = gather(l2ball(1:other_result.iteration));

other_result.r_primal               = r_primal(1:i);
other_result.r_dual                 = r_dual(1:i);
% other_result.gamma                  = gamma(1:i);
