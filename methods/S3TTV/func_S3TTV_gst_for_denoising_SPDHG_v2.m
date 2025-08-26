%% Spatio-Spectral Structure Tensor Total Variation for Hyperspectral Image Denoising and Destriping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f(U,S,T) = |P(D(Dl(U)))|_* + L1ball(S) + L1ball(T) + 
%               L2ball(U+S+T) + box constraint(U) + Dv(T)=0
%
% f1(U,S,T) = 0
% f2(U,S,T) = box constraint(U) + L1ball(S) + L1ball(T)
% f3(U,S,T) = |P(D(Dl(U)))|_* + L2ball(U+S+T) + Dv(T)=0
%
% A = (PDDl O O; I I I; O O Dv)
%
% Algorithm is based on SPDHG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HSI_restored, removed_noise, other_result] ...
     = func_S3TTV_gst_for_denoising_SPDHG_v2(HSI_clean, HSI_noisy, params)
fprintf('** Running func_S3TTV_gst_for_denoising_SPDHG_v2 **\n');
HSI_clean = single(HSI_clean);
HSI_noisy  = single(HSI_noisy);
HSI_noisy_gpu = gpuArray(single(HSI_noisy));
[n1, n2, n3] = size(HSI_noisy);

alpha       = gpuArray(single(params.alpha));
beta        = gpuArray(single(params.beta));
epsilon     = gpuArray(single(params.epsilon));
blocksize   = gpuArray(single(params.blocksize));
maxiter     = gpuArray(single(params.maxiter));
stopcri     = gpuArray(single(params.stopcri));

prob_patch  = gpuArray(single(params.prob_patch));
eta         = gpuArray(single(params.eta));
clip_c      = log(2);
seed        = 'default';

%% Setting params
dispiter    = unique([1:9, 10:10:maxiter]);
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
T = zeros([n1, n2, n3], 'single', 'gpuArray');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dual variables
% Y1: term of S3TTV
% Y2: term of box constraint
% Y3: term of sparse noise 
% Y4: term of stripe noise 
% Y5: stripe noise
%
% Y1 = Pk(D(Dl(U)))
% Y2 = U + S + T
% Y3 = Dv(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y1 = zeros([n1, n2, n3, 2, blocksize(1), blocksize(2)], 'single', 'gpuArray');
Y2 = zeros([n1, n2, n3], 'single', 'gpuArray');
Y3 = zeros([n1, n2, n3], 'single', 'gpuArray');


%% Setting operators
% Difference operators
D       = @(z) cat(4, z([2:end, 1],:,:) - z, z(:,[2:end, 1],:) - z);
Dt      = @(z) z([end,1:end-1],:,:,1) - z(:,:,:,1) + z(:,[end,1:end-1],:,2) - z(:,:,:,2);
Dv      = @(z) z([2:end, 1],:,:) - z;
Dvt     = @(z) z([end,1:end-1],:,:) - z(:,:,:);
Dl      = @(z) z(:, :, [2:end, 1], :) - z;
Dlt     = @(z) z(:,:,[end,1:end-1],:) - z(:,:,:,:);

% Expansion operators
P = @(z) func_PeriodicExpansion(z, blocksize);
Pt = @(z) func_PeriodicExpansionTrans(z);


%% Setting stepsize parameters for P-PDS
sq_opnorm_D = 8;
sq_opnorm_Dl = 4;
sq_opnorm_Dv = 4;
sq_opnorm_P = (n1*n2)^2;

% sigma_YL = gpuArray(single(1/(sq_opnorm_P * sq_opnorm_Dl * sq_opnorm_D + 1)));
sigma_YL = gpuArray(single(1/(sq_opnorm_Dl * sq_opnorm_D)));
sigma_Y2 = gpuArray(single(1/3));
sigma_Y3 = gpuArray(single(1/sq_opnorm_Dv));

tau = gpuArray(single(1/(2 + n1*n2*p)));


%% main loop (P-PDS)
converge_rate_U = zeros([1, maxiter], "single");
converge_rate_S = zeros([1, maxiter], "single");
converge_rate_T = zeros([1, maxiter], "single");
converge_rate_N = zeros([1, maxiter], "single");
move_mpsnr = zeros([1, maxiter], "single");
move_mssim = zeros([1, maxiter], "single");
running_time = zeros([1, maxiter], "single");
l2ball = zeros([1, maxiter], "single");

fprintf('~~~ P-PDS STARTS ~~~\n');

for i = 1:maxiter
    tic;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating U
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U_tmp   = U - tau*(Dlt(Dt(Pt(Y1))) + Y2);
    U_next  = ProjBox(U_tmp, 0, 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating S
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_tmp   = S - tau*(Y2);
    S_next  = ProjFastL1Ball(S_tmp, alpha);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating T
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T_tmp   = T - tau*(Y2 + Dvt(Y3));
    T_next  = ProjFastL1Ball(T_tmp, beta);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---（新）Y1更新：確率的ミニバッチ + SPDHG extrapolation ---
    % 1) ミニバッチ選択（b1*b2 の中から m 個）
    idx = randperm(L, m);
    sel_mask = false(b1, b2); 
    sel_mask(idx) = true;

    % 選んだシフト部分だけ更新、それ以外はY1_new = Y1_tmp = Y1;
    Y1_tmp = Y1 + sigma_YL * (P(D(Dl(U_next))).* reshape(sel_mask, [1, 1, 1, 1, b1, b2]));
    Y1_new = Y1_tmp - sigma_YL * Prox_S3TTV_patch(Y1_tmp/sigma_YL, 1/sigma_YL, sel_mask);
    Y1_next = Y1_new + ((1 / p) * (Y1_new - Y1).* reshape(sel_mask, [1, 1, 1, 1, b1, b2]));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y2_tmp  = Y2 + sigma_Y2*(U_next + S_next + T_next);
    Y2_new  = Y2_tmp - sigma_Y2*ProjL2ball(Y2_tmp/sigma_Y2, HSI_noisy_gpu, epsilon);
    Y2_next = 2*Y2_new - Y2;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y3_next = Y3 + 2*sigma_Y3*Dv(T_next);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculating error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = HSI_noisy - U - S - T;
    N_next = HSI_noisy - U_next - S_next - T_next;

    converge_rate_U(i) = norm(U_next(:) - U(:),2)/norm(U(:),2);
    converge_rate_S(i) = norm(S_next(:) - S(:),2)/norm(S(:),2);
    converge_rate_T(i) = norm(T_next(:) - T(:),2)/norm(T(:),2);
    converge_rate_N(i) = norm(N_next(:) - N(:),2)/norm(N(:),2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating stepsizes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r_primal = norm(U_next(:) - U(:),2) + norm(S_next(:) - S(:),2) + norm(T_next(:) - T(:),2);
    r_dual = norm(Y1_next(:) - Y1(:),2) + norm(Y2_next(:) - Y2(:),2) ...
        + norm(Y3_next(:) - Y3(:),2);

    % gamma   = eta * (r_primal - r_dual);
    % gamma   = max(-clip_c, min(clip_c, gamma));  % clip
    % gamma   = exp(gamma);

    gamma = 1;

    tau      = tau / gamma;
    sigma_YL = sigma_YL * gamma;
    sigma_Y2 = sigma_Y2 * gamma;
    sigma_Y3 = sigma_Y3 * gamma;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating all variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U   = U_next;
    S   = S_next;
    T   = T_next;
    
    Y1  = Y1_next;
    Y2  = Y2_next;
    Y3  = Y3_next;

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
removed_noise.stripe_noise          = gather(T);
removed_noise.gaussian_noise        = gather(HSI_noisy_gpu - U - S - T);
removed_noise.all_noise             = gather(HSI_noisy_gpu - U);

other_result.converge_rate_U        = gather(converge_rate_U(1:other_result.iteration));
other_result.converge_rate_S        = gather(converge_rate_S(1:other_result.iteration));
other_result.converge_rate_T        = gather(converge_rate_T(1:other_result.iteration));
other_result.converge_rate_N        = gather(converge_rate_N(1:other_result.iteration));

other_result.move_mpsnr             = gather(move_mpsnr(1:other_result.iteration));
other_result.move_mssim             = gather(move_mssim(1:other_result.iteration));
other_result.running_time           = gather(running_time(1:other_result.iteration));

other_result.l2ball                 = gather(l2ball(1:other_result.iteration));
