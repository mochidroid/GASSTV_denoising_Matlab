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
% Algorithm is based on Naganuma's P-PDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HSI_restored, removed_noise, other_result] ...
     = func_S3TTV_gs_for_denoising_NPPDS(HSI_clean, HSI_noisy, params)
fprintf('** Running func_S3TTV_gs_for_denoising_NPPDS **\n');
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

%% Setting params
dispiter    = unique([1:10, 1000:1000:maxiter]);
dispband    = round(n3/2);


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
% Y2: term of l2ball
% Y3: term of stripe noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y1 = zeros([n1, n2, n3, 2, blocksize(1), blocksize(2)], 'single', 'gpuArray');
Y2 = zeros([n1, n2, n3], 'single', 'gpuArray');
% Y3 = zeros([n1, n2, n3], 'single', 'gpuArray');


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
gamma1_U    = gpuArray(single(1./(prod(blocksize)^2 * 8*4 + 1)));
gamma1_S    = gpuArray(single(1));
% gamma1_T    = gpuArray(single(1/(4 + 1)));
% gamma2      = gpuArray(single(1/3));
gamma2      = gpuArray(single(1/2));


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
    U_tmp   = U - gamma1_U*(Dlt(Dt(Pt(Y1))) + Y2);
    U_next  = ProjBox(U_tmp, 0, 1);
    U_res = 2*U_next - U;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating S
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_tmp   = S - gamma1_S*(Y2);
    S_next  = ProjFastL1Ball(S_tmp, alpha);
    S_res = 2*S_next - S;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating T
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % T_tmp   = T - gamma1_T*(Y2 + Dvt(Y3));
    % T_next  = ProjFastL1Ball(T_tmp, beta);
    % T_res = 2*T_next - T;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y1_tmp  = Y1 + gamma2*(P(D(Dl(U_res))));
    Y1_next = Y1_tmp - gamma2*Prox_S3TTV(Y1_tmp/gamma2, 1/gamma2, blocksize);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Y2_tmp  = Y2 + gamma2*(U_res + S_res + T_res);
    Y2_tmp  = Y2 + gamma2*(U_res + S_res);
    Y2_next = Y2_tmp - gamma2*ProjL2ball(Y2_tmp/gamma2, HSI_noisy, epsilon);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Y3_next = Y3 + gamma2.*Dv(T_res);

    
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
    % Updating all variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U   = U_next;
    S   = S_next;
    % T   = T_next;
    
    Y1  = Y1_next;
    Y2  = Y2_next;
    % Y3  = Y3_next;

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
