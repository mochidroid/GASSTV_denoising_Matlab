%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f(U,S,T) = lambda*|D(Ds(U))|_1 + Σ_{k=1}^{K}|P_k(U)|_* + L1ball(S) + L1ball(T) + 
%               L2ball(U+S+T) + box constraint(U) + Dv(T)=0
%
% f1(U,S,T) = 0
% f2(U,S,T) = box constraint(U) + L1ball(S) + L1ball(T)
% f3(U,S,T) = lambda*|D(Ds(U))|_1 + |P(U)|_*,N + L2ball(U+S+T) + Dv(T)=0
%
% A = (lambda*DDs O O; P O O; I I I; O O Dv)
%
% Algorithm is based on Pock's P-PDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [HSI_restored, removed_noise, other_result] ...
     = func_SSTV_LRMR_gs_for_denoising_by_PsPPDS(HSI_clean, HSI_noisy, params)
fprintf('** Running func_SSTV_LRMR_gs_for_denoising_by_PsPPDS **\n');
HSI_clean = single(HSI_clean);
HSI_noisy  = single(HSI_noisy);
HSI_noisy_gpu = gpuArray(single(HSI_noisy));
[n1, n2, n3] = size(HSI_noisy);

epsilon     = gpuArray(single(params.epsilon));
alpha       = gpuArray(single(params.alpha));
% beta        = gpuArray(single(params.beta));
lambda      = gpuArray(single(params.lambda));
blocksize   = gpuArray(single(params.blocksize));
shiftstep   = gpuArray(single(params.shiftstep));
maxiter     = gpuArray(single(params.maxiter));
stopcri     = gpuArray(single(params.stopcri));
tmp_save_name   = params.tmp_save_name;

%% Setting params
disprate    = gpuArray(single(100));
saverate    = gpuArray(single(1000));

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
% Y1: term of SSTV
% Y2: term of l2ball
% Y3: term of stripe noise
%
% Y1 = D(Ds(U))
% Y2 = P(U)
% Y3 = U + S + T
% Y4 = Dv(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y1 = zeros([n1, n2, n3, 2], 'single', 'gpuArray');
Y2 = zeros([n1, n2, n3, 1, blocksize(1)/shiftstep(1), blocksize(2)/shiftstep(2)], 'single', 'gpuArray');
Y3 = zeros([n1, n2, n3], 'single', 'gpuArray');
% Y4 = zeros([n1, n2, n3], 'single', 'gpuArray');

%% Setting operator and stepsize for Pock's P-PDS
% Difference operators with circulant boundary
D       = @(z) cat(4, z([2:end, 1],:,:) - z, z(:,[2:end, 1],:) - z);
Dt      = @(z) z([end,1:end-1],:,:,1) - z(:,:,:,1) + z(:,[end,1:end-1],:,2) - z(:,:,:,2);
% Dv      = @(z) z([2:end, 1],:,:) - z;
% Dvt     = @(z) z([end,1:end-1],:,:) - z(:,:,:);
Ds      = @(z) z(:, :, [2:end, 1], :) - z;
Dst     = @(z) z(:,:,[end,1:end-1],:) - z(:,:,:,:);

% Stepsize parameters
gamma1_U    = gpuArray(single(1./(2*2 + 2*2 + (prod(blocksize)/prod(shiftstep)) + 1)));
gamma1_S    = gpuArray(single(1));
% gamma1_T    = gpuArray(single(1/(2 + 1)));
gamma2_Y1   = gpuArray(single(1/(2*2)));
gamma2_Y2   = gpuArray(single(1));
gamma2_Y3   = gpuArray(single(1/2));
% gamma2_Y4   = gpuArray(single(1/2));

% Expansion operators
P = @(z) func_PeriodicExpansion(z, blocksize, shiftstep);
Pt = @(z) func_PeriodicExpansionTrans(z);


%% main loop (DP-PDS)
fprintf('~~~ P-PDS STARTS ~~~\n');

converge_rate_U = zeros([1, maxiter], 'single');
converge_rate_S = zeros([1, maxiter], 'single');
% converge_rate_T = zeros([1, maxiter], 'single');
converge_rate_N = zeros([1, maxiter], 'single');
move_mpsnr = zeros([1, maxiter], 'single');
move_mssim = zeros([1, maxiter], 'single');
running_time = zeros([1, maxiter], 'single');
l2ball = zeros([1, maxiter], 'single');


for i = 1:maxiter
    tic;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating U
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U_tmp   = U - gamma1_U.*(Dst(Dt(Y1)) + Pt(Y2) + Y3);
    U_next  = ProjBox(U_tmp, 0, 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating S
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_tmp   = S - gamma1_S.*Y3;
    S_next  = ProjFastL1Ball(S_tmp, alpha);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating T
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % T_tmp   = T - gamma1_T.*(Y3 + Dvt(Y4));
    % T_next  = ProjFastL1Ball(T_tmp, beta);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y1_tmp  = Y1 + gamma2_Y1.*(D(Ds(2*U_next - U)));
    Y1_next = Y1_tmp - gamma2_Y1.*Prox_l1norm(Y1_tmp./gamma2_Y1, lambda./gamma2_Y1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y2_tmp  = Y2 + gamma2_Y2*(P(2*U_next - U));    
    Y2_next = Y2_tmp - gamma2_Y2*Prox_S3TTV(Y2_tmp./gamma2_Y2, 1./gamma2_Y2, blocksize);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Y3_tmp  = Y3 + gamma2_Y3.*(2*(U_next + S_next + T_next) - (U + S + T));    
    % Y3_next = Y3_tmp - gamma2_Y3.*ProjL2ball(Y3_tmp./gamma2_Y3, HSI_noisy_gpu, epsilon);
    Y3_tmp  = Y3 + gamma2_Y3.*(2*(U_next + S_next) - (U + S));    
    Y3_next = Y3_tmp - gamma2_Y3.*ProjL2ball(Y3_tmp./gamma2_Y3, HSI_noisy_gpu, epsilon);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Y4_next = Y4 + gamma2_Y4.*Dv(2*T_next - T);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculating error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % N = HSI_noisy - U - S - T;
    % N_next = HSI_noisy - U_next - S_next - T_next;
    N = HSI_noisy - U - S;
    N_next = HSI_noisy - U_next - S_next;

    move_U = norm(U_next(:) - U(:),2)/norm(U(:),2);
    move_S = norm(S_next(:) - S(:),2)/norm(S(:),2);
    % move_T = norm(T_next(:) - T(:),2)/norm(T(:),2);
    move_N = norm(N_next(:) - N(:),2)/norm(N(:),2);
    
    converge_rate_U(i) = move_U;   
    converge_rate_S(i) = move_S;
    % converge_rate_T(i) = move_T;
    converge_rate_N(i) = move_N;
    
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

 
    % Saving results per iter
    running_time(i) = toc;

    move_mpsnr(i) = calc_MPSNR(gather(U), HSI_clean);
    move_mssim(i) = calc_MSSIM(gather(U), HSI_clean);

    l2ball(i) = norm(gather(N(:)), 2);

    
    if i>=2 && converge_rate_U(i) < stopcri
        break
    end
    if (mod(i, disprate) == 0)
        fprintf('Iter: %d, Error: %0.6f, MPSNR: %#.4g, MSSIM: %#.4g, Time: %0.2f.\n', ...
            i, move_U, move_mpsnr(i), move_mssim(i), sum(running_time));
        figure(1);
            subplot(1,2,1)
            semilogy(converge_rate_U(1:i));
            title('converge rate U TV')
            
            subplot(1,2,2)
            semilogy(l2ball(1:i));
            title('l2ball')
            drawnow
    end
    if (mod(i, saverate) == 0)
        U_curr = gather(U);
        S_curr = gather(S);
        % T_curr = gather(T);
        Y1_curr = gather(Y1);
        Y2_curr = gather(Y2);
        Y3_curr = gather(Y3);
        % Y4_curr = gather(Y4);
        i_curr = gather(i);

        save(tmp_save_name, ...
            'U_curr', 'S_curr', 'Y1_curr', 'Y2_curr', 'Y3_curr', 'i_curr', ...
            'converge_rate_U', 'converge_rate_S', 'converge_rate_N', ...
            'running_time', 'move_mpsnr', 'move_mssim', 'l2ball', 'params', ...
            '-v7.3', '-nocompression');
    end
end

fprintf('Iter: %d, Error: %0.6f, MPSNR: %#.4g, MSSIM: %#.4g, Time: %0.2f.\n', ...
    i, move_U, move_mpsnr(i), move_mssim(i), sum(running_time));

fprintf('~~~ P-PDS ENDS ~~~\n');

%% Organizing results for output
HSI_restored                        = gather(U);
other_result.iteration              = gather(i);
removed_noise.sparse_noise          = gather(S);
% removed_noise.stripe_noise          = gather(T);
% removed_noise.gaussian_noise        = HSI_noisy - HSI_restored ...
%     - removed_noise.sparse_noise - removed_noise.stripe_noise;
removed_noise.gaussian_noise        = HSI_noisy - HSI_restored ...
    - removed_noise.sparse_noise;
removed_noise.all_noise             = HSI_noisy - HSI_restored;

other_result.converge_rate_U        = gather(converge_rate_U(1:other_result.iteration));
other_result.converge_rate_S        = gather(converge_rate_S(1:other_result.iteration));
% other_result.converge_rate_T        = gather(converge_rate_T(1:other_result.iteration));
other_result.converge_rate_N        = gather(converge_rate_N(1:other_result.iteration));

other_result.move_mpsnr             = gather(move_mpsnr(1:other_result.iteration));
other_result.move_mssim             = gather(move_mssim(1:other_result.iteration));
other_result.running_time           = gather(running_time(1:other_result.iteration));

other_result.l2ball                 = gather(l2ball(1:other_result.iteration));


%% Plotting result
figure;
subplot(1,2,1)
semilogy(other_result.converge_rate_U);
title('converge rate U TV')

subplot(1,2,2)
semilogy(other_result.l2ball);
title('l2ball')
drawnow