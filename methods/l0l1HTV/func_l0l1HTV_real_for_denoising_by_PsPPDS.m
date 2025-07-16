%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f(U,S,T) = |D(Ds(U))|_1 + l_{1,0}(D(U)) + L1ball(T) + 
%               L2ball(U+S+T) + box constraint(U) + Dv(T)=0
%
% f1(U,S,T) = 0
% f2(U,S,T) = box constraint(U) + L1ball(S) + L1ball(T)
% f3(U,S,T) = |D(Ds(U))|_1 + l_{1,0}(D(U)) + L2ball(U+S+T) + Dv(T)=0
%
% A = (DDs O O; D O O; I I I; O O Dv)
%
% Algorithm is based on Pock's P-PDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [HSI_restored, removed_noise, other_result] ...
     = func_l0l1HTV_real_for_denoising_by_PsPPDS(HSI_noisy, params)
fprintf('** Running func_l0l1HTV_real_for_denoising_by_PsPPDS **\n');
HSI_noisy  = single(HSI_noisy);
HSI_noisy_gpu = gpuArray(single(HSI_noisy));
[n1, n2, n3] = size(HSI_noisy);

epsilon             = gpuArray(single(params.epsilon));
eps_L10ball         = gpuArray(single(params.eps_L10ball));
alpha               = gpuArray(single(params.alpha));
beta                = gpuArray(single(params.beta));
stepsize_reduction  = gpuArray(single(params.stepsize_reduction));
boundary            = params.boundary;
maxiter             = gpuArray(single(params.maxiter));
stopcri             = gpuArray(single(params.stopcri));

%% Setting params
disprate    = gpuArray(single(1000));

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
% Y1: term of SSTV
% Y2: term of l_{1,0} norm
% Y3: term of l2ball
% Y4: term of stripe noise
%
% Y1 = D(Ds(U))
% Y2 = l_{1,0}(D(U))
% Y3 = U + S + T
% Y4 = Dv(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y1 = zeros([n1, n2, n3, 2], 'single', 'gpuArray');
Y2 = zeros([n1, n2, n3, 2], 'single', 'gpuArray');
Y3 = zeros([n1, n2, n3], 'single', 'gpuArray');
Y4 = zeros([n1, n2, n3], 'single', 'gpuArray');

%% Setting operator and stepsize for Pock's P-PDS

switch boundary
    case 'Neumann' % Difference operators with Neumann boundary
        D       = @(z) cat(4, z([2:end, end],:,:) - z, z(:,[2:end, end],:) - z);
        Dt      = @(z) cat(1, -z(1, :, :, 1), -z(2:end-1, :, :, 1) + z(1:end-2, :, :, 1), z(end-1, :, :, 1)) ...
                + cat(2, -z(:, 1, :, 2), -z(:, 2:end-1, :, 2) + z(:, 1:end-2, :, 2), z(:, end-1, :, 2));
        Dv      = @(z) z([2:end, end],:,:) - z;
        Dvt     = @(z) cat(1, -z(1, :, :), -z(2:(n1-1), :, :) + z(1:(n1-2), :, :), z(n1-1, :, :));
        Ds      = @(z) z(:, :, [2:end, end], :) - z;
        Dst     = @(z) cat(3, -z(:, :, 1), -z(:, :, 2:n3-1) + z(:, :, 1:n3-2), z(:, :, n3-1));
        one = ones(n1, n2, n3);

        precon1_Dv = cat(1, one(1, :, :), 2*one(2:(n1-1), :, :), one(n1-1, :, :));
        precon1_Dh = cat(2, one(:, 1, :), 2*one(:, 2:(n2-1), :), one(:, n2-1, :));
        precon1_Db = cat(3, one(:, :, 1), 2*one(:, :, 2:(n3-1)), one(:, :, n3-1));
        
        zero_para = 100;
        
        precon2_Dv = cat(1, 2*ones(n1-1, n2), zero_para*ones(1, n2)).*ones(n1, n2, n3);
        precon2_Dh = cat(2, 2*ones(n1, n2-1), zero_para*ones(n1, 1)).*ones(n1, n2, n3);
        precon2_Db = cat(3, 2*ones(n1, n2, n3-1), zero_para*ones(n1, n2, 1));

        gamma1_U    = gpuArray(single(1./(precon1_Dv.*precon1_Db + precon1_Dh.*precon1_Db + ...
                        precon1_Dv + precon1_Dh + 1)));
        % gamma1_S    = gpuArray(single(1));
        gamma1_T    = gpuArray(single(1./(precon1_Dv + 1)));
        gamma2_Y1_1 = gpuArray(single(1/(precon2_Dv.*precon2_Db)));
        gamma2_Y1_2 = gpuArray(single(1/(precon2_Dh.*precon2_Db)));
        gamma2_Y1 = cat(4, gamma2_Y1_1, gamma2_Y1_2);
        gamma2_Y2_1 = gpuArray(single(1/precon2_Dv));
        gamma2_Y2_2 = gpuArray(single(1/precon2_Dh));
        gamma2_Y2   = cat(4, gamma2_Y2_1, gamma2_Y2_2);
        gamma2_Y3   = gpuArray(single(1));
        % gamma2_Y4   = gpuArray(single(1./(precon2_Dv)));

        %memory release
        clear precon1_Dv precon1_Dh precon1_Db
        clear precon2_Dv precon1_Dh precon1_Db
        clear zero_para one
        clear gamma2_Y1_1 gamma2_Y1_2 gamma2_Y2_1 gamma2_Y2_2

    case 'circulant' % Difference operators with circulant boundary
        D       = @(z) cat(4, z([2:end, 1],:,:) - z, z(:,[2:end, 1],:) - z);
        Dt      = @(z) z([end,1:end-1],:,:,1) - z(:,:,:,1) + z(:,[end,1:end-1],:,2) - z(:,:,:,2);
        Dv      = @(z) z([2:end, 1],:,:) - z;
        Dvt     = @(z) z([end,1:end-1],:,:) - z(:,:,:);
        Ds      = @(z) z(:, :, [2:end, 1], :) - z;
        Dst     = @(z) z(:,:,[end,1:end-1],:) - z(:,:,:,:);

        gamma1_U    = gpuArray(single(1./(2*2 + 2*2 + 2 + 2 + 1)));
        gamma1_S    = gpuArray(single(1));
        gamma1_T    = gpuArray(single(1/(2 + 1)));
        gamma2_Y1   = gpuArray(single(1/(2*2)));
        gamma2_Y2   = gpuArray(single(1/2));
        gamma2_Y3   = gpuArray(single(1/3));
        gamma2_Y4   = gpuArray(single(1/2));

end

%% main loop (P-PDS)
fprintf('~~~ P-PDS STARTS ~~~\n');

converge_rate_U = zeros([1, maxiter], 'single');
converge_rate_S = zeros([1, maxiter], 'single');
converge_rate_T = zeros([1, maxiter], 'single');
converge_rate_N = zeros([1, maxiter], 'single');
move_function_value = zeros([1, maxiter], 'single');
running_time = zeros([1, maxiter], 'single');
l2ball = zeros([1, maxiter], 'single');


for i = 1:maxiter
    tic;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating U
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U_tmp   = U - gamma1_U.*(Dst(Dt(Y1)) + Dt(Y2) + Y3);
    U_next  = ProjBox(U_tmp, 0, 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating S
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_tmp   = S - gamma1_S.*Y3;
    S_next  = ProjFastL1Ball(S_tmp, alpha);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating T
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T_tmp   = T - gamma1_T.*(Y3 + Dvt(Y4));
    T_next  = ProjFastL1Ball(T_tmp, beta);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y1_tmp  = Y1 + gamma2_Y1.*(D(Ds(2*U_next - U)));
    Y1_next = Y1_tmp - gamma2_Y1.*Prox_l1norm(Y1_tmp./gamma2_Y1, 1/gamma2_Y1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y2_tmp  = Y2 + gamma2_Y2.*(D(2*U_next - U));
    Y2_next = Y2_tmp - gamma2_Y2.*ProjL10ball(Y2_tmp./gamma2_Y2, eps_L10ball);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y3_tmp  = Y3 + gamma2_Y3.*(2*(U_next + S_next + T_next) - (U + S + T));    
    Y3_next = Y3_tmp - gamma2_Y3.*ProjL2ball(Y3_tmp./gamma2_Y3, HSI_noisy_gpu, epsilon);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y4_next = Y4 + gamma2_Y4.*Dv(2*T_next - T);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculating error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = HSI_noisy - U - S - T;
    N_next = HSI_noisy - U_next - S_next - T_next;

    move_U = norm(U_next(:) - U(:),2)/norm(U(:),2);
    move_S = norm(S_next(:) - S(:),2)/norm(S(:),2);
    move_T = norm(T_next(:) - T(:),2)/norm(T(:),2);
    move_N = norm(N_next(:) - N(:),2)/norm(N(:),2);
    
    converge_rate_U(i) = move_U;   
    converge_rate_S(i) = move_S;
    converge_rate_T(i) = move_T;
    converge_rate_N(i) = move_N;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating all variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U   = U_next;
    S   = S_next;
    T   = T_next;
    
    Y1  = Y1_next;
    Y2  = Y2_next;
    Y3  = Y3_next;
    Y4  = Y4_next;

    gamma1_U    = gamma1_U  * stepsize_reduction;
    gamma1_S    = gamma1_S  * stepsize_reduction;
    gamma1_T    = gamma1_T  * stepsize_reduction;
    gamma2_Y1   = gamma2_Y1 * stepsize_reduction;
    gamma2_Y2   = gamma2_Y2 * stepsize_reduction;
    gamma2_Y3   = gamma2_Y3 * stepsize_reduction;
    gamma2_Y4   = gamma2_Y4 * stepsize_reduction;
 
    % Saving results per iter
    running_time(i) = toc;

    move_function_value(i) = sum(abs(D(Ds(U))), "all");
    l2ball(i) = norm(gather(N(:)), 2);

    
    if i>=2 && converge_rate_U(i) < stopcri
        break
    end
    if (mod(i, disprate) == 0)
        fprintf('Iter: %d, Error: %0.6f, FV: %#.4g, Time: %0.2f.\n', ...
            i, move_U, move_function_value(i), sum(running_time));
    end
end

fprintf('Iter: %d, Error: %0.6f, FV: %#.4g, Time: %0.2f.\n', ...
            i, move_U, move_function_value(i), sum(running_time));

fprintf('~~~ P-PDS ENDS ~~~\n');

%% Organizing results for output
HSI_restored                        = gather(U);
other_result.iteration              = gather(i);
removed_noise.sparse_noise          = gather(S);
removed_noise.stripe_noise          = gather(T);
removed_noise.gaussian_noise    = HSI_noisy - HSI_restored - ...
                                    removed_noise.sparse_noise - removed_noise.stripe_noise;
removed_noise.all_noise             = HSI_noisy - HSI_restored;

other_result.converge_rate_U        = gather(converge_rate_U(1:other_result.iteration));
other_result.converge_rate_S        = gather(converge_rate_S(1:other_result.iteration));
other_result.converge_rate_T        = gather(converge_rate_T(1:other_result.iteration));
other_result.converge_rate_N        = gather(converge_rate_N(1:other_result.iteration));

other_result.move_function_value    = gather(move_function_value(1:other_result.iteration));
other_result.running_time           = gather(running_time(1:other_result.iteration));

other_result.l2ball                 = gather(l2ball(1:other_result.iteration));


%% Plotting result
figure;
subplot(1,3,1)
semilogy(other_result.converge_rate_U);
title('converge rate U TV')

subplot(1,3,2)
semilogy(other_result.move_function_value);
title('function value')

subplot(1,3,3)
semilogy(other_result.l2ball);
title('l2ball')
drawnow