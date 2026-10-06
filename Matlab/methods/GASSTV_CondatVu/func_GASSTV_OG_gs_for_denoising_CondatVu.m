%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f(U,S,T) = |D(Ds(U))|_1 + \lambda1|WspDsp(U)|_{1or2} + \lambda2|WsDsU|_{1or2} +
%             L1ball(S) + L1ball(T) + L2ball(U+S+T) + box constraint(U) + Dv(T)=0
%
% f1(U,S,T) = 0
% f2(U,S,T) = L2ball(U+S+T)
% f3(U,S,T) = |D(Ds(U))|_1 + \lambda1|WspDsp(U)|_{1or2} + \lambda2|WsDsU|_{1or2} +
%              box constraint(U) + L1ball(S) + L1ball(T) + Dv(T)=0
%
% A = (DDs O O; WsDs O O; I O O; O I O; O O I; O O Dv)
%
% Algorithm is based on Naganuma's P-PDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f(U,S) = |D(Ds(U))|_1 + \lambda1|WspDsp(U)|_{1or2} + \lambda2|WsDsU|_{1or2} +
%             L2ball(U+S) + box constraint(U) + L1ball(S)
%
% f1(U,S) = 0
% f2(U,S) = L2ball(U+S)
% f3(U,S) = |D(Ds(U))|_1 + \lambda1|WspDsp(U)|_{1or2} + \lambda2|WsDsU|_{1or2} +
%          box constraint(U) + L1ball(S)
%
% A = (DDs O; WspDsp O; WsDs O; I O; O I)
%
% Algorithm is based on Naganuma's P-PDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [HSI_restored, removed_noise, output] ...
     = func_GASSTV_OG_gs_for_denoising_CondatVu(HSI_clean, HSI_noisy, params)
fprintf('** Running func_GASSTV_OG_gs_for_denoising_CondatVu **\n');
HSI_clean = single(HSI_clean);
HSI_noisy  = single(HSI_noisy);
HSI_noisy_gpu = gpuArray(single(HSI_noisy));
[n1, n2, n3] = size(HSI_noisy);

epsilon         = gpuArray(single(params.epsilon));
alpha           = gpuArray(single(params.alpha));
% beta            = gpuArray(single(params.beta));
% sigma_s         = single(params.sigma_s);
lambda_rho_sp   = gpuArray(single(params.lambda_rho_sp));
lambda2         = gpuArray(single(params.lambda2));
maxiter         = gpuArray(single(params.maxiter));
stopcri         = gpuArray(single(params.stopcri));

% Spatial graph
sigma_sp        = params.sigma_sp;

% Spectral graph
sigma_l         = params.sigma_l;
num_segments    = single(params.num_segments);
order_filt      = single(params.order_filt);

k_lap = single(10);


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
% Y1: term of SSTV
% Y2: term of spatial graph reg
% Y3: term of flatness of stripe noise
% Y4: term of data-fidelity
%
% Y1 = (D(Ds(U)))
% Y2 = Wsp.*Dsp(U)
% Y3 = Dv(T) 
% Y4 = U + S + T 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y1 = zeros([n1, n2, n3, 2], 'single', 'gpuArray');
Y2 = zeros([n1, n2, n3, 4], 'single', 'gpuArray');
% Y3 = zeros([n1, n2, n3], 'single', 'gpuArray');
Y4 = zeros([n1, n2, n3], 'single', 'gpuArray');


%% Setting operators
% Difference operators with Neumann boundary
Dsp     = @(z) D4_Neumann_GPU(z);
Dspt    = @(z) D4t_Neumann_GPU(z);
D       = @(z) cat(4, z([2:end, end],:,:) - z, z(:,[2:end, end],:) - z);
Dt      = @(z) cat(1, -z(1, :, :, 1), -z(2:end-1, :, :, 1) + z(1:end-2, :, :, 1), z(end-1, :, :, 1)) ...
                + cat(2, -z(:, 1, :, 2), -z(:, 2:end-1, :, 2) + z(:, 1:end-2, :, 2), z(:, end-1, :, 2));
% Dv      = @(z) z([2:end, end],:,:) - z;
% Dvt     = @(z) cat(1, -z(1, :, :), -z(2:(n1-1), :, :) + z(1:(n1-2), :, :), z(n1-1, :, :));
Dl      = @(z) z(:, :, [2:end, end], :) - z;
Dlt     = @(z) cat(3, -z(:, :, 1), -z(:, :, 2:n3-1) + z(:, :, 1:n3-2), z(:, :, n3-1));


%% Constructng graphs
% (A) Spatial Graph Weights
Wsp = Create_SpatialGraphWeight(HSI_clean, sigma_sp);
lambda_sp = sum(abs(Wsp.*Dsp(HSI_clean)), "all") * lambda_rho_sp;


% (B) Spectral Graph Laplacians
% returns L_delta: [K x K x S], B: [K x K x S], lam_max: [S x 1]
[L_delta_cpu, ~, lam_max_vec, info_l] = ...
    Create_SpectralGraphLaplacian(HSI_clean, num_segments, sigma_l, k_lap, order_filt);

L_delta = gpuArray(single(L_delta_cpu)); % Upload Laplacians to GPU
segID_gpu = gpuArray(int32(info_l.labels)); % [n1 x n2] Segment labels
S_seg = size(L_delta, 3);


lam_max_vec(~isfinite(lam_max_vec)) = 0;
max_eig_L = max(single(lam_max_vec(:)));


%% Setting stepsize parameters
% Calculate Max Lipschitz Constant for GLR term
% Lipschitz_GLR = || omega * D_lam^T * L * D_lam ||_op
%          <= omega * ||D_lam||^2 * ||L||
% ||D_lam||^2 <= 4, ||L|| = max(eigenvalues)
Lipschitz_GLR = 4 * lambda2 * max_eig_L;

spopnorm_Dvh = 8;
spopnorm_Dl = 4;
spopnorm_Wsp = 1;
spopnorm_Dsp = 16;
% spopnorm_Dv = 4;


gamma1_U    = gpuArray(single(1/((Lipschitz_GLR/2) + ...
                spopnorm_Dvh*spopnorm_Dl + spopnorm_Wsp*spopnorm_Dsp + 1)));
gamma1_S    = gpuArray(single(1));
% gamma1_T    = gpuArray(single(1/(spopnorm_Dv + 1)));
% gamma2      = gpuArray(single(1/3));
gamma2      = gpuArray(single(1/2));


%% main loop (P-PDS)
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
    grad_GLR = compute_GLR_gradient(U, L_delta, segID_gpu, lambda2);
    U_tmp   = U - gamma1_U*(grad_GLR + Dlt(Dt(Y1)) + Dspt(Wsp.*Y2) + Y4);
    U_next  = ProjBox(U_tmp, 0, 1);
    U_res   = 2*U_next - U; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating S
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_tmp   = S - gamma1_S*Y4;
    S_next  = ProjFastL1Ball(S_tmp, alpha);
    S_res   = 2*S_next - S;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating T
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % T_tmp   = T - gamma1_T*(Dvt(Y3) + Y4);
    % T_next  = ProjFastL1Ball(T_tmp, beta);
    % T_res   = 2*T_next - T;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y1_tmp  = Y1 + gamma2*D(Dl(U_res));
    Y1_next = Y1_tmp - gamma2*ProxL1norm(Y1_tmp/gamma2, 1/gamma2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y2_tmp  = Y2 + gamma2*Wsp.*(Dsp(U_res));
    Y2_next = Y2_tmp - gamma2*ProjFastL1Ball(Y2_tmp, lambda_sp);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Y3_next = Y3 + gamma2*Dv(T_res);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Y4_tmp  = Y4 + gamma2*(U_res + S_res + T_res); 
    Y4_tmp  = Y4 + gamma2*(U_res + S_res);
    Y4_next = Y4_tmp - gamma2*ProjL2ball(Y4_tmp/gamma2, HSI_noisy_gpu, epsilon);
    

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
    Y4  = Y4_next;


 
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

fprintf("~~~ P-PDS ENDS ~~~\n");

%% Organizing results for output
HSI_restored                    = gather(U);
output.iter                     = gather(i);
removed_noise.all_noise         = HSI_noisy - HSI_restored;
removed_noise.sparse_noise      = gather(S);
% removed_noise.stripe_noise      = gather(T);
% removed_noise.gaussian_noise    = HSI_noisy - HSI_restored - ...
%                                     removed_noise.sparse_noise - removed_noise.stripe_noise;
removed_noise.gaussian_noise    = HSI_noisy - HSI_restored - ...
                                    removed_noise.sparse_noise;

output.converge_rate_U        = gather(converge_rate_U(1:output.iter));
output.converge_rate_S        = gather(converge_rate_S(1:output.iter));
% output.converge_rate_T        = gather(converge_rate_T(1:output.iter));
output.converge_rate_N        = gather(converge_rate_N(1:output.iter));

output.move_mpsnr             = gather(move_mpsnr(1:output.iter));
output.move_mssim             = gather(move_mssim(1:output.iter));
output.running_time           = gather(running_time(1:output.iter));

output.l2ball                 = gather(l2ball(1:output.iter));

output.Wsp                    = gather(Wsp);
output.L_delta                = gather(L_delta);

end


%% Caluculating gradient of GLR
function grad = compute_GLR_gradient(U, L_delta, segID_gpu, omega)
    % 1. Compute spectral differences: V = D_lam(U)  [n1 x n2 x K]
    [n1, n2, n3] = size(U);
    K = n3-1;
    Dl_GLR  = @(z) z(:, :, 2:end) - z(:,:,1:end-1);
    Dlt_GLR  = @(z) cat(3, -z(:, :, 1), -z(:, :, 2:n3-1) + z(:, :, 1:n3-2), z(:, :, n3-1));

    V = Dl_GLR(U);
    
    % 2. Segment-wise Matrix Multiplication: W = L * V
    % We compute this by flattening spatial dimensions and iterating segments.
    % L_delta: [K x K x S]
    
    W = zeros(n1*n2, K, 'single', 'gpuArray');
    V_flat = reshape(V, n1*n2, K);
    segID_flat = reshape(segID_gpu, n1*n2, 1);
    
    S_seg = size(L_delta, 3);
    
    % Loop over segments (S is small, ~50)
    for s = 1:S_seg
        % Indices of pixels in segment s
        idx = (segID_flat == s);
        if ~any(idx), continue; end
        
        % Extract vectors: [N_s x K]
        v_s = V_flat(idx, :);
        
        % Multiply by Laplacian L^{(s)} [K x K]
        % Result [N_s x K] = [N_s x K] * [K x K]^T
        % Since L is symmetric, transpose is same, but dimensionally:
        % v_s is row vectors. w_s = v_s * L
        L_s = L_delta(:, :, s);
        w_s = v_s * L_s; % MATLAB * is matrix mul. (1xK)*(KxK) = (1xK)
        
        % Store back
        W(idx, :) = w_s;
    end
    
    % Reshape back to image
    W = reshape(W, n1, n2, K);
    
    % 3. Apply Adjoint Difference and Weight: grad = omega * D_lam^T(W)
    grad = omega * Dlt_GLR(W);
end