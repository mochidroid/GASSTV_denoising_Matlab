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
     = func_GASSTV_OraGuide_gs_for_denoising_by_NsPPDS(HSI_clean, HSI_noisy, params)
fprintf('** Running func_GASSTV_OraGuide_gs_for_denoising_by_NsPPDS **\n');
HSI_clean = single(HSI_clean);
HSI_noisy  = single(HSI_noisy);
HSI_noisy_gpu = gpuArray(single(HSI_noisy));
[n1, n2, n3] = size(HSI_noisy);

epsilon         = gpuArray(single(params.epsilon));
alpha           = gpuArray(single(params.alpha));
% beta            = gpuArray(single(params.beta));
sigma_sp        = gpuArray(single(params.sigma_sp));
sigma_s         = single(params.sigma_s);
lambda1         = gpuArray(single(params.lambda1));
lambda2         = gpuArray(single(params.lambda2));
num_segments    = single(params.num_segments);
maxiter         = gpuArray(single(params.maxiter));
stopcri         = gpuArray(single(params.stopcri));

L               = params.L;

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
% Y3: term of spectral graph reg
% Y4: term of box constraint
% Y5: term of sparse noise 
% Y6: term of stripe noise 
% Y7: term of flatness of stripe noise
%
% Y1 = (D(Ds(U)))
% Y2 = Wsp.*Dsp(U)
% Y3 = Ws.*Ds(U)
% Y4 = U
% Y5 = S
% Y6 = T
% Y7 = Dv(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y1 = zeros([n1, n2, n3, 2], 'single', 'gpuArray');
Y2 = zeros([n1, n2, n3, 4], 'single', 'gpuArray');
Y3 = zeros([n1, n2, n3], 'single', 'gpuArray');
Y4 = zeros([n1, n2, n3], 'single', 'gpuArray');
Y5 = zeros([n1, n2, n3], 'single', 'gpuArray');
% Y6 = zeros([n1, n2, n3], 'single', 'gpuArray');
% Y7 = zeros([n1, n2, n3], 'single', 'gpuArray');


%% Setting stepsize parameters
gamma1_U    = gpuArray(single(1/(8*4 + 16 + 4 + 1)));
% gamma1_U    = gpuArray(single(1/(8*4 + 4 + 1)));
gamma1_S    = gpuArray(single(1));
% gamma1_T    = gpuArray(single(1/(4 + 1)));
% gamma2      = gpuArray(single(1/3));
gamma2      = gpuArray(single(1/2));


%% Setting operators
% Difference operators with Neumann boundary
Dsp     = @(z) D4_Neumann_GPU(z);
Dspt    = @(z) D4t_Neumann_GPU(z);
D       = @(z) cat(4, z([2:end, end],:,:) - z, z(:,[2:end, end],:) - z);
Dt      = @(z) cat(1, -z(1, :, :, 1), -z(2:end-1, :, :, 1) + z(1:end-2, :, :, 1), z(end-1, :, :, 1)) ...
                + cat(2, -z(:, 1, :, 2), -z(:, 2:end-1, :, 2) + z(:, 1:end-2, :, 2), z(:, end-1, :, 2));
% Dv      = @(z) z([2:end, end],:,:) - z;
% Dvt     = @(z) cat(1, -z(1, :, :), -z(2:(n1-1), :, :) + z(1:(n1-2), :, :), z(n1-1, :, :));
Ds      = @(z) z(:, :, [2:end, end], :) - z;
Dst     = @(z) cat(3, -z(:, :, 1), -z(:, :, 2:n3-1) + z(:, :, 1:n3-2), z(:, :, n3-1));


% Graph based weight matrix
% Wsp = Create_GraphWeight_OnlyLumi(HSI_noisy, sigma_sp);
Wsp = Create_GraphWeight_OnlyLumi(HSI_clean, sigma_sp);
Ws = Create_GraphWeight_Spectral(HSI_clean, sigma_s, num_segments);

switch L
    case 'L1' % p = 1, prox of L1 norm
        Prox_Y2 = @(z) ProxL1norm(z/gamma2, lambda1/gamma2);
        Prox_Y3 = @(z) ProxL1norm(z/gamma2, lambda2/gamma2);
    case 'L2' % p = 2, prox of L2 norm
        Prox_Y2 = @(z) ProxL2norm(z/gamma2, lambda1/gamma2);
        Prox_Y3 = @(z) ProxL2norm(z/gamma2, lambda2/gamma2);
end


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
    % Updating Primal Variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U_tmp   = U - gamma1_U*(Dst(Dt(Y1)) + Dspt(Wsp.*Y2) + Dst(Ws.*Y3) + Y4);
    S_tmp   = S - gamma1_S*Y5;
    % T_tmp   = T - gamma1_T.*(Y6 + Dvt(Y7));

    % Primal_sum = U_tmp + S_tmp + T_tmp;
    Primal_sum = U_tmp + S_tmp;
    Primal_sum = ProjL2ball(Primal_sum, HSI_noisy_gpu, epsilon) - Primal_sum;
    % 
    % U_next = U_tmp + Primal_sum/3;
    % S_next = S_tmp + Primal_sum/3;
    % T_next = T_tmp + Primal_sum/3;
    U_next = U_tmp + Primal_sum/2;
    S_next = S_tmp + Primal_sum/2;

    U_res = 2*U_next - U;
    S_res = 2*S_next - S;
    % T_res = 2*T_next - T;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y1_tmp  = Y1 + gamma2*D(Ds(U_res));
    Y1_next = Y1_tmp - gamma2*ProxL1norm(Y1_tmp/gamma2, 1/gamma2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y2_tmp  = Y2 + gamma2*Wsp.*(Dsp(U_res));
    Y2_next = Y2_tmp - gamma2*Prox_Y2(Y2_tmp);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y3_tmp  = Y3 + gamma2*Ws.*(Ds(U_res));
    Y3_next = Y3_tmp - gamma2*Prox_Y3(Y3_tmp);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y4_tmp  = Y4 + gamma2*U_res;
    Y4_next = Y4_tmp - gamma2*ProjBox(Y4_tmp/gamma2, 0, 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y5_tmp  = Y5 + gamma2*S_res;
    Y5_next = Y5_tmp - gamma2*ProjFastL1Ball(Y5_tmp/gamma2, alpha);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Y6_tmp  = Y6 + gamma2*T_res;
    % Y6_next = Y6_tmp - gamma2*ProjFastL1Ball(Y6_tmp/gamma2, beta);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y7
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Y7_next = Y7 + gamma2*Dv(T_res);

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
    Y3  = Y3_next;
    Y4  = Y4_next;
    Y5  = Y5_next;
    % Y6  = Y6_next;
    % Y7  = Y7_next;

 
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
        semilogy(move_function_value(1:i));
        title("Function value")
        
        subplot(2,3,6)
        semilogy(l2ball(1:i));
        title("L2ball")
        drawnow
    end
end

fprintf("Iter: %d, Error: %0.6f, FV: %#.4g, MPSNR: %#.4g, MSSIM: %#.4g, Time: %0.2f.\n", ...
    i, converge_rate_U(i), move_function_value(i), move_mpsnr(i), move_mssim(i), sum(running_time));

fprintf("~~~ P-PDS ENDS ~~~\n");

%% Organizing results for output
HSI_restored                    = gather(U);
output.iter                     = gather(i);
removed_noise.all_noise         = HSI_noisy - HSI_restored;
removed_noise.sparse_noise      = gather(S);
% removed_noise.stripe_noise      = gather(T);
% removed_noise.gaussian_noise    = HSI_noisy - HSI_restored - ...
%                                     removed_noise.sparse_noise - removed_noise.stripe_noise;
removed_noise.gaussian_noise    = HSI_noisy - HSI_restored - removed_noise.sparse_noise;

output.converge_rate_U        = gather(converge_rate_U(1:output.iter));
output.converge_rate_S        = gather(converge_rate_S(1:output.iter));
% output.converge_rate_T        = gather(converge_rate_T(1:output.iter));
output.converge_rate_N        = gather(converge_rate_N(1:output.iter));

output.move_mpsnr             = gather(move_mpsnr(1:output.iter));
output.move_mssim             = gather(move_mssim(1:output.iter));
output.running_time           = gather(running_time(1:output.iter));

output.l2ball                 = gather(l2ball(1:output.iter));

output.Wsp                    = gather(Wsp);
output.Ws                     = gather(Ws);