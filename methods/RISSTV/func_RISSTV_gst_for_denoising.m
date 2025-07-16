%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f(U,W1,W2,S,T) = \lambda|W1|_{1,2} + |W2|_{1,2} + 
%           (L*(W1) = D(U)) + (L*(W2) = D(Ds(U))) +
%           L2ball(U+S+T) + box constraint(U) + L1ball(S) + L1ball(T) + (Dv(T)=0)
%
% f1(U,W1,W2,S,T) = 0
% f2(U,W1,W2,S,T) = L2ball(U+S+T) + \lambda|W1|_{1,2} + |W2|_{1,2}              
% f3(U,W1,W2,S,T) = (L*(W1) = D(U)) + (L*(W2) = D(Ds(U))) + 
%                   box constraint(U) + L1ball(S) + L1ball(T) + (Dv(T)=0)
%                       
%
% A = (D -L* O O O; DDs O -L* O O; I O O O O; O O O I O; O O O O I; O O O O Dv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [HSI_restored, removed_noise, other_result] ...
     = func_RISSTV_gst_for_denoising(HSI_clean, HSI_noisy, params)
fprintf("** Running func_RISSTV_gst_for_denoising **\n");
HSI_clean = single(HSI_clean);
HSI_noisy  = single(HSI_noisy);
HSI_noisy_gpu = gpuArray(single(HSI_noisy));
[n1, n2, n3] = size(HSI_noisy);

epsilon     = gpuArray(single(params.epsilon));
alpha       = gpuArray(single(params.alpha));
beta        = gpuArray(single(params.beta));
lambda      = gpuArray(single(params.lambda));
maxiter     = gpuArray(single(params.maxiter));
stopcri     = gpuArray(single(params.stopcri));

%% Setting params
dispiter    = unique([1:10, 1000:1000:maxiter]);
dispband    = round(n3/2);


%% Setting operator
% Difference operators with Neumann boundary
D       = @(z) cat(4, z([2:end, end],:,:) - z, z(:,[2:end, end],:) - z);
Dt      = @(z) cat(1, -z(1, :, :, 1), -z(2:end-1, :, :, 1) + z(1:end-2, :, :, 1), z(end-1, :, :, 1)) ...
        + cat(2, -z(:, 1, :, 2), -z(:, 2:end-1, :, 2) + z(:, 1:end-2, :, 2), z(:, end-1, :, 2));

Dl      = @(z) z(:, :, [2:end, end]) - z;
Dlt     = @(z) cat(3, -z(:, :, 1), -z(:, :, 2:end-1) + z(:, :, 1:end-2), z(:, :, end-1));

Dv      = @(z) z([2:end, end], :, :) - z;
Dvt     = @(z) cat(1, -z(1, :, :), -z(2:(end-1), :, :) + z(1:(end-2), :, :), z(end-1, :, :));


%% Initializing primal and dual variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% primal variables
% U: clean HSI
% W1: interpolation (first-order difference)
% W2: interpolation (second-order difference)
% S: sparse noise(salt-and-pepper noise)
% T: stripe noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U  = zeros([n1, n2, n3], "single", "gpuArray");
W1 = zeros([n1, n2, n3, 2, 3], "single", "gpuArray");
W2 = zeros([n1, n2, n3, 2, 3], "single", "gpuArray");
S  = zeros([n1, n2, n3], "single", "gpuArray");
T  = zeros([n1, n2, n3], "single", "gpuArray");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dual variables
% Y1: term of connection with Du and L*w
% Y2: term of connection with DDsu and L*w 
% Y3: term of box constraint
% Y4: term of sparse noise 
% Y5: term of stripe noise 
% Y6: stripe noise
%
% Y1 = D(U) - L*(W1)
% Y2 = D(Ds(U)) - L*(W2) 
% Y3 = U
% Y4 = S
% Y5 = T
% Y6 = Dv(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y1 = zeros([n1, n2, n3, 2], "single", "gpuArray");
Y2 = zeros([n1, n2, n3, 2], "single", "gpuArray");
Y3 = zeros([n1, n2, n3], "single", "gpuArray");
Y4 = zeros([n1, n2, n3], "single", "gpuArray");
Y5 = zeros([n1, n2, n3], "single", "gpuArray");
Y6 = zeros([n1, n2, n3], "single", "gpuArray");


%% Setting step size for P-PDS
gamma1_U    = gpuArray(single(1/(8 + 8*4 + 1)));
gamma1_W1   = gpuArray(single(1/4));
gamma1_W2   = gpuArray(single(1/4));
gamma1_S    = gpuArray(single(1));
gamma1_T    = gpuArray(single(1/4));
gamma2      = gpuArray(single(1/5));


%% main loop (P-PDS)
fprintf("~~~ P-PDS STARTS ~~~\n");

converge_rate_U = zeros([1, maxiter], "single");
converge_rate_W1 = zeros([1, maxiter], "single");
converge_rate_W2 = zeros([1, maxiter], "single");
converge_rate_S = zeros([1, maxiter], "single");
converge_rate_T = zeros([1, maxiter], "single");
converge_rate_N = zeros([1, maxiter], "single");
move_function_value = zeros([1, maxiter], "single");
move_AomegaU = zeros([1, maxiter], "single");
move_mpsnr = zeros([1, maxiter], "single");
move_mssim = zeros([1, maxiter], "single");
running_time = zeros([1, maxiter], "single");
l2ball = zeros([1, maxiter], "single");


for i = 1:maxiter
    tic;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Primal Variables except W1 W2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U_tmp   = U - gamma1_U*(Dt(Y1) + Dlt(Dt(Y2)) + Y3);
    S_tmp   = S - gamma1_S*Y4;
    T_tmp   = T - gamma1_T*(Y5 + Dvt(Y6));

    Primal_sum = U_tmp + S_tmp + T_tmp;
    Primal_sum = ProjL2ball(Primal_sum, HSI_noisy_gpu, epsilon) - Primal_sum;

    U_next = U_tmp + Primal_sum/3;
    S_next = S_tmp + Primal_sum/3;
    T_next = T_tmp + Primal_sum/3;

    U_res = 2*U_next - U;
    S_res = 2*S_next - S;
    T_res = 2*T_next - T;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating W1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W1_tmp   = W1 + gamma1_W1*(opL(Y1));
    W1_next  = Prox_l12norm_d4(W1_tmp, lambda*gamma1_W1);
    W1_res   = 2*W1_next - W1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating W2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W2_tmp   = W2 + gamma1_W2*(opL(Y2));
    W2_next  = Prox_l12norm_d4(W2_tmp, gamma1_W2);
    W2_res   = 2*W2_next - W2;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y1_next  = Y1 + gamma2*(D(U_res)-opLadj(W1_res));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y2_next  = Y2 + gamma2*(D(Dl(U_res))-opLadj(W2_res));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y3_tmp  = Y3 + gamma2*U_res;
    Y3_next = Y3_tmp - gamma2*ProjBox(Y3_tmp/gamma2, 0, 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y4_tmp  = Y4 + gamma2*S_res;
    Y4_next = Y4_tmp - gamma2*ProjFastL1Ball(Y4_tmp/gamma2, alpha);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y5_tmp  = Y5 + gamma2*T_res;
    Y5_next = Y5_tmp - gamma2*ProjFastL1Ball(Y5_tmp/gamma2, beta);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y6_next = Y6 + gamma2*Dv(T_res);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculating error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = HSI_noisy_gpu - U - S - T;
    N_next = HSI_noisy_gpu - U_next - S_next - T_next;

    converge_rate_U(i) = norm(U_next(:) - U(:),2)/norm(U(:),2);
    converge_rate_W1(i) = norm(W1_next(:) - W1(:),2)/norm(W1(:),2);
    converge_rate_W2(i) = norm(W2_next(:) - W2(:),2)/norm(W2(:),2);
    converge_rate_S(i) = norm(S_next(:) - S(:),2)/norm(S(:),2);
    converge_rate_T(i) = norm(T_next(:) - T(:),2)/norm(T(:),2);
    converge_rate_N(i) = norm(N_next(:) - N(:),2)/norm(N(:),2);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating all variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U   = U_next;
    W1  = W1_next;
    W2  = W2_next;
    S   = S_next;
    T   = T_next;
    
    Y1  = Y1_next;
    Y2  = Y2_next;
    Y3  = Y3_next;
    Y4  = Y4_next;
    Y5  = Y5_next;
    Y6  = Y6_next;

 
    % Saving results per iter
    running_time(i) = toc;

    move_function_value(i) = lambda*sum(sqrt(sum(W1.^2,4)), "all") + sum(sqrt(sum(W2.^2,4)), "all");
    move_AomegaU(i) = lambda*sum(abs(D(U)), "all") + sum(abs(D(Dl(U))), "all");
    move_mpsnr(i) = calc_MPSNR(gather(U), HSI_clean);
    move_mssim(i) = calc_MSSIM(gather(U), HSI_clean);

    l2ball(i) = norm(gather(N(:)), 2);

    if i>=2 && converge_rate_U(i) < stopcri
        break
    end

    % Displaying progress 
    if ismember(i, dispiter)
        fprintf("Iter: %d, Error: %0.6f, FV: %#.4g, MPSNR: %#.4g, MSSIM: %#.4g, Time: %0.2f.\n", ...
            i, converge_rate_U(i), move_function_value(i), move_mpsnr(i), move_mssim(i), sum(running_time));
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
HSI_restored                        = gather(U);
other_result.interpolation1         = gather(W1);
other_result.interpolation2         = gather(W2);
other_result.iteration              = gather(i);
removed_noise.sparse_noise          = gather(S);
removed_noise.stripe_noise          = gather(T);
removed_noise.gaussian_noise        = gather(HSI_noisy_gpu - U - S - T);
removed_noise.all_noise             = gather(HSI_noisy_gpu - U);

other_result.converge_rate_U        = gather(converge_rate_U(1:other_result.iteration));
other_result.converge_rate_W1       = gather(converge_rate_W1(1:other_result.iteration));
other_result.converge_rate_W2       = gather(converge_rate_W2(1:other_result.iteration));
other_result.converge_rate_S        = gather(converge_rate_S(1:other_result.iteration));
other_result.converge_rate_T        = gather(converge_rate_T(1:other_result.iteration));
other_result.converge_rate_N        = gather(converge_rate_N(1:other_result.iteration));

other_result.move_function_value    = gather(move_function_value(1:other_result.iteration));
other_result.move_AomegaU           = gather(move_AomegaU(1:other_result.iteration));
other_result.move_mpsnr             = gather(move_mpsnr(1:other_result.iteration));
other_result.move_mssim             = gather(move_mssim(1:other_result.iteration));
other_result.running_time           = gather(running_time(1:other_result.iteration));

other_result.l2ball                 = gather(l2ball(1:other_result.iteration));

end