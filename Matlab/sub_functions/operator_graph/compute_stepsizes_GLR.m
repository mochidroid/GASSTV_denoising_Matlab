function [gamma1, gamma2] = compute_stepsizes_GLR(L_delta, lambda2)
% Compute step sizes for GASSTV+GLR (unified gamma1 / gamma2)
% L_delta : spectral graph Laplacian (n3-1 x n3-1)
% lambda2 : weight of GLR term

    % ---- safety factors ----
    safety_prim = 0.8;   % primal side safety
    safety_dual = 0.9;   % dual side safety

    % ---- operator norm upper-bounds ----
    Ds_sq     = 4;   % ||Ds||^2 (spectral diff)
    DvhDs_sq  = 16;  % ||D_vh * Ds||^2 (SSTV)
    Dsp_sq    = 16;  % ||D_sp||^2 (4-neighbor spatial)
    I_sq      = 1;   % identity
    Dv_sq     = 4;   % ||Dv||^2 (vertical diff)

    % ---- max eigenvalue of L_delta ----
    try
        lam_max = gather(eigs(gather(L_delta), 1, 'lm'));
    catch
        lam_max = max(abs(eig(gather(L_delta))));
    end
    lam_max = max(single(lam_max), 0);

    % ---- Lipschitz constant of GLR term ----
    L_glr = lambda2 * Ds_sq * lam_max;
    L_total = max(L_glr, single(1e-6)); % avoid 0

    % ---- unified gamma1 ----
    gamma1 = gpuArray(single(safety_prim / L_total));

    % ---- unified gamma2 ----
    A_norm_sq = [DvhDs_sq; Dsp_sq; I_sq; I_sq; I_sq; Dv_sq];  
    gamma2_candidates = 1 ./ (gather(gamma1) * A_norm_sq);
    gamma2 = gpuArray(single(safety_dual * min(gamma2_candidates)));

    fprintf('[stepsizes] lam_max(LΔ)=%.3g, L_glr=%.3g, gamma1=%.3g, gamma2=%.3g\n', ...
        double(lam_max), double(gather(L_glr)), double(gather(gamma1)), double(gather(gamma2)));
end
