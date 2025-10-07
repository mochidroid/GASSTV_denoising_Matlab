function [HSI_restored, removed_noise, other_result] ...
     = select_func_S3TTV_for_denoising_SPDHG(HSI_clean, HSI_noisy, params, deg)

addpath('./methods/S3TTV');

% Selecting SSTV function based on noise conditions
if deg.sparse_rate == 0 && deg.stripe_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_S3TTV_g_for_denoising_SPDHG(HSI_clean, HSI_noisy, params);

elseif deg.stripe_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_S3TTV_gs_for_denoising_SPDHG_v3(HSI_clean, HSI_noisy, params);

elseif deg.sparse_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_S3TTV_gt_for_denoising_SPDHG(HSI_clean, HSI_noisy, params);

else
    [HSI_restored, removed_noise, other_result] = ...
        func_S3TTV_gst_for_denoising_SPDHG_v4(HSI_clean, HSI_noisy, params);
end