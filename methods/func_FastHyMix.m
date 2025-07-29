function [HSI_restored, removed_noise, other_result] = func_FastHyMix(~, HSI_noisy, params, ~)
addpath('./methods/HSI-MixedNoiseRemoval-FastHyMix-main/scripts');
k_subspace = params.k_subspace;

img_noisy = HSI_noisy;


dir_prev = pwd;

cd('./methods/HSI-MixedNoiseRemoval-FastHyMix-main/')
[img_FastHyMix, M, noise_std_Gaussion] = FastHyMix(img_noisy,  k_subspace);

cd(dir_prev)

%% Organizing results for output
HSI_restored = normalize01(img_FastHyMix);
% max(HSI_restored, [], "all")
% min(HSI_restored, [], "all")

removed_noise.all_noise = HSI_noisy - HSI_restored;

other_result.M = M;
other_result.noise_std_Gaussion = noise_std_Gaussion;
end

