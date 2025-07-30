function [HSI_restored, removed_noise, other_result] = func_FastHyMix(~, HSI_noisy, params, ~)
addpath('./methods/HSI-MixedNoiseRemoval-FastHyMix-main/scripts');
k_subspace = params.k_subspace;

img_noisy = HSI_noisy;

% Saving current directory
dir_prev = pwd;

cd('./methods/HSI-MixedNoiseRemoval-FastHyMix-main/')
try 
    [img_FastHyMix, M, noise_std_Gaussion] = FastHyMix(img_noisy,  k_subspace);
    
    % Organizing results for output
    HSI_restored = normalize01(img_FastHyMix);
    % max(HSI_restored, [], "all")
    % min(HSI_restored, [], "all")
    
    removed_noise.all_noise = HSI_noisy - HSI_restored;
    
    other_result.M = M;
    other_result.noise_std_Gaussion = noise_std_Gaussion;
catch
    HSI_restored = img_noisy;
    fprintf('Error\n');
    removed_noise = NaN;
    other_result = NaN;
end
    

% Returning current directory
cd(dir_prev)

end

