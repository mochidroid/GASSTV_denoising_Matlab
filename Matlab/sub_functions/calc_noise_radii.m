function params = calc_noise_radii(params, deg, hsi, HSI_noisy)
% params.rho_radius から epsilon (ガウシアン), alpha (スパース), beta (ストライプ) の半径を計算する.
% deadline ノイズはスパースノイズの項に含める (exp_denoising_gstd_GeoSSTV.m, 4a66b09 と同じ式).
%   cr_dl : deadline の被覆率 (deadline の平均幅を 2 として 1 - exp(-2 * deadline_rate))
%   mu    : 観測の (0, 1) 内の画素値の平均 (deadline で 0 になった画素の元の値の目安)
r_sp = deg.sparse_rate;
dl_mean_width = 2;
cr_dl = 1 - exp(-dl_mean_width * deg.deadline_rate);

mask_valid = (HSI_noisy > 0) & (HSI_noisy < 1);
mu = mean(double(HSI_noisy(mask_valid)), "all");

params.epsilon = params.rho_radius * deg.gaussian_sigma * sqrt(hsi.N * (1 - r_sp) * (1 - cr_dl));
params.alpha = params.rho_radius * hsi.N * (0.5 * r_sp + mu * cr_dl);
params.beta = params.rho_radius * hsi.N * (1 - r_sp) * (1 - cr_dl) ...
    * deg.stripe_rate * deg.stripe_intensity / 2;
end
