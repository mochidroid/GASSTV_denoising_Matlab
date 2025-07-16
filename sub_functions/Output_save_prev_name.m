function [save_folder_name] = Output_save_prev_name(deg, image)

is_sparse       = deg.is_sparse;
is_stripe       = deg.is_stripe;


gaussian_sigma = deg.gaussian_sigma;

save_folder_name = "./result_TGRS/denoising_" + image + "/" + ...
    "g" + num2str(gaussian_sigma);

if is_sparse == 1
    sparse_rate = deg.sparse_rate;

    save_folder_name = save_folder_name + "_ps" + num2str(sparse_rate);
end

if is_stripe == 1
    stripe_rate = deg.stripe_rate;
    stripe_intensity = deg.stripe_intensity;

    save_folder_name = save_folder_name + ...
        "_pt" + num2str(stripe_rate) + "_tint" + num2str(stripe_intensity);
end


% rho = params.rho;
% save_folder_name = save_folder_name + "/r" + num2str(rho);

