close all;
%% Settings
num_scale = 5;
scale_width = 6;

types = 'parula';
% types = 'turbo';
% types = 'hot';
cmap = colormap(types);

num_ver_pixel = 256;
num_hor_pixel = 26;
max_colorbar = 1; % max colorbar
max_sign = 1; % max signal

is_save = 1;

%% Making colorbar matrix
colorbar_ind = reshape(num_ver_pixel:-1:1, [num_ver_pixel, 1]).*ones(num_ver_pixel, num_hor_pixel);
mycolorbar = ind2rgb(colorbar_ind, cmap);

% Enclosing in black
mycolorbar(1, :, :) = 0;
mycolorbar(end, :, :) = 0;
mycolorbar(:, 1, :) = 0;
mycolorbar(:, end, :) = 0;

%% Making scale of colorbar
pixel_start = floor((1 - max_sign/max_colorbar)*num_ver_pixel);

for i = 0:1:num_scale
    v_pix = floor(i*(num_ver_pixel - pixel_start)/(num_scale + 1) + pixel_start);

    if v_pix > 0
        mycolorbar(v_pix, end-scale_width:end, :) = 0;
    else
    end
end

%% Showing colorbar
figure;
imshow(mycolorbar)

%% Saving colorbar
if is_save
    save_folder_name = "./sub_functions/Save_image";
    mkdir(save_folder_name);
    imwrite(mycolorbar, save_folder_name + "/colorbar_" + types + "_b" + num_scale + ".png")
end