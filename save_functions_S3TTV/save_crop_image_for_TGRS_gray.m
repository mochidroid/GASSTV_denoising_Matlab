close all
addpath(genpath('sub_functions'))

%% Selecting case of experiment
expt_case = 'Case1';
% expt_case = 'Case2';
% expt_case = 'Case3';
% expt_case = 'Case6';

%% Setting parameters
param_multiple = 1;
% param_multiple = 1.5;

param_diff_multiple = 8;
% colormap = turbo;
% colormap = jet;
% colormap = parula;
colormap = hot;

is_load_image           = 1;
is_show_evaluation      = 1;
is_show_all             = 1;
is_show_restored_image  = 0;
is_show_diff            = 0;
is_save_image           = 0;
is_save_diff            = 0;


switch expt_case
    case 'Case1' %Jasper Ridge g0.05 ps0.05
        % start_pos = [52, 85];
        start_pos = [60, 79];
        % start_pos = [36, 8];
        crop_size = [20, 20];
        expansion_rate = 2;
        embed_tblr = 'bl';
        % embed_tblr = 'br';

        % save_band = 73;
        save_band = 131;

    case 'Case2' % Jasper Ridge g0.1 ps0.05
        % start_pos = [38, 53];
        start_pos = [60, 79];
        crop_size = [20, 20];
        expansion_rate = 2;
        embed_tblr = 'bl';

        % save_band = 43;
        save_band = 139;
        % save_band = 52;
        % save_band = 59;
        % save_band = 102;

    case 'Case3' %PaviaU g0.05 pt0.05
        % start_pos = [95, 13];
        % start_pos = [78, 3];
        start_pos = [15, 60];
        crop_size = [20, 20];
        expansion_rate = 2;
        % embed_tblr = 'br';
        embed_tblr = 'bl';

        % param_diff_multiple = 10;

        % save_band = 41; %PaviaU
        save_band = 58;

        % save_band = 46;

    case 'Case6' % g0.1 ps0.05 pt0.05
        % start_pos = [15, 54];
        start_pos = [13, 56];
        crop_size = [20, 20];
        expansion_rate = 2;
        embed_tblr = 'bl';

        save_band = 59; %Jasper Ridge1
        % save_band = 129; %Jasper Ridge2
end

%% Loading result
if is_load_image
    switch expt_case
        case 'Case1'
            deg.gaussian_sigma      = 0.05;
            deg.sparse_rate         = 0.05;
            
            image = 'JasperRidge';
    
            save_folder_name = "./result/denoising_" + image + "/" + ...
                "g" + num2str(deg.gaussian_sigma) + ...
                "_ps" + num2str(deg.sparse_rate);

        case 'Case2'
            deg.gaussian_sigma      = 0.1;
            deg.sparse_rate         = 0.05;
            
            image = 'JasperRidge';
    
            save_folder_name = "./result/denoising_" + image + "/" + ...
                "g" + num2str(deg.gaussian_sigma) + ...
                "_ps" + num2str(deg.sparse_rate);
    
        case 'Case3'
            deg.gaussian_sigma      = 0.05;
            deg.stripe_rate         = 0.05;
            deg.stripe_intensity    = 0.5;
            
            image = 'PaviaU120';
    
            save_folder_name = "./result/denoising_" + image + "/" + ...
                "g" + num2str(deg.gaussian_sigma) + ...
                "_pt" + num2str(deg.stripe_rate) + ...
                "_tint" + num2str(deg.stripe_intensity);
    
        case 'Case6'
            deg.gaussian_sigma      = 0.1;
            deg.sparse_rate         = 0.05;
            deg.stripe_rate         = 0.05;
            deg.stripe_intensity    = 0.5;
            
            image = 'JasperRidge';
    
            save_folder_name = "./result/denoising_" + image + "/" + ...
                "g" + num2str(deg.gaussian_sigma) + ...
                "_ps" + num2str(deg.sparse_rate) + ...
                "_pt" + num2str(deg.stripe_rate) + ...
                "_tint" + num2str(deg.stripe_intensity);
    end
    
    load(save_folder_name + "/S3TTV_cir_PsPPDS_" + image + "_bl10_st1_r0.95_stop1e-5")
    
    load(save_folder_name + "/SSTV_cir_PsPPDS_" + image + "_r0.95_stop1e-5")
    load(save_folder_name + "/HSSTV_L1_cir_PsPPDS_" + image + "_o0.05_r0.95_stop1e-5")
    load(save_folder_name + "/HSSTV_L12_cir_PsPPDS_" + image + "_o0.05_r0.95_stop1e-5")
    load(save_folder_name + "/l0l1HTV_cir_PsPPDS_" + image + "_sr0.999_r0.95_maxiter20000")
    load(save_folder_name + "/TPTV_" + image)
    load(save_folder_name + "/STV_cir_PsPPDS_" + image + "_bl10_st1_r0.95_stop1e-5")
    load(save_folder_name + "/SSST_cir_PsPPDS_" + image + "_bl10_st1_r0.95_stop1e-5")
    load(save_folder_name + "/LRTDTV_" + image + "_stop1e-5")
end

%% Showing evaluation

if is_show_evaluation
    % Calucalting evaluation of noisy image
    mpsnr_noisy  = MPSNR(HSI_noisy, HSI_clean);
    mssim_noisy  = MSSIM(HSI_noisy, HSI_clean);
    ergas_noisy  = ERGAS(HSI_noisy, HSI_clean, 1); % GSD ratio = 1 for recovery problem;
    sam_noisy    = SAM(HSI_noisy, HSI_clean);

    fprintf('~~~ RESULTS ~~~\n');
    fprintf('\t\t\t MPSNR\t MSSIM\t ERGAS\t SAM\n')
    
    fprintf('noisy:\t\t %#0.4g\t %#0.4g\t %#0.4g\t %#0.4g\n', ...
        mpsnr_noisy, mssim_noisy, ergas_noisy, sam_noisy);
    fprintf('SSTV:\t\t %#0.4g\t %#0.4g\t %#0.4g\t %#0.4g\n', ...
        mpsnr_SSTV, mssim_SSTV, ergas_SSTV, sam_SSTV);
    fprintf('HSSTV L1:\t %#0.4g\t %#0.4g\t %#0.4g\t %#0.4g\n', ...
        mpsnr_HSSTV_L1, mssim_HSSTV_L1, ergas_HSSTV_L1, sam_HSSTV_L1);
    fprintf('HSSTV_L12:\t %#0.4g\t %#0.4g\t %#0.4g\t %#0.4g\n', ...
        mpsnr_HSSTV_L12, mssim_HSSTV_L12, ergas_HSSTV_L12, sam_HSSTV_L12);
    fprintf('l0l1HTV:\t %#0.4g\t %#0.4g\t %#0.4g\t %#0.4g\n', ...
        mpsnr_l0l1HTV, mssim_l0l1HTV, ergas_l0l1HTV, sam_l0l1HTV);
    fprintf('TPTV:\t\t %#0.4g\t %#0.4g\t %#0.4g\t %#0.4g\n', ...
        mpsnr_TPTV, mssim_TPTV, ergas_TPTV, sam_TPTV);
    fprintf('STV:\t\t %#0.4g\t %#0.4g\t %#0.4g\t %#0.4g\n', ...
        mpsnr_STV, mssim_STV, ergas_STV, sam_STV);
    fprintf('SSST:\t\t %#0.4g\t %#0.4g\t %#0.4g\t %#0.4g\n', ...
        mpsnr_SSST, mssim_SSST, ergas_SSST, sam_SSST);
    fprintf('LRTDTV:\t\t %#0.4g\t %#0.4g\t %#0.4g\t %#0.4g\n', ...
        mpsnr_LRTDTV, mssim_LRTDTV, ergas_LRTDTV, sam_LRTDTV);
    fprintf('S3TTV:\t\t %#0.4g\t %#0.4g\t %#0.4g\t %#0.4g\n', ...
        mpsnr_S3TTV, mssim_S3TTV, ergas_S3TTV, sam_S3TTV);
end

%% Cropping images
image_clean = HSI_clean(:,:,save_band)*param_multiple;

image_noisy = HSI_noisy(:,:,save_band)*param_multiple;
diff_image_noisy = abs(image_noisy-image_clean);
diff_image_crop_noisy = Crop_Embed_image(diff_image_noisy, ...
        start_pos, crop_size, expansion_rate, embed_tblr);

image_SSTV = SSTV_restored_HSI(:,:,save_band)*param_multiple;
diff_image_SSTV = abs(image_SSTV-image_clean)*param_diff_multiple;
diff_image_crop_SSTV = Crop_Embed_image(diff_image_SSTV, ...
    start_pos, crop_size, expansion_rate, embed_tblr);

image_HSSTV_L1 = HSSTV_L1_restored_HSI(:,:,save_band)*param_multiple;
diff_image_HSSTV_L1 = abs(image_HSSTV_L1-image_clean)*param_diff_multiple;
diff_image_crop_HSSTV_L1 = Crop_Embed_image(diff_image_HSSTV_L1, ...
    start_pos, crop_size, expansion_rate, embed_tblr);

image_HSSTV_L12 = HSSTV_L12_restored_HSI(:,:,save_band)*param_multiple;
diff_image_HSSTV_L12 = abs(image_HSSTV_L12-image_clean)*param_diff_multiple;
diff_image_crop_HSSTV_L12 = Crop_Embed_image(diff_image_HSSTV_L12, ...
    start_pos, crop_size, expansion_rate, embed_tblr);

image_l0l1HTV = l0l1HTV_restored_HSI(:,:,save_band)*param_multiple;
diff_image_l0l1HTV = abs(image_l0l1HTV-image_clean)*param_diff_multiple;
diff_image_crop_l0l1HTV = Crop_Embed_image(diff_image_l0l1HTV, ...
    start_pos, crop_size, expansion_rate, embed_tblr);

image_TPTV = TPTV_restored_HSI(:,:,save_band)*param_multiple;
diff_image_TPTV = abs(image_TPTV-image_clean)*param_diff_multiple;
diff_image_crop_TPTV = Crop_Embed_image(diff_image_TPTV, ...
    start_pos, crop_size, expansion_rate, embed_tblr);

image_STV = STV_restored_HSI(:,:,save_band)*param_multiple;
diff_image_STV = abs(image_STV-image_clean)*param_diff_multiple;
diff_image_crop_STV = Crop_Embed_image(diff_image_STV, ...
    start_pos, crop_size, expansion_rate, embed_tblr);

image_SSST = SSST_restored_HSI(:,:,save_band)*param_multiple;
diff_image_SSST = abs(image_SSST-image_clean)*param_diff_multiple;
diff_image_crop_SSST = Crop_Embed_image(diff_image_SSST, ...
    start_pos, crop_size, expansion_rate, embed_tblr);

image_LRTDTV = LRTDTV_restored_HSI(:,:,save_band)*param_multiple;
diff_image_LRTDTV = abs(image_LRTDTV-image_clean)*param_diff_multiple;
diff_image_crop_LRTDTV = Crop_Embed_image(diff_image_LRTDTV, ...
    start_pos, crop_size, expansion_rate, embed_tblr);

image_S3TTV = S3TTV_restored_HSI(:,:,save_band)*param_multiple;
diff_image_S3TTV = abs(image_S3TTV-image_clean)*param_diff_multiple;
diff_image_crop_S3TTV = Crop_Embed_image(diff_image_S3TTV, ...
    start_pos, crop_size, expansion_rate, embed_tblr);


cat_restored_image = cat(2, image_clean, image_noisy, ...
    image_SSTV, image_HSSTV_L1, image_HSSTV_L12, ...
    image_l0l1HTV, image_TPTV, image_STV, ...
    image_SSST, image_LRTDTV, image_S3TTV);

cat_diff = cat(2, zeros([hsi.n1, hsi.n2]), diff_image_noisy, ...
    diff_image_SSTV, diff_image_HSSTV_L1, diff_image_HSSTV_L12, ...
    diff_image_l0l1HTV, diff_image_TPTV, diff_image_STV, ...
    diff_image_SSST, diff_image_LRTDTV, diff_image_S3TTV);

cat_diff_crop = cat(2, zeros([hsi.n1, hsi.n2]), diff_image_crop_noisy, ...
    diff_image_crop_SSTV, diff_image_crop_HSSTV_L1, diff_image_crop_HSSTV_L12, ...
    diff_image_crop_l0l1HTV, diff_image_crop_TPTV, diff_image_crop_STV, ....
    diff_image_crop_SSST, diff_image_crop_LRTDTV, diff_image_crop_S3TTV);


%% Showing images
if is_show_all
    cat_image = cat(1, cat_restored_image, cat_diff, cat_diff_crop);
    
    figure;
    imshow(cat_image);
end

if is_show_restored_image
    figure;
    imshow(cat_restored_image);
end

if is_show_diff
    figure;
    % imshow(cat_diff, colormap=turbo)
    % imshow(cat_diff, colormap=jet)
    imshow(cat_diff, colormap=parula)
end


%% Adding arrow
% imshow(diff_image_TPTV, 'colormap', colormap, 'Border','tight')
% hold on
% 
% % 矢印の始点と終点の座標を指定
% startX = 60; % 矢印の始点の x 座標
% startY = 53; % 矢印の始点の y 座標
% lenX = 0;   % 矢印の終点の x 座標
% lenY = 20;   % 矢印の終点の y 座標
% 
% % 矢印を描画
% % quiver(startY, startX, endY - startY, endX - startX, 'w', 'LineWidth', 2);
% scaleFactor = 1;  % 矢印のスケールファクター（サイズを調整するための係数）
% maxHeadSize = 20;
% quiver(startY, startX, lenY, lenX, scaleFactor, 'w', 'LineWidth', 4, 'MaxHeadSize', maxHeadSize);
% 
% % quiver(startX, startY, endX - startX, endY - startY, scaleFactor, 'Color', 'red', 'LineWidth', 2);
% 
% hold off; % 描画の終了
% 
% SaveFigPDF(TPTV_figure,"test.pdf")
% 
% % この下の行のコメントアウトを解除するとより余白を削ったPDFを作成します。
% % set(gca, 'LooseInset', get(gca, 'TightInset'));

%% Saving images
if is_save_image
    save_restored_image_name = "./result_image_for_TGRS/" + ...
        expt_case + "_" + image + "_restored_image" + ...
        "_b" + num2str(save_band) + "_m" + num2str(param_multiple);
    mkdir(save_restored_image_name);

    imwrite(image_clean, save_restored_image_name + "/image_clean.png");
    imwrite(image_noisy, save_restored_image_name + "/image_noisy.png");
    imwrite(image_SSTV, save_restored_image_name + "/image_SSTV.png");
    imwrite(image_HSSTV_L1, save_restored_image_name + "/image_HSSTV_L1.png");
    imwrite(image_HSSTV_L12, save_restored_image_name + "/image_HSSTV_L12.png");
    imwrite(image_l0l1HTV, save_restored_image_name + "/image_l0l1HTV.png");
    imwrite(image_TPTV, save_restored_image_name + "/image_TPTV.png");
    imwrite(image_STV, save_restored_image_name + "/image_STV.png");
    imwrite(image_SSST, save_restored_image_name + "/image_SSST.png");
    imwrite(image_LRTDTV, save_restored_image_name + "/image_LRTDTV.png");
    imwrite(image_S3TTV, save_restored_image_name + "/image_S3TTV.png");
end

%% Saving images
if is_save_diff
    save_diff_image_name = "./result_image_for_TGRS/" + ...
        expt_case + "_" + image + "_diff_image" + ...
        "_b" + num2str(save_band) + "_m" + num2str(param_diff_multiple);
    mkdir(save_diff_image_name);

    imwrite(diff_image_crop_noisy, save_diff_image_name + "/diff_image_noisy.png");
    imwrite(diff_image_crop_SSTV, save_diff_image_name + "/diff_image_SSTV.png");
    imwrite(diff_image_crop_HSSTV_L1, save_diff_image_name + "/diff_image_HSSTV_L1.png");
    imwrite(diff_image_crop_HSSTV_L12, save_diff_image_name + "/diff_image_HSSTV_L12.png");
    imwrite(diff_image_crop_l0l1HTV, save_diff_image_name + "/diff_image_l0l1HTV.png");
    imwrite(diff_image_crop_LRTDTV, save_diff_image_name + "/diff_image_LRTDTV.png");
    imwrite(diff_image_crop_TPTV, save_diff_image_name + "/diff_image_TPTV.png");
    imwrite(diff_image_crop_STV, save_diff_image_name + "/diff_image_STV.png");
    imwrite(diff_image_crop_SSST, save_diff_image_name + "/diff_image_SSST.png");
    imwrite(diff_image_crop_S3TTV, save_diff_image_name + "/diff_image_S3TTV.png");
end


%%
my_map = hot;
I = randn(256, 256);
I = I/max(I, [], "all");
I_ind = uint8(255*I);
I_map = ind2rgb(I_ind, my_map);