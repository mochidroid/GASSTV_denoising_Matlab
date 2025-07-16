function  [HSI_noisy, hsi] = Load_real_HSI(set_image)
%% Loading HSI
switch set_image
    case 'IndianPines'
        load('H:/マイドライブ/MATLAB_Shere/HSIData/IndianPine/Indian_pines.mat')
        HSI_noisy = normalize01(indian_pines(1:end-25, 26:end, :));
        
    case 'Suwannee'
        load('H:/マイドライブ/MATLAB_Shere/HSIData/Suwannee_original.mat')
        start_pos = [401, 220];
        U_tmp = u_org(start_pos(1):start_pos(1)+99, start_pos(2):start_pos(2)+99, :);
        HSI_noisy = normalize01_by_band(U_tmp);


    case 'Urban'
        load('H:/マイドライブ/MATLAB_Shere/HSIData/Urban/Urban_R162.mat');
        U_tmp = reshape(Y, [162, nCol, nRow]);
        HSI_noisy = normalize01(permute(U_tmp, [2,3,1]));

    case 'IndianPines_v1'
        load('H:/マイドライブ/MATLAB_Shere/HSIData/IndianPine/Indian_pines.mat')
        HSI_noisy = normalize01(indian_pines);

    case 'IndianPines120'
        load('H:/マイドライブ/MATLAB_Shere/HSIData/IndianPine/Indian_pines.mat')
        HSI_noisy = normalize01(indian_pines(1:end-25, 26:end, :));

    case 'Swannee_v1'
        load('H:/マイドライブ/MATLAB_Shere/HSIData/Suwannee_original.mat')
        U_tmp = u_org(641:740, 116:215, :);
        HSI_noisy = normalize01_by_band(U_tmp);

    case 'Swannee_251_550'
        load('H:/マイドライブ/MATLAB_Shere/HSIData/Suwannee_original.mat')
        U_tmp = u_org(251:550, :, :);
        HSI_noisy = normalize01_by_band(U_tmp);

    case 'Swannee_401_201'
        load('H:/マイドライブ/MATLAB_Shere/HSIData/Suwannee_original.mat')
        U_tmp = u_org(401:500, 201:300, :);
        HSI_noisy = normalize01_by_band(U_tmp);

    case 'Swannee_261_221'
        load('H:/マイドライブ/MATLAB_Shere/HSIData/Suwannee_original.mat')
        start_pos = [261, 221];
        U_tmp = u_org(start_pos(1):start_pos(1)+99, start_pos(2):start_pos(2)+99, :);
        HSI_noisy = normalize01_by_band(U_tmp);

    case 'Swannee_401_220'
        load('H:/マイドライブ/MATLAB_Shere/HSIData/Suwannee_original.mat')
        start_pos = [401, 220];
        U_tmp = u_org(start_pos(1):start_pos(1)+99, start_pos(2):start_pos(2)+99, :);
        HSI_noisy = normalize01_by_band(U_tmp);
end

hsi.sizeof = size(HSI_noisy);
[hsi.n1,hsi.n2,hsi.n3] = size(HSI_noisy);
hsi.N = prod(hsi.sizeof);