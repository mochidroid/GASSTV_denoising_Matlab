function  [HSI_clean, hsi] = Load_HSI(set_image)
rng('default')
%% Loading HSI
switch set_image
    case 'JasperRidge'
        load('H:/マイドライブ/MATLAB_Share/HSIData/JasperRidge/jasperRidge2_R198.mat');
        U_tmp = reshape(Y, [198, 100, 100]);
        HSI_clean = normalize01(permute(U_tmp, [2,3,1]));
        % U_tmp = permute(U_tmp, [2,3,1]);
        % HSI_clean = normalize01(U_tmp(:,:,[1:102, 110:143, 147:end]));

    case 'JasperRidge64'
        load('H:/マイドライブ/MATLAB_Share/HSIData/JasperRidge/jasperRidge2_R198.mat');
        U_tmp = reshape(Y, [198, 100, 100]);
        U_tmp = permute(U_tmp, [2,3,1]);
        hsi.start_pos = [1, 37, 1];
        hsi_size = [64, 64, size(U_tmp, 3)-hsi.start_pos(3)+1];
        hsi.end_pos = hsi.start_pos + hsi_size - 1;
        HSI_clean = normalize01( ...
            U_tmp(hsi.start_pos(1):hsi.end_pos(1), ...
            hsi.start_pos(2):hsi.end_pos(2), ...
            [1:102, 110:143, 147:end]) ...
            );


    case 'PaviaU'
        load('H:/マイドライブ/MATLAB_Share/HSIData/PaviaU/PaviaU.mat')
        hsi.start_pos = [170, 200, 5];
        hsi_size = [140, 140, size(paviaU, 3)-hsi.start_pos(3)+1];
        hsi.end_pos = hsi.start_pos + hsi_size - 1;
        HSI_clean = normalize01(paviaU(hsi.start_pos(1):hsi.end_pos(1), hsi.start_pos(2):hsi.end_pos(2), hsi.start_pos(3):hsi.end_pos(3)));

    case 'PaviaU120'
        load('H:/マイドライブ/MATLAB_Share/HSIData/PaviaU/PaviaU.mat')
        hsi.start_pos = [170, 210, 5];
        hsi_size = [120, 120, size(paviaU, 3)-hsi.start_pos(3)+1];
        hsi.end_pos = hsi.start_pos + hsi_size - 1;
        HSI_clean = normalize01(paviaU(hsi.start_pos(1):hsi.end_pos(1), hsi.start_pos(2):hsi.end_pos(2), hsi.start_pos(3):hsi.end_pos(3)));
        
    case 'PaviaU64'
        load('H:/マイドライブ/MATLAB_Share/HSIData/PaviaU/PaviaU.mat')
        hsi.start_pos = [211, 211, 5];
        hsi_size = [64, 64, size(paviaU, 3)-hsi.start_pos(3)+1];
        hsi.end_pos = hsi.start_pos + hsi_size - 1;
        HSI_clean = normalize01(paviaU(hsi.start_pos(1):hsi.end_pos(1), hsi.start_pos(2):hsi.end_pos(2), hsi.start_pos(3):hsi.end_pos(3)));


    case 'WashingtonDC'
        load('H:/マイドライブ/MATLAB_Share/HSIData/WashingtonDC/WashingtonDC_image.mat')
        hsi.start_pos = [657, 140, 1];
        hsi_size = [100, 100, size(WashingtonDC, 3)];
        hsi.end_pos = hsi.start_pos + hsi_size - 1;
        HSI_clean = normalize01(WashingtonDC(hsi.start_pos(1):hsi.end_pos(1), ...
            hsi.start_pos(2):hsi.end_pos(2), hsi.start_pos(3):hsi.end_pos(3)));

    case 'Beltsville'
        load('H:/マイドライブ/MATLAB_Share/HSIData/Beltsville.mat');
        image = "Beltsville";
        org_image = u_org;
        start_pos = [145, 157, 1];
        HSI_size = [100, 100, size(org_image, 3)];
        end_pos = start_pos + HSI_size - 1;
        org_HSI_tmp = org_image(start_pos(1):end_pos(1), start_pos(2):end_pos(2), start_pos(3):end_pos(3));
        HSI_clean = normalize01(org_HSI_tmp);

    case 'MoffettField128'
        load('H:/マイドライブ/MATLAB_Share/HSIData/MoffettField.mat');
        image = "MoffettField";
        org_image = I_REF;
        start_pos = [21, 11, 1];
        HSI_size = [128, 128, size(org_image, 3)];
        end_pos = start_pos + HSI_size - 1;
        org_HSI_tmp = org_image(start_pos(1):end_pos(1), start_pos(2):end_pos(2), start_pos(3):end_pos(3));
        HSI_clean = normalize01(org_HSI_tmp);

    case 'MoffettField64'
        load('H:/マイドライブ/MATLAB_Share/HSIData/MoffettField.mat');
        image = "MoffettField";
        org_image = I_REF;
        % start_pos = [60, 70, 1];
        % start_pos = [43, 10, 1];
        start_pos = [11, 120, 1];
        HSI_size = [64, 64, size(org_image, 3)];
        end_pos = start_pos + HSI_size - 1;

        org_HSI_tmp = org_image(start_pos(1):end_pos(1), start_pos(2):end_pos(2), start_pos(3):end_pos(3));
        HSI_clean = normalize01(org_HSI_tmp);

    case 'Salinas'
        load('H:/マイドライブ/MATLAB_Share/HSIData/Salinas_corrected.mat');
        org_image = salinas_corrected;
        start_pos = [221, 71, 7];
        HSI_size = [100, 100, size(org_image, 3)-6];
        end_pos = start_pos + HSI_size - 1;
        org_HSI_tmp = org_image(start_pos(1):end_pos(1), start_pos(2):end_pos(2), start_pos(3):end_pos(3));
        HSI_clean = normalize01(org_HSI_tmp);
   
end


hsi.sizeof = size(HSI_clean);
[hsi.n1,hsi.n2,hsi.n3] = size(HSI_clean);
hsi.N = prod(hsi.sizeof);
