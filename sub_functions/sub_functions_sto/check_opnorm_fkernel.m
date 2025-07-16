
%% filterの準備

% blocksize = [3,3]; % block size of ASTV
% fkernel = sqrt(fspecial('gaussian', blocksize, 0.5));

% size_image = [11 11];

% clear
% close all;

addpath(genpath('sub_functions'))
% addpath('function_STSSTV')
rng('default')

fprintf('******* initium *******\n');

%% Generating observation
%%%%%%%%%%%%%%%%%%%%% User settings of experiment %%%%%%%%%%%%%%%%%%%%%%%%%%%%
image = 'Beltsville'; % image set = {'Beltsville','Pavia', 'MoffettField'};
deg.sigma       = 0.05; % standard derivation of Gaussian noise
deg.sp_rate     = 0.05; % ratio of salt-and-pepper noise
deg.dl_rate     = 0; % ratio of dead line noise

% Loading HSI
image_file_name_tmp = image + "_MixedNoise_" + ...
    "g" + num2str(deg.sigma) + "_sp" + num2str(deg.sp_rate) + "_dl" + num2str(deg.dl_rate) + ...
    "_size_64_64";
image_file_name = "./images/" + image_file_name_tmp + ".mat";
load(image_file_name)

% deg, gaussian_noise, hsi, noisy_HSI, org_HSI, sparse_noise, start_pos, u_org

% [hsi.n1, hsi.n2, hsi.n3] = size(noised_HSI);
% hsi.siseof = size(noised_HSI);
% hsi.N = hsi.n1 * hsi.n2 * hsi.n3;

%% Setting parameters
%%%%%%%%%%%%%%%%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%
blocksize = [8,8]; % block size of STSSTV
shiftstep = [1,1]; % [1,1] means full overlap
kernel = 'Uniform'; % weighting for each block
% params.kernel = 'Gaussian';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n1,n2,n3] = size(noisy_HSI);

%% Setting operator
D = @(z) cat(4, z([2:end, 1],:,:) - z, z(:,[2:end, 1],:) - z); % difference operator with circular boundary
Dt = @(z) z([end,1:end-1],:,:,1) - z(:,:,:,1) + z(:,[end,1:end-1],:,2) - z(:,:,:,2);
Db = @(z) z(:,:,[2:end, 1],:) - z; % difference operator with circular boundary
Dbt = @(z) z(:,:,[end,1:end-1],:) - z(:,:,:,:);
isTF = 0; % isTF = 1 makes P a tight frame
P = @(z) func_PeriodicExpansionGrad(z, blocksize, shiftstep, isTF);
Pt = @(z) func_PeriodicExpansionTransGrad(z, isTF);
switch kernel
    case 'Gaussian'
        fkernel = sqrt(fspecial('gaussian', blocksize, 0.5));
        fkernel = repmat(fkernel,[fix(n1/blocksize(1))+1, fix(n2/blocksize(2))+1, n3, 2, blocksize(1)/shiftstep(1), blocksize(2)/shiftstep(2)]);
        fkernel = fkernel(1:n1,1:n2,:,:,:,:);
        W = @(z) z.*fkernel;
        Wt = @(z) z.*fkernel;
    otherwise
        fkernel = ones(n1,n2,n3,2,blocksize(1)/shiftstep(1),blocksize(2)/shiftstep(2)) / prod(blocksize);
        W = @(z) z/prod(blocksize);
        Wt = @(z) z/prod(blocksize);
end


%% Setting step size for DP-PDS
W_opnorm = sqrt(max(fft2(circshift(fkernel, [2, 2])), [], 'all'));
P_opnorm = sqrt(prod(blocksize./shiftstep));
D_opnorm = 2*sqrt(2);

%% fftで求めた作用素ノルム

% fkernel_tmp = zeros(size_image);
% fkernel_tmp(1:(floor(blocksize(1)/2) + 1), 1:(floor(blocksize(2)/2) + 1)) = ...
%     fkernel((floor(blocksize(1)/2) + 1):end, (floor(blocksize(2)/2) + 1):end);
% 
% fkernel_tmp(1:(floor(blocksize(1)/2) + 1), (size_image(2) - floor(blocksize(2)/2) + 1):size_image(2)) = ...
%     fkernel((floor(blocksize(1)/2) + 1):end, 1:floor(blocksize(2)/2));
% 
% fkernel_tmp((size_image(1) - floor(blocksize(1)/2) + 1):size_image(1), 1:(floor(blocksize(2)/2) + 1)) = ...
%     fkernel(1:floor(blocksize(1)/2), (floor(blocksize(2)/2) + 1):end);
% 
% fkernel_tmp((size_image(1) - floor(blocksize(1)/2) + 1):size_image(1), (size_image(2) - floor(blocksize(2)/2) + 1):size_image(2)) = ...
%     fkernel(1:floor(blocksize(1)/2), 1:floor(blocksize(2)/2));
% % fkernel_tmp(1:blocksize(1), 1:blocksize(2)) = fkernel;
% 
% % ones_fft = fft(fft(ones(size_image), [], 1)'/sqrt(size_image(1)), [], 2)'/sqrt(size_image(2));
% % fkernel_fft = fft(fft(fkernel_tmp.*ones_fft, [], 1)/sqrt(size_image(1)), [], 2)/sqrt(size_image(2));
% % norm1 = sqrt(max(fkernel_fft.*conj(fkernel_fft), [], "all"));
% 
% % ones_fft = fft2(ones(size_image));
% % fkernel_fft = ifft2(fkernel_tmp.*ones_fft);
% % norm1 = sqrt(max(fkernel_fft.*conj(fkernel_fft), [], "all"));
% 
% norm1 = sqrt(max(fft2(circshift(fkernel, [2, 2])), [], "all"));

%% svdで求めた作用素ノルム

size_image = size(noisy_HSI);
func_fkernel = @(x, tflag) reshape(imfilter(reshape(x, size_image), gather(fkernel), 'circular'), [prod(size_image) 1]);
s = svds(func_fkernel, [prod(size_image) prod(size_image)], 1, 'largest', 'Tolerance', 1e-3, 'MaxIterations', 100);
norm2 = max(s);


%%
% X = rand(size_image);
% X_imfilt = imfilter(X, fkernel, 'circular');
% % fkernel_fft = fft(fft(fkernel_tmp, [], 1)/sqrt(size_image(1)), [], 2)/sqrt(size_image(2));
% % X_fftfilt = fkernel_fft.*fft(fft(X, [], 1)/sqrt(size_image(1)), [], 2)/sqrt(size_image(2));
% % X_fftfilt = fft(fft(X_fftfilt, [], 1)'/sqrt(size_image(1)), [], 2)'/sqrt(size_image(2));
% fkernel_fft = fft2(fkernel_tmp);
% X_fftfilt = fkernel_fft.*fft2(X);
% X_fftfilt = ifft2(X_fftfilt);
% diff_X = sum(abs(X_imfilt - X_fftfilt), "all");

%%
% disp('*******************************************')
% disp(append('norm1 : ', num2str(norm1)));
% disp(append('norm2 : ', num2str(norm2)));
% disp(append('diff_X : ', num2str(diff_X)));
% disp('*******************************************')


