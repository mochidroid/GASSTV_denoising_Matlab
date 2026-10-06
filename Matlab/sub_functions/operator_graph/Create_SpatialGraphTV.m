function [Wsp, eta_sp] = Create_SpatialGraphTV(X, sigma_sp)
% Creating graph based weight matrix
guide_image = mean(X, 3);
[n1,n2,n3] = size(X);

diff_1_0 = cat(1, guide_image(2:n1, :), Inf(1, n2));
diff_0_1 = cat(2, guide_image(:, 2:n2), Inf(n1, 1));

diff_1_1 = cat(2, diff_1_0(:, 2:n2), Inf(n1, 1));
diff_1_2 = cat(1, Inf(1, n2), diff_0_1(1:n1-1, :));

diff_mat = [diff_1_0(:), diff_0_1(:), diff_1_1(:), diff_1_2(:)];
grad_mat = guide_image(:).*ones(1, 4) - diff_mat;

% dist_mat = [1, 1, 2, 2];

% W_mat_tmp = exp(-(grad_mat.^2)./(sigma_x.^2)/2).*exp(-dist_mat./(sigma_l.^2)/2);
W_mat_tmp = exp(-(grad_mat.^2)/(sigma_sp^2)/2);
% W_mat_tmp = exp(-abs(grad_mat)/(sigma_x));
W_mat = repmat(W_mat_tmp, n3, 1);

Wsp = reshape(W_mat, [n1, n2, n3, 4]);
Wsp = single(Wsp);

eta_sp = single(norm(reshape(W_mat_tmp.*grad_mat, [n1*n2*4, 1]), 1) * n3);

% implay(cat(2, W(:,:,:,1), W(:,:,:,2), W(:,:,:,3), W(:,:,:,4)))