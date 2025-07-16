%% create Graph_Based spatial weighted difference operator
% guide_image = mean(Noi, 3);
guide_image = mean(HSI_noisy, 3);
[n1,n2,n3] = size(HSI_noisy);

diff_1_0 = cat(1, guide_image(2:n1, :), Inf(1, n2));
diff_0_1 = cat(2, guide_image(:, 2:n2), Inf(n1, 1));

diff_1_1 = cat(2, diff_1_0(:, 2:n2), Inf(n1, 1));
diff_1_2 = cat(1, Inf(1, n2), diff_0_1(1:n1-1, :));

diff_mat = [diff_1_0(:), diff_0_1(:), diff_1_1(:), diff_1_2(:)];
grad_mat = guide_image(:).*ones(1, 4) - diff_mat;

dist_mat = [1, 1, 2, 2];

W_mat_tmp = exp(-(grad_mat.^2)./(sigma_x.^2)/2).*exp(-dist_mat./(sigma_l.^2)/2);
W_mat = repmat(W_mat_tmp, n3, 1);

end_index = [1, n1, n1+1, n1-1];

Dg_1_0 = spdiags([-W_mat(:,1), W_mat(:,1)], [0 end_index(1)],  n1*n2*n3, n1*n2*n3+1);
Dg_0_1 = spdiags([-W_mat(:,2), W_mat(:,2)], [0 end_index(2)],  n1*n2*n3, n1*n2*n3+1);
Dg_1_1 = spdiags([-W_mat(:,3), W_mat(:,3)], [0 end_index(3)],  n1*n2*n3, n1*n2*n3+1);
Dg_1_2 = spdiags([-W_mat(:,4), W_mat(:,4)], [0 end_index(4)],  n1*n2*n3, n1*n2*n3+1);

Dg_1_0 = Dg_1_0(:,1:n1*n2*n3);
Dg_0_1 = Dg_0_1(:,1:n1*n2*n3);
Dg_1_1 = Dg_1_1(:,1:n1*n2*n3);
Dg_1_2 = Dg_1_2(:,1:n1*n2*n3);

Dg = cat(1, Dg_1_0, Dg_0_1, Dg_1_1, Dg_1_2);

%% create spectral difference operator
Db_op = spdiags([-ones(n1*n2*n3,1) ones(n1*n2*n3,1)], [0, n1*n2], n1*n2*n3, n1*n2*n3);
Db_op(n1*n2*(n3-1)+1 : n1*n2*n3, :) = 0;

DgDb = Dg * Db_op;
DgDbt = DgDb';

cube2vec = @(x) reshape(x, [n1*n2*n3,1]);
vec2cube = @(x) reshape(x, [n1,n2,n3]);