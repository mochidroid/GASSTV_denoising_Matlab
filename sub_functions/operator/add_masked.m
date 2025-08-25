function Y = add_masked(Y, X, sel_mask, blocksize)
% Y <- Y + X を、sel_mask のブロック領域だけで行う
% Y, X は [n1, n2, n3, 2, s1, s2]、sel_mask は [nb1, nb2]
[n1,n2,~,~,s1,s2] = size(Y);
b1 = blocksize(1); b2 = blocksize(2);
nb1 = n1/b1; nb2 = n2/b2;

% マスクを [n1,n2] に拡張（各ブロック内は 1）
tile_mask = kron(sel_mask, ones(b1,b2,'logical'));
tile_mask = reshape(tile_mask, [n1, n2]);

% 6次元にブロードキャスト
tile_mask_6d = reshape(tile_mask, [n1, n2, 1, 1, 1, 1]);
tile_mask_6d = repmat(tile_mask_6d, [1,1, size(Y,3), size(Y,4), s1, s2]);

Y = Y + X .* tile_mask_6d;


z(1+b1*(i-1):b1*i, 1+b2*(j-1):b2*j, :, :, k, l);