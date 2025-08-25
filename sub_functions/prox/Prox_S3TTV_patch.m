function[z] = Prox_S3TTV_patch(z, gamma, sel_mask)
% z: [v, h, c, d, s1, s2]
% sel_mask: [nb1, nb2] の logical。true の (i,j) ブロックのみ prox をかける
% 省略時は全ブロック更新

[v, h, c, d, b1, b2] = size(z);
nb1 = v / b1; nb2 = h / b2;

for k = 1:b1
    for l = 1:b2
        if ~sel_mask(k, l)
            z(:,:,:,:,k,l) = 0;
            continue; 
        end
        for i = 1:nb1
            for j = 1:nb2                
                block_tensor = z(1+b1*(i-1):b1*i, 1+b2*(j-1):b2*j, :, :, k, l);
                M = reshape(block_tensor, [b1*b2, c*d]);
                % SVD soft-threshold
                [U_, S_, V_] = svd(M, 'econ'); 
                Sthre = diag(max(0, diag(S_) - gamma));
                z(1+b1*(i-1):b1*i, 1+b2*(j-1):b2*j, :, :, k, l) ...
                    = reshape(U_ * Sthre * V_', [b1, b2, c, d]);
            end
        end
    end
end