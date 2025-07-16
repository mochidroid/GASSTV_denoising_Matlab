function z = Dg_operator(z, sigma_x)
[v, h, c, d, s1, s2] = size(z);

num_block_v = fix(v/blocksize(1));
num_block_h = fix(h/blocksize(2));

for i = 1:num_block_v
    for j = 1:num_block_h
        for k = 1:s1
            for l = 1:s2
                if (i == num_block_v) && (k-1-v+(num_block_v*blocksize(1)) > 0) || ...
                        (j == num_block_h) && (l-1-h+(num_block_v*blocksize(2)) > 0)
                    z(1+blocksize(1)*(i-1):blocksize(1)*i,1+blocksize(2)*(j-1):blocksize(2)*j,:,:,k,l)...
                        = zeros([blocksize, c, d]);
                end
            end
        end
    end
end
end