function D4z = D4_Neumann_GPU(z)
    [n1, n2, n3] = size(z);
    D4z = zeros([n1, n2, n3, 4], 'single', 'gpuArray');
    D4z(:,:,:,[1,2]) = cat(4, z([2:end, end],:,:) - z, z(:,[2:end, end],:) - z);
    D4z(1:n1-1, 1:n2-1, :, 3) = z(2:n1, 2:n2, :) - z(1:n1-1, 1:n2-1, :); % direction of down right
    D4z(2:n1, 1:n2-1, :, 4) = z(1:n1-1, 2:n2, :) - z(2:n1, 1:n2-1, :); % direction of up right
end