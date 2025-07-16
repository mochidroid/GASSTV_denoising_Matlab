function result = D4_Neumann(z)
    result = cat(4, z([2:end, end],:,:) - z, z(:,[2:end, end],:) - z, ...
                    z([1, 1:end-1], :, :) - z, z(:, [1, 1:end-1], :) - z);
end