function result = Dvhb_Neumann(z)
    result = cat(4, z([2:end, end],:,:) - z, z(:,[2:end, end],:) - z, z(:, :, [2:end, end]) - z);
end