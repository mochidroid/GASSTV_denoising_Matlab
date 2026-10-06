function result = Dvh_circulant(z)
    result = cat(4, z([2:end, 1],:,:) - z, z(:,[2:end, 1],:) - z);
end