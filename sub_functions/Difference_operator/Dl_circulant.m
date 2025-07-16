function result = Db_circulant(z)
    result = z(:, :, [2:end, 1], :) - z;
end