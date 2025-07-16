function result = Dlt_Neumann(z)
    result = cat(3, -z(:, :, 1), -z(:, :, 2:end-1) + z(:, :, 1:end-2), z(:, :, end-1));
end