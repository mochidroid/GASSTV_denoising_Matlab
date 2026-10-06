function result = Dvt_Neumann(z)
    result = cat(1, -z(1, :, :), -z(2:end-1, :, :) + z(1:end-2, :, :), z(end-1, :, :));
end