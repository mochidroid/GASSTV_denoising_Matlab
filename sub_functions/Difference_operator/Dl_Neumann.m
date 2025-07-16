function result = Dl_Neumann(z)
    result = z(:, :, [2:end, end]) - z;
end