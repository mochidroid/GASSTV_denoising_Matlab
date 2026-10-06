function result = Dh_Neumann(z)
    n2 = size(z, 2);
    result = z(:, [2:n2, n2], :) - z;
end