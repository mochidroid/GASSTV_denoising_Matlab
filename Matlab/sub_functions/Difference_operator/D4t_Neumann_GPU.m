function D4tz = D4t_Neumann_GPU(z)
    D4tz = cat(1, -z(1, :, :, 1), -z(2:end-1, :, :, 1) + z(1:end-2, :, :, 1), z(end-1, :, :, 1)) ...
        + cat(2, -z(:, 1, :, 2), -z(:, 2:end-1, :, 2) + z(:, 1:end-2, :, 2), z(:, end-1, :, 2));

    z_lt = z(:, :, :, 3); % direction of down right
    z_lt(end, :, :) = 0;
    z_lt(:, end, :) = 0;
    D4tz = D4tz + circshift(z_lt, [1 1 0]) - z_lt;

    z_lb = z(:, :, :, 4); % direction of up right
    z_lb(1, :, :) = 0;
    z_lb(:, end, :) = 0;
    D4tz = D4tz + circshift(z_lb, [-1 1 0]) - z_lb;
end