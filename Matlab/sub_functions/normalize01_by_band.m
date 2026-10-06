function normalized01_x = normalize01_by_band(x)
n3 = size(x, 3);
normalized01_x = zeros(size(x));


for i = 1:n3
    x_max = max(x(:,:,i), [], 'all');
    x_min = min(x(:,:,i), [], 'all');

    normalized01_x(:,:,i) = (x(:,:,i) - x_min) / (x_max - x_min);
end