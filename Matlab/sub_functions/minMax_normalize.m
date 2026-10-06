function minMax_normalized_x = minMax_normalize(x, min_value, max_value)
    x_max = max(x, [], 'all');
    x_min = min(x, [], 'all');

    minMax_normalized_x = (x - x_min) / (x_max - x_min) * (max_value - min_value) + min_value;
end