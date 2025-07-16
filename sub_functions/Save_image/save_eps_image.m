function save_eps_image(image, save_file_name)

fig = figure;
imshow(image)
axis off
set(gca, 'Position', [0 0 1 1]); % 余白をなくす
exportgraphics(gca, save_file_name, 'ContentType', 'vector');

close(fig)