% Generating an image for plot, where the region of interest is enlarged and
% embedded in the output image 
function[output_image, output_crop_image] = Crop_Sep_image(...
    u, start_pos, crop_size)
end_pos = start_pos + crop_size - 1;

imsize = size(u);

%% Cropping image
output_image = u;
%for i=1:numel(croptblr)
output_crop_image = u(start_pos(1):end_pos(1), start_pos(2):end_pos(2),:);

%% Enclosing cropped area
output_image(start_pos(1)-1:end_pos(1)+1, start_pos(2)-1, :) = 1; % left
output_image(start_pos(1)-1:end_pos(1)+1, end_pos(2)+1, :) = 1; % right
output_image(start_pos(1)-1, start_pos(2)-1:end_pos(2)+1, :) = 1; % top
output_image(end_pos(1)+1, start_pos(2)-1:end_pos(2)+1, :) = 1; % bottom

end