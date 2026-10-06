% Generating an image for plot, where the region of interest is enlarged and
% embedded in the output image 
function[output_image] = Crop_Embed_image(...
    u, start_pos, crop_size, expansion_rate, embed_tblr)
end_pos = start_pos + crop_size - 1;

imsize = size(u);
if numel(imsize) == 2
    n3 = 1; % gray image
else
    n3 = 3; % RGB image
end

%% Cropping image
output_image = u;
%for i=1:numel(croptblr)
u_crop = u(start_pos(1):end_pos(1), start_pos(2):end_pos(2),:);
u_crop = imresize(u_crop, expansion_rate, 'nearest');

%% Generating embedded image
u_embed = ones(crop_size(1)*expansion_rate+4, crop_size(2)*expansion_rate+4, n3);
u_embed(3:end-2, 3:end-2, :) = u_crop;
embed_size = size(u_embed);

%% Enclosing cropped area
output_image(start_pos(1)-1:end_pos(1)+1, start_pos(2)-1, :) = 1; % left
output_image(start_pos(1)-1:end_pos(1)+1, end_pos(2)+1, :) = 1; % right
output_image(start_pos(1)-1, start_pos(2)-1:end_pos(2)+1, :) = 1; % top
output_image(end_pos(1)+1, start_pos(2)-1:end_pos(2)+1, :) = 1; % bottom

%% Enbedding cropped image
switch embed_tblr % selecting position: t=top, b=bottom, l=left, r=right
    case 'tl'
        output_image(1:embed_size(1),1:embed_size(2),:) = u_embed;
    case 'tr'
        output_image(1:embed_size(1),end-embed_size(2)+1:end,:) = u_embed;
    case 'bl'
        output_image(end-embed_size(1)+1:end,1:embed_size(2),:) = u_embed;
    case 'br'
        output_image(end-embed_size(1)+1:end,end-embed_size(2)+1:end,:) = u_embed;
end
%end