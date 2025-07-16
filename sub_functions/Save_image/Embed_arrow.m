function [output_image] = Embed_arrow(u, head_pos, length, handle_width, head_width, dir_tblr)
head_length = (handle_width+1)/2 + head_width;
width = handle_width + head_width*2;
output_image = u;

imsize = size(u);
if numel(imsize) == 2
    n3 = 1; % gray image
else
    n3 = 3; % RGB image
end

%% Generating arrow
arrow_image = zeros([length, width, n3]); %initial

% Generating arrow head
for i = 1:head_length
    arrow_image(i, ...
        head_length - i + 1 : head_length + i - 1, :) = 1;
end

% Generating arrow handle
arrow_image(head_length+1:end, ...
    head_width+1:head_width+handle_width, :) = 1;

%% Embedding arrow
Embed_arrow_image = zeros(size(u));

switch dir_tblr
    case 't' % top
        Embed_arrow_image(head_pos(1):head_pos(1)+length-1, ...
            head_pos(2)-head_width:head_pos(2)+handle_width+head_width-1, :) ...
            = arrow_image;

    case 'b' % bottom
        arrow_image = flip(arrow_image, 1);
        Embed_arrow_image(head_pos(1)-length+1:head_pos(1), ...
            head_pos(2)-head_width:head_pos(2)+handle_width+head_width-1, :) ...
            = arrow_image;

    case 'l' % left
        arrow_image = rot90(arrow_image, 1);
        Embed_arrow_image(head_pos(1)-head_width:head_pos(1)+handle_width+head_width-1, ...
            head_pos(2):head_pos(2)+length-1, :) ...
            = arrow_image;

    case 'r' % right
        arrow_image = rot90(arrow_image, -1);
        Embed_arrow_image(head_pos(1)-head_width:head_pos(1)+handle_width+head_width-1, ...
            head_pos(2)-length+1:head_pos(2), :) ...
            = arrow_image;

end

output_image(Embed_arrow_image == 1) = 1;