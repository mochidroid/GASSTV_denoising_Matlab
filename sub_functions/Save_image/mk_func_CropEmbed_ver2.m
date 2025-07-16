% make image for plot, where the region of interest is enlarged and
% embedded in the output image 
function[out] = mk_func_CropEmbed_ver2(u, cropind, cropsize, croprate, croptblr, varargin)

imsize = size(u);
if numel(imsize) == 2
    c = 1;
else
    c = 3;
end



out = u;
%for i=1:numel(croptblr)
ucrop = u(cropind(1,1):cropind(1,1)+cropsize(1)-1, cropind(1,2):cropind(1,2)+cropsize(2)-1,:);
ucrop = imresize(ucrop, croprate, 'nearest');
uembed = ones(cropsize(1)*croprate+4, cropsize(2)*croprate+4, c);


if ~isempty(varargin)
    switch varargin{1}
        case 'r'
            uembed(:,:,2:3) = 0;
        case 'y'
            uembed(:,:,3) = 0;
            out(cropind(1,1)-1:cropind(1,1)+cropsize(1), cropind(1,2)-1:cropsize(2)+1:cropind(1,2)+cropsize(2),1:2) = 1;
            out(cropind(1,1)-1:cropind(1,1)+cropsize(1), cropind(1,2)-1:cropsize(2)+1:cropind(1,2)+cropsize(2),3) = 0;
            out(cropind(1,1)-1:cropsize(1)+1:cropind(1,1)+cropsize(1), cropind(1,2)-1:cropind(1,2)+cropsize(2),1:2) = 1;
            out(cropind(1,1)-1:cropsize(1)+1:cropind(1,1)+cropsize(1), cropind(1,2)-1:cropind(1,2)+cropsize(2),3) = 0;
        case 'k'
            uembed = uembed*0;
        otherwise
    end
end
uembed(3:end-2, 3:end-2, :) = ucrop;
embedsize = size(uembed);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% naganuma shuusei
% 周りの色を変える
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
srrounding_color = 1;
out(cropind(1):(cropind(1)+cropsize(1)), cropind(2), :) = srrounding_color;
out(cropind(1):(cropind(1)+cropsize(1)), cropind(2)+cropsize(2), :) = srrounding_color;
out(cropind(1), cropind(2):(cropind(2)+cropsize(2)), :) = srrounding_color;
out(cropind(1)+cropsize(1), cropind(2):(cropind(2)+cropsize(2)), :) = srrounding_color;

switch croptblr%{i}
    case 'tl'
        out(1:embedsize(1),1:embedsize(2),:) = uembed;
    case 'tr'
        out(1:embedsize(1),end-embedsize(2)+1:end,:) = uembed;
    case 'bl'
        out(end-embedsize(1)+1:end,1:embedsize(2),:) = uembed;
    case 'br'
        out(end-embedsize(1)+1:end,end-embedsize(2)+1:end,:) = uembed;
end
%end