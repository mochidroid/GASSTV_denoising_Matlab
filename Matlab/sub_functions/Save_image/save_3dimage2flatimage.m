function [rpt_gray_frames] = save_3dimage2flatimage(gray_frames)
% set gray_frames to 3D signals
n1 = size(gray_frames, 1);
n2 = size(gray_frames, 2);
n3 = size(gray_frames, 3);
% projection transforming -> pt
pt_n1 = n1 + n3 - 1;
pt_n2 = n2 + n3 - 1;
pt_gray_frames = ones(pt_n1, pt_n2);
for i = n3:-1:1
    pt_gray_frames((pt_n1 - n1 - i + 2):(pt_n1 - i + 1), i:(i + n2 - 1)) = gray_frames(:, :, i);
end
imshow(pt_gray_frames);
% ER_3rd = 1/2;
ER_3rd = 1/3;
num_frame = floor(n1*ER_3rd);
step = floor(n3/n1/ER_3rd);
resize_n1 = n1 + num_frame;
resize_n2 = n2 + num_frame;
rpt_gray_frames = ones(resize_n1, resize_n2);
for i = num_frame:-1:1
%     if mod(i, 15) == 1
    if mod(i, 1) == 0
        rpt_gray_frames((resize_n1 - n1 - i + 2):(resize_n1 - i + 1), i:(i + n2 - 1)) = gray_frames(:, :, step*(i - 1) + 1);
    end
end
% imshow(rpt_gray_frames);