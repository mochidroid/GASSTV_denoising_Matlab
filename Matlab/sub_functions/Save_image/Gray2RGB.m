function [u_rgb] = Gray2RGB(u, cmap)
u_ind = uint8(255*u);
u_rgb = ind2rgb(u_ind, cmap);