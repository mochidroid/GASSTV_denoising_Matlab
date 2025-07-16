function [result_spectrum] = Extract_spectrum(HSI, x, y)

n3 = size(HSI, 3);

result_spectrum = reshape(HSI(x,y,:), [n3, 1]);