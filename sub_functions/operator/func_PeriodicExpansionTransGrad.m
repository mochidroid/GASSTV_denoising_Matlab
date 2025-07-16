function[u] = func_PeriodicExpansionTransGrad(PDu, isTF)

[v, h, c, d, s1, s2] = size(PDu);
u = zeros(v,h,c,d);
for i = 1:s1
    for j = 1:s2
        u = u + circshift(PDu(:,:,:,:,i,j), [-i+1, -j+1, 0]);
    end
end

if isTF == 1 % to make it parseval tight frame
    u = u/sqrt(s1*s2);
end