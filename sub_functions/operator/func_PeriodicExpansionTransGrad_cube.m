function[u] = func_PeriodicExpansionTransGrad_cube(PDu, isTF)

[v, h, c, d, s1, s2, s3] = size(PDu);
u = zeros(v,h,c,d);
for i = 1:s1
    for j = 1:s2
        for k = 1:s3
            u = u + circshift(PDu(:,:,:,:,i,j,k), [-i+1, -j+1, -k+1]);
        end
    end
end

if isTF == 1 % to make it parseval tight frame
    u = u/sqrt(s1*s2);
end