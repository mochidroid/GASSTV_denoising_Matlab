function[PDu] = func_PeriodicExpansionGrad_cube(Du, blocksize, shiftstep, isTF)

[v, h, c, d] = size(Du);
PDu = zeros(v, h, c, d, blocksize(1)/shiftstep(1), blocksize(2)/shiftstep(2), blocksize(3)/shiftstep(3));

for i = 1:blocksize(1)/shiftstep(1)
    for j = 1:blocksize(2)/shiftstep(2)
        for k = 1:blocksize(3)/shiftstep(3)
            PDu(:,:,:,:,i,j,k) = circshift(Du, [i-1, j-1, k-1]);
        end
    end
end

if isTF == 1 % to make it parseval tight frame
    PDu = PDu/sqrt(prod(blocksize)/prod(shiftstep));
end