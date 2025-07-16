addpath(genpath('sub_functions'))

n1 = 50;
n2 = 50;
n3 = 10;

hsi.n1 = n1;
hsi.n2 = n2;
hsi.n3 = n3;
hsi.sizeof = [n1, n2, n3];

Phi = @(z) dct_hsi_gpu(z, hsi.sizeof);
Phit = @(z) idct_hsi_gpu(z, hsi.sizeof);

% [n1, n2, n3, n4, blocksize(1)/shiftstep(1), blocksize(2)/shiftstep(2)]
s = svds(@(x,tflag) Afun(x, tflag, Phi, Phit, hsi), [n1*n2*n3, n1*n2*n3], 1);
% s = svds(@(x,tflag) ...
%       Afun(x, tflag, P, Pt, size), [n1*n2*n3*n4, n1*n2*n3*n4*blocksize(1)/shiftstep(1)*blocksize(2)/shiftstep(2)], 1);

fprintf("op_norm = %f\n", s);

function y = Afun(x,tflag, Phi, Phit, hsi)
n1 = hsi.n1;
n2 = hsi.n2;
n3 = hsi.n3;

if strcmp(tflag,'notransp')
    z = reshape(x, n1, n2, n3);
    Y = Phi(z);
    y = reshape(Y, [n1*n2*n3, 1]);
else
    z = reshape(x, n1, n2, n3);
    Y = Phit(z);
    y = reshape(Y, [n1*n2*n3, 1]);
end
end
