addpath(genpath('sub_functions'))

blocksize = [10, 10];
% shiftstep = [1, 1];
shiftstep = [3, 3];


n1 = 50;
n2 = 50;
n3 = 10;
n4 = 2;
n5 = fix(blocksize(1) / shiftstep(1));
n6 = fix(blocksize(2) / shiftstep(2));

hsi.n1 = n1;
hsi.n2 = n2;
hsi.n3 = n3;
hsi.n4 = n4;
hsi.n5 = n5;
hsi.n6 = n6;
hsi.blocksize = blocksize;
hsi.shiftstep = shiftstep;

W = @(z) W_operator(z, blocksize);

% [n1, n2, n3, n4, blocksize(1)/shiftstep(1), blocksize(2)/shiftstep(2)]
s = svds(@(x,tflag) Afun(x, tflag, W, hsi), [n1*n2*n3*n4*n5*n6, n1*n2*n3*n4*n5*n6], 1);
% s = svds(@(x,tflag) Afun(x, tflag, P, Pt, size), [n1*n2*n3*n4, n1*n2*n3*n4*blocksize(1)/shiftstep(1)*blocksize(2)/shiftstep(2)], 1);

function y = Afun(x,tflag, W, hsi)
n1 = hsi.n1;
n2 = hsi.n2;
n3 = hsi.n3;
n4 = hsi.n4;
n5 = hsi.n5;
n6 = hsi.n6;
blocksize = hsi.blocksize;
if strcmp(tflag,'notransp')
    z = reshape(x, n1, n2, n3, n4, n5, n6);
    Y = W(z);
    y = reshape(Y, [n1*n2*n3*n4*n5*n6, 1]);
else
    z = reshape(x, n1, n2, n3, n4, n5, n6);
    Y = W(z);
    y = reshape(Y, [n1*n2*n3*n4*n5*n6, 1]);
end
end
