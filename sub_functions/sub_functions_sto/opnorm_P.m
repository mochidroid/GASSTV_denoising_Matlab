addpath(genpath('sub_functions'))

blocksize = [3, 3];
shiftstep = [1, 1];
% shiftstep = [5, 5];


n1 = 50;
n2 = 50;
n3 = 10;
n4 = 2;

size.n1 = n1;
size.n2 = n2;
size.n3 = n3;
size.n4 = n4;
size.blocksize = blocksize;
size.shiftstep = shiftstep;

P = @(z) func_PeriodicExpansion(z, blocksize, shiftstep);
Pt = @(z) func_PeriodicExpansionTrans(z);

% [n1, n2, n3, n4, blocksize(1)/shiftstep(1), blocksize(2)/shiftstep(2)]
s = svds(@(x,tflag) Afun(x, tflag, P, Pt, size), [n1*n2*n3*n4*blocksize(1)/shiftstep(1)*blocksize(2)/shiftstep(2), n1*n2*n3*n4], 1);
% s = svds(@(x,tflag) Afun(x, tflag, P, Pt, size), [n1*n2*n3*n4, n1*n2*n3*n4*blocksize(1)/shiftstep(1)*blocksize(2)/shiftstep(2)], 1);

function y = Afun(x,tflag, P, Pt, size)
n1 = size.n1;
n2 = size.n2;
n3 = size.n3;
n4 = size.n4;
blocksize = size.blocksize;
shiftstep = size.shiftstep;
if strcmp(tflag,'notransp')
    z = reshape(x, n1, n2, n3, n4);
    Y = P(z);
    y = reshape(Y, [n1*n2*n3*n4*blocksize(1)/shiftstep(1)*blocksize(2)/shiftstep(2), 1]);
else
    z = reshape(x, n1, n2, n3, n4, blocksize(1)/shiftstep(1), blocksize(2)/shiftstep(2));
    Y = Pt(z);
    y = reshape(Y, [n1*n2*n3*n4, 1]);
end
end
