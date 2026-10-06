% Author: Shunsuke Ono (ono@isl.titech.ac.jp)
% Last version: May. 1, 2017

function[out] = ProxL1norm(x, gamma)

x(abs(x) <= gamma) = 0;
x(abs(x) > gamma) = x(abs(x) > gamma) - sign(x(abs(x) > gamma)).*gamma;
out = x;
