% Author: Shunsuke Ono (ono@c.titech.ac.jp)

function[u] = ProxDCweightedL1(u, gamma, w)

u(:,:,1,:) = sign(u(:,:,1,:)).*max(abs(u(:,:,1,:))-w*gamma, 0);
u(:,:,2:end,:) = sign(u(:,:,2:end,:)).*max(abs(u(:,:,2:end,:))-gamma, 0);