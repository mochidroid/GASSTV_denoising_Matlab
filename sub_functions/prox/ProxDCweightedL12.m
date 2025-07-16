% Author: Shunsuke Ono (ono@c.titech.ac.jp)

function[u] = ProxDCweightedL12(u, gamma, w)

u(:,:,1,:) = sign(u(:,:,1,:)).*max(abs(u(:,:,1,:))-w*gamma, 0);

u(:,:,2:end,:) = sign(u(:,:,2:end,:)).*max(abs(u(:,:,2:end,:))-gamma, 0);
thresh = ((sqrt(sum(u(:,:,2:end,:).^2, 4))).^(-1))*gamma;
thresh(thresh>1) = 1;
thresh = ones(size(thresh)) - thresh;
u(:,:,2:end,:) = repmat(thresh,1,1,size(u,3)-1,1).*u(:,:,2:end,:);