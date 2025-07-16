function result = Prox_l12norm_d4(X, gamma)
    T = max(1 - gamma./sqrt(sum(X.*X, 4)), 0);
    result = T.*X;

end