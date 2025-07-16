% inf_{Y} 1/2||Y-V||^{2} + 1/(2*gamma)||X-Y||^{2}

function result = prox_l2_norm(X, V, gamma)
    result = (gamma.*V + X)./(1 + gamma);
end
