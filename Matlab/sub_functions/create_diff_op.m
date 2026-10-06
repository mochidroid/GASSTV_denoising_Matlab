%% create difference operator
Dv_op = spdiags([-ones(hsi.n1,1) ones(hsi.n1,1)], [0,1], hsi.n1, hsi.n1);
Dv_op(hsi.n1, :) = 0;
Dv_op = kron(speye(hsi.n2 * hsi.n3), Dv_op);

Dh_op = spdiags([-ones(hsi.n1*hsi.n2,1) ones(hsi.n1*hsi.n2,1)], [0, hsi.n1], hsi.n1*hsi.n2, hsi.n1*hsi.n2);
Dh_op(hsi.n1*(hsi.n2-1) + 1 : hsi.n1*hsi.n2, :) = 0;
Dh_op = kron(speye(hsi.n3), Dh_op);

Db_op = spdiags([-ones(hsi.N,1) ones(hsi.N,1)], [0, hsi.n1*hsi.n2], hsi.N, hsi.N);
Db_op(hsi.n1*hsi.n2*(hsi.n3-1)+1 : hsi.N, :) = 0;