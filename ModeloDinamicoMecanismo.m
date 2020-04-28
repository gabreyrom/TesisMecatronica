syms q1 q2 My Mx m1 q1_p q2_p q1_bip q2_bip

jacobiano = jacobian([((Mx*cos(q2)/(sin(q1-q2))))*sin(q1),My - ((Mx*cos(q2)/(sin(q1-q2)))*cos(q1))],[q1,q2]);
jacobiano_tras = transpose(jacobiano);
M = simplify(m1 * jacobiano_tras * jacobiano);

q = [q1;q2];
q_p = [q1_p;q2_p];
q_bip = [q1_bip;q2_bip];

M_punto = simplify(diff(M,q1)*q_p(1,1)+diff(M,q2)*q_p(2,1));

%M_punto = simplify(diff(M,q1)+diff(M,q2));

aux = simplify(M*q_p);

C_aux = simplify(jacobian([aux(1),aux(2)],[q1,q2]));

C = simplify(C_aux - 0.5*transpose(C_aux));

%C = M_punto * q_p - diff(transpose(q_p) * M * q_p ,q1,q2);

inv_Jac = simplify(inv(jacobiano));
inv_Jac_punto = simplify(diff(inv_Jac,q1)*q_p(1,1)+diff(inv_Jac,q2)*q_p(2,1));
inv_Jac_tras = transpose(inv_Jac);

MatrizC_bar = inv_Jac_tras * inv_Jac_punto + inv_Jac_tras*C*inv_Jac;

M_bar = inv_Jac_tras * M * inv_Jac;

M_bar = simplify(M_bar);

MatrizC_bar = simplify(MatrizC_bar);