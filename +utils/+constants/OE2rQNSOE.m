function rQNSOE = OE2rQNSOE(chief_OE, deputy_OE)
% Input: chief_OE = [a_c, e_c, i_c, Omega_c, omega_c, M_c]
%        deputy_OE = [a_d, e_d, i_d, Omega_d, omega_d, M_d]

a_c     = chief_OE(1);
e_c     = chief_OE(2);
i_c     = chief_OE(3);
Omega_c = chief_OE(4);
omega_c = chief_OE(5);
M_c     = chief_OE(6);

a_d     = deputy_OE(1);
e_d     = deputy_OE(2);
i_d     = deputy_OE(3);
Omega_d = deputy_OE(4);
omega_d = deputy_OE(5);
M_d     = deputy_OE(6);

rQNSOE = [
    (a_d - a_c)/a_c;
    (M_d + omega_d) - (M_c + omega_c) + (Omega_d - Omega_c)*cos(i_c);
    e_d*cos(omega_d) - e_c*cos(omega_c);
    e_d*sin(omega_d) - e_c*sin(omega_c);
    i_d - i_c;
    (Omega_d - Omega_c)*sin(i_c)
];
end
