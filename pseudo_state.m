function b = pseudo_state(chief_OE, deputy_rSOE0, deputy_rSOEd, tau, mu, J2, Re)
% Pull inputs
a_c        = chief_OE(1);
e_c        = chief_OE(2);
i_c        = chief_OE(3);
Omega_c    = chief_OE(4);
omega_c    = chief_OE(5);
nu_c       = chief_OE(6);
E_c = 2*atan2(tan(nu_c/2), sqrt((1+e_c)/(1-e_c)));
M_c = E_c - e_c*sin(E_c);

da         = deputy_rSOE0(1);
dM         = deputy_rSOE0(2);
de         = deputy_rSOE0(3);
domega     = deputy_rSOE0(4);
di         = deputy_rSOE0(5);
dOmega     = deputy_rSOE0(6);

da_des     = deputy_rSOEd(1);
dM_des     = deputy_rSOEd(2);
de_des     = deputy_rSOEd(3);
domega_des = deputy_rSOEd(4);
di_des     = deputy_rSOEd(5);
dOmega_des = deputy_rSOEd(6);

% Compute initial deputy singular absolute OEs
a_d = a_c+da;
e_d = e_c+de;
i_d = i_c+di;
Omega_d = Omega_c+dOmega;
omega_d = omega_c+domega;
M_d = M_c+dM;
E_d = newton_raphson(M_d, e_d);
nu_d = 2 * atan2(sqrt(1+e_d)*sin(E_d/2), sqrt(1-e_d)*cos(E_d/2));
deputy_OE_0 = [a_d, e_d, i_d, Omega_d, omega_d, nu_d];

% Compute desired deputy orbit singular absolute OEs
a_d_des = a_c+da_des;
e_d_des = e_c+de_des;
i_d_des = i_c+di_des;
Omega_d_des = Omega_c+dOmega_des;
omega_d_des = omega_c+domega_des;
M_d_des = M_c+dM_des;
E_d_des = newton_raphson(M_d_des, e_d_des);
nu_d_des = 2 * atan2(sqrt(1+e_d_des)*sin(E_d_des/2), sqrt(1-e_d_des)*cos(E_d_des/2));
deputy_OE_des = [a_d_des, e_d_des, i_d_des, Omega_d_des, omega_d_des, nu_d_des];

% Compute initial and desired Deputy QNS ROEs
deputy_rQNSOE_0 = OE2rQNSOE(chief_OE, deputy_OE_0);
deputy_rQNSOE_des = OE2rQNSOE(chief_OE, deputy_OE_des);

% Compute Phi and pseudo state
Phi = stm_qns_J2(a_c, e_c, i_c, omega_c, mu, J2, Re, tau);
b = deputy_rQNSOE_des - Phi*deputy_rQNSOE_0;
end


function E = newton_raphson(M, e, epsilon)
    if nargin < 3
        epsilon = 1e-10;
    end
    E = M;
    max_iter = 1e5;
    for i = 1:max_iter
        f_E = E - e * sin(E) - M;
        f_prime_E = 1 - e * cos(E);
        increment = f_E / f_prime_E;
        E = E - increment;
        if abs(increment) <= epsilon
            break;
        end
    end
end