function Phi = stm_qns_J2(a, e, inc, omega, mu, J2, Re, tau)

  % mean motion
  n = sqrt(mu/(a^3));

  % eccentricity-dependent terms
  eta   = sqrt(1 - e^2);
  kappa = (3*J2*Re^2*sqrt(mu))/(4*a^(7/2)*eta^4);
  E     = 1 + eta;
  F     = 4 + 3*eta;
  G     = 1/eta^2;

  % inclination-dependent terms
  P = 3*cos(inc)^2 - 1;
  Q = 5*cos(inc)^2 - 1;
  S = sin(2*inc);
  T = sin(inc)^2;

  % combined factors
  EP  = E*P;
  FGP = F*G*P;
  GQ = G*Q;

  % exi, exf, eyi, eyf
  omega_dot = kappa*Q;
  omega_f = omega + omega_dot*tau;
  exi = e*cos(omega);
  eyi = e*sin(omega);
  exf = e*cos(omega_f);
  eyf = e*sin(omega_f);


  % assemble the STM
  Phi = [
    1,                                  0,                  0,                    0,                   0,                0;
   -(1.5*n + 1.5*kappa*EP)*tau,        1,     kappa*exi*FGP*tau,     kappa*eyi*FGP*tau,    -kappa*F*S*tau,    0;
    3.5*kappa*eyf*Q*tau,                 0, cos(omega_dot*tau)-4*kappa*exi*eyf*GQ*tau, -sin(omega_dot*tau)-4*kappa*eyi*eyf*GQ*tau,  5*kappa*eyf*S*tau, 0;
   -3.5*kappa*exf*Q*tau,                 0, sin(omega_dot*tau)+4*kappa*exi*exf*GQ*tau,  cos(omega_dot*tau)+4*kappa*eyi*exf*GQ*tau, -5*kappa*exf*S*tau, 0;
    0,                                  0,                  0,                    0,                   1,                0;
    3.5*kappa*S*tau,                    0, -4*kappa*exi*G*S*tau,   -4*kappa*eyi*G*S*tau,  2*kappa*T*tau,      1
  ];
end
