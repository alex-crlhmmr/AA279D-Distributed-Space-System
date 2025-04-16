function [r, v] = kepler_to_ijk(orb, mu)
a = orb(1);
e = orb(2);
i = deg2rad(orb(3));
O = deg2rad(orb(4));
w = deg2rad(orb(5));
nu = deg2rad(orb(6));

% Constants
E = 2*atan2(sqrt((1-e)/(1+e))*tan(nu/2), 1);
n = sqrt(mu/a^3);

% Rotation matrices
Rz_Omega = [cos(-O) sin(-O) 0; -sin(-O) cos(-O) 0; 0 0 1];
Rx_i = [1 0 0; 0 cos(-i) sin(-i); 0 -sin(-i) cos(-i)];
Rz_w = [cos(-w) sin(-w) 0; -sin(-w) cos(-w) 0; 0 0 1];
R = Rz_Omega * Rx_i * Rz_w;

% Calculate r and v
r_pqw = [a*(cos(E)-e); a*sqrt(1-e^2)*sin(E); 0];
r = R * r_pqw;
v_pqw = ((a*n)/(1-e*cos(E))) * [-sin(E); sqrt(1-e^2)*cos(E); 0];
v = R * v_pqw;

end