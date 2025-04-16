function [a, e, i, O, w, nu] = ECItoKepler(r, v)
mu = 398600.4418;
h = cross(r, v);
i = atan2(sqrt(h(1)^2 + h(2)^2), h(3));
O = atan2(h(1), -h(2));
P = norm(h)^2 / mu;
a = (2/norm(r) - norm(v)^2 / mu)^(-1);
n = sqrt(mu/a^3);
e = sqrt(1 - P/a);
E = atan2(dot(r, v) / (a^(2)*n), 1-norm(r)/a);
nu = 2*atan2(sqrt((1+e)/(1-e))*tan(E/2), 1);
u = atan2(r(3)/sin(i), r(1)*cos(O) + r(2)*sin(O));
w = u-nu;

i  = mod(rad2deg(i), 360);
O  = mod(rad2deg(O), 360);
w  = mod(rad2deg(w), 360);
nu = mod(rad2deg(nu), 360);
end
