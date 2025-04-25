function [f_vals, states] = linear_ecc_mapping(quasi_elements, alpha0, f_change, step_size)
% Solve for constants
a = alpha0(1);         
e = alpha0(2);
w = alpha0(5);
i = alpha0(3);
mu = 398600.4418;    
n = sqrt(mu/a^3);

% Calculate values needed to get constants
eta = sqrt(1-e^2);
ex = e*cos(w);
ey = e*sin(w);

% F vals
f_initial = alpha0(6);
f_final = f_initial + f_change;
f_vals = f_initial:step_size:f_final;

% Calculate t matrix
E = 2*atan2(sqrt((1 - e)./(1 + e)).*tan(f_vals./2), 1);
E = unwrap(E);
M = E - e.*sin(E);               
M0 = M(1);                       
ts = (M - M0) ./ n;    

% Initialize states
states = zeros(length(f_vals), 6);

for j = 1:length(f_vals)
% assume scalars k, eta, n, t, ex, ey, u, i are already defined
f = f_vals(j);
k = 1 + e*cos(f);
kp = -e*sin(f);
u = f + w;
t = ts(j);

% Compute B matrix entries (x component shown as example)
bx1 = 1 / k + (3 / 2) * kp * n / eta^3 * t;
bx2 = -kp / eta^3;
bx3 = (1 / eta^3) * (ex * (k - 1) / (1 + eta) - cos(u));
bx4 = (1 / eta^3) * (ey * (k - 1) / (1 + eta) - sin(u));
bx6 = kp / eta^3 * cot(i);
by1 = -(3 / 2) * k * n / eta^3 * t;
by2 = k / eta^3;
by3 = (1 / eta^2) * ((1 + 1 / k) * sin(u) + ey / k + k / eta * (ey / (1 + eta)));
by4 = -(1 / eta^2) * ((1 + 1 / k) * cos(u) + ex / k + k / eta * (ex / (1 + eta)));
by6 = (1 / k - k / eta^3) * cot(i);
bz5 = (1 / k) * sin(u);
bz6 = -(1 / k) * cos(u);
% Velocities
bxd1 = kp / 2 + (3 / 2) * k^2 * (1 - k) * n / eta^3 * t;
bxd2 = k^2 / eta^3 * (k - 1);
bxd3 = k^2 / eta^3 * (eta * sin(u) + ey * (k - 1) / (1 + eta));
bxd4 = -k^2 / eta^3 * (eta * cos(u) + ex * (k - 1) / (1 + eta));
bxd6 = -k^2 / eta^3 * (k - 1) * cot(i);
byd1 = -(3 / 2) * k * (1 + k * kp * n / eta^3 * t);
byd2 = k^2 / eta^3 * kp;
byd3 = (1 + k^2 / eta^3) * cos(u) + ex / eta^2 * k  * (1 + k / eta * (1 - k) / (1 + eta));
byd4 = (1 + k^2 / eta^3) * sin(u) + ey / eta^2 * k  * (1 + k / eta * (1 - k) / (1 + eta));
byd6 = -(1 + k^2 / eta^3) * kp * cot(i);
bzd5 = cos(u)+ex;
bzd6 = sin(u)+ey;
% Assemble B matrix
B = [bx1, bx2, bx3, bx4, 0, bx6;
     by1, by2, by3, by4, 0, by6;
     0,    0,    0,    0, bz5, bz6;
     bxd1, bxd2, bxd3, bxd4, 0, bxd6;
     byd1, byd2, byd3, byd4, 0, byd6;
     0,    0,    0,    0, bzd5, bzd6];
M = [ ...
    a*eta^2*eye(3),    zeros(3); ...
    zeros(3),        (a*n)/eta*eye(3) ...
];

Phi = M * B;

state = Phi*quasi_elements;
states(j,:) = state;


end


end
