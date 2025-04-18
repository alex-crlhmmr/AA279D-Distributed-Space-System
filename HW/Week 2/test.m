%% Reset
clc; clear; close all;
%%
% Initial Keplerian orbital elements [a, e, i, RAAN, omega, nu]
alpha0 = [6771; 0.0005; deg2rad(51.64); deg2rad(257); 0; deg2rad(30)]; % Chief
alpha1 = [6771; 0.0006; deg2rad(51.69); deg2rad(257.5); deg2rad(0.5); deg2rad(25)]; % Deputy

% Calculate initial state [ECI position and velocity] of chief
mu = 398600.4418;
[r0, v0] = kepler_to_ijk(alpha0, mu);

% Calculate initial state [RTN position and velocity] of deputy
[r1, v1] = kepler_to_ijk(alpha1, mu);
[rRTN, vRTN] = ECI2RTN(r0, v0, r1, v1);

% Initialize state
state0 = zeros(12, 1);
state0(1:3) = r0;
state0(4:6) = v0;
state0(7:9) = rRTN;
state0(10:12) = vRTN;

% Time span
orbit_period = 2*pi*sqrt(alpha0(1)^3 / mu);
t0 = 0;
t_end =  5 * orbit_period; % Propagate for 2 orbits
tspan = t0:0.1:t_end;

% ODE options with precise tolerances
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% PROPOGATE
[t, state] = ode113(@state_dot, tspan, state0, options);

% PLOT
x1 = state(:, 7);
y1 = state(:, 8);
z1 = state(:, 9);
vx1 = state(:, 10);
vy1 = state(:, 11);
vz1 = state(:, 12);

figure(1)
subplot(3,1,1)
plot(t/3600, x1)
ylabel('R (km)')
subplot(3,1,2)
plot(t/3600, y1)
ylabel('T (km)')
subplot(3,1,3)
plot(t/3600, z1)
ylabel('N (km)')
xlabel('Time (s)')



figure(2)
subplot(3,1,1)
plot(t/3600, vx1)
ylabel('R (km/s)')
subplot(3,1,2)
plot(t/3600, vy1)
ylabel('T (km/s)')
subplot(3,1,3)
plot(t/3600, vz1)
ylabel('N (km/s)')
xlabel('Time (hr)')



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


function [rRTN, vRTN] = ECI2RTN(r_chief_ECI, v_chief_ECI, r_deputy_ECI, v_deputy_ECI)
r0 = r_chief_ECI;
v0 = v_chief_ECI;
r1 = r_deputy_ECI;
v1 = v_deputy_ECI;

% ECI TO CHIEF RTN ROTATION MATRIX
R0_hat = r0 / norm(r0);
N0 = cross(r0, v0);
N0_hat = N0 / norm(N0);
T0 = cross(N0_hat, R0_hat);
T0_hat = T0 / norm(T0);
R_ECI2RTN = [R0_hat'; T0_hat'; N0_hat'];

% CALCULATE DEPUTY RTN COORDS
rRTN = R_ECI2RTN * (r1-r0);
vRTN = R_ECI2RTN * (v1-v0);
end

function [statedot] = state_dot(t, state)
statedot = zeros(12, 1);
mu = 398600.4418;

% CHIEF ECI INTEGRATOR --------------------------
% Unpack state
x0 = state(1);
y0 = state(2);
z0 = state(3);
vx0 = state(4);
vy0 = state(5);
vz0 = state(6);
r0_vec = [x0; y0; z0];
v0_vec = [vx0; vy0; vz0];
r0 = norm(r0_vec);

% Calculate state derivative
statedot(1:3) = state(4:6);     
accel_kepler = -mu / r0^3 * r0_vec;
perturbations = [0; 0; 0];
statedot(4:6) = accel_kepler + perturbations; 

% DEPUTY RTN INTEGRATOR
% Unpack state
x1 = state(7);
y1 = state(8);
z1 = state(9);
x1dot = state(10);
y1dot = state(11);
z1dot = state(12);

% Positional derivative = velocity
statedot(7:9) = state(10:12);

% Compute chief theta derivatives and r derivatives
h = cross(r0_vec, v0_vec);
theta0dot = norm(h) / r0^2;
r0dot = dot(r0_vec, v0_vec) / r0;
theta0dot2 = (-2*r0dot*theta0dot) / r0;

% Calculate theta/derivates and r/derivatives from chief state
statedot(10) = 2*theta0dot*y1dot + theta0dot2*y1 + theta0dot^2*x1 - (mu*(r0+x1)) / ((r0+x1)^2 + y1^2 + z1^2)^1.5 + mu/r0^2;
statedot(11) = -2*theta0dot*x1dot - theta0dot2*x1 + theta0dot^2*y1 - (mu*y1) / ((r0+x1)^2 + y1^2 + z1^2)^1.5;
statedot(12) = -(mu*z1) / ((r0+x1)^2 + y1^2 + z1^2)^1.5;

end