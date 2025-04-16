% Initial Keplerian orbital elements [a, e, i, RAAN, omega, nu]
alpha0 = [6771; 0.0005; deg2rad(51.64); deg2rad(257); 0; deg2rad(30)]; % Chief
alpha1 = [6771; 0.0006; deg2rad(51.69); deg2rad(257.5); deg2rad(0.5); deg2rad(25)]; % Deputy

% Calculate initial state [ECI position and velocity] of chief
mu = 398600.4418;
[r0, v0] = kepler_to_ijk(alpha0, mu);

% Calculate initial RTN position and velocity of deputy
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
t_end = 1 * orbit_period; % Propagate for 2 orbits
tspan = t0:0.1:t_end;

% ODE options with precise tolerances
options = odeset('RelTol', 1e-20, 'AbsTol', 1e-20);

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
plot(x1, y1)
xlabel('R (km)')
ylabel('T (km)')
axis equal

figure(2)
plot(vx1, vy1)
xlabel('R (km/s)')
ylabel('T (km/s)')
axis equal
