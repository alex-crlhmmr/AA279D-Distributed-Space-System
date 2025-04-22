% Initial Keplerian orbital elements [a, e, i, RAAN, omega, nu]
alpha0 = [6771; 0.0005; deg2rad(51.64); deg2rad(257); deg2rad(0); deg2rad(30)]; % Chief
alpha1 = [6771; 0.0006; deg2rad(51.64+0.5); deg2rad(257.5); deg2rad(0.5); deg2rad(25)]; % Deputy

% Calculate initial state [ECI position and velocity] of chief
mu = 398600.4418;
[r0, v0] = utils.OE2ECI(alpha0, mu);

% Calculate initial state [RTN position and velocity] of deputy
[r1, v1] = utils.OE2ECI(alpha1, mu);
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
subplot(2,2,1)
plot(x1, y1)
xlabel('R (km)')
ylabel('T (km)')
axis equal
subplot(2,2,2)
plot(x1, z1)
xlabel('R (km)')
ylabel('N (km)')
axis equal
subplot(2,2,3)
plot(y1, z1)
xlabel('T (km)')
ylabel('N (km)')
axis equal
subplot(2,2,4)
plot3(x1, y1, z1)
xlabel('R (km)')
ylabel('T (km)')
zlabel('N (km)')
axis equal


figure(3)
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
