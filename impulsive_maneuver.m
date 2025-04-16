%% INTEGRATE OVER 2 PERIODS AS IN (B) ------------------------------------
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


%% CALCULATE DELTA A AND APPLY VELOCITY CHANGE ---------------------------
r_chief_ECI = state(end, 1:3)';
v_chief_ECI = state(end, 4:6)';
r_deputy_RTN = state(end, 7:9)';
v_deputy_RTN = state(end, 10:12)';
[r_deputy_ECI, v_deputy_ECI] = RTN2ECI(r_chief_ECI, v_chief_ECI, r_deputy_RTN, v_deputy_RTN);

[a_chief, e, i, O, w, nu] = ECItoKepler(r_chief_ECI, v_chief_ECI);
[a_deputy, e, i, O, w, nu] = ECItoKepler(r_deputy_ECI, v_deputy_ECI);
da = a_chief - a_deputy;
dvT_deputy = (mu*da) / (2*a_deputy^2*norm(v_deputy_ECI));

state0_new = state(end, :);
state0_new(11) = v_deputy_RTN(2) + dvT_deputy;

%% RECALCULATE NEW ORBIT
% Time span
orbit_period = 2*pi*sqrt(alpha0(1)^3 / mu);
t0 = 0;
t_end =  5 * orbit_period; % Propagate for 2 orbits
tspan2 = t_end:0.1:2*t_end;

% ODE options with precise tolerances
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% PROPOGATE
[t2, state2] = ode113(@state_dot, tspan2, state0_new, options);

states = [state; state2];
ts = [t; t2];

% PLOT
x1 = states(:, 7);
y1 = states(:, 8);
z1 = states(:, 9);
vx1 = states(:, 10);
vy1 = states(:, 11);
vz1 = states(:, 12);

figure(1)
subplot(3,1,1)
plot(ts/3600, x1)
ylabel('R (km)')
subplot(3,1,2)
plot(ts/3600, y1)
ylabel('T (km)')
subplot(3,1,3)
plot(ts/3600, z1)
ylabel('N (km)')
xlabel('Time (s)')

figure(2)
subplot(3,1,1)
plot(ts/3600, vx1)
ylabel('R (km/s)')
subplot(3,1,2)
plot(ts/3600, vy1)
ylabel('T (km/s)')
subplot(3,1,3)
plot(ts/3600, vz1)
ylabel('N (km/s)')
xlabel('Time (hr)')
%}