%% Reset
clc; clear; close all;
%% (a)

% Orbited body
body = 'earth';
const = utils.getConstants({body});

% Chief initial orbital elements
a_0 = 6771000 ;         % semi-major axis [m]
e_0 = 0.0005;           % eccentricity [-]
i_0 = deg2rad(51.64);    % inclination [rad]
W_0 = deg2rad(257);   % RAAN [rad]
w_0 = deg2rad(0);    % argument of perigee [rad]
f_0 = 30;                   % true anomaly [rad]

% Chief initial state in OE & ECI
initial_state_0_OE = [a_0, e_0, i_0, W_0, w_0, f_0];
initial_state_0_ECI = utils.OE2ECI(initial_state_0_OE, const, body);
r_initial_0_ECI = initial_state_0_ECI(1:3);
v_initial_0_ECI = initial_state_0_ECI(4:6);

% Deputy initial orbital elements (ISS)
a_1 = a_0;                   % semi-major axis [m]
e_1 = e_0+0.0001;            % eccentricity [-]
i_1 = i_0 + deg2rad(0.5);   % inclination [rad]
W_1 = W_0 + deg2rad(0.5);   % RAAN [rad]
w_1 = w_0 + deg2rad(0.5);   % argument of perigee [rad]
f_1 = f_0 - deg2rad(5);      % true anomaly [rad]

% Deputy initial state in OE & ECI
initial_state_1_OE = [a_1, e_1, i_1, W_1, w_1, f_1];
initial_state_1_ECI = utils.OE2ECI(initial_state_1_OE, const, body);
r_initial_1_ECI = initial_state_1_ECI(1:3);
v_initial_1_ECI = initial_state_1_ECI(4:6);

% Deputy initial relative state to chief in ECI
initial_relative_state_ECI = initial_state_1_ECI-initial_state_0_ECI;
initial_delta_r_ECI = initial_relative_state_ECI(1:3);
initial_delta_v_ECI = initial_relative_state_ECI(4:6);



%% (b) - proposal 1

function vec_RTN = ECI2RTN(vec_ECI, r_sat_ECI, v_sat_ECI)
    % Compute the radial unit vector (R-hat)
    r_hat = r_sat_ECI / norm(r_sat_ECI);

    % Compute the angular momentum vector and then the normal unit vector (N-hat)
    h_vec = cross(r_sat_ECI, v_sat_ECI);
    n_hat = h_vec / norm(h_vec);

    % Compute the transverse unit vector (T-hat) via the cross product of N-hat and R-hat
    t_hat = cross(n_hat, r_hat);

    % Form the transformation (rotation) matrix from ECI to RTN
    R_ECI2RTN = [r_hat.'; t_hat.'; n_hat.'];

    % Transform the input ECI vector to the RTN frame.
    vec_RTN = R_ECI2RTN * vec_ECI;
end


function statedot = getStatedot(t, state, const, body)
    mu = const.(body).mu;
    statedot = zeros(12, 1);

    % Chiefin ECI 
    r0 = state(1:3);       
    v0 = state(4:6);       
    r0_norm = norm(r0);

    % Deputy relative to chief in chief's RTN ===
    rho = state(7:9);      
    drho = state(10:12);   

    x = rho(1); y = rho(2); z = rho(3);
    dx = drho(1); dy = drho(2); dz = drho(3);

    % Chief motion
    a0 = -mu / r0_norm^3 * r0;
    statedot(1:3) = v0;
    statedot(4:6) = a0;

    % RTN frame angular rates 
    h_vec = cross(r0, v0);
    h = norm(h_vec);
    theta_dot = h / r0_norm^2;
    r0_dot = dot(r0, v0) / r0_norm;
    theta_ddot = -2 * r0_dot * theta_dot / r0_norm;

    % Relative acceleration in rotating RTN frame 
    
    denom = ((r0_norm + x)^2 + y^2 + z^2)^(3/2);

    ax = 2*theta_dot*dy + theta_ddot*y + theta_dot^2*x - mu * ((r0_norm + x)/denom - 1/r0_norm^2);
    ay =  -2*theta_dot*dx - theta_ddot*x+ theta_dot^2*y - mu*y/denom;
    az = -mu*z/denom;
    
    % Relative motion 
    statedot(7:9) = drho;
    statedot(10:12) = [ax; ay; az];
end



% Initial states
initial_delta_r_RTN = ECI2RTN(initial_delta_r_ECI,r_initial_0_ECI, v_initial_0_ECI);
initial_delta_v_RTN = ECI2RTN(initial_delta_v_ECI,r_initial_0_ECI, v_initial_0_ECI);
initial_state = [r_initial_0_ECI, v_initial_0_ECI, initial_delta_r_RTN,initial_delta_v_RTN];

% Orbital period
mu = const.(body).mu;
T = 2 * pi * sqrt(a_0^3 / mu); % Orbital period [s]

% Simulation parameters
tstart = 0;
tint = 10;
tend = 5*T;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% Run the propagator
odefun = @(t, state) getStatedot(t,state,const,body);                               
[t, state_out] = ode113(odefun, (tstart:tint:tend)', initial_state, options);

% Extract relative position and velocity of deputy
relative_pos_RTN = state_out(:, 7:9).';   
relative_vel_RTN = state_out(:,10:12).'; 

% Plot relative position in RTN
figure;
subplot(3,1,1); plot(t, relative_pos_RTN(1,:), 'r'); ylabel('R [m]'); grid on;
title('Relative Position in Chief RTN Frame');
subplot(3,1,2); plot(t, relative_pos_RTN(2,:), 'g'); ylabel('T [m]'); grid on;
subplot(3,1,3); plot(t, relative_pos_RTN(3,:), 'b'); ylabel('N [m]'); xlabel('Time [s]'); grid on;

% Plot relative velocity in RTN
figure;
subplot(3,1,1); plot(t, relative_vel_RTN(1,:), 'r'); ylabel('Ṙ [m/s]'); grid on;
title('Relative Velocity in Chief RTN Frame');
subplot(3,1,2); plot(t, relative_vel_RTN(2,:), 'g'); ylabel('Ṫ [m/s]'); grid on;
subplot(3,1,3); plot(t, relative_vel_RTN(3,:), 'b'); ylabel('Ṅ [m/s]'); xlabel('Time [s]'); grid on;


% Animation

N = length(t);
r_chief_ECI = state_out(:,1:3);
v_chief_ECI = state_out(:,4:6);
delta_r_RTN = state_out(:,7:9);

r_deputy_ECI = zeros(N, 3);

for k = 1:N
    rc = r_chief_ECI(k, :)';  % 3x1
    vc = v_chief_ECI(k, :)';  % 3x1

    r_hat = rc / norm(rc);
    h_vec = cross(rc, vc);
    n_hat = h_vec / norm(h_vec);
    t_hat = cross(n_hat, r_hat);

    R_RTN2ECI = [r_hat, t_hat, n_hat];  

    rho_RTN = delta_r_RTN(k, :)';  

    r_deputy_ECI(k, :) = (rc + R_RTN2ECI * rho_RTN)';  
end


R_earth = const.(body).R;
[xe, ye, ze] = sphere(100);

figure;
hold on;
axis equal;
grid on;
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('Chief and Deputy Orbits in ECI Frame (Animated)');
view(35, 25);

surf(R_earth*xe, R_earth*ye, R_earth*ze, ...
    'FaceColor', 'cyan', 'EdgeAlpha', 0.2, 'FaceAlpha', 0.3);

plot3(r_chief_ECI(:,1), r_chief_ECI(:,2), r_chief_ECI(:,3), 'r:', 'LineWidth', 0.5);
plot3(r_deputy_ECI(:,1), r_deputy_ECI(:,2), r_deputy_ECI(:,3), 'b:', 'LineWidth', 0.5);

h_chief = plot3(NaN, NaN, NaN, 'ro', 'MarkerFaceColor', 'r');
h_deputy = plot3(NaN, NaN, NaN, 'bo', 'MarkerFaceColor', 'b');

legend('Earth', 'Chief Trajectory', 'Deputy Trajectory', 'Chief', 'Deputy');

for k = 1:10:N
    rc = r_chief_ECI(k,:);
    rd = r_deputy_ECI(k,:);
    
    set(h_chief, 'XData', rc(1), 'YData', rc(2), 'ZData', rc(3));
    set(h_deputy, 'XData', rd(1), 'YData', rd(2), 'ZData', rd(3));
    
    drawnow;
end

%% (b) - proposal 2

% Initial Keplerian orbital elements [a, e, i, RAAN, omega, nu]
alpha0 = [6771; 0.0005; deg2rad(51.64); deg2rad(257); 0; deg2rad(30)]; % Chief
alpha1 = [6771; 0.0006; deg2rad(52.14); deg2rad(257.5); deg2rad(0.5); deg2rad(25)]; % Deputy

% Calculate initial state [ECI position and velocity] of chief
mu = 398600.4418;
[r0, v0] = kepler_to_ijk(alpha0, mu);

% Calculate initial state [RTN position and velocity] of deputy
[r1, v1] = kepler_to_ijk(alpha1, mu);
[rRTN, vRTN] = ECI2RTN2(r0, v0, r1, v1);

% Initialize state
state0 = zeros(12, 1);
state0(1:3) = r0;
state0(4:6) = v0;
state0(7:9) = rRTN;
state0(10:12) = vRTN;

% Time span
orbit_period = 2*pi*sqrt(alpha0(1)^3 / mu);
t0 = 0;
t_end =  5 * orbit_period; 
tspan = t0:0.1:t_end;

% ODE options with precise tolerances
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% PROPOGATE
[t_b, state_b] = ode113(@state_dot, tspan, state0, options);

% PLOT
x1 = state_b(:, 7);
y1 = state_b(:, 8);
z1 = state_b(:, 9);
vx1 = state_b(:, 10);
vy1 = state_b(:, 11);
vz1 = state_b(:, 12);

figure(1)
subplot(3,1,1)
plot(t_b/3600, x1)
ylabel('R [km]')
subplot(3,1,2)
plot(t_b/3600, y1)
ylabel('T [km]')
subplot(3,1,3)
plot(t_b/3600, z1)
ylabel('N [km]')
xlabel('Time (hr)')



figure(2)
subplot(3,1,1)
plot(t_b/3600, vx1)
ylabel('R [km/s]')
subplot(3,1,2)
plot(t_b/3600, vy1)
ylabel('T [km/s]')
subplot(3,1,3)
plot(t_b/3600, vz1)
ylabel('N [km/s]')
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


function [rRTN, vRTN] = ECI2RTN2(r_chief_ECI, v_chief_ECI, r_deputy_ECI, v_deputy_ECI)
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


%% (c)

% Run the unperturbed FODE propagator for chief in ECI
perturbated = false;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);                               
[t_c, state_ECI_0_c] = ode113(odefun, (tstart:tint:tend)', initial_state_0_ECI, options);


% Run the unperturbed FODE for deputy in ECI
perturbated = false;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);                               
[t_c, state_ECI_1_c] = ode113(odefun, (tstart:tint:tend)', initial_state_1_ECI, options);


relative_state_ECI = state_ECI_1_c - state_ECI_0_c;

% Preallocate RTN vectors over time
N = length(t);
relative_pos_RTN = zeros(3, N);
relative_vel_RTN = zeros(3, N);

for k = 1:N
    % Extract chief and deputy ECI position and velocity
    r_chief = state_ECI_0_c(k, 1:3).';
    v_chief = state_ECI_0_c(k, 4:6).';
    
    r_deputy = state_ECI_1_c(k, 1:3).';
    v_deputy = state_ECI_1_c(k, 4:6).';
    
    % Compute relative position and velocity in ECI
    dr_ECI = r_deputy - r_chief;
    dv_ECI = v_deputy - v_chief;
    
    % Transform to chief's RTN frame
    dr_RTN = ECI2RTN(dr_ECI, r_chief, v_chief);
    dv_RTN = ECI2RTN(dv_ECI, r_chief, v_chief);
    
    % Relative psoition and velocity of deputy to chief in chief's RTN
    relative_pos_RTN(:, k) = dr_RTN;
    relative_vel_RTN(:, k) = dv_RTN;
end

% Position in RTN
figure;
subplot(3,1,1); plot(t_c/3600, relative_pos_RTN(1,:)/1e3, 'r'); ylabel('R [km]'); grid on;
subplot(3,1,2); plot(t_c/3600, relative_pos_RTN(2,:)/1e3, 'g'); ylabel('T [km]'); grid on;
subplot(3,1,3); plot(t_c/3600, relative_pos_RTN(3,:)/1e3, 'b'); ylabel('N [km]'); xlabel('Time [hr]'); grid on;

% Velocity in RTN
figure;
subplot(3,1,1); plot(t_c/3600, relative_vel_RTN(1,:)/1e3, 'r'); ylabel('Ṙ [km/s]'); grid on;
subplot(3,1,2); plot(t_c/3600, relative_vel_RTN(2,:)/1e3, 'g'); ylabel('Ṫ [km/s]'); grid on;
subplot(3,1,3); plot(t_c/3600, relative_vel_RTN(3,:)/1e3, 'b'); ylabel('Ṅ [km/s]'); xlabel('Time [hr]'); grid on;


% Animation of chief and deputy orbiting Earth in ECI

% Earth setup
R_earth = const.(body).R;
[xe, ye, ze] = sphere(100);

% Extract satellite positions
r_chief_ECI = state_ECI_0_c(:,1:3);
r_deputy_ECI = state_ECI_1_c(:,1:3);
N = length(t);

figure;
hold on;
axis equal;
grid on;
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('Chief and Deputy Orbits in ECI Frame (Animated)');
view(35, 25);

% Plot Earth
surf(R_earth*xe, R_earth*ye, R_earth*ze, ...
    'FaceColor', 'cyan', 'EdgeAlpha', 0.2, 'FaceAlpha', 0.3);

% Plot full trajectories (transparent as reference)
plot3(r_chief_ECI(:,1), r_chief_ECI(:,2), r_chief_ECI(:,3), 'r:', 'LineWidth', 0.5);
plot3(r_deputy_ECI(:,1), r_deputy_ECI(:,2), r_deputy_ECI(:,3), 'b:', 'LineWidth', 0.5);

% Initialize satellite markers
h_chief = plot3(NaN, NaN, NaN, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Chief');
h_deputy = plot3(NaN, NaN, NaN, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Deputy');

legend;

% Animate
for k = 1:10:N
    rc = r_chief_ECI(k, :);
    rd = r_deputy_ECI(k, :);
    
    % Update satellite positions
    set(h_chief, 'XData', rc(1), 'YData', rc(2), 'ZData', rc(3));
    set(h_deputy, 'XData', rd(1), 'YData', rd(2), 'ZData', rd(3));
    
    drawnow;
end

% Convert to km
x1_km = x1 ;
y1_km = y1 ;
z1_km = z1 ;

rR_km = relative_pos_RTN(1,:) / 1e3;
rT_km = relative_pos_RTN(2,:) / 1e3;
rN_km = relative_pos_RTN(3,:) / 1e3;

% Time in hours
t_hr = t_b / 3600;  % from section (b) — same length as state_b
t_c_hr = t_c / 3600; % from section (c)

% If time arrays differ, interpolate section (b) onto t_c
if ~isequal(t_hr, t_c_hr)
    x1_km = interp1(t_hr, x1_km, t_c_hr, 'spline');
    y1_km = interp1(t_hr, y1_km, t_c_hr, 'spline');
    z1_km = interp1(t_hr, z1_km, t_c_hr, 'spline');
    t_hr = t_c_hr;  % Use t_c as master time
end

figure;

subplot(3,1,1);
plot(t_hr, x1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_hr, rR_km, 'r--', 'LineWidth', 1.2);
ylabel('R [km]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

subplot(3,1,2);
plot(t_hr, y1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_hr, rT_km, 'r--', 'LineWidth', 1.2);
ylabel('T [km]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

subplot(3,1,3);
plot(t_hr, z1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_hr, rN_km, 'r--', 'LineWidth', 1.2);
ylabel('N [km]');
xlabel('Time [hr]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

% Convert to km/s
vx1_km = vx1;
vy1_km = vy1;
vz1_km = vz1;

vR_km = relative_vel_RTN(1,:) / 1e3;
vT_km = relative_vel_RTN(2,:) / 1e3;
vN_km = relative_vel_RTN(3,:) / 1e3;

% Use same time as before: t_hr == t_c_hr
if ~isequal(t_b, t_c)
    vx1_km = interp1(t_b / 3600, vx1_km, t_c_hr, 'spline');
    vy1_km = interp1(t_b / 3600, vy1_km, t_c_hr, 'spline');
    vz1_km = interp1(t_b / 3600, vz1_km, t_c_hr, 'spline');
end

figure;

subplot(3,1,1);
plot(t_c_hr, vx1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_c_hr, vR_km, 'r--', 'LineWidth', 1.2);
ylabel('Ṙ [km/s]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

subplot(3,1,2);
plot(t_c_hr, vy1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_c_hr, vT_km, 'r--', 'LineWidth', 1.2);
ylabel('Ṫ [km/s]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

subplot(3,1,3);
plot(t_c_hr, vz1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_c_hr, vN_km, 'r--', 'LineWidth', 1.2);
ylabel('Ṅ [km/s]');
xlabel('Time [hr]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;


%% (d)

% Initial Keplerian orbital elements [a, e, i, RAAN, omega, nu]
alpha0 = [6771; 0.0005; deg2rad(51.64); deg2rad(257); 0; deg2rad(30)]; % Chief
alpha1 = [6771+30; 0.0006; deg2rad(52.14); deg2rad(257.5); deg2rad(0.5); deg2rad(25)]; % Deputy

% Calculate initial state [ECI position and velocity] of chief
mu = 398600.4418;
[r0, v0] = kepler_to_ijk(alpha0, mu);

% Calculate initial state [RTN position and velocity] of deputy
[r1, v1] = kepler_to_ijk(alpha1, mu);
[rRTN, vRTN] = ECI2RTN2(r0, v0, r1, v1);

% Initialize state
state0 = zeros(12, 1);
state0(1:3) = r0;
state0(4:6) = v0;
state0(7:9) = rRTN;
state0(10:12) = vRTN;

% Time span
orbit_period = 2*pi*sqrt(alpha0(1)^3 / mu);
t0 = 0;
t_end =  5 * orbit_period; 
tspan = t0:0.1:t_end;

% ODE options with precise tolerances
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% PROPOGATE
[t_b, state_b] = ode113(@state_dot, tspan, state0, options);


% Deputy initial orbital elements (ISS)
a_1 = a_0;                   % semi-major axis [m]
e_1 = e_0+0.0001;            % eccentricity [-]
i_1 = i_0 + deg2rad(0.5);   % inclination [rad]
W_1 = W_0 + deg2rad(0.5);   % RAAN [rad]
w_1 = w_0 + deg2rad(0.5);   % argument of perigee [rad]
f_1 = f_0 - deg2rad(5);      % true anomaly [rad]

% Deputy initial state in OE & ECI
initial_state_1_OE = [a_1+30000, e_1, i_1, W_1, w_1, f_1];
initial_state_1_ECI = utils.OE2ECI(initial_state_1_OE, const, body);
r_initial_1_ECI = initial_state_1_ECI(1:3);
v_initial_1_ECI = initial_state_1_ECI(4:6);


% Run the unperturbed FODE propagator for chief in ECI
perturbated = false;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);                               
[t_c, state_ECI_0_c] = ode113(odefun, (tstart:tint:tend)', initial_state_0_ECI, options);


% Run the unperturbed FODE for deputy in ECI
perturbated = false;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);                               
[t_c, state_ECI_1_c] = ode113(odefun, (tstart:tint:tend)', initial_state_1_ECI, options);


relative_state_ECI = state_ECI_1_c - state_ECI_0_c;

% Preallocate RTN vectors over time
N = length(t);
relative_pos_RTN = zeros(3, N);
relative_vel_RTN = zeros(3, N);

for k = 1:N
    % Extract chief and deputy ECI position and velocity
    r_chief = state_ECI_0_c(k, 1:3).';
    v_chief = state_ECI_0_c(k, 4:6).';
    
    r_deputy = state_ECI_1_c(k, 1:3).';
    v_deputy = state_ECI_1_c(k, 4:6).';
    
    % Compute relative position and velocity in ECI
    dr_ECI = r_deputy - r_chief;
    dv_ECI = v_deputy - v_chief;
    
    % Transform to chief's RTN frame
    dr_RTN = ECI2RTN(dr_ECI, r_chief, v_chief);
    dv_RTN = ECI2RTN(dv_ECI, r_chief, v_chief);
    
    % Relative psoition and velocity of deputy to chief in chief's RTN
    relative_pos_RTN(:, k) = dr_RTN;
    relative_vel_RTN(:, k) = dv_RTN;
end



% Convert to km
x1_km = x1 ;
y1_km = y1 ;
z1_km = z1 ;

rR_km = relative_pos_RTN(1,:) / 1e3;
rT_km = relative_pos_RTN(2,:) / 1e3;
rN_km = relative_pos_RTN(3,:) / 1e3;

% Time in hours
t_hr = t_b / 3600;  % from section (b) — same length as state_b
t_c_hr = t_c / 3600; % from section (c)

% If time arrays differ, interpolate section (b) onto t_c
if ~isequal(t_hr, t_c_hr)
    x1_km = interp1(t_hr, x1_km, t_c_hr, 'spline');
    y1_km = interp1(t_hr, y1_km, t_c_hr, 'spline');
    z1_km = interp1(t_hr, z1_km, t_c_hr, 'spline');
    t_hr = t_c_hr;  % Use t_c as master time
end

figure;

subplot(3,1,1);
plot(t_hr, x1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_hr, rR_km, 'r--', 'LineWidth', 1.2);
ylabel('R [km]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

subplot(3,1,2);
plot(t_hr, y1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_hr, rT_km, 'r--', 'LineWidth', 1.2);
ylabel('T [km]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

subplot(3,1,3);
plot(t_hr, z1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_hr, rN_km, 'r--', 'LineWidth', 1.2);
ylabel('N [km]');
xlabel('Time [hr]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

% Convert to km/s
vx1_km = vx1;
vy1_km = vy1;
vz1_km = vz1;

vR_km = relative_vel_RTN(1,:) / 1e3;
vT_km = relative_vel_RTN(2,:) / 1e3;
vN_km = relative_vel_RTN(3,:) / 1e3;

% Use same time as before: t_hr == t_c_hr
if ~isequal(t_b, t_c)
    vx1_km = interp1(t_b / 3600, vx1_km, t_c_hr, 'spline');
    vy1_km = interp1(t_b / 3600, vy1_km, t_c_hr, 'spline');
    vz1_km = interp1(t_b / 3600, vz1_km, t_c_hr, 'spline');
end

figure;

subplot(3,1,1);
plot(t_c_hr, vx1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_c_hr, vR_km, 'r--', 'LineWidth', 1.2);
ylabel('Ṙ [km/s]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

subplot(3,1,2);
plot(t_c_hr, vy1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_c_hr, vT_km, 'r--', 'LineWidth', 1.2);
ylabel('Ṫ [km/s]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

subplot(3,1,3);
plot(t_c_hr, vz1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_c_hr, vN_km, 'r--', 'LineWidth', 1.2);
ylabel('Ṅ [km/s]');
xlabel('Time [hr]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;


%% (f)

% Orbited body
body = 'earth';
const = utils.getConstants({body});

% Chief initial orbital elements
a_0 = 6771000 ;         % semi-major axis [m]
e_0 = 0.0005;           % eccentricity [-]
i_0 = deg2rad(51.64);    % inclination [rad]
W_0 = deg2rad(257);   % RAAN [rad]
w_0 = deg2rad(0);    % argument of perigee [rad]
f_0 = 30;                   % true anomaly [rad]

% Chief initial state in OE & ECI
initial_state_0_OE = [a_0, e_0, i_0, W_0, w_0, f_0];
initial_state_0_ECI = utils.OE2ECI(initial_state_0_OE, const, body);

% Deputy initial orbital elements (ISS)
a_1 = a_0+50000;                   % semi-major axis [m]
e_1 = e_0+0.0001;            % eccentricity [-]
i_1 = i_0 + deg2rad(0.5);   % inclination [rad]
W_1 = W_0 + deg2rad(0.5);   % RAAN [rad]
w_1 = w_0 + deg2rad(0.5);   % argument of perigee [rad]
f_1 = f_0 - deg2rad(5);      % true anomaly [rad]

% Deputy initial state in OE & ECI
initial_state_1_OE = [a_1, e_1, i_1, W_1, w_1, f_1];
initial_state_1_ECI = utils.OE2ECI(initial_state_1_OE, const, body);

% Deputy initial relative state to chief in ECI
initial_relative_state_ECI = initial_state_1_ECI-initial_state_0_ECI;
initial_delta_r_ECI = initial_relative_state_ECI(1:3);
initial_delta_v_ECI = initial_relative_state_ECI(4:6);

% Orbital period
mu = const.(body).mu;
T = 2 * pi * sqrt(a_0^3 / mu); % Orbital period [s]

% Simulation parameters
tstart = 0;
tint = 10;
tend = 2*T;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% Run the unperturbed FODE propagator for chief in ECI
perturbated = false;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);                               
[t, state_ECI_0] = ode113(odefun, (tstart:tint:tend)', initial_state_0_ECI, options);


% Run the unperturbed FODE for deputy in ECI
perturbated = false;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);                               
[t, state_ECI_1] = ode113(odefun, (tstart:tint:tend)', initial_state_1_ECI, options);


relative_state_ECI = state_ECI_1 - state_ECI_0;

% Preallocate RTN vectors over time
N = length(t);
relative_pos_RTN = zeros(3, N);
relative_vel_RTN = zeros(3, N);

for k = 1:N
    % Extract chief and deputy ECI position and velocity
    r_chief = state_ECI_0(k, 1:3).';
    v_chief = state_ECI_0(k, 4:6).';
    
    r_deputy = state_ECI_1(k, 1:3).';
    v_deputy = state_ECI_1(k, 4:6).';
    
    % Compute relative position and velocity in ECI
    dr_ECI = r_deputy - r_chief;
    dv_ECI = v_deputy - v_chief;
    
    % Transform to chief's RTN frame
    dr_RTN = ECI2RTN(dr_ECI, r_chief, v_chief);
    dv_RTN = ECI2RTN(dv_ECI, r_chief, v_chief);
    
    % Relative psoition and velocity of deputy to chief in chief's RTN
    relative_pos_RTN(:, k) = dr_RTN;
    relative_vel_RTN(:, k) = dv_RTN;
end

% Position in RTN
figure;
subplot(3,1,1); plot(t/3600, relative_pos_RTN(1,:)/1e3, 'r'); ylabel('R [km]'); grid on;
subplot(3,1,2); plot(t/3600, relative_pos_RTN(2,:)/1e3, 'g'); ylabel('T [km]'); grid on;
subplot(3,1,3); plot(t/3600, relative_pos_RTN(3,:)/1e3, 'b'); ylabel('N [km]'); xlabel('Time [hr]'); grid on;

% Velocity in RTN
figure;
subplot(3,1,1); plot(t/3600, relative_vel_RTN(1,:)/1e3, 'r'); ylabel('Ṙ [km/s]'); grid on;
subplot(3,1,2); plot(t/3600, relative_vel_RTN(2,:)/1e3, 'g'); ylabel('Ṫ [km/s]'); grid on;
subplot(3,1,3); plot(t/3600, relative_vel_RTN(3,:)/1e3, 'b'); ylabel('Ṅ [km/s]'); xlabel('Time [hr]'); grid on;


% Animation of chief and deputy orbiting Earth in ECI

% Earth setup
R_earth = const.(body).R;
[xe, ye, ze] = sphere(100);

% Extract satellite positions
r_chief_ECI = state_ECI_0_c(:,1:3);
r_deputy_ECI = state_ECI_1_c(:,1:3);
N = length(t);

figure;
hold on;
axis equal;
grid on;
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('Chief and Deputy Orbits in ECI Frame (Animated)');
view(35, 25);

% Plot Earth
surf(R_earth*xe, R_earth*ye, R_earth*ze, ...
    'FaceColor', 'cyan', 'EdgeAlpha', 0.2, 'FaceAlpha', 0.3);

% Plot full trajectories (transparent as reference)
plot3(r_chief_ECI(:,1), r_chief_ECI(:,2), r_chief_ECI(:,3), 'r:', 'LineWidth', 0.5);
plot3(r_deputy_ECI(:,1), r_deputy_ECI(:,2), r_deputy_ECI(:,3), 'b:', 'LineWidth', 0.5);

% Initialize satellite markers
h_chief = plot3(NaN, NaN, NaN, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Chief');
h_deputy = plot3(NaN, NaN, NaN, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Deputy');

legend;

% Animate
for k = 1:10:N
    rc = r_chief_ECI(k, :);
    rd = r_deputy_ECI(k, :);
    
    % Update satellite positions
    set(h_chief, 'XData', rc(1), 'YData', rc(2), 'ZData', rc(3));
    set(h_deputy, 'XData', rd(1), 'YData', rd(2), 'ZData', rd(3));
    
    drawnow;
end

% Convert to km
x1_km = x1 ;
y1_km = y1 ;
z1_km = z1 ;

rR_km = relative_pos_RTN(1,:) / 1e3;
rT_km = relative_pos_RTN(2,:) / 1e3;
rN_km = relative_pos_RTN(3,:) / 1e3;

% Time in hours
t_hr = t_b / 3600;  % from section (b) — same length as state_b
t_c_hr = t_c / 3600; % from section (c)

% If time arrays differ, interpolate section (b) onto t_c
if ~isequal(t_hr, t_c_hr)
    x1_km = interp1(t_hr, x1_km, t_c_hr, 'spline');
    y1_km = interp1(t_hr, y1_km, t_c_hr, 'spline');
    z1_km = interp1(t_hr, z1_km, t_c_hr, 'spline');
    t_hr = t_c_hr;  % Use t_c as master time
end

figure;

subplot(3,1,1);
plot(t_hr, x1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_hr, rR_km, 'r--', 'LineWidth', 1.2);
ylabel('R [km]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

subplot(3,1,2);
plot(t_hr, y1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_hr, rT_km, 'r--', 'LineWidth', 1.2);
ylabel('T [km]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

subplot(3,1,3);
plot(t_hr, z1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_hr, rN_km, 'r--', 'LineWidth', 1.2);
ylabel('N [km]');
xlabel('Time [hr]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

% Convert to km/s
vx1_km = vx1;
vy1_km = vy1;
vz1_km = vz1;

vR_km = relative_vel_RTN(1,:) / 1e3;
vT_km = relative_vel_RTN(2,:) / 1e3;
vN_km = relative_vel_RTN(3,:) / 1e3;

% Use same time as before: t_hr == t_c_hr
if ~isequal(t_b, t_c)
    vx1_km = interp1(t_b / 3600, vx1_km, t_c_hr, 'spline');
    vy1_km = interp1(t_b / 3600, vy1_km, t_c_hr, 'spline');
    vz1_km = interp1(t_b / 3600, vz1_km, t_c_hr, 'spline');
end

figure;

subplot(3,1,1);
plot(t_c_hr, vx1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_c_hr, vR_km, 'r--', 'LineWidth', 1.2);
ylabel('Ṙ [km/s]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

subplot(3,1,2);
plot(t_c_hr, vy1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_c_hr, vT_km, 'r--', 'LineWidth', 1.2);
ylabel('Ṫ [km/s]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;

subplot(3,1,3);
plot(t_c_hr, vz1_km, 'b', 'LineWidth', 1.2); hold on;
plot(t_c_hr, vN_km, 'r--', 'LineWidth', 1.2);
ylabel('Ṅ [km/s]');
xlabel('Time [hr]');
legend('Nonlinear (b)', 'FODE diff (c)');
grid on;