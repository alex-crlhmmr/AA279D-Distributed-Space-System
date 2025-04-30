%% Reset

clc; clear; close all;
%% Q1

% Orbited body
body = 'moon';
const = utils.getConstants({body});

% Initial orbital elements
h_0 = 200e3;              % orbit height [m]
R_0 = const.(body).R;     % body radius [m]
a_0 = R_0 + h_0;          % semi-major axis [m]
e_0 = 0.01;               % eccentricity [-]
i_0 = deg2rad(45);        % inclination [rad]
W_0 = deg2rad(20);         % RAAN [rad]
w_0 = deg2rad(20);        % argument of perigee [rad]
f_0 = 0;                  % true anomaly [rad]

%% Q2

% Initial state in OE
state_OE_0 = [a_0, e_0, i_0, W_0, w_0, f_0];

% Initial state in ECI
state_ECI_0 = utils.OE2ECI(state_OE_0, const, body)

%% Q3

% Orbital period
mu = const.(body).mu;
T = 2 * pi * sqrt(a_0^3 / mu); % Orbital period [s]

% Simulation parameters
tstart = 0;
tint = 10;
tend = 20*T;
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Run the propagator
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body);                               
[t, state] = ode113(odefun, (tstart:tint:tend)', state_ECI_0, options);

% Extract position and velocity vectors
r = state(:, 1:3);
v = state(:, 4:6);

% Plot position over time
figure;
subplot(2,1,1);
plot(t/3600, vecnorm(r, 2, 2) / 1e3); % distance in km
xlabel('Time [h]');
ylabel('Radius [km]');
title('Position magnitude over time');
grid on;

% Plot velocity over time
subplot(2,1,2);
plot(t/3600, vecnorm(v, 2, 2) / 1e3); % speed in km/s
xlabel('Time [h]');
ylabel('Velocity [km/s]');
title('Velocity magnitude over time');
grid on;

% Plot 3D orbit in ECI
utils.plotOrbit3D_ECI(r, const, body);

%% Q4

% Run the unpertubed propagator
perturbated = false;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);   
[t, state] = ode113(odefun, (tstart:tint:tend)', state_ECI_0, options);

% Extract position and velocity vectors
r_unpertubed = state(:, 1:3);
v_unpertubed = state(:, 4:6);

function [r_ECI_all, v_ECI_all] = propagateKeplerian(state_OE_0, t, const, body)
    % Extract constants
    mu = const.(body).mu;
    
    % Unpack initial orbital elements
    a = state_OE_0(1);
    e = state_OE_0(2);
    i = state_OE_0(3);
    Omega = state_OE_0(4);
    omega = state_OE_0(5);
    f0 = state_OE_0(6);
    
    % Compute initial eccentric anomaly
    E0 = 2 * atan( sqrt((1 - e) / (1 + e)) * tan(f0 / 2) );
    if E0 < 0
        E0 = E0 + 2*pi; % Ensure positive
    end

    % Mean motion
    n = sqrt(mu / a^3);
    
    % Initial mean anomaly
    M0 = E0 - e * sin(E0);
    
    % Prepare outputs
    r_ECI_all = zeros(length(t), 3);
    v_ECI_all = zeros(length(t), 3);
    
    % Loop over each time step
    for idx = 1:length(t)
        % Current mean anomaly
        M = M0 + n * (t(idx) - t(1));

        % Solve Kepler's Equation for E
        E = solveKepler(M, e);
        
        % True anomaly
        f = 2 * atan2( sqrt(1 + e) * sin(E/2), sqrt(1 - e) * cos(E/2) );
        
        % Update state vector in perifocal frame
        state_OE = [a; e; i; Omega; omega; f];
        state_perifocal = utils.OE2Perifocal(state_OE, const, body);
        
        % Perifocal to ECI
        Q = perifocal2ECI(i, Omega, omega);
        
        % Separate position and velocity
        r_perifocal = state_perifocal(1:3);
        v_perifocal = state_perifocal(4:6);
        
        % Rotate to ECI frame
        r_ECI = Q * r_perifocal;
        v_ECI = Q * v_perifocal;
        
        % Store
        r_ECI_all(idx, :) = r_ECI';
        v_ECI_all(idx, :) = v_ECI';
    end
end

function E = solveKepler(M, e)
    % Solve Kepler's equation: M = E - e*sin(E)
    % Initial guess
    E = M;
    
    % Iterative solution (Newton-Raphson)
    tol = 1e-10;
    ratio = 1;
    while abs(ratio) > tol
        ratio = (E - e * sin(E) - M) / (1 - e * cos(E));
        E = E - ratio;
    end
end

function Q = perifocal2ECI(i, Omega, omega)
    % Rotation matrix from perifocal to ECI frame
    R3_Omega = [cos(Omega), -sin(Omega), 0;
                sin(Omega),  cos(Omega), 0;
                         0,           0, 1];
                     
    R1_i = [1,        0,         0;
            0, cos(i), -sin(i);
            0, sin(i),  cos(i)];
        
    R3_omega = [cos(omega), -sin(omega), 0;
                sin(omega),  cos(omega), 0;
                          0,           0, 1];
                      
    Q = (R3_Omega * R1_i * R3_omega)';
end

% Analytical propagation (Keplerian, unperturbed)
[r_analytical, v_analytical] = propagateKeplerian(state_OE_0, t, const, body);


% Errors in ECI frame
error_r = r - r_analytical;
error_v = v - v_analytical;


% Preallocate
error_r_RTN = zeros(size(error_r));
error_v_RTN = zeros(size(error_v));

for idx = 1:length(t)
    r_vec = r(idx, :)';
    v_vec = v(idx, :)';
    
    R_hat = r_vec / norm(r_vec);
    T_hat = cross([0; 0; 1], R_hat); % Approximate transverse
    if norm(T_hat) < 1e-8
        T_hat = v_vec / norm(v_vec); % use velocity if near poles
    else
        T_hat = T_hat / norm(T_hat);
    end
    N_hat = cross(R_hat, T_hat);
    
    % Build RTN frame
    RTN = [R_hat, T_hat, N_hat];
    
    % Express error vectors in RTN
    error_r_RTN(idx, :) = RTN' * error_r(idx, :)';
    error_v_RTN(idx, :) = RTN' * error_v(idx, :)';
end


figure;
subplot(2,1,1);
plot(t/3600, abs(error_r_RTN) / 1e3); % km
xlabel('Time [h]');
ylabel('Position error [km]');
legend('Radial', 'Transverse', 'Normal');
title('Position error in RTN frame');
grid on;

subplot(2,1,2);
plot(t/3600, abs(error_v_RTN) / 1e3); % km/s
xlabel('Time [h]');
ylabel('Velocity error [km/s]');
legend('Radial', 'Transverse', 'Normal');
title('Velocity error in RTN frame');
grid on;




%% Q5 - Unperturbed

% Run the unpertubed propagator
perturbated = false;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);   
[t, state] = ode113(odefun, (tstart:tint:tend)', state_ECI_0, options);

% Extract position and velocity vectors
r = state(:, 1:3);
v = state(:, 4:6);

% Preallocate OEs array
numSteps = length(t);
elements = zeros(numSteps, 6);
h_mag = zeros(numSteps, 1);
e_mag = zeros(numSteps, 1);
energy = zeros(numSteps, 1);

% Compute OEs, angular momentum, eccentricty vector,specific mechanical energy at each time step
for idx = 1:numSteps

    % Position and velocity
    r_vec = r(idx, :);
    v_vec = v(idx, :);

    % OEs
    [a, e, i, W, w, f] = utils.ECI2OE(r_vec, v_vec, const, body);
    elements(idx, :) = [a, e, i, W, w, f];

    % Angular momentum vector and magnitude
    h_vec = cross(r_vec, v_vec);
    h_mag(idx) = norm(h_vec);

    % Eccentricity vector and magnitude
    e_vec = cross(v_vec, h_vec) / mu - r_vec / norm(r_vec);
    e_mag(idx) = norm(e_vec);
    
    % Specific mechanical energy
    energy(idx) = norm(v_vec)^2 / 2 - mu / norm(r_vec);
end

% Unpack OEs
a = elements(:, 1);
e = elements(:, 2);
i = elements(:, 3);
W = elements(:, 4);
w = elements(:, 5);
f = unwrap(elements(:, 6));

% Plot OES
figure;

subplot(3,2,1);
plot(t/3600, a/1e3);
xlabel('Time [h]'); ylabel('Semi-major axis [km]'); grid on;

subplot(3,2,2);
plot(t/3600, e);
xlabel('Time [h]'); ylabel('Eccentricity'); grid on;

subplot(3,2,3);
plot(t/3600, i);
xlabel('Time [h]'); ylabel('Inclination [deg]'); grid on;

subplot(3,2,4);
plot(t/3600, W);
xlabel('Time [h]'); ylabel('RAAN [deg]'); grid on;

subplot(3,2,5);
plot(t/3600, w);
xlabel('Time [h]'); ylabel('Argument of periapsis [deg]'); grid on;

subplot(3,2,6);
plot(t/3600, f);
xlabel('Time [h]'); ylabel('True anomaly [deg]'); grid on;

sgtitle('Orbital elements over time');


% Plot the angular momentum, eccentricty vector and specific mechanical energy over time
figure;

subplot(3,1,1);
plot(t/3600, h_mag);
xlabel('Time [h]'); ylabel('Angular momentum [m^2/s]'); grid on;

subplot(3,1,2);
plot(t/3600, e_mag);
xlabel('Time [h]'); ylabel('Eccentricity vector magnitude'); grid on;

subplot(3,1,3);
plot(t/3600, energy/1e6);
xlabel('Time [h]'); ylabel('Specific mechanical energy [MJ/kg]'); grid on;

sgtitle('Angular Momentum, Eccentricity Vector, Specific Mechanical Energy');

%% Q5 - Perturbed

% Run the perturbed propagator
perturbated = true;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);   
[t, state] = ode113(odefun, (tstart:tint:tend)', state_ECI_0, options);

% Extract position and velocity vectors
r = state(:, 1:3);
v = state(:, 4:6);

% Preallocate OEs array
numSteps = length(t);
elements = zeros(numSteps, 6);
h_mag = zeros(numSteps, 1);
e_mag = zeros(numSteps, 1);
energy = zeros(numSteps, 1);

% Compute OEs, angular momentum, eccentricty vector,specific mechanical energy at each time step
for idx = 1:numSteps

    % Position and velocity
    r_vec = r(idx, :);
    v_vec = v(idx, :);

    % OEs
    [a, e, i, W, w, f] = utils.ECI2OE(r_vec, v_vec, const, body);
    elements(idx, :) = [a, e, i, W, w, f];

    % Angular momentum vector and magnitude
    h_vec = cross(r_vec, v_vec);
    h_mag(idx) = norm(h_vec);

    % Eccentricity vector and magnitude
    e_vec = cross(v_vec, h_vec) / mu - r_vec / norm(r_vec);
    e_mag(idx) = norm(e_vec);
    
    % Specific mechanical energy
    energy(idx) = norm(v_vec)^2 / 2 - mu / norm(r_vec);
end

% Unpack OEs
a = elements(:, 1);
e = elements(:, 2);
i = elements(:, 3);
W = elements(:, 4);
w = elements(:, 5);
f = unwrap(elements(:, 6));

% Plot OES
figure;

subplot(3,2,1);
plot(t/3600, a/1e3);
xlabel('Time [h]'); ylabel('Semi-major axis [km]'); grid on;

subplot(3,2,2);
plot(t/3600, e);
xlabel('Time [h]'); ylabel('Eccentricity'); grid on;

subplot(3,2,3);
plot(t/3600, i);
xlabel('Time [h]'); ylabel('Inclination [deg]'); grid on;

subplot(3,2,4);
plot(t/3600, W);
xlabel('Time [h]'); ylabel('RAAN [deg]'); grid on;

subplot(3,2,5);
plot(t/3600, w);
xlabel('Time [h]'); ylabel('Argument of periapsis [deg]'); grid on;

subplot(3,2,6);
plot(t/3600, f);
xlabel('Time [h]'); ylabel('True anomaly [deg]'); grid on;

sgtitle('Orbital elements over time');


% Plot the angular momentum, eccentricty vector and specific mechanical energy over time
figure;

subplot(3,1,1);
plot(t/3600, h_mag);
xlabel('Time [h]'); ylabel('Angular momentum [m^2/s]'); grid on;

subplot(3,1,2);
plot(t/3600, e_mag);
xlabel('Time [h]'); ylabel('Eccentricity vector magnitude'); grid on;

subplot(3,1,3);
plot(t/3600, energy/1e6);
xlabel('Time [h]'); ylabel('Specific mechanical energy [MJ/kg]'); grid on;

sgtitle('Angular Momentum, Eccentricity Vector, Specific Mechanical Energy');


%% Q6 - Superimposing FODE & Keplerian OE propagators

% Initial state in OE
state_OE_0 = [a_0, e_0, i_0, W_0, w_0, f_0];

% Initial state in ECI
state_ECI_0 = utils.OE2ECI(state_OE_0, const, body);

% Run the perturbed FODE propagator
perturbated = true;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t, state, const, body, perturbated);
[t, state] = ode113(odefun, (tstart:tint:tend)', state_ECI_0, options);

% Extract position and velocity vectors
r = state(:, 1:3);
v = state(:, 4:6);

% Compute orbital elements from ECI (osculating elements)
numSteps = length(t);
elements_FODE = zeros(numSteps, 6);

for idx = 1:numSteps
    [a, e, i, W, w, f] = utils.ECI2OE(r(idx, :), v(idx, :), const, body);
    elements_FODE(idx, :) = [a, e, i, W, w, f];
end

% Unpack
a_FODE = elements_FODE(:, 1);
e_FODE = elements_FODE(:, 2);
i_FODE = elements_FODE(:, 3);
W_FODE = elements_FODE(:, 4);
w_FODE = elements_FODE(:, 5);
f_FODE = unwrap(elements_FODE(:, 6));


% Run the perturbed Keplerian OE propagator
perturbated = true;
odefun = @(t, state) propagators.KeplerianPropagator.getStatedot(t, state, const, body, perturbated);
[t, state_OE] = ode113(odefun, (tstart:tint:tend)', state_OE_0, options);

% Unpack elements
a_Kep = state_OE(:, 1);
e_Kep = state_OE(:, 2);
i_Kep = state_OE(:, 3);
W_Kep = state_OE(:, 4);
w_Kep = state_OE(:, 5);
f_Kep = state_OE(:, 6);


figure;

subplot(3,2,1);
plot(t/3600, a_FODE/1e3, 'b', t/3600, a_Kep/1e3, 'r--');
xlabel('Time [h]'); ylabel('Semi-major axis [km]'); grid on;
legend('FODE Propagator', 'Keplerian OE Propagator');

subplot(3,2,2);
plot(t/3600, e_FODE, 'b', t/3600, e_Kep, 'r--');
xlabel('Time [h]'); ylabel('Eccentricity'); grid on;

subplot(3,2,3);
plot(t/3600, rad2deg(i_FODE), 'b', t/3600, rad2deg(i_Kep), 'r--');
xlabel('Time [h]'); ylabel('Inclination [deg]'); grid on;

subplot(3,2,4);
plot(t/3600, rad2deg(W_FODE), 'b', t/3600, rad2deg(W_Kep), 'r--');
xlabel('Time [h]'); ylabel('RAAN [deg]'); grid on;

subplot(3,2,5);
plot(t/3600, rad2deg(w_FODE), 'b', t/3600, rad2deg(w_Kep), 'r--');
xlabel('Time [h]'); ylabel('Argument of Periapsis [deg]'); grid on;

subplot(3,2,6);
plot(t/3600, rad2deg(f_FODE), 'b', t/3600, rad2deg(f_Kep), 'r--');
xlabel('Time [h]'); ylabel('True Anomaly [deg]'); grid on;

sgtitle('Comparison: FODE vs. Keplerian OE propagators');


%% Q6 - Superimposing FODE & Mean J2 effect OE propagators

function statedot = getMeanStatedot(t, state, const, body)
    % Extract constants
    mu = const.(body).mu;
    J2 = const.(body).J2;
    Re = const.(body).R;

    % Extract mean orbital elements
    a = state(1);
    e = state(2);
    i = state(3);
    W = state(4); % RAAN
    w = state(5); % argument of perigee
    M = state(6); % mean anomaly

    % Derived quantities
    n = sqrt(mu / a^3); % mean motion
    p = a * (1 - e^2);

    % Mean rates (secular drift)
    dadt = 0;
    dedt = 0;
    didt = 0;
    dWdt = - (3/2) * J2 * (Re^2 / p^2) * n * cos(i);
    dwdt = (3/4) * J2 * (Re^2 / p^2) * n * (5 * cos(i)^2 - 1);
    dMdt = n + (3/4) * J2 * (Re^2 / p^2) * n * sqrt(1 - e^2) * (3 * cos(i)^2 - 1);

    % Pack derivatives
    statedot = [dadt; dedt; didt; dWdt; dwdt; dMdt];
end

% Initial state in OE,ECI,mean
state_OE_0 = [a_0, e_0, i_0, W_0, w_0, f_0];
state_ECI_0 = utils.OE2ECI(state_OE_0, const, body);
E_0 = 2 * atan( sqrt((1 - e_0) / (1 + e_0)) * tan(f_0 / 2) );
E_0 = mod(E_0, 2*pi);
M_0 = E_0 - e_0 * sin(E_0);
state_mean_0 = [a_0, e_0, i_0, W_0, w_0, M_0]; 



% Run the perturbed mean OE propagator
odefun = @(t, state) getMeanStatedot(t, state, const, body);
[t, state_mean] = ode113(odefun, (tstart:tint:tend)', state_mean_0, options);
a_mean = state_mean(:, 1);
e_mean = state_mean(:, 2);
i_mean = state_mean(:, 3);
W_mean = state_mean(:, 4);
w_mean = state_mean(:, 5);
M_mean = state_mean(:, 6);

% Run the perturbed FODE propagator
perturbated = true;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t, state, const, body, perturbated);
[t, state] = ode113(odefun, (tstart:tint:tend)', state_ECI_0, options);
r = state(:, 1:3);
v = state(:, 4:6);
numSteps = length(t);
elements_FODE = zeros(numSteps, 6);
for idx = 1:numSteps
    [a, e, i, W, w, f] = utils.ECI2OE(r(idx, :), v(idx, :), const, body);
    elements_FODE(idx, :) = [a, e, i, W, w, f];
end
a_FODE = elements_FODE(:, 1);
e_FODE = elements_FODE(:, 2);
i_FODE = elements_FODE(:, 3);
W_FODE = elements_FODE(:, 4);
w_FODE = elements_FODE(:, 5);
f_FODE = unwrap(elements_FODE(:, 6));

% Get M_FODE
E_FODE = 2 * atan( sqrt((1 - e_FODE) ./ (1 + e_FODE)) .* tan(f_FODE / 2) );
E_FODE = mod(E_FODE, 2*pi);
M_FODE = E_FODE - e_FODE .* sin(E_FODE);
M_FODE_unwrapped = unwrap(M_FODE);


figure;

subplot(3,2,1);
plot(t/3600, a_FODE/1e3, 'b', t/3600, a_mean/1e3, 'r--');
xlabel('Time [h]'); ylabel('Semi-major axis [km]'); grid on;
legend('FODE (Osculating)', 'Mean Theory');

subplot(3,2,2);
plot(t/3600, e_FODE, 'b', t/3600, e_mean, 'r--');
xlabel('Time [h]'); ylabel('Eccentricity'); grid on;

subplot(3,2,3);
plot(t/3600, rad2deg(i_FODE), 'b', t/3600, rad2deg(i_mean), 'r--');
xlabel('Time [h]'); ylabel('Inclination [deg]'); grid on;

subplot(3,2,4);
plot(t/3600, rad2deg(W_FODE), 'b', t/3600, rad2deg(W_mean), 'r--');
xlabel('Time [h]'); ylabel('RAAN [deg]'); grid on;

subplot(3,2,5);
plot(t/3600, rad2deg(w_FODE), 'b', t/3600, rad2deg(w_mean), 'r--');
xlabel('Time [h]'); ylabel('Argument of periapsis [deg]'); grid on;

subplot(3,2,6);
plot(t/3600, rad2deg(M_FODE_unwrapped), 'b', t/3600, rad2deg(M_mean), 'r--');
xlabel('Time [h]'); ylabel('Mean anomaly [deg]'); grid on;

sgtitle('Comparison FODE vs Mean OE (J2) propagators');