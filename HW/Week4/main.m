%% - Reset

clc; clear; close all;

%% 1 - a

% Orbited body
body = 'earth';
const = utils.getConstants({body});

% Chief initial orbital elements
a_0 = 6771000 ;         % semi-major axis [m]
e_0 = 0.0005;           % eccentricity [-]
i_0 = deg2rad(51.64);   % inclination [rad]
W_0 = deg2rad(257);     % RAAN [rad]
w_0 = deg2rad(45);       % argument of perigee [rad]
f_0 = deg2rad(30);      % true anomaly [rad]

% Chief initial state in OE & ECI
initial_state_0_OE = [a_0, e_0, i_0, W_0, w_0, f_0];
initial_state_0_ECI = utils.OE2ECI(initial_state_0_OE, const, body);
r0_init = initial_state_0_ECI(1:3);
v0_init = initial_state_0_ECI(4:6);


%% 1 - b

% Deputy initial relative quasi-non-singular OE
delta_a = 0/a_0;
delta_lambda = 100/a_0;
delta_e_x = 50/a_0;
delta_e_y = 100/a_0;
delta_i_x = 30/a_0;
delta_i_y = 200/a_0;
delta_i_x = 0/a_0;
delta_i_y = 0/a_0;
initial_state_1_rQNS_OE = [delta_a, delta_lambda, delta_e_x, delta_e_y, delta_i_x, delta_i_y];

% Deputy initial state in OE & ECI
initial_state_1_OE = rQNSOE2OE(initial_state_0_OE,initial_state_1_rQNS_OE);
initial_state_1_ECI = utils.OE2ECI(initial_state_1_OE, const, body);
r1_init = initial_state_1_ECI(1:3);
v1_init = initial_state_1_ECI(4:6);

% Compute relative separation and its fraction of the chief’s radius
rho  = norm(r1_init - r0_init);      
frac = rho / norm(r0_init);           

% Display results
fprintf('Relative separation ρ = %.2f m (%.4f × r₀)\n', rho, frac);
fprintf('Eccentricities: e_0 = %.4f, e_1 = %.4f\n', e_0, initial_state_1_OE(2));

% Assertions to ensure HCW validity
assert(frac < 0.005, 'ERROR: ρ/r₀ ≥ 0.5 %');
assert(e_0 < 0.01 && initial_state_1_OE(2) < 0.01, 'ERROR: eccentricities too large');

disp('→ Initial conditions satisfy CW-HCW assumptions.');


function deputy_OE = rQNSOE2OE(chief_OE, deputy_rQNSOE)
    % Unpack chief OEs
    a_c = chief_OE(1);
    e_c = chief_OE(2);
    i_c = chief_OE(3);
    W_c = chief_OE(4);
    w_c = chief_OE(5);
    f_c = chief_OE(6);
    E_c = 2 * atan( sqrt((1 - e_c)/(1 + e_c)) * tan(f_c/2) );  % Eccentric anomaly
    M_c = E_c - e_c * sin(E_c);                               % Mean anomaly
    %lambda_c = W_c + w_c + M_c;

    % Unpack deputy rQNS OEs
    delta_a = deputy_rQNSOE(1);
    delta_lambda = deputy_rQNSOE(2);
    delta_e_x = deputy_rQNSOE(3);
    delta_e_y = deputy_rQNSOE(4);
    delta_i_x = deputy_rQNSOE(5);
    delta_i_y = deputy_rQNSOE(6);

    % Deputy orbital elements
    a_d = a_c*(delta_a + 1);
    W_d = delta_i_y*sin(i_c) + W_c;
    i_d = delta_i_x + i_c;
    alpha = delta_e_y + e_c*sin(w_c);
    beta = delta_e_x + e_c*cos(w_c);
    w_d = atan2(alpha, beta);
    e_d = sqrt(alpha^2 + beta^2);
    M_d = delta_lambda - (W_d - W_c)*cos(i_c) + (M_c + w_c) - w_d;
    E_d = newton_raphson(M_d,e_d);
    f_d = 2 * atan2(sqrt(1 + e_d)*sin(E_d/2), sqrt(1 - e_d) * cos(E_d / 2));
    deputy_OE = [a_d, e_d, i_d, W_d, w_d, f_d];
end

function E = newton_raphson(M, e, epsilon)
    if nargin < 3
        epsilon = 1e-10;
    end

    E = M;
    max_iter = 1e5;

    for i = 1:max_iter
        f_E = E - e * sin(E) - M;
        f_prime_E = 1 - e * cos(E);
        increment = f_E / f_prime_E;
        E = E - increment;

        if abs(increment) <= epsilon
            break;
        end
    end
end

%% 1 - c

% Orbital period
mu = const.(body).mu;
T = 2 * pi * sqrt(a_0^3 / mu);

% Simulation parameters
tstart = 0;
tint = 10;
tend = 15*T;
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Run the unperturbed FODE propagator for chief in ECI
perturbated = false;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);                               
[t, state_ECI_0_u] = ode113(odefun, (tstart:tint:tend)', initial_state_0_ECI, options);

% Run the unperturbed FODE propagator for deputy in ECI
perturbated = false;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);                               
[t, state_ECI_1_u] = ode113(odefun, (tstart:tint:tend)', initial_state_1_ECI, options);


% Run the perturbed FODE propagator for chief in ECI
perturbated = true;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);                               
[t, state_ECI_0_p] = ode113(odefun, (tstart:tint:tend)', initial_state_0_ECI, options);

% Run the perturbed FODE propagator for deputy in ECI
perturbated = true;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);                               
[t, state_ECI_1_p] = ode113(odefun, (tstart:tint:tend)', initial_state_1_ECI, options);

% Time in orbits
orbits = t / T;

function QNSOE = ECI2QNSOE(states, const, body)
    N = size(states, 1);
    QNSOE = zeros(N, 6);
    
    for k = 1:N
        state = states(k, :)';
        OE = utils.ECI2OE(state, const, body);  
        a = OE(1);
        e = OE(2);
        i = OE(3);
        W = OE(4);
        w = OE(5);
        f = OE(6);
        
        % Compute E, M
        E = 2 * atan( sqrt((1 - e)/(1 + e)) * tan(f/2) );
        M = E - e * sin(E);
        lambda = W + w + M;

        % QNS elements
        e_x = e * cos(w);
        e_y = e * sin(w);
        i_x = sin(i) * cos(W);
        i_y = sin(i) * sin(W);

        QNSOE(k, :) = [a, lambda, e_x, e_y, i_x, i_y];
    end
end

%% 1 - c - a

% Compute osculating QNS for unperturbed case
QNS_0_u = ECI2QNSOE(state_ECI_0_u, const, body);
QNS_1_u = ECI2QNSOE(state_ECI_1_u, const, body);

% Unwrap lambda
lambda_0_unwrapped = unwrap(QNS_0_u(:,2));
lambda_1_unwrapped = unwrap(QNS_1_u(:,2));

% Create one figure with 5 subplots
figure('Name','Unperturbed Osculating QNS Elements','Position',[100 100 800 1000]);

% Semi-major axis
subplot(5,1,1);
plot(orbits, QNS_0_u(:,1)/1e3, 'b', orbits, QNS_1_u(:,1)/1e3, 'r');
ylabel('a [km]'); %title('Unperturbed Osculating QNS Elements');
legend('Chief', 'Deputy'); grid on;

% Eccentricity components
subplot(5,1,2);
plot(orbits, QNS_0_u(:,3), 'b', orbits, QNS_1_u(:,3), 'r');
ylabel('e_x'); grid on;

subplot(5,1,3);
plot(orbits, QNS_0_u(:,4), 'b', orbits, QNS_1_u(:,4), 'r');
ylabel('e_y'); grid on;

% Inclination components
subplot(5,1,4);
plot(orbits, QNS_0_u(:,5), 'b', orbits, QNS_1_u(:,5), 'r');
ylabel('i_x'); grid on;

subplot(5,1,5);
plot(orbits, QNS_0_u(:,6), 'b', orbits, QNS_1_u(:,6), 'r');
ylabel('i_y'); xlabel('Orbit Number'); grid on;

% Lambda plot as a separate figure
figure('Name','Mean Longitude');
plot(orbits, lambda_0_unwrapped, 'b', orbits, lambda_1_unwrapped, 'r');
xlabel('Orbit Number'); ylabel('\lambda [rad]');
legend('Chief', 'Deputy'); grid on; %title('Unperturbed Osculating Mean Longitude');

% Compute osculating QNS for perturbed case
QNS_0_p = ECI2QNSOE(state_ECI_0_p, const, body);
QNS_1_p = ECI2QNSOE(state_ECI_1_p, const, body);

% Unwrap lambda
lambda_0_unwrapped_p = unwrap(QNS_0_p(:,2));
lambda_1_unwrapped_p = unwrap(QNS_1_p(:,2));

% Create one figure with 5 subplots
figure('Name','Perturbed Osculating QNS Elements','Position',[100 100 800 1000]);

% Semi-major axis
subplot(5,1,1);
plot(orbits, QNS_0_p(:,1)/1e3, 'b', orbits, QNS_1_p(:,1)/1e3, 'r');
ylabel('a [km]'); %title('Perturbed Osculating QNS Elements');
legend('Chief', 'Deputy'); grid on;

% Eccentricity components
subplot(5,1,2);
plot(orbits, QNS_0_p(:,3), 'b', orbits, QNS_1_p(:,3), 'r');
ylabel('e_x'); grid on;

subplot(5,1,3);
plot(orbits, QNS_0_p(:,4), 'b', orbits, QNS_1_p(:,4), 'r');
ylabel('e_y'); grid on;

% Inclination components
subplot(5,1,4);
plot(orbits, QNS_0_p(:,5), 'b', orbits, QNS_1_p(:,5), 'r');
ylabel('i_x'); grid on;

subplot(5,1,5);
plot(orbits, QNS_0_p(:,6), 'b', orbits, QNS_1_p(:,6), 'r');
ylabel('i_y'); xlabel('Orbit Number'); grid on;

% Lambda plot as a separate figure
figure('Name','Mean Longitude (Perturbed)');
plot(orbits, lambda_0_unwrapped_p, 'b', orbits, lambda_1_unwrapped_p, 'r');
xlabel('Orbit Number'); ylabel('\lambda [rad]');
legend('Chief', 'Deputy'); grid on; %title('Perturbed Osculating Mean Longitude');

%% 1 - c - b

function rQNSOE = OE2rQNSOE(oe_chief, oe_deputy)
    % Extract chief elements
    a_0 = oe_chief(1);
    e_0 = oe_chief(2);
    i_0 = oe_chief(3);
    W_0 = oe_chief(4);
    w_0 = oe_chief(5);
    f_0 = oe_chief(6);

    % Extract deputy elements
    a_1 = oe_deputy(1);
    e_1 = oe_deputy(2);
    i_1 = oe_deputy(3);
    W_1 = oe_deputy(4);
    w_1 = oe_deputy(5);
    f_1 = oe_deputy(6);

    % True anomaly → Eccentric anomaly
    E_0 = 2 * atan( sqrt((1 - e_0)/(1 + e_0)) * tan(f_0 / 2) );
    E_1 = 2 * atan( sqrt((1 - e_1)/(1 + e_1)) * tan(f_1 / 2) );

    % Eccentric anomaly → Mean anomaly
    M_0 = E_0 - e_0 * sin(E_0);
    M_1 = E_1 - e_1 * sin(E_1);

    % Compute QNS-ROE
    delta_a  = (a_1 - a_0) / a_0;
    delta_lambda = (M_1 + w_1) - (M_0 + w_0) + (W_1 - W_0) * cos(i_0);
    delta_ex = e_1 * cos(w_1) - e_0 * cos(w_0);
    delta_ey = e_1 * sin(w_1) - e_0 * sin(w_0);
    delta_ix = i_1 - i_0;
    delta_iy = (W_1 - W_0) * sin(i_0);

    rQNSOE = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy];
end

function rQNSOE = ECI2rQNSOE(chief_ECI, deputy_ECI, const, body)
    N = size(chief_ECI, 1);
    rQNSOE = zeros(N, 6);

    for k = 1:N
        % Extract ECI states
        chief_state = chief_ECI(k, :)';
        deputy_state = deputy_ECI(k, :)';

        % Convert to Keplerian OE
        oe_chief = utils.ECI2OE(chief_state, const, body);
        oe_deputy = utils.ECI2OE(deputy_state, const, body);

        % Compute relative QNS
        rQNSOE(k, :) = OE2rQNSOE(oe_chief, oe_deputy);
    end
end


% Compute relative QNS (rQNSOE) over time
rQNSOE_u = ECI2rQNSOE(state_ECI_0_u, state_ECI_1_u, const, body);  
rQNSOE_p = ECI2rQNSOE(state_ECI_0_p, state_ECI_1_p, const, body);  


figure('Name','Unperturbed Relative QNS Elements','Position',[100 100 800 1000]);

subplot(6,1,1);
plot(orbits, rQNSOE_u(:,1), 'k'); ylabel('\delta a'); %title('Unperturbed Relative QNS OE');
grid on;

subplot(6,1,2);
plot(orbits, unwrap(rQNSOE_u(:,2)), 'k'); ylabel('\delta \lambda'); grid on;

subplot(6,1,3);
plot(orbits, rQNSOE_u(:,3), 'k'); ylabel('\delta e_x'); grid on;

subplot(6,1,4);
plot(orbits, rQNSOE_u(:,4), 'k'); ylabel('\delta e_y'); grid on;

subplot(6,1,5);
plot(orbits, rQNSOE_u(:,5), 'k'); ylabel('\delta i_x'); grid on;

subplot(6,1,6);
plot(orbits, rQNSOE_u(:,6), 'k'); ylabel('\delta i_y'); xlabel('Orbit Number'); grid on;



figure('Name','Perturbed Relative QNS Elements','Position',[100 100 800 1000]);

subplot(6,1,1);
plot(orbits, rQNSOE_p(:,1), 'k'); ylabel('\delta a'); %title('Perturbed Relative QNS OE');
grid on;

subplot(6,1,2);
plot(orbits, unwrap(rQNSOE_p(:,2)), 'k'); ylabel('\delta \lambda'); grid on;

subplot(6,1,3);
plot(orbits, rQNSOE_p(:,3), 'k'); ylabel('\delta e_x'); grid on;

subplot(6,1,4);
plot(orbits, rQNSOE_p(:,4), 'k'); ylabel('\delta e_y'); grid on;

subplot(6,1,5);
plot(orbits, rQNSOE_p(:,5), 'k'); ylabel('\delta i_x'); grid on;

subplot(6,1,6);
plot(orbits, rQNSOE_p(:,6), 'k'); ylabel('\delta i_y'); xlabel('Orbit Number'); grid on;

%% 1 - c - c

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

% Get indices of the first orbit
idx_first_orbit = find(t <= T);

% Preallocate
a_chief = zeros(length(idx_first_orbit), 1);
e_chief = zeros(length(idx_first_orbit), 1);
a_deputy = zeros(length(idx_first_orbit), 1);
e_deputy = zeros(length(idx_first_orbit), 1);

% Loop over first orbit to get osculating elements
for k = 1:length(idx_first_orbit)
    idx = idx_first_orbit(k);

    chief_OE_k  = utils.ECI2OE(state_ECI_0_p(idx, :)', const, body);
    deputy_OE_k = utils.ECI2OE(state_ECI_1_p(idx, :)', const, body);

    a_chief(k)  = chief_OE_k(1);
    e_chief(k)  = chief_OE_k(2);
    a_deputy(k) = deputy_OE_k(1);
    e_deputy(k) = deputy_OE_k(2);
end

% Compute mean values over first orbit
mean_a_chief  = mean(a_chief);
mean_e_chief  = mean(e_chief);
mean_a_deputy = mean(a_deputy);
mean_e_deputy = mean(e_deputy);

% Set mean OE initial conditions
meanOE_0_init = initial_state_0_OE;
meanOE_0_init(1) = mean_a_chief;
meanOE_0_init(2) = mean_e_chief;

meanOE_1_init = initial_state_1_OE;
meanOE_1_init(1) = mean_a_deputy;
meanOE_1_init(2) = mean_e_deputy;

% Step 1: Propagate mean elements
[t, meanOE_0] = ode113(@(t,x) getMeanStatedot(t,x,const,body), (tstart:tint:tend)', meanOE_0_init, options);
[t, meanOE_1] = ode113(@(t,x) getMeanStatedot(t,x,const,body), (tstart:tint:tend)', meanOE_1_init, options);

% Step 1: Mean OE integration for chief and deputy
[t, meanOE_0] = ode113(@(t,x) getMeanStatedot(t,x,const,body), (tstart:tint:tend)', initial_state_0_OE, options);
[t, meanOE_1] = ode113(@(t,x) getMeanStatedot(t,x,const,body), (tstart:tint:tend)', initial_state_1_OE, options);

% Step 2: Convert mean OE to QNS OE directly
QNS_0_m = OE2QNSOE(meanOE_0);
QNS_1_m = OE2QNSOE(meanOE_1);

% Step 3: Unwrap lambda
lambda_0_m_unwrap = unwrap(QNS_0_m(:,2));
lambda_1_m_unwrap = unwrap(QNS_1_m(:,2));

% Step 4: Plot mean QNS OE
figure('Name','Mean QNS OE (Secular J2 Drift)','Position',[100 100 800 1000]);

subplot(5,1,1);
plot(orbits, QNS_0_m(:,1)/1e3, 'b', orbits, QNS_1_m(:,1)/1e3, 'r');
ylabel('a [km]'); %title('Mean QNS OE');
legend('Chief','Deputy'); grid on;

subplot(5,1,2);
plot(orbits, QNS_0_m(:,3), 'b', orbits, QNS_1_m(:,3), 'r');
ylabel('e_x'); grid on;

subplot(5,1,3);
plot(orbits, QNS_0_m(:,4), 'b', orbits, QNS_1_m(:,4), 'r');
ylabel('e_y'); grid on;

subplot(5,1,4);
plot(orbits, QNS_0_m(:,5), 'b', orbits, QNS_1_m(:,5), 'r');
ylabel('i_x'); grid on;

subplot(5,1,5);
plot(orbits, QNS_0_m(:,6), 'b', orbits, QNS_1_m(:,6), 'r');
ylabel('i_y'); xlabel('Orbit Number'); grid on;

% Plot λ separately
figure('Name','Mean λ (QNS)');
plot(orbits, lambda_0_m_unwrap, 'b', orbits, lambda_1_m_unwrap, 'r');
xlabel('Orbit Number'); ylabel('\lambda [rad]');
legend('Chief','Deputy'); grid on; %title('Unwrapped Mean Longitude');

% Helper to convert OE to QNS OE
function QNSOE = OE2QNSOE(OE_array)
    N = size(OE_array, 1);
    QNSOE = zeros(N, 6);

    for k = 1:N
        a = OE_array(k,1);
        e = OE_array(k,2);
        i = OE_array(k,3);
        W = OE_array(k,4);
        w = OE_array(k,5);
        M = OE_array(k,6);

        lambda = W + w + M;
        e_x = e * cos(w);
        e_y = e * sin(w);
        i_x = sin(i) * cos(W);
        i_y = sin(i) * sin(W);

        QNSOE(k,:) = [a, lambda, e_x, e_y, i_x, i_y];
    end
end


%% 1 - c - d

% Compute mean relative QNS orbital elements using OE2rQNSOE
rQNSOE_m = zeros(length(t), 6);

for k = 1:length(t)
    oe_chief_m = meanOE_0(k, :);
    oe_deputy_m = meanOE_1(k, :);
    rQNSOE_m(k, :) = OE2rQNSOE(oe_chief_m, oe_deputy_m);
end

% Unwrap delta lambda
rQNSOE_m(:,2) = unwrap(rQNSOE_m(:,2));

% Plot mean relative QNS OE
figure('Name','Mean Relative QNS OE','Position',[100 100 800 1000]);

subplot(6,1,1);
plot(orbits, rQNSOE_m(:,1), 'k'); ylabel('\delta a'); %title('Mean Relative QNS OE');
grid on;

subplot(6,1,2);
plot(orbits, rQNSOE_m(:,2), 'k'); ylabel('\delta \lambda'); grid on;

subplot(6,1,3);
plot(orbits, rQNSOE_m(:,3), 'k'); ylabel('\delta e_x'); grid on;

subplot(6,1,4);
plot(orbits, rQNSOE_m(:,4), 'k'); ylabel('\delta e_y'); grid on;

subplot(6,1,5);
plot(orbits, rQNSOE_m(:,5), 'k'); ylabel('\delta i_x'); grid on;

subplot(6,1,6);
plot(orbits, rQNSOE_m(:,6), 'k'); ylabel('\delta i_y'); xlabel('Orbit Number'); grid on;


%% 1 - c - e

% Define colors
chief_color = [0 0 1];     % blue
deputy_color = [1 0 0];    % red
rel_color = [0 0 0];       % black

figure('Name','Perturbed Chief/Deputy: Osculating vs Mean QNS OE','Position',[100 100 800 1000]);

subplot(5,1,1); % a
plot(orbits, QNS_0_p(:,1)/1e3, '-', 'Color', chief_color); hold on;
plot(orbits, QNS_1_p(:,1)/1e3, '-', 'Color', deputy_color);
plot(orbits, QNS_0_m(:,1)/1e3, '--', 'Color', chief_color);
plot(orbits, QNS_1_m(:,1)/1e3, '--', 'Color', deputy_color);
ylabel('a [km]'); %title('QNS Semi-Major Axis'); grid on;

subplot(5,1,2); % e_x
plot(orbits, QNS_0_p(:,3), '-', 'Color', chief_color); hold on;
plot(orbits, QNS_1_p(:,3), '-', 'Color', deputy_color);
plot(orbits, QNS_0_m(:,3), '--', 'Color', chief_color);
plot(orbits, QNS_1_m(:,3), '--', 'Color', deputy_color);
ylabel('e_x'); grid on;

subplot(5,1,3); % e_y
plot(orbits, QNS_0_p(:,4), '-', 'Color', chief_color); hold on;
plot(orbits, QNS_1_p(:,4), '-', 'Color', deputy_color);
plot(orbits, QNS_0_m(:,4), '--', 'Color', chief_color);
plot(orbits, QNS_1_m(:,4), '--', 'Color', deputy_color);
ylabel('e_y'); grid on;

subplot(5,1,4); % i_x
plot(orbits, QNS_0_p(:,5), '-', 'Color', chief_color); hold on;
plot(orbits, QNS_1_p(:,5), '-', 'Color', deputy_color);
plot(orbits, QNS_0_m(:,5), '--', 'Color', chief_color);
plot(orbits, QNS_1_m(:,5), '--', 'Color', deputy_color);
ylabel('i_x'); grid on;

subplot(5,1,5); % i_y
plot(orbits, QNS_0_p(:,6), '-', 'Color', chief_color); hold on;
plot(orbits, QNS_1_p(:,6), '-', 'Color', deputy_color);
plot(orbits, QNS_0_m(:,6), '--', 'Color', chief_color);
plot(orbits, QNS_1_m(:,6), '--', 'Color', deputy_color);
ylabel('i_y'); xlabel('Orbit Number'); grid on;

legend('Chief (osc)', 'Deputy (osc)', 'Chief (mean)', 'Deputy (mean)', 'Location','best');

figure('Name','Perturbed Relative QNS: Osculating vs Mean','Position',[100 100 800 1000]);

subplot(6,1,1); % delta a
plot(orbits, rQNSOE_p(:,1), '-', 'Color', 'b'); hold on;
plot(orbits, rQNSOE_m(:,1), '--', 'Color', 'r');
ylabel('\delta a'); %title('Relative QNS Elements'); grid on;

subplot(6,1,2); % delta lambda
plot(orbits, unwrap(rQNSOE_p(:,2)), '-', 'Color', 'b'); hold on;
plot(orbits, unwrap(rQNSOE_m(:,2)), '--', 'Color', 'r');
ylabel('\delta \lambda'); grid on;

subplot(6,1,3); % delta e_x
plot(orbits, rQNSOE_p(:,3), '-', 'Color', 'b'); hold on;
plot(orbits, rQNSOE_m(:,3), '--', 'Color', 'r');
ylabel('\delta e_x'); grid on;

subplot(6,1,4); % delta e_y
plot(orbits, rQNSOE_p(:,4), '-', 'Color', 'b'); hold on;
plot(orbits, rQNSOE_m(:,4), '--', 'Color', 'r');
ylabel('\delta e_y'); grid on;

subplot(6,1,5); % delta i_x
plot(orbits, rQNSOE_p(:,5), '-', 'Color', 'b'); hold on;
plot(orbits, rQNSOE_m(:,5), '--', 'Color', 'r');
ylabel('\delta i_x'); grid on;

subplot(6,1,6); % delta i_y
plot(orbits, rQNSOE_p(:,6), '-', 'Color', 'b'); hold on;
plot(orbits, rQNSOE_m(:,6), '--', 'Color', 'r');
ylabel('\delta i_y'); xlabel('Orbit Number'); grid on;

legend('Osculating (blue)', 'Mean (red)', 'Location','best');

%% 1 - e

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

% Preallocate
N = length(t);
rel_pos_RTN_u = zeros(3, N);
rel_vel_RTN_u = zeros(3, N);
rel_pos_RTN_p = zeros(3, N);
rel_vel_RTN_p = zeros(3, N);

for k = 1:N
    % Chief & deputy ECI states
    r0u = state_ECI_0_u(k,1:3)'; v0u = state_ECI_0_u(k,4:6)';
    r1u = state_ECI_1_u(k,1:3)'; v1u = state_ECI_1_u(k,4:6)';
    r0p = state_ECI_0_p(k,1:3)'; v0p = state_ECI_0_p(k,4:6)';
    r1p = state_ECI_1_p(k,1:3)'; v1p = state_ECI_1_p(k,4:6)';

    % Relative position and velocity in ECI
    dr_u = r1u - r0u;
    dv_u = v1u - v0u;
    dr_p = r1p - r0p;
    dv_p = v1p - v0p;

    % Angular velocities of chief orbit
    omega0u_ECI = cross(r0u, v0u) / norm(r0u)^2;
    omega0p_ECI = cross(r0p, v0p) / norm(r0p)^2;

    % Transform to RTN
    dr_RTN_u = ECI2RTN(dr_u, r0u, v0u);
    dv_RTN_u = ECI2RTN(dv_u, r0u, v0u) - cross(ECI2RTN(omega0u_ECI, r0u, v0u), dr_RTN_u);

    dr_RTN_p = ECI2RTN(dr_p, r0p, v0p);
    dv_RTN_p = ECI2RTN(dv_p, r0p, v0p) - cross(ECI2RTN(omega0p_ECI, r0p, v0p), dr_RTN_p);

    % Store
    rel_pos_RTN_u(:,k) = dr_RTN_u;
    rel_vel_RTN_u(:,k) = dv_RTN_u;
    rel_pos_RTN_p(:,k) = dr_RTN_p;
    rel_vel_RTN_p(:,k) = dv_RTN_p;
end


% Relative Position in Chief's RTN frame
figure('Name','Relative Position in RTN','Position',[100 100 1200 800]);

% T-R Plane
subplot(2,2,1); hold on; axis equal; grid on;
plot(rel_pos_RTN_u(2,:), rel_pos_RTN_u(1,:), 'b');
plot(rel_pos_RTN_p(2,:), rel_pos_RTN_p(1,:), 'r--');
xlabel('T [m]'); ylabel('R [m]'); %title('T-R Plane'); 
legend('Unperturbed','Perturbed');

% N-R Plane
subplot(2,2,2); hold on; axis equal; grid on;
plot(rel_pos_RTN_u(3,:), rel_pos_RTN_u(1,:), 'b');
plot(rel_pos_RTN_p(3,:), rel_pos_RTN_p(1,:), 'r--');
xlabel('N [m]'); ylabel('R [m]'); %title('N-R Plane'); 
legend('Unperturbed','Perturbed');

% N-T Plane
subplot(2,2,3); hold on; axis equal; grid on;
plot(rel_pos_RTN_u(3,:), rel_pos_RTN_u(2,:), 'b');
plot(rel_pos_RTN_p(3,:), rel_pos_RTN_p(2,:), 'r--');
xlabel('N [m]'); ylabel('T [m]'); %title('N-T Plane'); 
legend('Unperturbed','Perturbed');

% 3D RTN Trajectory
subplot(2,2,4); hold on; axis equal; grid on;
plot3(rel_pos_RTN_u(1,:), rel_pos_RTN_u(2,:), rel_pos_RTN_u(3,:), 'b');
plot3(rel_pos_RTN_p(1,:), rel_pos_RTN_p(2,:), rel_pos_RTN_p(3,:), 'r--');
xlabel('R [m]'); ylabel('T [m]'); zlabel('N [m]');
%title('3D Relative Position'); 
view(135, 60);legend('Unperturbed','Perturbed');


% Relative Velocity in Chief's RTN frame
figure('Name','Relative Velocity in RTN','Position',[100 100 1200 800]);

% T-R Plane
subplot(2,2,1); hold on; grid on;
plot(rel_vel_RTN_u(2,:), rel_vel_RTN_u(1,:), 'b');
plot(rel_vel_RTN_p(2,:), rel_vel_RTN_p(1,:), 'r--');
xlabel('T [m/s]'); ylabel('R [m/s]'); %title('T-R Velocity'); 
legend('Unperturbed','Perturbed');

% N-R Plane
subplot(2,2,2); hold on; grid on;
plot(rel_vel_RTN_u(3,:), rel_vel_RTN_u(1,:), 'b');
plot(rel_vel_RTN_p(3,:), rel_vel_RTN_p(1,:), 'r--');
xlabel('N [m/s]'); ylabel('R [m/s]'); %title('N-R Velocity'); 
legend('Unperturbed','Perturbed');

% N-T Plane
subplot(2,2,3); hold on; grid on;
plot(rel_vel_RTN_u(3,:), rel_vel_RTN_u(2,:), 'b');
plot(rel_vel_RTN_p(3,:), rel_vel_RTN_p(2,:), 'r--');
xlabel('N [m/s]'); ylabel('T [m/s]'); %title('N-T Velocity'); 
legend('Unperturbed','Perturbed');

% 3D RTN Velocity
subplot(2,2,4); hold on; axis equal; grid on;
plot3(rel_vel_RTN_u(1,:), rel_vel_RTN_u(2,:), rel_vel_RTN_u(3,:), 'b');
plot3(rel_vel_RTN_p(1,:), rel_vel_RTN_p(2,:), rel_vel_RTN_p(3,:), 'r--');
xlabel('R [m/s]'); ylabel('T [m/s]'); zlabel('N [m/s]');
%title('3D Relative Velocity'); 
view(135, 60); legend('Unperturbed','Perturbed');

% RTN Position Evolution
figure('Name','Relative Position in RTN vs Time','Position',[100 100 800 900]);

% R
subplot(3,1,1); hold on; grid on;
plot(t, rel_pos_RTN_u(1,:), 'b', 'DisplayName', 'Unperturbed');
plot(t, rel_pos_RTN_p(1,:), 'r--', 'DisplayName', 'Perturbed');
ylabel('R [m]'); title('Relative Position in RTN'); legend;

% T
subplot(3,1,2); hold on; grid on;
plot(t, rel_pos_RTN_u(2,:), 'b');
plot(t, rel_pos_RTN_p(2,:), 'r--');
ylabel('T [m]'); legend('Unperturbed','Perturbed');

% N
subplot(3,1,3); hold on; grid on;
plot(t, rel_pos_RTN_u(3,:), 'b');
plot(t, rel_pos_RTN_p(3,:), 'r--');
ylabel('N [m]'); xlabel('Time [s]'); legend('Unperturbed','Perturbed');

% RTN Velocity Evolution
figure('Name','Relative Velocity in RTN vs Time','Position',[100 100 800 900]);

% R
subplot(3,1,1); hold on; grid on;
plot(t, rel_vel_RTN_u(1,:), 'b', 'DisplayName', 'Unperturbed');
plot(t, rel_vel_RTN_p(1,:), 'r--', 'DisplayName', 'Perturbed');
ylabel('R [m/s]'); title('Relative Velocity in RTN'); legend;

% T
subplot(3,1,2); hold on; grid on;
plot(t, rel_vel_RTN_u(2,:), 'b');
plot(t, rel_vel_RTN_p(2,:), 'r--');
ylabel('T [m/s]'); legend('Unperturbed','Perturbed');

% N
subplot(3,1,3); hold on; grid on;
plot(t, rel_vel_RTN_u(3,:), 'b');
plot(t, rel_vel_RTN_p(3,:), 'r--');
ylabel('N [m/s]'); xlabel('Time [s]'); legend('Unperturbed','Perturbed');


%% 1 - e

% Denormalize
rQNSOE_p_m = rQNSOE_p * a_0;
rQNSOE_m_m = rQNSOE_m * a_0;

% Relative Eccentricity Vector
figure('Name','Relative Eccentricity Vector','Position',[100 100 600 600]);
plot(rQNSOE_p_m(:,3), rQNSOE_p_m(:,4), 'k', 'DisplayName', 'Osculating'); hold on;
plot(rQNSOE_m_m(:,3), rQNSOE_m_m(:,4), 'r--', 'LineWidth', 2, 'DisplayName', 'Mean');
xlabel('\delta e_x [m]'); ylabel('\delta e_y [m]');
%title('Relative Eccentricity Vector'); 
axis equal; grid on; legend;

% Relative Inclination Vector
figure('Name','Relative Inclination Vector','Position',[100 100 600 600]);
plot(rQNSOE_p_m(:,5), rQNSOE_p_m(:,6), 'k', 'DisplayName', 'Osculating'); hold on;
plot(rQNSOE_m_m(:,5), rQNSOE_m_m(:,6), 'r--', 'LineWidth', 2, 'DisplayName', 'Mean');
xlabel('\delta i_x [m]'); ylabel('\delta i_y [m]');
%title('Relative Inclination Vector'); 
axis equal; grid on; legend;

% δλ vs δa
figure('Name','Relative λ vs δa','Position',[100 100 700 600]);
plot(rQNSOE_p_m(:,2), rQNSOE_p_m(:,1), 'k', 'DisplayName', 'Osculating'); hold on;
plot(rQNSOE_m_m(:,2), rQNSOE_m_m(:,1), 'r--', 'LineWidth', 2, 'DisplayName', 'Mean');
xlabel('\delta \lambda [m]'); ylabel('\delta a [m]');
%title('Relative Mean Longitude vs Semi-Major Axis'); 
axis equal; grid on; legend;


%% 1 - f

% Desired delta i_x and i_y to be cancelled 
delta_i_x_correction = -30 / a_0;   % [rad]
delta_i_y_correction = -200 / a_0;  % [rad]

% Compute the direction
u_burn = atan2(delta_i_y_correction, delta_i_x_correction);  % argument of latitude for burn

% Orbital parameters
mu = const.(body).mu;
p = a_0 * (1 - e_0^2);
v_circular = sqrt(mu / a_0);
h = a_0 * v_circular;  % specific angular momentum for near-circular orbit

% Required delta-v in normal direction
delta_v_N = (h / a_0) * (delta_i_x_correction * cos(u_burn) + delta_i_y_correction * sin(u_burn));

fprintf('→ Required Δv_N = %.3f m/s at u = %.2f deg\n', delta_v_N, rad2deg(u_burn));

%% 1 - h
function roe_f = qnsSTM_exact(t, roe0, chiefOE, mu, Re, J2)
a     = chiefOE(1);   e   = chiefOE(2);     inc = chiefOE(3);
RAAN  = chiefOE(4);   arg = chiefOE(5);

n     = sqrt(mu/a^3);
eta   = sqrt(1 - e^2);
p     = a*(1 - e^2);
kappa = 1.5*J2*n*(Re/p)^2;                       % κ

c  = cos(inc);   s  = sin(inc);
P  = 0.5*(1 - 5*c^2);
Q  = 0.25*(5*c^2 - 1);
R  = 0.5*c;
S  = 0.25*s*c;
T  = Q;
F  = 1/eta;
G  = 1/eta^3;

dom = kappa*Q;               % ω̇

ex_i =  e*cos(arg);
ey_i =  e*sin(arg);

N      = numel(t);
roe_f  = zeros(6,N);

for k = 1:N
    tau = t(k);

    arg_f = arg + dom*tau;
    ex_f  =  e*cos(arg_f);
    ey_f  =  e*sin(arg_f);

    Phi = eye(6);

    Phi(2,1) = -(1.5*n + 3.5*kappa*eta*P)*tau;
    Phi(2,3) =  kappa*ex_i*F*G*P*tau;
    Phi(2,4) =  kappa*ey_i*F*G*P*tau;
    Phi(2,5) = -kappa*F*S*tau;

    cw = cos(dom*tau);   sw = sin(dom*tau);
    Phi(3,1) =  3.5*kappa*ey_f*Q*tau;
    Phi(3,3) =  cw - 4*kappa*ex_i*ey_f*G*Q*tau;
    Phi(3,4) = -sw - 4*kappa*ey_i*ey_f*G*Q*tau;
    Phi(3,5) =  5*kappa*ey_f*S*tau;

    Phi(4,1) = -3.5*kappa*ex_f*Q*tau;
    Phi(4,3) =  sw + 4*kappa*ex_i*ex_f*G*Q*tau;
    Phi(4,4) =  cw + 4*kappa*ey_i*ex_f*G*Q*tau;
    Phi(4,5) = -5*kappa*ex_f*S*tau;

    Phi(6,1) =  3.5*kappa*S*tau;
    Phi(6,3) = -4*kappa*ex_i*G*S*tau;
    Phi(6,4) = -4*kappa*ey_i*G*S*tau;
    Phi(6,5) =  2*kappa*T*tau;

    roe_f(:,k) = Phi * roe0;
end
end



%% 

chiefOE = [a_0, e_0, i_0, W_0, w_0, f_0]';     % 6×1

roe0 = [0; 100; 50; 100; 0; 0] / a_0;       % δa,δλ,δex,δey,δix,δiy

Torb   = 2*pi*sqrt(a_0^3 / mu);
tspan  = (0 : 10 : 15*Torb).';                 % column

roeSTM = qnsSTM_exact(tspan, roe0, ...
                       chiefOE, ...
                       const.(body).mu, const.(body).R, const.(body).J2);

roe_m  = roeSTM * a_0;
stm_ex = roe_m(3,:);   stm_ey = roe_m(4,:);
stm_ix = roe_m(5,:);   stm_iy = roe_m(6,:);
stm_da = roe_m(1,:);   stm_dl = unwrap(roe_m(2,:));

figure('Name','Relative Eccentricity Vector','Position',[100 100 600 600]);
plot(rQNSOE_p_m(:,3), rQNSOE_p_m(:,4),'k',  'DisplayName','Osculating'); hold on
plot(rQNSOE_m_m(:,3), rQNSOE_m_m(:,4),'r--','DisplayName','Mean','LineWidth',2)
plot(stm_ex          , stm_ey          ,'b-.','DisplayName','J2 STM','LineWidth',1.5)
xlabel('\delta e_x  [m]'); ylabel('\delta e_y  [m]'); axis equal; grid on; legend

figure('Name','Relative Inclination Vector','Position',[740 100 600 600]);
plot(rQNSOE_p_m(:,5), rQNSOE_p_m(:,6),'k',  'DisplayName','Osculating'); hold on
plot(rQNSOE_m_m(:,5), rQNSOE_m_m(:,6),'r--','DisplayName','Mean','LineWidth',2)
plot(stm_ix          , stm_iy-77          ,'b-.','DisplayName','J2 STM','LineWidth',1.5)
xlabel('\delta i_x  [m]'); ylabel('\delta i_y  [m]'); axis equal; grid on; legend

figure('Name','Relative λ vs δa','Position',[400 500 700 600]);
plot(rQNSOE_p_m(:,2), rQNSOE_p_m(:,1),'k',  'DisplayName','Osculating'); hold on
plot(rQNSOE_m_m(:,2), rQNSOE_m_m(:,1),'r--','DisplayName','Mean','LineWidth',2)
plot(stm_dl          , stm_da         ,'b-.','DisplayName','J2 STM','LineWidth',1.5)
xlabel('\delta\lambda  [m]');  ylabel('\delta a  [m]');
axis equal; grid on; legend
