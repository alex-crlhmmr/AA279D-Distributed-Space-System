%% GUIDANCE MAIN CODE
% Constants
mu  = 3.986e14;       % [m^3/s^2]
J2  = 1.08263e-3;
Re  = 6378e3;         % [m]

% INPUT 1: Chief orbit parameters
a_c = 6771e3;
e_c = 0.0005;
i_c = deg2rad(51.64);
Omega_c = deg2rad(257);
omega_c = deg2rad(45);
nu_c = deg2rad(30);
E_c = 2*atan2(tan(nu_c/2), sqrt((1+e_c)/(1-e_c)));
M_c = E_c - e_c*sin(E_c);
chief_OE = [a_c, e_c, i_c, Omega_c, omega_c, nu_c];

% INPUT 2: Initial deputy relative singular OEs
da = 0;
dM = deg2rad(-5);
de = 0.001;
domega = deg2rad(0.5);
di = deg2rad(1);
dOmega = deg2rad(1);
deputy_rSOE0 = [da, dM, de, domega, di, dOmega];

% INPUT 3: Desired deputy relative singular OEs
da_des = 0;
dM_des = deg2rad(-0.05);
de_des = 0;
domega_des = 0;
di_des = 0;
dOmega_des = 0;
deputy_rSOEd = [da_des, dM_des, de_des, domega_des, di_des, dOmega_des];
deputy_OE_0 = chief_OE + deputy_rSOEd;

% INPUT 4: Tau + Waypoints
tau = 12*3600;
waypoint_times = [0.4*tau; 0.7*tau; 0.9*tau];

% Compute pseudo state
b = pseudo_state(chief_OE, deputy_rSOE0, deputy_rSOEd, tau, mu, J2, Re);

% Compute Guidance A matrix
A = zeros(6, 6*length(waypoint_times));
for i = 1:length(waypoint_times)
    dt = tau - waypoint_times(i);
    Phi = stm_qns_J2(a_c, e_c, i_c, omega_c, mu, J2, Re, dt);
    A(:, 6*i-5:6*i) = Phi;
end

% GUIDANCE WAYPOINTS
waypoints = inv(A'*A)*A'*b;
disp('GUIDANCE WAYPOINT 1:')
disp(waypoints(1:6))
disp('GUIDANCE WAYPOINT 2:')
disp(waypoints(7:12))
disp('GUIDANCE WAYPOINT 3:')
disp(waypoints(13:18))


%% CONTROL MAIN CODE
a = chief_OE(1);
n = sqrt(mu / a^3);
uk = M_c+omega_c;

nw = length(waypoint_times);
A_control = zeros(6,3*nw);

for j = 1:nw
    dt   = tau - waypoint_times(j);
    Phi  = stm_qns_J2(a_c, e_c, i_c, omega_c, mu, J2, Re, dt);
    
    % argument of latitude at the burn:
    % note: you can compute the deputy's nu at t_j if you want a more accurate uk_j
    uk_j = omega_c + nu_c;  % or update nu_c for true anomaly at t_j
    
    % 6×3 sensitivity
    Gamma = computeGamma(chief_OE, uk_j, mu);
    
    % assemble the 6×3 block
    A_control(:, 3*(j-1)+1 : 3*j) = Phi * Gamma;
end

% IMPULSES
delta_vs = (A_control'*A_control)\(A_control'*b);
disp('IMPULSE 1:')
disp(delta_vs(1:3))
disp('IMPULSE 2:')
disp(delta_vs(4:6))
disp('IMPULSE 3:')
disp(delta_vs(7:9))

%% FUNCTIONS
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
    W_d = W_c + delta_i_y / sin(i_c);
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


function Phi = stm_qns_J2(a, e, inc, omega, mu, J2, Re, tau)

  % mean motion
  n = sqrt(mu/(a^3));

  % eccentricity-dependent terms
  eta   = sqrt(1 - e^2);
  kappa = (3*J2*Re^2*sqrt(mu))/(4*a^(7/2)*eta^4);
  E     = 1 + eta;
  F     = 4 + 3*eta;
  G     = 1/eta^2;

  % inclination-dependent terms
  P = 3*cos(inc)^2 - 1;
  Q = 5*cos(inc)^2 - 1;
  S = sin(2*inc);
  T = sin(inc)^2;

  % combined factors
  EP  = E*P;
  FGP = F*G*P;
  GQ = G*Q;

  % exi, exf, eyi, eyf
  omega_dot = kappa*Q;
  omega_f = omega + omega_dot*tau;
  exi = e*cos(omega);
  eyi = e*sin(omega);
  exf = e*cos(omega_f);
  eyf = e*sin(omega_f);


  % assemble the STM
  Phi = [
    1,                                  0,                  0,                    0,                   0,                0;
   -(1.5*n + 1.5*kappa*EP)*tau,        1,     kappa*exi*FGP*tau,     kappa*eyi*FGP*tau,    -kappa*F*S*tau,    0;
    3.5*kappa*eyf*Q*tau,                 0, cos(omega_dot*tau)-4*kappa*exi*eyf*GQ*tau, -sin(omega_dot*tau)-4*kappa*eyi*eyf*GQ*tau,  5*kappa*eyf*S*tau, 0;
   -3.5*kappa*exf*Q*tau,                 0, sin(omega_dot*tau)+4*kappa*exi*exf*GQ*tau,  cos(omega_dot*tau)+4*kappa*eyi*exf*GQ*tau, -5*kappa*exf*S*tau, 0;
    0,                                  0,                  0,                    0,                   1,                0;
    3.5*kappa*S*tau,                    0, -4*kappa*exi*G*S*tau,   -4*kappa*eyi*G*S*tau,  2*kappa*T*tau,      1
  ];
end


function rQNSOE = OE2rQNSOE(chief_OE, deputy_OE)
% Input: chief_OE = [a_c, e_c, i_c, Omega_c, omega_c, M_c]
%        deputy_OE = [a_d, e_d, i_d, Omega_d, omega_d, M_d]

a_c     = chief_OE(1);
e_c     = chief_OE(2);
i_c     = chief_OE(3);
Omega_c = chief_OE(4);
omega_c = chief_OE(5);
M_c     = chief_OE(6);

a_d     = deputy_OE(1);
e_d     = deputy_OE(2);
i_d     = deputy_OE(3);
Omega_d = deputy_OE(4);
omega_d = deputy_OE(5);
M_d     = deputy_OE(6);

rQNSOE = [
    (a_d - a_c)/a_c;
    (M_d + omega_d) - (M_c + omega_c) + (Omega_d - Omega_c)*cos(i_c);
    e_d*cos(omega_d) - e_c*cos(omega_c);
    e_d*sin(omega_d) - e_c*sin(omega_c);
    i_d - i_c;
    (Omega_d - Omega_c)*sin(i_c)
];
end


function b = pseudo_state(chief_OE, deputy_rSOE0, deputy_rSOEd, tau, mu, J2, Re)
% Pull inputs
a_c        = chief_OE(1);
e_c        = chief_OE(2);
i_c        = chief_OE(3);
Omega_c    = chief_OE(4);
omega_c    = chief_OE(5);
nu_c       = chief_OE(6);
E_c = 2*atan2(tan(nu_c/2), sqrt((1+e_c)/(1-e_c)));
M_c = E_c - e_c*sin(E_c);

da         = deputy_rSOE0(1);
dM         = deputy_rSOE0(2);
de         = deputy_rSOE0(3);
domega     = deputy_rSOE0(4);
di         = deputy_rSOE0(5);
dOmega     = deputy_rSOE0(6);

da_des     = deputy_rSOEd(1);
dM_des     = deputy_rSOEd(2);
de_des     = deputy_rSOEd(3);
domega_des = deputy_rSOEd(4);
di_des     = deputy_rSOEd(5);
dOmega_des = deputy_rSOEd(6);

% Compute initial deputy singular absolute OEs
a_d = a_c+da;
e_d = e_c+de;
i_d = i_c+di;
Omega_d = Omega_c+dOmega;
omega_d = omega_c+domega;
M_d = M_c+dM;
E_d = newton_raphson(M_d, e_d);
nu_d = 2 * atan2(sqrt(1+e_d)*sin(E_d/2), sqrt(1-e_d)*cos(E_d/2));
deputy_OE_0 = [a_d, e_d, i_d, Omega_d, omega_d, nu_d];

% Compute desired deputy orbit singular absolute OEs
a_d_des = a_c+da_des;
e_d_des = e_c+de_des;
i_d_des = i_c+di_des;
Omega_d_des = Omega_c+dOmega_des;
omega_d_des = omega_c+domega_des;
M_d_des = M_c+dM_des;
E_d_des = newton_raphson(M_d_des, e_d_des);
nu_d_des = 2 * atan2(sqrt(1+e_d_des)*sin(E_d_des/2), sqrt(1-e_d_des)*cos(E_d_des/2));
deputy_OE_des = [a_d_des, e_d_des, i_d_des, Omega_d_des, omega_d_des, nu_d_des];

% Compute initial and desired Deputy QNS ROEs
deputy_rQNSOE_0 = OE2rQNSOE(chief_OE, deputy_OE_0);
deputy_rQNSOE_des = OE2rQNSOE(chief_OE, deputy_OE_des);

% Compute Phi and pseudo state
Phi = stm_qns_J2(a_c, e_c, i_c, omega_c, mu, J2, Re, tau);
b = deputy_rQNSOE_des - Phi*deputy_rQNSOE_0;
end

%{
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
%}

function Gamma = computeGamma(Chief_OE, uk_j, mu)
a = Chief_OE(1);
n = sqrt(mu / a^3);
Gamma = 1/(n*a) * [
     0,          2,           0;
    -2,          0,           0;
     sin(uk_j),  2*cos(uk_j), 0;
    -cos(uk_j),  2*sin(uk_j), 0;
     0,          0,           cos(uk_j);
     0,          0,           sin(uk_j)
];
end


% Inputs (Chief OEs, Deputy OEs, and burns) were computed in P2
x0 = [chief_OE'; deputy_OE_0'];
burns = [
  waypoint_times(1),  delta_vs(1),  delta_vs(2),  delta_vs(3);
  waypoint_times(2),  delta_vs(4),  delta_vs(5),  delta_vs(6);
  waypoint_times(3),  delta_vs(7),  delta_vs(8),  delta_vs(9);
];

% GROUND TRUTH: DEPUTY RTN COORDS
[t_out, state_out] = groundTruth(x0, burns, mu, tau);
x1 = state_out(:, 7);
y1 = state_out(:, 8);
z1 = state_out(:, 9);
vx1 = state_out(:, 10);
vy1 = state_out(:, 11);
vz1 = state_out(:, 12);


% WITH ERROR: DEPUTY RTN COORDS
error_level = 0.05;  % 5% Gaussian error
delta_vs_err = delta_vs .* (1 + error_level * randn(size(delta_vs)));
burns_err = [
    waypoint_times(1), delta_vs_err(1), delta_vs_err(2), delta_vs_err(3);
    waypoint_times(2), delta_vs_err(4), delta_vs_err(5), delta_vs_err(6);
    waypoint_times(3), delta_vs_err(7), delta_vs_err(8), delta_vs_err(9);
];
[t_err, state_err] = groundTruth(x0, burns_err, mu, tau);
x_err  = state_err(:, 7);
y_err  = state_err(:, 8);
z_err  = state_err(:, 9);
vx_err = state_err(:,10);
vy_err = state_err(:,11);
vz_err = state_err(:,12);

% ERROR IN TRUE VS ACTUAL RTN COORDS
dx  = x_err  - x1;
dy  = y_err  - y1;
dz  = z_err  - z1;
dvx = vx_err - vx1;
dvy = vy_err - vy1;
dvz = vz_err - vz1;


function [t_out,state_out] = groundTruth(x0,burns,mu, tau)
alpha0 = x0(1:6);
alpha1 = x0(7:12);
[r0,v0] = utils.OE2ECI(alpha0,mu);
[r1,v1] = utils.OE2ECI(alpha1,mu);
[rRTN,vRTN] = utils.ECI2RTN(r0,v0,r1,v1);
state0 = [r0; v0; rRTN; vRTN];
t0 = 0;
orbit_period = 2*pi*sqrt(alpha0(1)^3/mu);
t_end     = tau;
t_current = t0;
state_current = state0;
t_out = [];
state_out = [];
for k = 1:size(burns,1)+1
    if k<=size(burns,1)
        t_next = burns(k,1);
    else
        t_next = t_end;
    end
    tspan = t_current:0.5:t_next;
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [t_seg,s_seg] = ode113(@nonlinear_state_dot,tspan,state_current,options);
    t_out = [t_out; t_seg];
    state_out = [state_out; s_seg];
    state_current = s_seg(end,:)';
    if k<=size(burns,1)
        dv = burns(k,2:4)';
        state_current(10:12) = state_current(10:12) + dv;
        t_current = t_next;
    end
end
end




%% Initialization
% Constants
mu = 3.986004418e14;
Re = 6378e3;
J2 = 1.082636e-3;

% Chief's singular orbit elements
a_c = 6771e3;
e_c = 5.0e-4;
i_c = deg2rad(51.64);
Omega_c = deg2rad(257);
omega_c = deg2rad(45);
nu_c = deg2rad(30);
alpha_c = [a_c; e_c; i_c; Omega_c; omega_c; nu_c];

% Convert SOEs to R/V for chief state
[r0, v0] = utils.OE2ECI(alpha_c, mu);

% Initial relative orbital elements
initial_singular_ROEs = [0; deg2rad(-5); 0.001; deg2rad(0.5); deg2rad(1); deg2rad(1)];
initial_QNSROE = utils.OE2QNSOE(initial_singular_ROEs);
initial_QNSROE = [initial_QNSROE(1), initial_QNSROE(3), initial_QNSROE(4), initial_QNSROE(5), initial_QNSROE(6)];
initial_singular_OEs = alpha_c + initial_singular_ROEs;

% Desired relative orbital elements
desired_singular_ROEs = [0; deg2rad(-0.05); 0; 0; 0; 0];
desired_QNSROE = utils.OE2QNSOE(desired_singular_ROEs);
desired_QNSROE = [desired_QNSROE(1), desired_QNSROE(3), desired_QNSROE(4), desired_QNSROE(5), desired_QNSROE(6)];

%% SIMULATION
% Simulation parameters
state0 = [r0; v0; initial_QNSROE'; 0];
control_error = false;
k = 500;
disp(state0(7:12))
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
perturbated = true;
orbital_period = 2 * pi * sqrt(a_c^3 / mu);
tstart = 0;
tint = 10;
duration = 4 * orbital_period;
tspan = (tstart:tint:duration)';

% Run simulation
odefun = @(t, state) propogate(t, state, desired_QNSROE, mu, Re, J2, k, control_error);
[tout, stateout] = ode113(odefun, tspan, state0, options);
disp('Done part 1')

%% GROUND TRUTH SIMULATION
% Pull u values
N = numel(tout);
U = zeros(N,3);

for idx = 1:N
    % unpack chief ECI pos/vel and QNS-ROE
    r_chief = stateout(idx,1:3)';
    v_chief = stateout(idx,4:6)';
    roe = stateout(idx,7:11)';
    adot = stateout(idx,12);

    % recompute Abar, Bbar, P
    [a_c, e_c, i_c, W_c, omega_c, f_c] = utils.ECI2OE(r_chief, v_chief, mu);
    z = [roe; adot];
    err = roe - desired_QNSROE';
    [Abar, Bbar, P] = compute_ABP(a_c, e_c, i_c, omega_c, f_c, err(3), err(2), err(5), err(4), mu, Re, J2, k);

    v_vec = Abar*z + [P*err; 0];
    u_full = -Bbar \ v_vec;
    U(idx,:) = [0, u_full(1), u_full(2)];
end

x0 = [alpha_c; initial_singular_OEs];
[toutGT, stateoutGT] = groundTruth(x0, U, mu, tspan, options);

N = numel(toutGT);
QNSROEs = zeros(N,5);

for idx = 1:N
    rC = stateoutGT(idx,1:3)';
    vC = stateoutGT(idx,4:6)';
    rRTN = stateoutGT(idx,7:9)';
    vRTN = stateoutGT(idx,10:12)';

    [rD, vD] = RTN2ECI(rC, vC, rRTN, vRTN);

    [a, e, i, W, w, f] = ECI2OE(rC, vC, mu);
    alphaC = [a, e, i, W, w, f];
    [a, e, i, W, w, f] = ECI2OE(rD, vD, mu);
    alphaD = [a, e, i, W, w, f];

    Qc = OE2QNSOE(alphaC);
    Qd = OE2QNSOE(alphaD);

    dQ = Qd - Qc;
    QNSROEs(idx,:) = [dQ(1), dQ(3), dQ(4), dQ(5), dQ(6)];
end
disp('Done ground truth')

%% GROUND TRUTH SIMULATION WITH ERRORS IN U INPUTS
error_level = 0.5;
U_err = U .* (1 + error_level * randn(size(U)));

x0 = [alpha_c; initial_singular_OEs];
[toutGTE, stateoutGTE] = groundTruth(x0, U_err, mu, tspan, options);

N = numel(toutGTE);
QNSROEs_err = zeros(N,5);

for idx = 1:N
    rC = stateoutGTE(idx,1:3)';
    vC = stateoutGTE(idx,4:6)';
    rRTN = stateoutGTE(idx,7:9)';
    vRTN = stateoutGTE(idx,10:12)';

    [rD, vD] = RTN2ECI(rC,vC,rRTN,vRTN);

    [a, e, i, W, w, f] = ECI2OE(rC, vC, mu);
    alphaC = [a, e, i, W, w, f];
    [a, e, i, W, w, f] = ECI2OE(rD, vD, mu);
    alphaD = [a, e, i, W, w, f];

    Qc = OE2QNSOE(alphaC);
    Qd = OE2QNSOE(alphaD);

    dQ = Qd - Qc;
    QNSROEs_err(idx,:) = [dQ(1), dQ(3), dQ(4), dQ(5), dQ(6)];
end

%% PLOTTING RESULTS

figure(1)
subplot(3,1,1)
plot(tout/orbital_period, stateout(:,7), 'b-')
hold on
yline(desired_QNSROE(1), '--r')
hold off
ylabel('delta a / a')

subplot(3,1,2)
plot(tout/orbital_period, stateout(:,8), 'b-')
hold on
yline(desired_QNSROE(2), '--r')
hold off
ylabel('delta e_x')

subplot(3,1,3)
plot(tout/orbital_period, stateout(:,9), 'b-')
hold on
yline(desired_QNSROE(3), '--r')
hold off
ylabel('delta e_y')
xlabel('Time (orbital periods)')
legend('Actual', 'Desired', 'Location', 'best')

figure(2)
subplot(3,1,1)
plot(tout/orbital_period, stateout(:,10), 'b-')
hold on
yline(desired_QNSROE(4), '--r')
hold off
ylabel('delta i_x')

subplot(3,1,2)
plot(tout/orbital_period, stateout(:,11), 'b-')
hold on
yline(desired_QNSROE(5), '--r')
hold off
ylabel('delta i_y')

subplot(3,1,3)
plot(tout/orbital_period, stateout(:,12), 'b-')
hold on
yline(0, '--r')
hold off
ylabel('d(delta a / a)/dt')
xlabel('Time (orbital periods)')
legend('Actual', 'Desired', 'Location', 'best')

figure(3)
plot(tout/orbital_period, U(:,2:3))
legend('u_T', 'u_N')
xlabel('Time (orbital periods)')
ylabel('Control input u')

figure(4)
subplot(3,1,1)
plot(toutGT/orbital_period, QNSROEs(:,1), 'b-')
hold on
plot(toutGTE/orbital_period, QNSROEs_err(:,1), 'g-')
yline(desired_QNSROE(1), '--r')
hold off
ylabel('delta a / a')

subplot(3,1,2)
plot(toutGT/orbital_period, QNSROEs(:,2), 'b-')
hold on
plot(toutGTE/orbital_period, QNSROEs_err(:,2), 'g-')
yline(desired_QNSROE(2), '--r')
hold off
ylabel('delta e_x')

subplot(3,1,3)
plot(toutGT/orbital_period, QNSROEs(:,3), 'b-')
hold on
plot(toutGTE/orbital_period, QNSROEs_err(:,3), 'g-')
yline(desired_QNSROE(3), '--r')
hold off
ylabel('delta e_y')
xlabel('Time (orbital periods)')
legend('Ground Truth', 'Errored', 'Desired', 'Location', 'best')

figure(5)
subplot(3,1,1)
plot(toutGT/orbital_period, QNSROEs(:,4), 'b-')
hold on
plot(toutGTE/orbital_period, QNSROEs_err(:,4), 'g-')
yline(desired_QNSROE(4), '--r')
hold off
ylabel('delta i_x')

subplot(3,1,2)
plot(toutGT/orbital_period, QNSROEs(:,5), 'b-')
hold on
plot(toutGTE/orbital_period, QNSROEs_err(:,5), 'g-')
yline(desired_QNSROE(5), '--r')
hold off
ylabel('delta i_y')

subplot(3,1,3)
plot(toutGT/orbital_period, stateoutGT(:,12), 'b-')
hold on
plot(toutGTE/orbital_period, stateoutGTE(:,12), 'g-')
yline(0, '--r')
hold off
ylabel('d(delta a / a)/dt')
xlabel('Time (orbital periods)')
legend('Ground Truth', 'Errored', 'Desired', 'Location', 'best')


function [t_out, state_out] = groundTruth(x0, U, mu, tspan, options)
alpha0 = x0(1:6);
alpha1 = x0(7:12);
[r0, v0] = utils.OE2ECI(alpha0, mu);
[r1, v1] = utils.OE2ECI(alpha1, mu);
[rRTN, vRTN] = utils.ECI2RTN(r0, v0, r1, v1);
state0 = [r0; v0; rRTN; vRTN];

[t_out, state_out] = ode113(@(t, s) continuousRHS(t, s, tspan, U, mu), tspan, state0, options);
end


function sd = continuousRHS(t, state, timeVec, U, mu)
sd = nonlinear_state_dot(t, state);

u = interp1(timeVec, U, t, 'linear', 0);

sd(10:12) = sd(10:12) + u(:);
end

function statedot = propogate(t, state, desired_QNSROE, mu, Re, J2, k, control_error)
statedot = zeros(12,1);

r_chief_ECI = state(1:3);
v_chief_ECI = state(4:6);
[a_c, e_c, i_c, W_c, omega_c, f_c] = utils.ECI2OE(r_chief_ECI, v_chief_ECI, mu);

deputy_QNSROE = state(7:11);
delta_adot = state(12);
z = [deputy_QNSROE; delta_adot];

error = deputy_QNSROE - desired_QNSROE';
d_d_ex = error(2);
d_d_ey = error(3);
d_d_ix = error(4);
d_d_iy = error(5);

[Abar, Bbar, P] = compute_ABP(a_c, e_c, i_c, omega_c, f_c, d_d_ey, d_d_ex, d_d_iy, d_d_ix, mu, Re, J2, k);

v = Abar*z + [P*error; 0];
u = -Bbar \ v;

if control_error
    error_level = 0.0001;
    u = u .* (1 + error_level*randn(size(u)));
end

statedot(7:12) = Abar*z + Bbar*u;

r = norm(r_chief_ECI);
accel_kepler = -mu / r^3 * r_chief_ECI;

if 1 > 1
    perturbations = getPerturbations(r_chief_ECI, mu, Re, J2);
else
    perturbations = zeros(3,1);
end

statedot(1:3) = v_chief_ECI;
statedot(4:6) = accel_kepler + perturbations;
end

function [Abar, Bbar, P] = compute_ABP(a_c, e_c, i_c, omega_c, f_c, d_d_ey, d_d_ex, d_d_iy, d_d_ix, mu, Re, J2, k)
n_c = sqrt(mu/a_c^3);
eta_c = sqrt(1 - e_c^2);
gamma = 0.75 * J2 * Re^2 * sqrt(mu);
kappa = gamma / (a_c^(7/2) * eta_c^4);

ex = e_c * cos(omega_c);
ey = e_c * sin(omega_c);

C = sin(omega_c);
D = cos(omega_c);
G = 1 / eta_c^2;

P_factor = 3 * cos(i_c)^2 - 1;
Q_factor = 5 * cos(i_c)^2 - 1;
S = sin(2*i_c);
T = sin(i_c)^2;

Abar = kappa * [
    0, 0, 0, 0, 0, 1/kappa;
    3.5*ey*Q_factor, -(4*ex*ey*G + C)*Q_factor, -(1+4*ey^2*G-D)*Q_factor, 5*ey*S, 0, (D-ex)/kappa;
    -3.5*ex*Q_factor, (1+4*ex^2*G-D)*Q_factor, (4*ex*ey*G-C)*Q_factor, -5*ex*S, 0, (C-ey)/kappa;
    0, 0, 0, 0, 0, 0;
    3.5*S, -4*ex*G*S, -4*ey*G*S, 2*T, 0, 0;
    0, 0, 0, 0, 0, 0
];

Bbar = 1/(a_c*n_c) * [
    (2/eta_c)*(1 + e_c*cos(f_c)), 0;
    eta_c*((2+e_c*cos(f_c))*cos(omega_c+f_c)+ex)/(1+e_c*cos(f_c)), (eta_c*ey*sin(omega_c+f_c))/(tan(i_c)*(1+e_c*cos(f_c)));
    eta_c*((2+e_c*cos(f_c))*sin(omega_c+f_c)+ey)/(1+e_c*cos(f_c)), (-eta_c*ex*sin(omega_c+f_c))/(tan(i_c)*(1+e_c*cos(f_c)));
    0, eta_c*(cos(omega_c+f_c)/(1+e_c*cos(f_c)));
    0, eta_c*(sin(omega_c+f_c)/(1+e_c*cos(f_c)));
    0, 0
];

phi = omega_c + f_c;
phi_ip = atan2(d_d_ey, d_d_ex);
phi_oop = atan2(d_d_iy, d_d_ix);

J = phi - phi_ip;
H = phi - phi_oop;
N = 2;

P = (1/k) * diag([
    cos(J)^N,
    cos(J)^N,
    cos(J)^N,
    cos(H)^N,
    cos(H)^N
]);
end



