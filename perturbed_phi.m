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


%% Ground Truth Simulation
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


%% GOUNRD TRUTH AND ERROR SIMULATION
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

% Plot 1: RTN Coords over time
figure(1)
subplot(3,1,1)
hold on
plot(t_out/3600,   x1,  'b')
plot(t_err/3600, x_err, '--r')
ylabel('R (m)')
legend('True','With error')

subplot(3,1,2)
hold on
plot(t_out/3600,   y1,  'b')
plot(t_err/3600, y_err, '--r')
ylabel('T (m)')

subplot(3,1,3)
hold on
plot(t_out/3600,   z1,  'b')
plot(t_err/3600, z_err, '--r')
ylabel('N (m)')
xlabel('Time (hours)')

% Plot 2: Error in RTN coords over time
figure(2)
subplot(3,1,1)
plot(t_out/3600,   dx,  'b')
ylabel('R_{error} (m)')

subplot(3,1,2)
plot(t_out/3600,   dy,  'b')
ylabel('T_{error} (m)')

subplot(3,1,3)
plot(t_out/3600,   dz,  'b')
ylabel('N_{error} (m)')
xlabel('Time (hours)')

% Plot 3: Delta v over time
dv = reshape(delta_vs, 3, [])';  % each row is [dV_R, dV_T, dV_N] for a burn

figure(3)
hold on
stem(waypoint_times/3600, dv(:,1), 'r', 'LineWidth', 1.5, 'Marker', 'o', 'DisplayName', '\delta v_R')
stem(waypoint_times/3600, dv(:,2), 'g', 'LineWidth', 1.5, 'Marker', 's', 'DisplayName', '\delta v_T')
stem(waypoint_times/3600, dv(:,3), 'b', 'LineWidth', 1.5, 'Marker', 'd', 'DisplayName', '\delta v_N')
xlabel('Time (hours)')
ylabel('\Delta v (m/s)')
legend('Location','best')
grid on


function [t_out,state_out] = groundTruth(x0,burns,mu, tau)
alpha0 = x0(1:6);
alpha1 = x0(7:12);
[r0,v0] = utils.OE2ECI(alpha0,mu);
[r1,v1] = utils.OE2ECI(alpha1,mu);
[rRTN,vRTN] = ECI2RTN(r0,v0,r1,v1);
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
