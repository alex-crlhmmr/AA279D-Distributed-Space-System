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
w_0 = deg2rad(0);       % argument of perigee [rad]
f_0 = deg2rad(30);      % true anomaly [rad]

% Chief initial state in OE & ECI
initial_state_0_OE = [a_0, e_0, i_0, W_0, w_0, f_0];
initial_state_0_ECI = utils.OE2ECI(initial_state_0_OE, const, body);
r0_init = initial_state_0_ECI(1:3);
v0_init = initial_state_0_ECI(4:6);

% Deputy initial orbital elements (ISS)
a_1 = a_0;                 % semi-major axis [m]
e_1 = e_0+0.0001;          % eccentricity [-]
i_1 = i_0+deg2rad(0.05);   % inclination [rad]
W_1 = W_0+deg2rad(0.05);   % RAAN [rad]
w_1 = w_0+deg2rad(0.05);   % argument of perigee [rad]
f_1 = f_0 - deg2rad(0.05); % true anomaly [rad]

% Deputy initial state in OE & ECI
initial_state_1_OE = [a_1, e_1, i_1, W_1, w_1, f_1];
initial_state_1_ECI = utils.OE2ECI(initial_state_1_OE, const, body);
r1_init = initial_state_1_ECI(1:3);
v1_init = initial_state_1_ECI(4:6);

% Compute relative separation and its fraction of the chief’s radius
rho  = norm(r1_init - r0_init);      
frac = rho / norm(r0_init);           

% Display results
fprintf('Relative separation ρ = %.2f m (%.4f × r₀)\n', rho, frac);
fprintf('Eccentricities: e_0 = %.4f, e_1 = %.4f\n', e_0, e_1);

% Assertions to ensure HCW validity
assert(frac < 0.005, 'ERROR: ρ/r₀ ≥ 0.5 %');
assert(e_0 < 0.01 && e_1 < 0.01, 'ERROR: eccentricities too large');

disp('→ Initial conditions satisfy CW-HCW assumptions.');



%% 1 - b

% Chief
fprintf('\n\nChief (ECI):\n');
fprintf('  r0 = [%.3f, %.3f, %.3f] km\n', r0_init/1e3);
fprintf('  v0 = [%.3f, %.3f, %.3f] km/s\n', v0_init/1e3);
oe_c = initial_state_0_OE;
fprintf('  OE_c = [a=%.3f km, e=%.6f, i=%.3f°, Ω=%.3f°, ω=%.3f°, f=%.3f°]\n\n',...
    oe_c(1)/1e3, oe_c(2), rad2deg(oe_c(3:6)));

% Deputy
fprintf('Deputy (ECI):\n');
fprintf('  r1 = [%.3f, %.3f, %.3f] km\n', r1_init/1e3);
fprintf('  v1 = [%.3f, %.3f, %.3f] km/s\n', v1_init/1e3);
oe_d = initial_state_1_OE;
fprintf('  OE_d = [a=%.3f km, e=%.6f, i=%.3f°, Ω=%.3f°, ω=%.3f°, f=%.3f°]\n\n',...
    oe_d(1)/1e3, oe_d(2), rad2deg(oe_d(3:6)));


% Deputy relative initital state in ECI
dr_init_ECI  = r1_init - r0_init;
dv_init_ECI = v1_init - v0_init;

% Deputy relative initital state in chief's RTN
R_ECI2RTN = utils.R2RTN(r0_init,v0_init);
dr_init_RTN = R_ECI2RTN * dr_init_ECI;
r0_init_norm = norm(r0_init);
r_hat   = r0_init / r0_init_norm;
h_vec   = cross(r0_init, v0_init);
omega_init = [0; 0; norm(h_vec)/r0_init_norm^2];
dv_init_RTN_uncorrected = R_ECI2RTN * dv_init_ECI;
dv_init_RTN = dv_init_RTN_uncorrected - cross(omega_init, dr_init_RTN);


% Deputy relative in chief's RTN
fprintf('Relative (RTN) at t0:\n');
fprintf('  Δr_RTN = [%.3f, %.3f, %.3f] m\n', dr_init_RTN);
fprintf('  Δv_RTN = [%.6f, %.6f, %.6f] m/s\n\n', dv_init_RTN);


% Initial state (absolute chief state in ECI and deputy relative state in chief's RTN)
initial_state = [r0_init; v0_init; dr_init_RTN; dv_init_RTN];

% Initial Quasi-Non-Singular ROE 

% Mean anomalies 
M_0 = f_0;
M_1 = f_1;

% Quasi non singular orbital elements
delta_a = (a_1 - a_0)/a_0;
delta_lambda = (M_1 + w_1) - (M_0 + w_0) + (W_1 - W_0)*cos(i_0);
delta_ex = e_1*cos(w_1) - e_0*cos(w_0);
delta_ey = e_1*sin(w_1) - e_0*sin(w_0);
delta_ix = i_1 - i_0;
delta_iy = (W_1 - W_0)*sin(i_0);
qns_init = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy];

% Print QNS–ROE
fprintf('\nQuasi-Non-Singular ROE (using your _0/_1 suffixes):\n');
fprintf('δa=%+8.3e\n', delta_a);
fprintf('δλ=%+8.3e rad\n', delta_lambda);
fprintf('δε_x=%+8.3e\n', delta_ex);
fprintf('δε_y=%+8.3e\n', delta_ey);
fprintf('δi_x=%+8.3e rad\n', delta_ix);
fprintf('  δi_y   = %+8.3e rad\n\n', delta_iy);

%% 1 - c

% Initial RTN relative state at t0:
x0 = dr_init_RTN(1);
y0 = dr_init_RTN(2);
z0 = dr_init_RTN(3);
xdot0 = dv_init_RTN(1);
ydot0 = dv_init_RTN(2);
zdot0 = dv_init_RTN(3);

% Chief mean motion n
mu = const.(body).mu;             
r0 = norm(r0_init);              
n = sqrt(mu/r0^3);  

% Constants c1–c6
c1 = - (3*x0 + 2*ydot0/n);
c2 = xdot0/n;
c3 = 4*x0 + 2*ydot0/n;
%c3 = 0;
c4 = y0 - 2*xdot0/n;
c5 = z0;
c6 = zdot0/n;

% Print them
fprintf('\nHCW integration constants [c1 … c6]:\n');
fprintf('c1 = %+8.3e\n', c1);
fprintf('c2 = %+8.3e\n', c2);
fprintf('c3 = %+8.3e\n', c3);
fprintf('c4 = %+8.3e\n', c4);
fprintf('c5 = %+8.3e\n', c5);
fprintf('c6 = %+8.3e\n\n', c6);

%% 1 - d

% 15 chief orbital periods
T_orbit = 2*pi/n;
t_final = 15 * T_orbit;
N_steps = 20000;
t = linspace(0, t_final, N_steps);

% preallocate arrays
x = zeros(size(t));
y = zeros(size(t));
z = zeros(size(t));
xd = zeros(size(t));
yd = zeros(size(t));
zd = zeros(size(t));

for k = 1:length(t)
    tau = t(k);
    sn = sin(n*tau);
    cn = cos(n*tau);
    
    % position
    x(k) = c3 + c1*cn + c2*sn;
    y(k) = c4 - 2*c1*sn + 2*c2*cn -(3/2)*n*c3*tau;
    z(k) = c5*cn + c6*sn;
    
    % velocity (time derivatives)
    xd(k) = -c1*n*sn + c2*n*cn;
    yd(k) = -2*c1*n*cn - 2*c2*n*sn - (3/2)*n*c3;
    zd(k) = -c5*n*sn + c6*n*cn;
end



figure('Name','HCW Position Projections + 3D','Color','w','Units','normalized','Position',[.1 .1 .8 .6]);

% T–R Plane
subplot(2,2,1);
plot(x, y, 'LineWidth',1.2);
xlabel('T (m)'); ylabel('R (m)');
%title('T–R Plane');
axis equal; grid on;

% N–R Plane
subplot(2,2,2);
plot(z, x, 'LineWidth',1.2);
xlabel('N (m)'); ylabel('R (m)');
%title('N–R Plane');
axis equal; grid on;

% T–N Plane
subplot(2,2,3);
plot(y, z, 'LineWidth',1.2);
xlabel('T (m)'); ylabel('N (m)');
%title('T–N Plane');
axis equal; grid on;

% 3D Position Plot 
subplot(2,2,4);
plot3(x, y, z, 'LineWidth',1.2);
xlabel('R (m)'); ylabel('T (m)'); zlabel('N (m)');
%title('3D Relative Position');
axis equal; grid on;


figure('Name','HCW Velocity Projections + 3D','Color','w','Units','normalized','Position',[.1 .1 .8 .6]);

% Ṫ–Ṙ Plane 
subplot(2,2,1);
plot(yd, xd, 'LineWidth',1.2); ̇
xlabel('Ṫ (m/s)', 'Interpreter', 'latex');
ylabel('Ṙ (m/s)', 'Interpreter', 'latex');
%title('Ṫ–Ṙ Plane', 'Interpreter', 'latex');
axis equal; grid on;

% Ṅ–Ṙ Plane 
subplot(2,2,2);
plot(zd, xd, 'LineWidth',1.2);̇
xlabel('Ṅ (m/s)', 'Interpreter', 'latex');
ylabel('Ṙ (m/s)', 'Interpreter', 'latex');
%title('Ṅ–Ṙ Plane', 'Interpreter', 'latex');
axis equal; grid on;

% Ṫ–Ṅ Plane 
subplot(2,2,3);
plot(yd, zd, 'LineWidth',1.2);  ̇
xlabel('Ṫ (m/s)', 'Interpreter', 'latex');
ylabel('Ṅ (m/s)', 'Interpreter', 'latex');
%title('Ṫ–Ṅ Plane', 'Interpreter', 'latex');
axis equal; grid on;

% 3D Relative Velocity Plot *
subplot(2,2,4);
plot3(xd, yd, zd, 'LineWidth',1.2);
xlabel('Ṙ (m/s)', 'Interpreter', 'latex');
ylabel('Ṫ (m/s)', 'Interpreter', 'latex');
zlabel('Ṅ (m/s)', 'Interpreter', 'latex');
%title('3D Relative Velocity', 'Interpreter', 'latex');
axis equal; grid on;





%% 1 - d - bis

figure('Name','RTN Position over Time','Color','w','Units','normalized','Position',[.2 .2 .6 .6]);

% R vs t
subplot(3,1,1);
plot(t, x, 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('R (m)');
title('Radial (R) Component vs Time');
axis tight; grid on;

% T vs t
subplot(3,1,2);
plot(t, y, 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('T (m)');
title('Along-track (T) Component vs Time');
axis tight; grid on;

% N vs t
subplot(3,1,3);
plot(t, z, 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('N (m)');
title('Cross-track (N) Component vs Time');
axis tight; grid on;


figure('Name','RTN Velocity over Time','Color','w','Units','normalized','Position',[.2 .2 .6 .6]);

% Ṙ vs t
subplot(3,1,1);
plot(t, xd, 'LineWidth',1.2);
xlabel('Time (s)');
ylabel('Ṙ (m/s)', 'Interpreter', 'latex');
title('Radial Velocity (Ṙ) Component vs Time', 'Interpreter', 'latex');
axis tight; grid on;

% Ṫ vs t
subplot(3,1,2);
plot(t, yd, 'LineWidth',1.2);
xlabel('Time (s)');
ylabel('Ṫ (m/s)', 'Interpreter', 'latex');
title('Along-track Velocity (Ṫ) Component vs Time', 'Interpreter', 'latex');
axis tight; grid on;

% Ṅ vs t
subplot(3,1,3);
plot(t, zd, 'LineWidth',1.2);
xlabel('Time (s)');
ylabel('Ṅ (m/s)', 'Interpreter', 'latex');
title('Cross-track Velocity (Ṅ) Component vs Time', 'Interpreter', 'latex');
axis tight; grid on;
