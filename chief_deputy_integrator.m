% Initial Keplerian orbital elements [a, e, i, RAAN, omega, nu]
alpha0 = [6771; 0.1005; deg2rad(51.64); deg2rad(257); deg2rad(0); deg2rad(30)]; % Chief
alpha1 = [6771; 0.1006; deg2rad(51.69); deg2rad(257.05); deg2rad(0.05); deg2rad(29.95)]; % Deputy

% Calculate initial state [ECI position and velocity] of chief
mu = 398600.4418;
[r0, v0] = utils.OE2ECI(alpha0, mu);

% Calculate initial state [RTN position and velocity] of deputy
[r1, v1] = utils.OE2ECI(alpha1, mu);
[rRTN, vRTN] = ECI2RTN(r0, v0, r1, v1);

norm(rRTN)

% Initialize state
state0 = zeros(12, 1);
state0(1:3) = r0;
state0(4:6) = v0;
state0(7:9) = rRTN;
state0(10:12) = vRTN;

% Time span
orbit_period = 2*pi*sqrt(alpha0(1)^3 / mu);
t0 = 0;
t_end =  15 * orbit_period; % Propagate for 15 orbits
tspan = t0:0.5:t_end;

% ODE options with precise tolerances
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% PROPOGATE
[t, state] = ode113(@nonlinear_state_dot, tspan, state0, options);

% PLOT
x1 = state(:, 7);
y1 = state(:, 8);
z1 = state(:, 9);
vx1 = state(:, 10);
vy1 = state(:, 11);
vz1 = state(:, 12);


%% YA SOLUTION
% Calculate initial state [ECI position and velocity] of chief
mu = 398600.4418;
[r0, v0] = utils.OE2ECI(alpha0, mu);

% Calculate initial state [RTN position and velocity] of deputy
[r1, v1] = utils.OE2ECI(alpha1, mu);

% Solve for constants
a = alpha0(1);         
e = alpha0(2);
f = alpha0(6);
mu = 398600.4418;    
n = sqrt(mu/a^3);
h = norm(cross(r0,v0));

% Calculate values needed to get constants
fdot = h/norm(r0)^2;
eta = sqrt(1-e^2);
k = 1 + e*cos(f);
c = k*cos(f);
s = k*sin(f);

% Initial state
[rRTN, vRTN] = ECI2RTN(r0, v0, r1, v1);
xbar0 = [ ...
   rRTN(1)/a;     
   vRTN(1)/(fdot*a); 
   rRTN(2)/a;  
   vRTN(2)/(fdot*a);    
   rRTN(3)/a;         
   vRTN(3)/(fdot*a); 
];


% Get constants
phi_inv = (1/eta^2) * [
    -3*s*(k+e^2)/k^2,            c - 2*e,     0,  -s*(k+1)/k,    0,               0;
    -3*(e + c/k),                -s,          0,  -(c*(k+1)/k + e), 0,            0;
     3*k - eta^2,                e*s,         0,   k^2,            0,               0;
    -3*e*s*(k+1)/k^2,           -2 + e*c,     eta^2, -e*s*(k+1)/k,  0,               0;
     0,                          0,           0,   0,              eta^2*cos(f),   -eta^2*sin(f);
     0,                          0,           0,   0,              eta^2*sin(f),    eta^2*cos(f)
];

Constants = phi_inv * xbar0;

% Propogate
step_size = (15*2*pi) / (length(tspan)-1);
f_change = 15*2*pi;
[f_vals, states] = get_YA_states(alpha0, Constants, step_size, f_change, n, mu, h);

x = states(:, 1);
vx = states(:, 2);
y = states(:, 3);
vy = states(:, 4);
z = states(:, 5);
vz = states(:, 6);

% Remove normalization
x = x*a;
vx = vx*a*fdot;
y = y*a;
vy = vy*a*fdot;
z = z*a;
vz = vz*a*fdot;



%% QUASI SOLUTION
a_0 = alpha0(1);
e_0 = alpha0(2);
i_0 = alpha0(3);
W_0 = alpha0(4);
w_0 = alpha0(5);
f_0 = alpha0(6);

a_1 = alpha1(1);
e_1 = alpha1(2);
i_1 = alpha1(3);
W_1 = alpha1(4);
w_1 = alpha1(5);
f_1 = alpha1(6);

E_0 = 2*atan2(sqrt((1 - e_0)/(1 + e_0))*tan(f_0/2), 1);
M_0 = E_0 - e_0*sin(E_0);

E_1 = 2*atan2(sqrt((1 - e_1)/(1 + e_1))*tan(f_1/2), 1);
M_1 = E_1 - e_1*sin(E_1);

delta_a = (a_1 - a_0)/a_0;
delta_lambda = (M_1 + w_1) - (M_0 + w_0) + (W_1 - W_0)*cos(i_0);
delta_ex = e_1*cos(w_1) - e_0*cos(w_0);
delta_ey = e_1*sin(w_1) - e_0*sin(w_0);
delta_ix = i_1 - i_0;
delta_iy = (W_1 - W_0)*sin(i_0);
qns_init = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy];
% Print QNS–ROE
fprintf('\nQuasi-Non-Singular ROE (using your _0/_1 suffixes):\n');
fprintf('  δa     = %+8.3e\n', delta_a);
fprintf('  δλ     = %+8.3e rad\n', delta_lambda);
fprintf('  δε_x   = %+8.3e\n', delta_ex);
fprintf('  δε_y   = %+8.3e\n', delta_ey);
fprintf('  δi_x   = %+8.3e rad\n', delta_ix);
fprintf('  δi_y   = %+8.3e rad\n\n', delta_iy);

quasi_elements = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy]';
alpha0 = [6771; 0.1005; deg2rad(51.64); deg2rad(257); deg2rad(0); deg2rad(30)]; % Chief

% [a, ex, ey, u, ix, iy

[f_vals2, states] = linear_ecc_mapping(quasi_elements, alpha0, 15*pi*2, step_size);

x2 = states(:, 1);
y2 = states(:, 2);
z2 = states(:, 3);
vx2 = states(:, 4);
vy2 = states(:, 5);
vz2 = states(:, 6);


%% ERROR
figure(4)
hold on
erYA = [x1, y1, z1] - [x, y, z];
erQNS = [x1, y1, z1] - [x2, y2, z2];
xYA = erYA(:,1);
yYA = erYA(:,2);
zYA = erYA(:,3);

xQN = erQNS(:,1);
yQN = erQNS(:,2);
zQN = erQNS(:,3);

%% YA only
figure(4)
subplot(2,2,1)
plot(xYA, yYA)
xlabel('R (km)')
ylabel('T (km)')
axis equal

subplot(2,2,2)
plot(xYA, zYA)
xlabel('R (km)')
ylabel('N (km)')
axis equal

subplot(2,2,3)
plot(yYA, zYA)
xlabel('T (km)')
ylabel('N (km)')
axis equal

subplot(2,2,4)
plot3(xYA, yYA, zYA)
xlabel('R (km)')
ylabel('T (km)')
zlabel('N (km)')
axis equal
view(3)
axis vis3d


figure(5)
subplot(2,2,1)
plot(xQN, yQN)
xlabel('R (km)')
ylabel('T (km)')
axis equal

subplot(2,2,2)
plot(xQN, zQN)
xlabel('R (km)')
ylabel('N (km)')
axis equal

subplot(2,2,3)
plot(yQN, zQN)
xlabel('T (km)')
ylabel('N (km)')
axis equal

subplot(2,2,4)
plot3(xQN, yQN, zQN)
xlabel('R (km)')
ylabel('T (km)')
zlabel('N (km)')
axis equal
view(3)
axis vis3d




%% PLOTS
figure(1)
subplot(3,1,1)
hold on
plot(t/orbit_period, x1)
plot(f_vals / (2*pi), x);
plot(f_vals2 / (2*pi), x2);
ylabel('R (km)')

subplot(3,1,2)
hold on
plot(t/orbit_period, y1)
plot(f_vals / (2*pi), y);
plot(f_vals2 / (2*pi), y2);
ylabel('T (km)')

subplot(3,1,3)
hold on
plot(t/orbit_period, z1)
plot(f_vals / (2*pi), z);
plot(f_vals2 / (2*pi), z2);
ylabel('N (km)')
xlabel('Orbital Periods')
legend('Nonlinear', 'QNS','Location','best')


figure(2)
subplot(2,2,1)
hold on
plot(x1, y1)
plot(x, y)
plot(x2, y2)
xlabel('R (km)')
ylabel('T (km)')
axis equal

subplot(2,2,2)
hold on
plot(x1, z1)
plot(x, z)
plot(x2, z2)
xlabel('R (km)')
ylabel('N (km)')
axis equal

subplot(2,2,3)
hold on
plot(y1, z1)
plot(y, z)
plot(y2, z2)
xlabel('T (km)')
ylabel('N (km)')
axis equal
hold on

subplot(2,2,4)
hold on
plot3(x1, y1, z1)
plot3(x, y, z)
plot3(x2, y2, z2)
xlabel('R (km)')
ylabel('T (km)')
zlabel('N (km)')
axis equal
view(3)   
axis vis3d
legend('Nonlinear', 'YA', 'QNS','Location','best')

figure(3)
subplot(3,1,1)
hold on
plot(t/orbit_period, vx1)
plot(f_vals / (2*pi), vx)
plot(f_vals2 / (2*pi), vx2)
ylabel('V_R (km/s)')

subplot(3,1,2)
hold on
plot(t/orbit_period, vy1)
plot(f_vals / (2*pi), vy)
plot(f_vals2 / (2*pi), vy2)
ylabel('V_T (km/s)')

subplot(3,1,3)
hold on
plot(t/orbit_period, vz1)
plot(f_vals / (2*pi), vz)
plot(f_vals2 / (2*pi), vz2)
ylabel('V_N (km/s)')
xlabel('Orbital Periods')
legend('Nonlinear', 'YA', 'QNS','Location','best')



figure(6)
subplot(2,2,1)
hold on
plot(vx1, vy1)
plot(vx, vy)
plot(vx2, vy2)
xlabel('V_R (km/s)')
ylabel('V_T (km/s)')
axis equal

subplot(2,2,2)
hold on
plot(vx1, vz1)
plot(vx, vz)
plot(vx2, vz2)
xlabel('V_R (km/s)')
ylabel('V_N (km/s)')
axis equal

subplot(2,2,3)
hold on
plot(vy1, vz1)
plot(vy, vz)
plot(vy2, vz2)
xlabel('V_T (km/s)')
ylabel('V_N (km/s)')
axis equal
hold on

subplot(2,2,4)
hold on
plot3(vx1, vy1, vz1)
plot3(vx, vy, vz)
plot3(vx2, vy2, vz2)
xlabel('V_R (km/s)')
ylabel('V_T (km/s)')
zlabel('V_N (km/s)')
axis equal
view(3)   
axis vis3d
legend('Nonlinear', 'YA', 'QNS', 'Location','best')