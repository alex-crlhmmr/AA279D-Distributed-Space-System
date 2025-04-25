alpha0 = [6771; 0.1005; deg2rad(51.64); deg2rad(257); deg2rad(0); deg2rad(30)]; % Chief
alpha1 = [6771; 0.1006; deg2rad(51.69); deg2rad(257.05); deg2rad(0.05); deg2rad(29.95)]; % Deputy

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

Constants = phi_inv * xbar0

% Propogate
step_size = 1*(pi/180);
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

figure(1)
subplot(3,1,1)
hold on
plot(f_vals/(2*pi), x)
ylabel('R (km)')
subplot(3,1,2)
hold on
plot(f_vals/(2*pi), y)
ylabel('T (km)')
subplot(3,1,3)
hold on
plot(f_vals/(2*pi), z)
ylabel('N (km)')
xlabel('Orbital Periods')

figure(2)
subplot(3,1,1)
hold on
plot(f_vals/(2*pi), vx)
ylabel('V_R (km/s)')
subplot(3,1,2)
hold on
plot(f_vals/(2*pi), vy)
ylabel('V_T (km/s)')
subplot(3,1,3)
hold on
plot(f_vals/(2*pi), vz)
ylabel('V_N (km/s)')
xlabel('Orbital Periods')

figure(3)
hold on
subplot(2,2,1)
hold on
plot(x, y)
xlabel('R (km)')
ylabel('T (km)')
subplot(2,2,2)
hold on
plot(x, z)
xlabel('R (km)')
ylabel('N (km)')
subplot(2,2,3)
hold on
plot(y, z)
xlabel('T (km)')
ylabel('N (km)')
subplot(2,2,4)
hold on
plot3(x, y, z)
xlabel('R (km)')
ylabel('T (km)')
zlabel('N (km)')



%% QUASI
[a_0,e_0,i_0,W_0,w_0,f_0] = deal( ...
    6771, ...                      % a₁
    0.1005, ...                    % e₁
    deg2rad(51.64), ...            % i₁
    deg2rad(257), ...              % Ω₁
    deg2rad(0), ...                % ω₁
    deg2rad(30) ...                % f₁
);
[a_1,e_1,i_1,W_1,w_1,f_1] = deal( ...
    6771, ...
    0.1006, ...
    deg2rad(51.69), ...
    deg2rad(257.05), ...
    deg2rad(0.05), ...
    deg2rad(29.95) ...
);


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

quasi_elements = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy]'
alpha0 = [6771; 0.1005; deg2rad(51.64); deg2rad(257); deg2rad(0); deg2rad(30)]; % Chief

% [a, ex, ey, u, ix, iy

[f_vals, states] = linear_ecc_mapping(quasi_elements, alpha0, 15*pi*2, 1*(pi/180));

x = states(:, 1);
y = states(:, 2);
z = states(:, 3);
vx = states(:, 4);
vy = states(:, 5);
vz = states(:, 6);

figure(1)
subplot(3,1,1)
plot(f_vals/(2*pi), x)
ylabel('R (km)')
legend('YA','QNS','Location','best')

subplot(3,1,2)
plot(f_vals/(2*pi), y)
ylabel('T (km)')
subplot(3,1,3)
plot(f_vals/(2*pi), z)
ylabel('N (km)')
xlabel('Orbital Periods')

figure(2)
subplot(3,1,1)
plot(f_vals/(2*pi), vx)
ylabel('V_R (km/s)')
subplot(3,1,2)
plot(f_vals/(2*pi), vy)
ylabel('V_T (km/s)')
subplot(3,1,3)
plot(f_vals/(2*pi), vz)
ylabel('V_N (km/s)')
xlabel('Orbital Periods')
legend('YA','QNS','Location','best')



figure(3)
subplot(2,2,1)
plot(x, y)
xlabel('R (km)')
ylabel('T (km)')
legend('YA','QNS','Location','best')

axis equal
subplot(2,2,2)
plot(x, z)
xlabel('R (km)')
ylabel('N (km)')
axis equal
subplot(2,2,3)
plot(y, z)
xlabel('T (km)')
ylabel('N (km)')
axis equal
subplot(2,2,4)
plot3(x, y, z)
xlabel('R (km)')
ylabel('T (km)')
zlabel('N (km)')
axis equal




