function [f_vals, states] = get_YA_states(alpha0, constants, step_size, f_change, n, mu, h) 
e = alpha0(2);
f_initial = alpha0(6);
f_final = f_initial + f_change;
f_vals = f_initial:step_size:f_final;
states = zeros(length(f_vals), 6);

E = 2*atan2(sqrt((1 - e)./(1 + e)).*tan(f_vals./2), 1);
E = unwrap(E);
M = E - e.*sin(E);                % continuous mean anomaly
M0 = M(1);                        % initial mean anomaly
t = (M - M0) ./ n;                % time *since* f0
Is = (mu^2/h^3) * t;              % integral I(f) = (μ²/h³)(t−t0)

for i = 1:length(f_vals)
    f = f_vals(i);
    I = Is(i);
    k  = 1 + e*cos(f);
    s  = k*sin(f);
    c  = k*cos(f);
    ds = k*cos(f) - e*sin(f)^2;    
    dc = -k*sin(f) - e*cos(f)^2;  
    
    phi = [ ...
    s,                  c,                 2 - 3*e*s*I,                0,       0,       0; ...
    ds,                 dc,                -3*e*(ds*I + s/k^2),        0,       0,       0; ...
    c*(1+1/k),         -s*(1+1/k),         -3*k^2*I,                   1,       0,       0; ...
   -2*s,               (e - 2*c),          -3*(1 - 2*e*s*I),           0,       0,       0; ...
    0,                  0,                  0,                         0,  cos(f),  sin(f); ...
    0,                  0,                  0,                         0, -sin(f),  cos(f)  ...
    ];
    
    state_i = phi*constants;
    states(i, :) = state_i';
    
end
end
