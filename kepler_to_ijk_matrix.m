function [rs, vs] = kepler_to_ijk_matrix(orb, mu)
% Extract orbital element matrices
num_points = size(orb,1);
as  = orb(:,1);
es  = orb(:,2);
is  = orb(:,3);
Os  = orb(:,4);
ws  = orb(:,5);
nus = orb(:,6);

% initialize rs / vs
rs = zeros(num_points, 3);
vs = zeros(num_points, 3);
Es = [];
ts = [];

for j = 1:num_points
    a = as(j);
    e = es(j);
    i = is(j);
    O = Os(j);
    w = ws(j);
    nu = nus(j);

    % Constants
    E = 2 * atan2( sqrt(1-e)*sin(nu/2), sqrt(1+e)*cos(nu/2) );
    n = sqrt(mu/a^3);
   
    % Rotation matrices
    Rz_Omega = [cos(-O) sin(-O) 0; -sin(-O) cos(-O) 0; 0 0 1];
    Rx_i = [1 0 0; 0 cos(-i) sin(-i); 0 -sin(-i) cos(-i)];
    Rz_w = [cos(-w) sin(-w) 0; -sin(-w) cos(-w) 0; 0 0 1];
    R = Rz_Omega * Rx_i * Rz_w;
    
    % Calculate r and v
    r_pqw = [a*(cos(E)-e); a*sqrt(1-e^2)*sin(E); 0];
    rs(j, :) = (R * r_pqw)';
    v_pqw = ((a*n)/(1-e*cos(E))) * [-sin(E); sqrt(1-e^2)*cos(E); 0];
    vs(j, :) = (R * v_pqw)';

    % Plot
    Es(end+1) = E;
    ts(end+1) = j;
end

end