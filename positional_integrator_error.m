%% Positional Integrator Error
mu = 4.9048695e12;
[r_numeric_RTN, v_numeric_RTN] = compute_errors_rv(ROUT, VOUT, mu);

figure(3)
plot3(r_numeric_RTN(:, 1), r_numeric_RTN(:, 2), r_numeric_RTN(:, 3))
title('Unperturbed Positional Error of Integrator in RTN frame')
xlabel('R (m)')
ylabel('T (m)')
zlabel('N (m)')

figure(4)
plot3(v_numeric_RTN(:, 1), v_numeric_RTN(:, 2), v_numeric_RTN(:, 3))
title('Unperturbed Velocity Error of Integrator in RTN frame')
xlabel('R (m/s)')
ylabel('T (m/s)')
zlabel('N (m/s)')


%% Function
function [r_numeric_RTN, v_numeric_RTN] = compute_errors_rv(r_numeric_ECI, v_numeric_ECI, mu)
% STEP 1: Convert ECI positions/velocities to orbital elements
[as, es, is, Os, ws, nus] = ECI_to_Kepler_matrix(r_numeric_ECI, v_numeric_ECI, mu);
orb = [as, es, is, Os, ws, nus];

% STEP 2: Compute keplerian orbital elements
num_points = size(r_numeric_ECI, 1);
asK  = orb(1,1)*ones(num_points,1);
esK  = orb(1,2)*ones(num_points,1);
isK  = orb(1,3)*ones(num_points,1);
OsK  = orb(1,4)*ones(num_points,1);
wsK  = orb(1,5)*ones(num_points,1);
nusK = atan2(r_numeric_ECI(:,2), r_numeric_ECI(:,1));
keplerian_elements = [asK, esK, isK, OsK, wsK, nusK];

% STEP 3: Convert keplerian orbital elements to ECI positions/velocities
[r_kepler_ECI, v_kepler_ECI] = kepler_to_ijk_matrix(keplerian_elements, mu);

% STEP 4: Convert Numerical points from ECI to RTN
r_numeric_RTN = zeros(num_points, 3);
v_numeric_RTN = zeros(num_points, 3);
for k = 1:num_points
    r_ref = r_kepler_ECI(k, :);
    v_ref = v_kepler_ECI(k, :);
    r_target = r_numeric_ECI(k, :);
    v_target = v_numeric_ECI(k, :);
    [rRTN, vRTN] = ECI2RTN(r_ref, v_ref, r_target, v_target);
    r_numeric_RTN(k, :) = rRTN;
    v_numeric_RTN(k, :) = vRTN;
end

end
