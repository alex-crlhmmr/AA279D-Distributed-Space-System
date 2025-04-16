function [rECI, vECI] = RTN2ECI(r_chief_ECI, v_chief_ECI, r_deputy_RTN, v_deputy_RTN)
r0 = r_chief_ECI;
v0 = v_chief_ECI;
r1 = r_deputy_RTN;
v1 = v_deputy_RTN;

% CHIEF RTN TO ECI ROTATION MATRIX
R0_hat = r0 / norm(r0);
N0 = cross(r0, v0);
N0_hat = N0 / norm(N0);
T0 = cross(N0_hat, R0_hat);
T0_hat = T0 / norm(T0);
R_ECI2RTN = [R0_hat'; T0_hat'; N0_hat']';

% CALCULATE DEPUTY RTN COORDS
rECI = R_ECI2RTN * (r1-r0);
vECI = R_ECI2RTN * (v1-v0);

end
