function [rRTN, vRTN] = ECI2RTN(r_chief_ECI, v_chief_ECI, r_deputy_ECI, v_deputy_ECI)
r0 = r_chief_ECI;
v0 = v_chief_ECI;
r1 = r_deputy_ECI;
v1 = v_deputy_ECI;

% ECI TO CHIEF RTN ROTATION MATRIX
R0_hat = r0 / norm(r0);
N0 = cross(r0, v0);
N0_hat = N0 / norm(N0);
T0 = cross(N0_hat, R0_hat);
T0_hat = T0 / norm(T0);
R_ECI2RTN = [R0_hat'; T0_hat'; N0_hat'];

% CALCULATE DEPUTY RTN COORDS
rRTN = R_ECI2RTN * (r1-r0);
vRTN = R_ECI2RTN * (v1-v0);

end
