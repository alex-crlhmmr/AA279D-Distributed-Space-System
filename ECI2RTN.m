function [rRTN, vRTN] = ECI2RTN(r_ref, v_ref, r_target, v_target)
% Inputs:
%   r_ref   - 3x1 position vector (ECI) for the reference (origin of the RTN frame)
%   v_ref   - 3x1 velocity vector (ECI) for the reference
%   r_target- 3x1 position vector (ECI) for the target point
%   v_target- 3x1 velocity vector (ECI) for the target point

% Outputs:
%   rRTN    - 3x1 relative position vector of the target in the RTN frame
%   vRTN    - 3x1 relative velocity vector of the target in the RTN frame

% Ensure the input vectors are columns.
if isrow(r_ref), r_ref = r_ref'; end
if isrow(v_ref), v_ref = v_ref'; end
if isrow(r_target), r_target = r_target'; end
if isrow(v_target), v_target = v_target'; end

% Compute the relative state differences in the ECI frame.
delta_r = r_target - r_ref;
delta_v = v_target - v_ref;

% Construct the RTN basis using the reference state:
% 1. Radial direction (R_hat): from the central body to the reference position.
R_hat = r_ref / norm(r_ref);

% 2. Normal direction (N_hat): perpendicular to the orbital plane.
N = cross(r_ref, v_ref);
N_hat = N / norm(N);

% 3. Transverse direction (T_hat): completes the right-handed system.
T_hat = cross(N_hat, R_hat);

% Form the rotation matrix.
R_matrix = [R_hat, T_hat, N_hat]';  % 3x3 matrix

% Transform the relative vectors from ECI to RTN using the rotation matrix.
rRTN = R_matrix * delta_r;
vRTN = R_matrix * delta_v;

end
