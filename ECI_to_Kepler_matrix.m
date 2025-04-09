function [as, es, is, Os, ws, nus] = ECI_to_Kepler_matrix(rin, vin, mu)
num_points = size(rin, 1);

as = zeros(num_points, 1);
es = zeros(num_points, 1);
is = zeros(num_points, 1);
Os = zeros(num_points, 1);
ws = zeros(num_points, 1);
nus = zeros(num_points, 1);


for k = 1:num_points
    r = rin(k, :);
    v = vin(k, :);
    h = cross(r, v);
    i = atan2(sqrt(h(1)^2 + h(2)^2), h(3));
    O = atan2(h(1), -h(2));
    P = norm(h)^2 / mu;
    a = (2/norm(r) - norm(v)^2 / mu)^(-1);
    n = sqrt(mu/a^3);
    e = sqrt(1 - P/a);
    E = atan2(dot(r, v) / (a^(2)*n), 1-norm(r)/a);
    nu = 2*atan2(sqrt((1+e)/(1-e))*tan(E/2), 1);
    u = atan2(r(3)/sin(i), r(1)*cos(O) + r(2)*sin(O));
    w = u-nu;
    i  = wrapTo2Pi(i);
    O  = wrapTo2Pi(O);
    w  = wrapTo2Pi(w);
    nu = wrapTo2Pi(nu);
    
    as(k) = a;
    es(k) = e;
    is(k) = i;
    Os(k) = O;
    ws(k) = w;
    nus(k) = nu;

end
end