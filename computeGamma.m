function Gamma = computeGamma(Chief_OE, uk_j, mu)
a = Chief_OE(1);
n = sqrt(mu / a^3);
Gamma = 1/(n*a) * [
     0,          2,           0;
    -2,          0,           0;
     sin(uk_j),  2*cos(uk_j), 0;
    -cos(uk_j),  2*sin(uk_j), 0;
     0,          0,           cos(uk_j);
     0,          0,           sin(uk_j)
];
end