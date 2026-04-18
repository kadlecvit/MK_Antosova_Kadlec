clc
clear

syms R u v

steps = 50;
format long g

% Gnomonic projection
x = R * tan(pi/2 - u) * cos(v);
y = R * tan(pi/2 - u) * sin(v);

% Partial derivatives
fu = diff(x, u);
fv = diff(x, v);
gu = diff(y, u);
gv = diff(y, v);

% Simplify
fu = simplify(fu, 'Steps', steps);
fv = simplify(fv, 'Steps', steps);
gu = simplify(gu, 'Steps', steps);
gv = simplify(gv, 'Steps', steps);

% Distortions
mp2 = (fu^2 + gu^2)/(R^2);
mr2 = (fv^2 + gv^2)/(R*cos(u))^2;
p = 2*(fu*fv + gu*gv)/(R*R*cos(u));

% Simplifications
mp2 = simplify(mp2, 'Steps', steps);
mr2 = simplify(mr2, 'Steps', steps);
p   = simplify(p, 'Steps', steps);

% Angle between projected meridian and parallel
num = gu*fv - fu*gv;
num = simplify(num, 'Steps', steps);
den = fu*fv + gu*gv;
den = simplify(den, 'Steps', steps);
fr = num/den;
fr = simplify(fr, 'Steps', steps);
omega = atan(fr);
omega = simplify(omega, 'Steps', steps);

% Area scale
P = num / (R^2*cos(u));
P = simplify(P, 'Steps', steps);

% Convergence
f = gu / fu;
f = simplify(f, "Steps", steps);
sigma = atan(f);
sigma = simplify(sigma, "Steps", steps);
conv = sigma - pi/2;

% Extreme azimuths
f = p / (mp2 - mr2);
f = simplify(f, "Steps", steps);
A = 0.5 * atan(f);

% Numerical computations (Point Q coordinates from v1)
un = 52.6226 * pi/180;
vn = 36 * pi/180;
Rn = 1;

% Coordinates
xn = double(subs(x, {u, v, R}, {un, vn, Rn}));
yn = double(subs(y, {u, v, R}, {un, vn, Rn}));

% Local linear scales
mpn = sqrt(double(subs(mp2, {u, v, R}, {un, vn, Rn})));
mrn = sqrt(double(subs(mr2, {u, v, R}, {un, vn, Rn})));
pn = double(subs(p, {u, v, R}, {un, vn, Rn}));

% Angle between meridian and parallel
omegan = double(subs(omega, {u, v, R}, {un, vn, Rn}));

% Area scale
Pn = double(subs(P, {u, v, R}, {un, vn, Rn}));

% Convergence
sigman = double(subs(sigma, {u, v, R}, {un, vn, Rn}));
convn = double(subs(conv, {u, v, R}, {un, vn, Rn}));

% Maximum angular distortion
mad = 2 * asin(abs(mpn - mrn) / (mpn + mrn));

% Extreme azimuths
A1n = double(subs(A, {u, v, R}, {un, vn, Rn}));
A2n = A1n + pi/2;

% Semi-axes of Tissot's indicatrix
an = sqrt(abs(mpn^2 * cos(A1n)^2 + mrn^2 * sin(A1n)^2 + pn * sin(A1n) * cos(A1n)));
bn = sqrt(abs(mpn^2 * cos(A2n)^2 + mrn^2 * sin(A2n)^2 + pn * sin(A2n) * cos(A2n)));

% DRAW
% Tissot's indicatrix
scale = 0.1;
rot = sigman - A1n;
t = 0:0.05:2*pi;
xe = scale * an * cos(t);
ye = scale * bn * sin(t);
xt = cos(rot) * xe - sin(rot) * ye + xn;
yt = sin(rot) * xe + cos(rot) * ye + yn;

% Graticule - normal aspect
umin = 45 * pi/180;
umax = 90 * pi/180;
vmin = -180 * pi/180;
vmax = 180 * pi/180;
Du = 10 * pi/180; 
Dv = Du;
du = pi/180;    
dv = du;

[XM, YM, XP, YP] = graticule(umin, umax, vmin, vmax, Du, Dv, du, dv, Rn, 90*pi/180, 0, 0, @gnom);

% Plotting
figure('Name', 'Tissot''s indicatrix');
hold on; axis equal; grid on;

plot(XM', YM', 'k', 'LineWidth', 0.5); 
plot(XP', YP', 'k', 'LineWidth', 0.5);
plot(xt, yt, 'r', 'LineWidth', 1.5);  

xlim([-2.5 2.5]); ylim([-2.5 2.5]);
xlabel('x'); ylabel('y');

set(gcf, 'Renderer', 'painters');