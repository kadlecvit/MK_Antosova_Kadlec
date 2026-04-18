clc
clear
format long g

% Common parameters
proj = @gnom;
s0 = 0;
R = 6380 * 1000;   
ax_lim = [-0.6e7 0.6e7 -0.6e7 0.6e7];

% Graticule steps
Du = 10 * pi/180;
Dv = Du;
du = 1 * pi/180;
dv = du;
steps = [Du, Dv, du, dv];

% Boundary points
us1 = 52.6226 * pi/180;  
uj1 = -us1;                
us2 = 10.8123 * pi/180;   
uj2 = -us2;                
ua  = 26.5651 * pi/180;    

% Continent files 
conts = {'eur.txt', 'afr.txt', 'afr_mad.txt', 'amer.txt', ...
         'austr.txt', 'antar.txt', 'greenl.txt', 'newzel1.txt', ...
         'newzel2.txt', 'tasm.txt'};

% Figure setup
figure('Name', 'Dodecahedron - 12 faces');

% Face 1 (k6)
% Input parameters
uk = -ua;
vk = 36 * pi/180;

% Graticule parameters
umin = -70 * pi/180;
umax = 30 * pi/180;
vmin = -20 * pi/180;
vmax = 90 * pi/180;

% Boundary points (A, F, G, H, B, A)
ub = [uj1, uj2, us2, uj2, uj1, uj1];
vb = [0, 0, 36, 72, 72, 0] * pi/180;

% Create globe face
subplot(3, 4, 1);
createGlobeFace([umin, umax, vmin, vmax], steps, R, uk, vk, s0, proj, conts, ub, vb)
axis(ax_lim);
title('Face 1 (K6)');

% Face 2 (k1)
% Input parameters
uk = ua;
vk = 0;

% Graticule parameters
umin = -30 * pi/180;
umax = 70 * pi/180;
vmin = -50 * pi/180;
vmax = 50 * pi/180;

% Boundary points (F, G, P, T, O, F)
ub = [uj2, us2, us1, us1, us2, uj2];
vb = [0, 36, 36, 324, 324, 0] * pi/180;

% Create globe face
subplot(3, 4, 2);
createGlobeFace([umin, umax, vmin, vmax], steps, R, uk, vk, s0, proj, conts, ub, vb)
axis(ax_lim);
title('Face 2 (K1)');

% Face 3 (k2)
% Input parameters
uk = ua;
vk = 72 * pi/180;

% Graticule parameters
umin = -30 * pi/180;
umax = 70 * pi/180;
vmin = 30 * pi/180;
vmax = 120 * pi/180;

% Boundary points (H, I, Q, P, G, H)
ub = [uj2, us2, us1, us1, us2, uj2];
vb = [72, 108, 108, 36, 36, 72] * pi/180;

% Create globe face
subplot(3, 4, 3);
createGlobeFace([umin, umax, vmin, vmax], steps, R, uk, vk, s0, proj, conts, ub, vb)
axis(ax_lim);
title('Face 3 (K2)');

% Face 4 (k3)
% Input parameters
uk = ua;
vk = 144 * pi/180;

% Graticule parameters
umin = -30 * pi/180;
umax = 70 * pi/180;
vmin = 90 * pi/180;
vmax = 200 * pi/180;

% Boundary points (J, K, R, Q, I, J)
ub = [uj2, us2, us1, us1, us2, uj2];
vb = [144, 180, 180, 108, 108, 144] * pi/180;

% Create globe face
subplot(3, 4, 4);
createGlobeFace([umin, umax, vmin, vmax], steps, R, uk, vk, s0, proj, conts, ub, vb)
axis(ax_lim);
title('Face 4 (K3)');

% Face 5 (k4)
% Input parameters
uk = ua;
vk = 216 * pi/180;

% Graticule parameters
umin = -30 * pi/180;
umax = 70 * pi/180;
vmin = 160 * pi/180;
vmax = 270 * pi/180;

% Boundary points (L, M, S, R, K, L)
ub = [uj2, us2, us1, us1, us2, uj2];
vb = [216, 252, 252, 180, 180, 216] * pi/180;

% Create globe face
subplot(3, 4, 5);
createGlobeFace([umin, umax, vmin, vmax], steps, R, uk, vk, s0, proj, conts, ub, vb)
axis(ax_lim);
title('Face 5 (K4)');

% Face 6 (k5)
% Input parameters
uk = ua;
vk = 288 * pi/180;

% Graticule parameters
umin = -30 * pi/180;
umax = 70 * pi/180;
vmin = 240 * pi/180;
vmax = 340 * pi/180;

% Boundary points (N, O, T, S, M, N)
ub = [uj2, us2, us1, us1, us2, uj2];
vb = [288, 324, 324, 252, 252, 288] * pi/180;

% Create globe face
subplot(3, 4, 6);
createGlobeFace([umin, umax, vmin, vmax], steps, R, uk, vk, s0, proj, conts, ub, vb)
axis(ax_lim);
title('Face 6 (K5)');

% Face 7 (k7)
% Input parameters
uk = -ua;
vk = 108 * pi/180;

% Graticule parameters
umin = -70 * pi/180;
umax = 30 * pi/180;
vmin = 50 * pi/180;
vmax = 160 * pi/180;

% Boundary points (B, C, J, I, H, B)
ub = [uj1, uj1, uj2, us2, uj2, uj1];
vb = [72, 144, 144, 108, 72, 72] * pi/180;

% Create globe face
subplot(3, 4, 7);
createGlobeFace([umin, umax, vmin, vmax], steps, R, uk, vk, s0, proj, conts, ub, vb)
axis(ax_lim);
title('Face 7 (K7)');

% Face 8 (k8)
% Input parameters
uk = -ua;
vk = 180 * pi/180;

% Graticule parameters
umin = -70 * pi/180;
umax = 30 * pi/180;
vmin = 120 * pi/180;
vmax = 230 * pi/180;

% Boundary points (C, D, L, K, J, C)
ub = [uj1, uj1, uj2, us2, uj2, uj1];
vb = [144, 216, 216, 180, 144, 144] * pi/180;

% Create globe face
subplot(3, 4, 8);
createGlobeFace([umin, umax, vmin, vmax], steps, R, uk, vk, s0, proj, conts, ub, vb)
axis(ax_lim);
title('Face 8 (K8)');

% Face 9 (k9)
% Input parameters
uk = -ua;
vk = 252 * pi/180;

% Graticule parameters
umin = -70 * pi/180;
umax = 30 * pi/180;
vmin = 200 * pi/180;
vmax = 300 * pi/180;

% Boundary points (D, E, N, M, L, D)
ub = [uj1, uj1, uj2, us2, uj2, uj1];
vb = [216, 288, 288, 252, 216, 216] * pi/180;

% Create globe face
subplot(3, 4, 9);
createGlobeFace([umin, umax, vmin, vmax], steps, R, uk, vk, s0, proj, conts, ub, vb)
axis(ax_lim);
title('Face 9 (K9)');

% Face 10 (k10)
% Input parameters
uk = -ua;
vk = 324 * pi/180;

% Graticule parameters
umin = -70 * pi/180;
umax = 30 * pi/180;
vmin = -100 * pi/180;
vmax = 20 * pi/180;
% Boundary points (E, A, F, O, N, E)
ub = [uj1, uj1, uj2, us2, uj2, uj1];
vb = [288, 0, 0, 324, 288, 288] * pi/180;
% Create globe face
subplot(3, 4, 10);
createGlobeFace([umin, umax, vmin, vmax], steps, R, uk, vk, s0, proj, conts, ub, vb)
axis(ax_lim);
title('Face 10 (K10)');

% Face 11 (north cap)
% Input parameters
uk = pi/2;
vk = 0;

% Graticule parameters
umin = 50 * pi/180;
umax = 90 * pi/180;
vmin = -180 * pi/180;
vmax = 180 * pi/180;

% Boundary points (P, Q, R, S, T, P)
ub = [us1, us1, us1, us1, us1, us1];
vb = [36, 108, 180, 252, 324, 36] * pi/180;

% Create globe face
subplot(3, 4, 11);
createGlobeFace([umin, umax, vmin, vmax], steps, R, uk, vk, s0, proj, conts, ub, vb)
axis(ax_lim);
title('Face 11 - North (K11)');

% Face 12 (south cap)
% Input parameters
uk = -pi/2;
vk = 0;

% Graticule parameters
umin = -90 * pi/180;
umax = -50 * pi/180;
vmin = -180 * pi/180;
vmax = 180 * pi/180;

% Boundary points (A, B, C, D, E, A)
ub = [uj1, uj1, uj1, uj1, uj1, uj1];
vb = [0, 72, 144, 216, 288, 0] * pi/180;

% Create globe face
subplot(3, 4, 12);
createGlobeFace([umin, umax, vmin, vmax], steps, R, uk, vk, s0, proj, conts, ub, vb)
axis(ax_lim);
title('Face 12 - South (K12)');

% Export to vector file
set(gcf, 'Renderer', 'painters');
