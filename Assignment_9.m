% Assignment 9 - Graydon Selanders
%% Question 1 - Plot orbit vs Time
mu = 398600.4418;
h = 69300; % km^2/s
e = 0.74;
RAAN = deg2rad(45);
inc = deg2rad(63.4); 
omega = deg2rad(270); % argument of perigee
theta = 0; % true anomaly

r_x = h^2/mu * 1/(1 + e*cos(theta)) * [cos(theta); sin(theta); 0];
v_x = mu/h * [-sin(theta); e + cos(theta); 0];

Q = [ cos(omega), sin(omega),    0; ... 
     -sin(omega), cos(omega),    0; ... 
         0,           0,         1] ... 
    * [  1,           0,         0; ... 
         0,        cos(inc),   sin(inc); ... 
         0,       -sin(inc),   cos(inc)] ... 
    *[cos(RAAN),  sin(RAAN),     0; ...
     -sin(RAAN),  cos(RAAN),     0; ... 
         0,           0,         1];

% Q eqn from slides is eci to pqw, so use transpose
r_X_ECI = Q.' * r_x;
v_X_ECI = Q.' * v_x;

% plot position vs time for one orbit
a = h^2/(mu*(1 - e^2));
T = 2*pi*sqrt(a^3/mu);
n = sqrt(mu/a^3);

n_pts = 720;
theta = deg2rad( linspace(0, 360 - 360/n_pts, n_pts));

r_x = h^2/mu .* (1./(1 + e.*cos(theta))) .* [cos(theta); sin(theta); zeros(1,n_pts)];
v_x = mu/h .* [-sin(theta); e + cos(theta); zeros(1, n_pts)];

r_ECI = Q.' * r_x;
v_ECI = Q.' * v_x;

E = 2*atan2( sqrt(1 - e).*sin(theta/2), ...
             sqrt(1 + e).*cos(theta/2) );
E = unwrap(E); % force 0..2*pi without jumps
M = E - e.*sin(E); % mean anomaly 
t = M/n; % time since periapsis 

% plot: xyz vs time
figure();
plot(t, r_ECI(1,:), 'LineWidth',1.3); hold on
plot(t, r_ECI(2,:), 'LineWidth',1.3);
plot(t, r_ECI(3,:), 'LineWidth',1.3); grid on; box on
xlabel('Time since periapsis [s]'); ylabel('ECI position [km]')
title(sprintf('ECI position vs time (a = %.0f km, e = %.2f, T = %.2f h)', a, e, T/3600))
legend('x_{ECI}','y_{ECI}','z_{ECI}','Location','best')


%% Question 2 - Calculate direction of the orientation vectors
%  Question 3 - Determine Yaw-Pitch-Roll angles for satellite body frame 
%               relative to the ECI frame at thetas from above

% use values from Q1
disp('---------- Question 2 & 3 ----------')
thetas = [0 90 180 270];
for th = thetas
    tx = [-sind(th); cosd(th); 0];
    rhat = [cosd(th); sind(th); 0];
    hhat = [0;0;1];
    ox = Q.' * tx; % v_par
    oy = -Q.' * hhat; % -ĥ
    oz = -Q.' * rhat; % nadir
    fprintf(['θ=%3d°  ox=[%.6f %.6f %.6f]\n' ...
             '        oy=[%.6f %.6f %.6f]\n' ...
             '        oz=[%.6f %.6f %.6f]\n'], ...
        th, ox, oy, oz);

    DCM = [ox.'; oy.'; oz.'];  
    alpha = atan2(DCM(1,2) , DCM(1,1)); % yaw
    beta  = asin(-DCM(1,3)); % pitch
    gamma = atan2(DCM(2,3) , DCM(3,3)); % roll

    fprintf('θ=%3d°  α=%.1f°  β=%.1f°  γ=%.1f°\n\n', ...
        th, rad2deg(alpha), rad2deg(beta), mod(rad2deg(gamma),360));
end

%% Question 4 - Sun-earth and sun-satellite characteristics
disp('---------- Question 4 ----------')
% A in word doc

% B sun-earth vector
% Julien date (from assignment 8)
Y = 2025; % year
Mon = 11; % month
D = 5; % day
H = 17.5; % hour of the day

JD = 367 * Y - floor(7*((Y + floor((Mon + 9) / 12)) / 4)) + floor((275 * Mon / 9)) + ...
    D + 1721013.5 + H/24;
n_j2000 = JD - 2451545.0;
epsilon = 23.439 - 3.56e-7 * n_j2000;

L = mod(280.459 + 0.98564736 * n_j2000, 360);
M = mod(357.529 + 0.98560023 * n_j2000, 360);

lambda = mod(L + 1.915*sind(M) + 0.02*sind(2*M), 360);

R_ES = (1.00014 - 0.01671*cosd(M) - 0.00014*cosd(2*M)) * 149597870.691;

r_es = R_ES .* [      cosd(lambda); ...
                sind(lambda)*cosd(epsilon); ...
                sind(lambda)*sind(epsilon)];
r_se = -r_es;
fprintf("Sun-Earth vector: [%.4d\n" + ...
        "                   %.4d\n" + ...
        "                   %.4d]\n\n", r_se);

% C - Sun-satellite vector
h = 52740; % km^2/s
e = 0;
RAAN = deg2rad(225);
inc = deg2rad(25); 
omega = deg2rad(0); % argument of perigee
theta = deg2rad(33); % true anomaly

r_x_es = h^2/mu * 1/(1 + e*cos(theta)) * [cos(theta); sin(theta); 0];

Q = [ cos(omega), sin(omega),    0; ... 
     -sin(omega), cos(omega),    0; ... 
         0,           0,         1] ... 
    * [  1,           0,         0; ... 
         0,        cos(inc),   sin(inc); ... 
         0,       -sin(inc),   cos(inc)] ... 
    *[cos(RAAN),  sin(RAAN),     0; ...
     -sin(RAAN),  cos(RAAN),     0; ... 
         0,           0,         1];

% Q eqn from slides is eci to pqw, so use transpose
r_ECI = Q.' * r_x_es;

% if r_sun-sat - r_eci = r_sun-earth, then r_sun-sat = r_sun-earth + r_eci
r_ss = r_se + r_ECI;
fprintf("Sun-Satellite vector: [%.4d\n" + ...
        "                       %.4d\n" + ...
        "                       %.4d]\n\n", r_ss);

%% Question 5 - Characterize sun-sat vector
disp('---------- Question 5 ----------')

v_x_es = mu/h * [-sin(theta); e + cos(theta); 0];

v_ECI = Q.' * v_x_es;

r_ss = r_se + r_ECI;
r_ss_hat = r_ss / norm(r_ss);

hvec = cross(r_ECI, v_ECI);
hhat = hvec/norm(hvec);

oy = -hhat;
oz = -r_ECI/norm(r_ECI); 
ox = cross(oy, oz);

DCM = [ox.'; oy.'; oz.']; 

s_orbit   = DCM * r_ss_hat;
s_xy  = [s_orbit(1); s_orbit(2); 0];
s_xy_hat = s_xy / norm(s_xy);

alpha = atan2d(s_xy(2), s_xy(1));
if alpha < 0
    alpha = alpha + 360; 
end
beta  = atan2d(s_xy(1), s_xy(2));
if beta  < 0
    beta  = beta  + 360; 
end

% Quadrant (A/B/C/D) as in the clarification figure
if     (s_xy(1) > 0) && (s_xy(2) > 0)
    quadrant = 'A  (x>0, y>0)';
elseif (s_xy(1) < 0) && (s_xy(2) > 0)
    quadrant = 'B  (x<0, y>0)';
elseif (s_xy(1) < 0) && (s_xy(2) < 0)
    quadrant = 'C  (x<0, y<0)';
else
    quadrant = 'D  (x>0, y<0)';
end

fprintf('Sun-satellite vector in orbit frame: \n');
fprintf('[%+.4f  %+.4f  %+.4f]\n\n', ...
        s_orbit(1), s_orbit(2), s_orbit(3));

fprintf("Projection onto body X–Y Plane (the sun-sensor Plane) and it's magnitude: \n");
fprintf('[%+.4f  %+.4f  %+.4f]\n\n', ...
        s_xy(1), s_xy(2), s_xy(3));
fprintf('%.4f\n\n', norm(s_xy));

fprintf('Unit direction in X–Y plane: \n');
fprintf('%+.4f  %+.4f  0]\n\n', ...
        s_xy_hat(1), s_xy_hat(2));

fprintf('Quadrant %s\n\n', quadrant);

fprintf('Angles in X–Y Plane:\n');
fprintf('From +X axis (CCW, 0–360°):   %.3f°\n', alpha);
fprintf('From +Y axis (CCW, 0–360°):   %.3f°\n', beta);



