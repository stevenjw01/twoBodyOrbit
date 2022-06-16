clear;
close all;
clc;

% --------- INITAL CONDITIONS ------------
rE = 6378; % Radius of Earth (km)
% Initial Conditions
mu = 3.986e5;
r0 = [15000,0,0];
v0 = [0,0,sqrt(mu/r0(1))];
%v0 = [0,0,sqrt(mu/r0(1))];
% --------- END INITAL CONDITIONS --------




% ----- CALCULATE ORBITAL ELEMENTS -------
r = norm(r0); % Calc distance
v = norm(v0); % Calc speed
vR = (r0(1)*v0(1) + r0(2)*v0(2) + r0(3)*v0(3))/r; %Radial velocity
hV = cross(r0,v0); % Vector specific angular momentum
h = norm(hV); % Magnitude of specific angular momentum

i = acosd(hV(3)/h); % Inclination
% Node line
K = [0,0,1];
NV = cross(K,hV); % Node line vector
N = norm(NV); % Magnitude of node line


omegaRA = acosd(NV(1)/N); % Right Ascension (RA) of ascending node

% Calculate eccentricity vector
eV = 1/mu*((v^2-mu/r)*r0-r*vR*v0);
e = norm(eV); % Magnitude of eccentricity vector

omega = acosd(dot(NV,eV)/(N*e)); % Argument of perigee

theta = acosd(dot(eV,r0)/(e*r)); % True Anomaly

rp = (h^2/mu)*(1/(1+e*cosd(0))); % Perigee radii
ra = (h^2/mu)*(1/(1+e*cosd(180))); % apogee radii


a = 1/2*(rp+ra); % Semimajor axis (km)
T = ((2*pi)/sqrt(mu))*a^(3/2); % Calculate period
% tMax = (2*pi)*(sqrt(a^3/mu));



% ----- END CALCULATE ORBITAL ELEMENTS -------


s0 = [r0;v0]; % State Vector



% Timespans
timeStep = 1;
timeSpan = 0:0.1:T; % (seconds) to end of period

% Output timeSPan and Solutions (sol)
% diffEq used to pass in variables
[~, sol] = ode45(@(t,s)diffEq(t,s,mu), timeSpan, s0); % Need 3 args

% Plot Solutions
figure;

hold on; grid on; grid minor;
% Create Earth
[X, Y, Z] = sphere(); 

surf(X*rE, Y*rE, Z*rE)
axis equal;
plot3(sol(:,1),sol(:,2),sol(:,3),'-b','LineWidth',2) % Plot in 3d


function sdot = diffEq(t,s,mu)

% r is first three elements, v is last three
rDQ = s(1:3);
vDQ = s(4:6);

sdot(1:3,1) = vDQ; %First three elements are velocity
sdot(4:6,1) = (-mu*rDQ)/(norm(rDQ))^3; % Actual two body EOM


end
