clear all;
close all;
clc;


rE = 6378; % Radius of Earth (km)
% Initial Conditions
mu = 3.986e5;
r0 = [18000,1,3];
v0 = [0,0,sqrt(mu/r0(1))];


s0 = [r0;v0]; % State Vector

tMax = (2*pi)*(sqrt(norm(r0)^3/mu));

% Timespans
timeStep = 1;
timeSpan = 0:1:50000; % (seconds) to end of period

% Output timeSPan and Solutions (sol)
% diffEq used to pass in variables
[~, sol] = ode45(@(t,s)diffEq(t,s,mu), timeSpan, s0); % Need 3 args

% Plot Solutions
figure;

hold on; grid on; grid minor;
% Create Earth
[X, Y, Z] = sphere(); 

surf(X*rE, Y*rE, Z*rE)
plot3(sol(:,1),sol(:,2),sol(:,3),'-b','LineWidth',2) % Plot in 3d


function sdot = diffEq(t,s,mu)

% r was first three elements, v was second three
r = s(1:3);
v = s(4:6);

sdot(1:3,1) = v; %First three elements are velocity
sdot(4:6,1) = (-mu*r)/(norm(r))^3; % Actual two body EOM


end