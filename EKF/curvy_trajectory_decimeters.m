close all;
clear all;
clc;
%load real_size_crop power_map linear_map poly_map
dt =.01; %
% Accelerations for the object
% first
a = zeros(1,1000);
    a(1,400:500)  =-(pi/2/101/dt );
    a(1,600:800)  = (pi/2/201/dt );
%% Straight
 % x = [3;3;4;4];

%% Curvy
 x = [3;3;6;3];
t = 0;
X = [];
T = [];
% for i=1:1500
    for i=1:1000
    F = [0 0 1 0;... % State transition matrix
        0 0  0 1;...
        0 0  0 a(i);...
        0 0 -a(i) 0];
    x = expm(F*dt)*x; %  State transition matrix (using matrix exponential) * the state
    t  = t + dt;
    X = [X x]; % Final state
    T = [T t];
end
figure;
%contour(power_map), hold on;
plot(X(1,:),X(2,:)) % plot the position of the receiver path
% xlim([0 2450])
% ylim([0 1900])
hold on

X1=round(X);
figure
%contour(power_map), hold on;
plot(X1(1,:),X1(2,:)) % round the position and plot the receiver path
% xlim([0 2450])
% ylim([0 1900])

% title('Walk Path')
xlabel('Dimension (pixel)')
ylabel('Dimension (pixel)')
legend('True Trajectory')
 % coefficients = polyfit(X(1,:), X(2,:), 1);

 save curvy_36x36_decimeters.mat
 %save linear_36x36_new.mat
