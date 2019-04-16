
clc; clear; close all;
tic
print_figures=1;
MCruns = 1;
set(0,'DefaultFigureWindowStyle','docked')


%% Sample walk path
% load curvy_36x36_decimeters.mat X X1 
randn('seed',9999)
% load walkpath1.mat traj traj1 
%x_dir=X1(1,:);
% x_dir(1)=1;
%y_dir=X1(2,:);
% y_dir(3)=1;
% y_dir(1:2)=1;
% % 
R=sqrt(7.5e-19);

%% basic conditions setting
x=50; % the size of the room
y=50;
z=30;

c=3*10^8;
FOV_LED=pi/2;
semi=60;
P_LED=0.0025;
Adet=0.01;
lambertian= -log(2)/log(cos(pi/180*semi));

%% The position of the LED source
LED_x1=[x/4];
LED_x2= [x-x/4];
LED_x3=[x/4];
LED_x4=[x-x/4];

LED_y1=[y/4];
LED_y2=[y-y/4];
LED_y3=[y-y/4];
LED_y4=[y/4];

LED_z1=[z];
LED_z2=[z];
LED_z3=[z];
LED_z4=[z];
h=z;
LED_pos1=[LED_x1;LED_y1; LED_z1];
LED_pos2=[LED_x2;LED_y2; LED_z2];
LED_pos3=[LED_x3;LED_y3; LED_z3];
LED_pos4=[LED_x4;LED_y4; LED_z4];

norm_rec=[0;0;1];
norm_wall_1=[0,1,0];
norm_wall_2=[0,-1,0];
norm_wall_3=[1,0,0];
norm_wall_4=[-1,0,0];

LED_number=length(LED_x1);
%% This part is to calculate the vectors of LEDs
% spherical coordinates
% theta=linspace(pi/2,pi,N)
angle_phi=10/180*pi;

% theta_3=[pi-angle_phi*2];
% theta_2=[pi-angle_phi];
theta_1=[pi];

% N=17; %
% phi_3=linspace(0,2*pi,N);
% 
% N=9; %
% phi_2=linspace(0+pi/4,2*pi+pi/4,N);

N=2;
phi_1=linspace(0,2*pi,N);
% theta=[pi-pi/4,pi-pi/5,pi-pi/6];
f=1;


for Mc = 1:MCruns  % for MC runs
Mc  % # of MC run
% 
% % Vector of Beam
% for j=1:length(phi_3)-1
%  
%       vector_beam(:,f)=[sin(theta_3)*cos(phi_3(j));
%                                sin(theta_3)*sin(phi_3(j));
%                                cos(theta_3)];
%     f=f+1;
% end
% 
% for j=1:length(phi_2)-1
%     vector_beam(:,f)=[sin(theta_2)*cos(phi_2(j));
%                              sin(theta_2)*sin(phi_2(j));
%                              cos(theta_2)];
%     f=f+1;
% end

for j=1:length(phi_1)-1
    vector_beam(:,f)=[sin(theta_1)*cos(phi_1(j));
                             sin(theta_1)*sin(phi_1(j));
                             cos(theta_1)];
    f=f+1;
end

[a, N_led]=size(vector_beam);

%%
P_wall_1=zeros(x,z);
P_wall_2=zeros(x,z);
P_wall_3=zeros(y,z);
P_wall_4=zeros(y,z);
P_F_1=0;
P_F_2=0;
P_F_3=0;
P_F_4=0;
%% incidence & irradiance angles matrix of floor and four walls
Incid_angle_floor=zeros(x,y);
Irrad_angle_floor=zeros(x,y);
Incid_angle_wall_1=zeros(x,z);
Irrad_angle_wall_1=zeros(x,z);
Incid_angle_wall_2=zeros(x,z);
Irrad_angle_wall_2=zeros(x,z);
Incid_angle_wall_3=zeros(y,z);
Irrad_angle_wall_3=zeros(y,z);
Incid_angle_wall_4=zeros(y,z);
Irrad_angle_wall_4=zeros(y,z);
Distance1=[];
Distance2=[];
Distance3=[];
Distance4=[];
%% calculate the floor incidence angles & irradiance angles
% LED1
 P_floor1=zeros(x,y); %% received power in the room

for j=1:y
    for i=1:x
        % for k=1:LED_number
        for m=1:N_led
       
        vector_ray=[i;j;0]-LED_pos1;
                vector_ray_total(i,j,:)=vector_ray;
                Incid_angle_floor(i,j)=angle(-vector_ray,norm_rec);
              % distance between transmitter and receiver
                d_t_r(i,j)=sqrt((i-LED_x1)^2+(j-LED_y1)^2+LED_z1^2);
                d_t_r_xy(i,j)=sqrt((i-LED_x1)^2+(j-LED_y1)^2);
                 d_z=z;
                c=0;d=.5;  %Zero Mean Normal Dist. Noise
%                     c=0;d=0.785;
%                     c=.5;d=-.5;
%                     c=0.25;d=-0.25;
%                     c=0.125;d=-0.125;
%                     e=c+(d)*randn(1,1);
                    e=1;
                Irrad_angle_floor(i,j)=angle(vector_ray,vector_beam(:,m));
                P1(i,j)=P_LED*Adet*((lambertian+1)*(cos(Irrad_angle_floor(i,j))^lambertian)*cos(Incid_angle_floor(i,j))/((2*pi)*(d_t_r(i,j)/10)^2));
                D_T_R1(i,j)=10*sqrt((P_LED*Adet*(lambertian+1)*(cos(Irrad_angle_floor(i,j))^lambertian)*cos(Incid_angle_floor(i,j))/(2*pi*P1(i,j))));
                D_T_R1_xy(i,j)=10*sqrt((P_LED*Adet*(lambertian+1)*(cos(Irrad_angle_floor(i,j))^lambertian)*cos(Incid_angle_floor(i,j))/(2*pi*P1(i,j)))-((0.1*h)^2));
                D_T_R1_xy_noisy(i,j)=10*sqrt(((P_LED*Adet*(lambertian+1)*(cos(Irrad_angle_floor(i,j))^lambertian)*cos(Incid_angle_floor(i,j))/(2*pi*(P1(i,j)+R*randn)))-((0.1*h)^2)));

                P_floor1(i,j)=P_floor1(i,j)+P1(i,j);
            
            end
          end
end

% LED2
P_floor2= zeros(x,y); %% received power in the room
for j=1:y
    for i=1:x
        %for k=1:LED_number
            for m=1:N_led
                vector_ray=[i;j;0]-LED_pos2;
                vector_ray_total(i,j,:)=vector_ray;
                Incid_angle_floor(i,j)=angle(-vector_ray,norm_rec);
                %Irrad_angle_floor(i,j,m)=angle(vector_ray,vector_beam(:,m))
                
                % distance between transmitter and receiver
                d_t_r(i,j)=sqrt((i-LED_x2)^2+(j-LED_y2)^2+LED_z2^2);
                d_t_r_xy(i,j)=sqrt((i-LED_x1)^2+(j-LED_y1)^2);
                d_z=z;
%                 c=0;d=.5;
%                     c=0;d=0.785;
%                         c=.5;d=-.5;
%                      c=0.25;d=-0.25;
%                     c=0.125;d=-0.125;
%                 e=c+(d-c)*randn(1,1);
           e=1;
                Irrad_angle_floor(i,j)=angle(vector_ray,vector_beam(:,m))*e;
                P2(i,j)=P_LED*Adet*((lambertian+1)/((2*pi)*(d_t_r(i,j)/10)^2))*cos(Irrad_angle_floor(i,j))^lambertian*cos(Incid_angle_floor(i,j));
                D_T_R2(i,j)=10*sqrt((P_LED*Adet*(lambertian+1)*(cos(Irrad_angle_floor(i,j))^lambertian)*cos(Incid_angle_floor(i,j))/(2*pi*P2(i,j))));
                D_T_R2_xy(i,j)=10*sqrt((P_LED*Adet*(lambertian+1)*(cos(Irrad_angle_floor(i,j))^lambertian)*cos(Incid_angle_floor(i,j))/(2*pi*P2(i,j)))-((0.1*h)^2));
                D_T_R2_xy_noisy(i,j)=10*sqrt((P_LED*Adet*(lambertian+1)*(cos(Irrad_angle_floor(i,j))^lambertian)*cos(Incid_angle_floor(i,j))/(2*pi*(P2(i,j)+R*randn)))-((0.1*h)^2));

                P_floor2(i,j)=P_floor2(i,j)+P2(i,j);
            end
        end
    %end
end


% LED3
P_floor3= zeros(x,y);  %% received power in the room
for j=1:y
    for i=1:x
       % for k=1:LED_number
            for m=1:N_led
                vector_ray=[i;j;0]-LED_pos3;
                vector_ray_total(i,j,:)=vector_ray;
                Incid_angle_floor(i,j)=angle(-vector_ray,norm_rec);
%              Irrad_angle_floor(i,j,m)=angle(vector_ray,vector_beam(:,m));
                % distance between transmitter and receiver
                d_t_r(i,j)=sqrt((i-LED_x3)^2+(j-LED_y3)^2+LED_z3^2);
                d_t_r_xy(i,j)=sqrt((i-LED_x1)^2+(j-LED_y1)^2);
                d_z=z;
             
%                     c=0;d=.5;
%                     c=0;d=0.785;
%                         c=.5;d=-.5;
%                      c=0.25;d=-0.25;
%                     c=0.125;d=-0.125;

%                 e=c+(d-c)*randn(1,1);
                e=1;
                Irrad_angle_floor(i,j)=angle(vector_ray,vector_beam(:,m))*e;
                P3(i,j)=P_LED*Adet*((lambertian+1)/((2*pi)*(d_t_r(i,j)/10)^2))*cos(Irrad_angle_floor(i,j))^lambertian*cos(Incid_angle_floor(i,j));
                D_T_R3(i,j)=10*sqrt((P_LED*Adet*(lambertian+1)*(cos(Irrad_angle_floor(i,j))^lambertian)*cos(Incid_angle_floor(i,j))/(2*pi*P3(i,j))));
                D_T_R3_xy(i,j)=10*sqrt((P_LED*Adet*(lambertian+1)*(cos(Irrad_angle_floor(i,j))^lambertian)*cos(Incid_angle_floor(i,j))/(2*pi*P3(i,j)))-((0.1*h)^2));
                D_T_R3_xy_noisy(i,j)=10*sqrt((P_LED*Adet*(lambertian+1)*(cos(Irrad_angle_floor(i,j))^lambertian)*cos(Incid_angle_floor(i,j))/(2*pi*(P3(i,j)+R*randn)))-((0.1*h)^2));

                P_floor3(i,j)=P_floor3(i,j)+P3(i,j);
            end
        end
    %end
end

% LED4
P_floor4= zeros(x,y);  %% received power in the room

for j=1:y
    for i=1:x
       % for k=1:LED_number
            for m=1:N_led
                vector_ray=[i;j;0]-LED_pos4;
                vector_ray_total(i,j,:)=vector_ray;
                Incid_angle_floor(i,j)=angle(-vector_ray,norm_rec);
%             Irrad_angle_floor(i,j,m)=angle(vector_ray,vector_beam(:,m));

                % distance between transmitter and receiver
                d_t_r(i,j)=sqrt((i-LED_x4)^2+(j-LED_y4)^2+LED_z4^2);
                d_t_r_xy(i,j)=sqrt((i-LED_x1)^2+(j-LED_y1)^2);
                d_z=z;
% %                      c=0;d=.5;
%                  c=0;d=0.785;
%                  c=.5;d=-.5;
%                  c=0.25;d=-0.25;
%                  c=0.125;d=-0.125;
%                     e=c+(d-c)*randn(1,1);
                    e=1;
                Irrad_angle_floor(i,j)=angle(vector_ray,vector_beam(:,m))*e;
                P4(i,j)=P_LED*Adet*((lambertian+1)/((2*pi)*(d_t_r(i,j)/10)^2))*cos(Irrad_angle_floor(i,j))^lambertian*cos(Incid_angle_floor(i,j));
                D_T_R4(i,j)=10*sqrt((P_LED*Adet*(lambertian+1)*(cos(Irrad_angle_floor(i,j))^lambertian)*cos(Incid_angle_floor(i,j))/(2*pi*P4(i,j))));
                D_T_R4_xy(i,j)=10*sqrt((P_LED*Adet*(lambertian+1)*(cos(Irrad_angle_floor(i,j))^lambertian)*cos(Incid_angle_floor(i,j))/(2*pi*P4(i,j))-((0.1*h)^2)));
                D_T_R4_xy_noisy(i,j)=10*sqrt((P_LED*Adet*(lambertian+1)*(cos(Irrad_angle_floor(i,j))^lambertian)*cos(Incid_angle_floor(i,j))/(2*pi*(P4(i,j)+R*randn))-((0.1*h)^2)));

                P_floor4(i,j)=P_floor4(i,j)+P4(i,j);
            end
        end
end
end
% 
% %% Measurement 1
% for j=1:length(y_dir)
% %      New1(y_dir(j),x_dir(j))=abs(P_floor1(y_dir(j),x_dir(j)));
%      Newer1(1,j)=D_T_R1_xy(x_dir(j),y_dir(j));
%      Newer1_noisy(1,j)=D_T_R1_xy_noisy(x_dir(j),y_dir(j));
% 
%      
% end
%  
% %% Measurement 2
% for j=1:length(y_dir)
% %     New2(y_dir(j),x_dir(j))=abs(P_floor2(y_dir(j),x_dir(j)));
%     Newer2(1,j)=D_T_R2_xy(x_dir(j),y_dir(j));
%     Newer2_noisy(1,j)=D_T_R2_xy_noisy(x_dir(j),y_dir(j));
% 
% end
% %% Measurement 3
% for j=1:length(y_dir)
% %     New3(y_dir(j),x_dir(j))=abs(P_floor3(y_dir(j),x_dir(j)));
%     Newer3(1,j)=D_T_R3_xy(x_dir(j),y_dir(j));
%     Newer3_noisy(1,j)=D_T_R3_xy_noisy(x_dir(j),y_dir(j));
% 
% end
% 
% %% Measurement 4
% for j=1:length(y_dir)
% %     New4(y_dir(j),x_dir(j))=abs(P_floor4(y_dir(j),x_dir(j)));
%     Newer4(1,j)=D_T_R4_xy(x_dir(j),y_dir(j));
%     Newer4_noisy(1,j)=D_T_R4_xy_noisy(x_dir(j),y_dir(j));
% 
%  end

if print_figures
figure;
surfc(P_floor1)
colorbar;
title('LED1 Power Dist.')

% hold on;
% contour(Newer1,50);
% colorbar;
% contour(New1,50);
% colorbar;
% title('Received Power Dist. from LED1')
% plot(x_dir,y_dir,'linewidth',2)

figure;
surfc(P_floor2)
colorbar;
title('LED2 Power Dist.')
% hold on;
% contour(Newer2,50);
% colorbar;
% contour(New2,50);
% colorbar;
% title('Received Power Dist. from LED2')
% plot(x_dir,y_dir,'linewidth',2)

figure;
surfc(P_floor3)
colorbar;
title('LED3 Power Dist.')
% hold on;
% contour(Newer3,50);
% colorbar;
% contour(New3,50);
% colorbar;
% title('Received Power Dist. from LED3')
% plot(x_dir,y_dir,'linewidth',2)

figure;
surfc(P_floor4)
colorbar;
title('LED4 Power Dist.')
% hold on;
% contour(Newer4,50);
% colorbar;
% contour(New4,50);
% colorbar;
% title('Received Power Dist. from LED4')
% plot(x_dir,y_dir,'linewidth',2)
figure
surfc(P_floor1+P_floor2+P_floor3+P_floor4)
figure
imagesc(P_floor1)
hold on
imagesc(P_floor2)
imagesc(P_floor3)
imagesc(P_floor4)
figure
contour(P_floor1+P_floor2+P_floor3+P_floor4)

end
save power_map_decimeters.mat P_floor1 P_floor2 P_floor3 P_floor4
toc