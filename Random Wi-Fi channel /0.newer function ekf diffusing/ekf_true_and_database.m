function [mean_rmse,CI]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor)

%%%
%close all;clear all;clc;

%%%%

print_figures=0;
MCruns=100;
EKF_RMSE=[];
EKF_RMSE_ROUND=[];
randn('seed',9999);
%load curvy_36x36_decimeters.mat X X1

x_dir=X1(1,:);
y_dir=X1(2,:);
%load sparse_kf P_floor r
% load power_map_decimeters P_floor

%load lin_comb_res matrixOut1 P_floor r
New1=P_floor;
New2=flip(fliplr(P_floor)); 
New3=fliplr(P_floor); 
New4=flipud(P_floor); 
%% Measurements
% *** I edited this section ***
%% Power Measurements
   
%Measurement 1
for j=1:length(y_dir)
    Newer1(1,j)=abs(New1(x_dir(j),y_dir(j)));
end
%Measurement 2
for j=1:length(y_dir)
    Newer2(1,j)=abs(New2(x_dir(j),y_dir(j)));
end
%Measurement 3
for j=1:length(y_dir)
    Newer3(1,j)=abs(New3(x_dir(j),y_dir(j)));
end
%Measurement 4
for j=1:length(y_dir)
    Newer4(1,j)=abs(New4(x_dir(j),y_dir(j)));
end
mse1st=[];
var1st=[];
bias1st=[];
signaltonoiseratio1=[];

%%
%load lin_comb_res matrixOut1

Map1=matrixOut1;
Map2=flip(fliplr(matrixOut1));
Map3=fliplr(matrixOut1);
Map4=flipud(matrixOut1);

%% Step 1
for MC=1:MCruns
    MC;
    dt =.01; %sampling interval
    N = length(Newer1);%no. of samples in simulation
    %state transition matrix [4x4]
    
    %qx =.3; % std
    % qy = 12; % noise on y position
    F = [0 0 1 0;
        0 0 0 1;
        0 0 0 0;
        0 0 0 0];
    [A,Q] = lti_disc(F,[],diag([qx qx 0 0]),dt); % Discretization for state transition matrix A and noise matrix Q
    
    %% Step 2: Initialize state and covariance
    P=[ 3 0 0 0; %Initial covariance [4x4]
        0 3 0 0;
        0 0 6 0;
        0 0 0 3];
    %add noise to measurements
    %Measurement Noise var.
    r1=vlc_noise; %variance
    
    %r2=200.01e-14;
    R=[r1 0 0 0; 0 r1 0 0;0 0 r1 0;0 0 0 r1]; % variance to kalman filter .^2
    %R=chol(R);
    w=chol(R)*zeros(4,length(Newer1));% measurement noise [4x50]
    %w=whiten(w,2);
    
    yt=[Newer1;Newer2;Newer3;Newer4]+w; % measurements from 4 Leds [4x50]
    %yt(yt<0)=0;
    
    %% Initialize and run EKF for comparison
    xe = zeros(4,N); % allocate memory
    xe(:,1) = [3 3 6 3]; % initial state for ekf
    % P = P0; % initial covariance
    XX_m=[];
    PP_M=[];
    XX_M1=[];
    YY_m=[];
    KK=[]; % Kalman gain
    VV=[];
    Hh=[];
    SS=[];
    PP=[];
    KK=[];
    XE=[xe(:,1) ];
    PPP=[];
    FF=[];
    HHH=[];
    for k=2:N
        % Prediction Step
        x_m =A*xe(:,k-1); %  Predicted state (state transition matrix times posterior state) [4x1]
        XX_m=[XX_m x_m];
        
        F= A*P*A';
        FF=[FF F];
        
        P_m =F + Q; % predicted covariance update part [4x4]
        PP_M=[PP_M P_m];
        
        %%
        if x_m(1)<=1.5
            x_m(1)=2;
        end
        
        if x_m(2)<=1.5
            x_m(2)=2;
        end
        if x_m(1)>=49.50
           %  x_m(1)=XE(1,k-1);
                x_m(1)=49;
        end
        
        if x_m(2)>=49.50
           %  x_m(2)=XE(2,k-1);
                  x_m(2)=49;
        end
        
        x_m1=round(x_m); % We are looking a discrete floor so need to find the discrete index [4x1]
        %% Numerical Jacobian Part
        % y_m is the observations obtained using the predicted positions [4x1]
        
        % *** I edited this section ***
        
        y_m = [abs(Map1(x_m1(1),x_m1(2)));
            abs(Map2(x_m1(1),x_m1(2)));
            abs(Map3(x_m1(1),x_m1(2)));
            abs(Map4(x_m1(1),x_m1(2)))];
        %%x-1
        yxplus1 = [abs(Map1(x_m1(1)+1,x_m1(2)));
            abs(Map2(x_m1(1)+1,x_m1(2)));
            abs(Map3(x_m1(1)+1,x_m1(2)));
            abs(Map4(x_m1(1)+1,x_m1(2)))];
        
        yyplus1 = [abs(Map1(x_m1(1),x_m1(2)+1));
            abs(Map2(x_m1(1),x_m1(2)+1));
            abs(Map3(x_m1(1),x_m1(2)+1));
            abs(Map4(x_m1(1),x_m1(2)+1))];
        
        yxminus1 = [abs(Map1(x_m1(1)-1,x_m1(2)));
            abs(Map2(x_m1(1)-1,x_m1(2)));
            abs(Map3(x_m1(1)-1,x_m1(2)));
            abs(Map4(x_m1(1)-1,x_m1(2)))];
        
        yyminus1 = [abs(Map1(x_m1(1),x_m1(2)-1));
            abs(Map2(x_m1(1),x_m1(2)-1));
            abs(Map3(x_m1(1),x_m1(2)-1));
            abs(Map4(x_m1(1),x_m1(2)-1));];
        
        
        % difference formula=predicted observations- the four directional power distributions
        % dx1=(y_m-yxplus1); %[4x1]
        dx1=(yxplus1(1)-yxminus1(1))/2;
        dx2=(yxplus1(2)-yxminus1(2))/2;
        dx3=(yxplus1(3)-yxminus1(3))/2;
        dx4=(yxplus1(4)-yxminus1(4))/2;
        
        
        dy1=(yyplus1(1)-yyminus1(1))/2;
        dy2=(yyplus1(2)-yyminus1(2))/2;
        dy3=(yyplus1(3)-yyminus1(3))/2;
        dy4=(yyplus1(4)-yyminus1(4))/2;
        %% The jacobian matrix
        H=[[dx1;dx2;dx3;dx4] [dy1;dy2;dy3;dy4] zeros(4,1) zeros(4,1) ]; %[4x4]
        %H=[[dx1; dx2] [dy1; dy2] zeros(2,1) zeros(2,1) ]; %[4x4]
        Hh=[Hh H];
        %% Measurement Update
        
        v=yt(:,k)-y_m; % measurement residual (innovation) [4x1]
        VV=[VV v];
        HH=H*P_m*H';
        HHH=[HHH HH];
        S=HH+(R); % Innovation covariance [4x4]
        SS=[SS S];
        K=(P_m*H')/S; % Gain [4x4]
        KK=[KK K];
        xe(:,k)=(x_m+K*v); %Updated state estimate [4x1]
        %XE=[XE xe(:,k)];
        
        
        %%
        
        P=P_m-(K*S*K'); % Updated state covariance [4x4]
        PPP(:,:,k)=P;
        
        
 rmse_ekf(MC)=sqrt(mean((X(1,:)-xe(1,:)).^2+(X(2,:)-xe(2,:)).^2));

    end
    %
      mean_rmse=mean(rmse_ekf);
    confi_int=rmse_ekf;
    SEM = std(confi_int)/sqrt(length(confi_int));               % Standard Error
    ts = tinv([0.01  0.99],length(confi_int)-1);      % T-Score
    CI = mean(confi_int) + ts*SEM;
    
end

end

%%%
