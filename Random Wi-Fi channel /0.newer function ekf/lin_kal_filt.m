function [kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal)
% clear all;
% close all;
% clc;
% print_figures=0;
% aviobj = VideoWriter('kal','MPEG-4');
% open(aviobj);

%%
% load tri_pos tri_res error_dist1
% load rand_paths tra
[m n v]=size(tri_res);
rng(seed_kal,'twister')

for i=1:m
    
    trajec=tra(i,:,:);
    
    trajectory=reshape(trajec, [n, v]);
    
    traj1=round(trajectory(1,:));
    traj1(traj1<=0)=1;
    traj1(traj1>=50)=50;
    
    traj2=round(trajectory(2,:));
    traj2(traj2<=0)=1;
    traj2(traj2>=50)=50;
    
    tri=tri_res(i,:,:);
    
    tri=reshape(tri, [n, v]);
    
    tri1=round(tri(1,:));
    tri1(tri1<=0)=1;
    tri1(tri1>=50)=50;
    
    tri2=round(tri(2,:));
    tri2(tri2<=0)=1;
    tri2(tri2>=50)=50;
    
    %% Define update equations (Coefficent matrices): A physics based model for where we expect the Quail to be [state transition (state + velocity)] + [input control (acceleration)]
    A = [1 0;0 1] ;
    ss=2;
    C = [1 0;0 1]; %
    os=2;
    
    %% define main variables
    Pos= [tri1(1); tri(1)]; %initized state
    Pos_estimate = Pos;  %x_estimate of initial location estimation
    Meas_noise_mag=0.5;
    Ex = Meas_noise_mag*eye(os); % Ex convert the measurement noise (stdv) into covariance matrix
    Ez =.5*eye(ss);% Ez convert the process noise (stdv) into covariance matrix
    P = Ex; % estimate of initial  position variance (covariance matrix)
    
    
    %% Measurements
    Meas_noise = Meas_noise_mag * randn(2,length(tri1));
    y = C * [tri1;tri2]+ Meas_noise;
    %% Do kalman filtering
    %initize estimation variables
    Pos_loc_estimate = []; %  position estimate
    Pos= [traj1(1); traj2(1)]; % re-initized state
    P_estimate = P;
    P_mag_estimate = [];
    predic_state = [];
    predic_var = [];
    for t = 1:length(traj1)
        % Predict next state of the quail with the last state and predicted motion.
        Pos_estimate = A * Pos_estimate;
        predic_state = [predic_state Pos_estimate] ;
        %predict next covariance
        P = A * P * A' + Ex;
        predic_var = [predic_var P] ;
        % predicted Ninja measurement covariance
        % Kalman Gain
        K = P*C'*inv(C*P*C'+Ez);
        % Update the state estimate.
        Pos_estimate = Pos_estimate + K * (y(:,t) -  C * Pos_estimate);
        
        % update covariance estimation.
        P =  (eye(2)-K*C)*P;
        %Store for plotting
        Pos_loc_estimate = [Pos_loc_estimate Pos_estimate];
        P_mag_estimate = [P_mag_estimate P];
        
    end
    rmse=sqrt(mean(((Pos_loc_estimate(1,:)-traj1).^2+(Pos_loc_estimate(2,:)-traj2).^2)));
    rmse1(i)=rmse;
    
    
%     if print_figures
%         figure(1)
%         drawnow
%         subplot 121
%         plot(Pos_loc_estimate(1,:),Pos_loc_estimate(2,:),':b',tri1,tri2,'--k',y(1,:),y(2,:),'or','linewidth',2,'markersize',8)
%         legend('Estimate','Trilateration','Meas')
%         title(['Trajectory #= ' num2str(i) ', RMSE of Kalman Filtering=' num2str(rmse) ' dm'])
%         subplot 122
%         plot(Pos_loc_estimate(1,:),Pos_loc_estimate(2,:),'-.r',traj1,traj2,':b','linewidth',2,'markersize',8)
%         legend('Estimate','True')
%         title(['Trajectory #= ' num2str(i) ', Estimate vs. True Traj. '])
%         fprintf('RMSE=%f\n',rmse1);
%         
%     end
    %     frame=getframe(figure(1));
    %     writeVideo(aviobj,frame);
    
    kal_res(i,:,:)=Pos_loc_estimate;
end
mean_rmse_kf=mean(rmse1);
end
%fprintf('Mean RMSE=%.3f\n',mean_rmse_kf);
%close(aviobj)
%save kal_filt_tri_result kal_res rmse1 mean_rmse_kf