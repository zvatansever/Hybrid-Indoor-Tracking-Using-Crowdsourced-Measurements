function [sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra)
% clear all;close all;clc;
% set(0,'DefaultFigureWindowStyle','docked')
% print_figures=0;
% aviobj = VideoWriter('fit_pro','MPEG-4');
% open(aviobj);

% % load power_map_decimeters P_floor1 P_floor2 P_floor3 P_floor4
%load power_map_decimeters P_floor1
%P_floor2 P_floor3 P_floor4
%load kal_filt_tri_result kal_res
%load rand_paths tra

[m n v]=size(kal_res);


rng(0,'twister');
r=sqrt(vlc_noise); %vlc_noise
% figure
% surf(P_floor1)
% pause
for k=1:m
    P_floor=P_floor1;
    w=r*randn(50,50);
    P_floor=P_floor+w;
    P_floor_all(k,:,:)=P_floor;
    
    
%     if print_figures
%         figure(1)
%         drawnow
%         
%         surf(P_floor)
%         title(['Trajectory #= ' num2str(k) ', Power Dist. '])
%     end
    trajec=tra(k,:,:);
    
    trajectory=reshape(trajec, [n, v]);
    
    traj1=round(trajectory(1,:));
    traj1(traj1<=0)=1;
    traj1(traj1>=50)=50;
    
    traj2=round(trajectory(2,:));
    traj2(traj2<=0)=1;
    traj2(traj2>=50)=50;
    
    kalman_est=kal_res(k,:,:);
    
    kalman_est=reshape(kalman_est, [n, v]);
    
    kal_est1=round(kalman_est(1,:));
    kal_est1(kal_est1<=0)=1;
    kal_est1(kal_est1>=50)=50;
    
    kal_est2=round(kalman_est(2,:));
    kal_est2(kal_est2<=0)=1;
    kal_est2(kal_est2>=50)=50;
    
    %% True results
    
    % Measurement 1
    for j=1:length(traj1)
        %      New1(y_dir(j),x_dir(j))=abs(P_floor1(y_dir(j),x_dir(j)));
        Newer1_true(j)=P_floor(traj1(j),traj2(j));
        % Newer1_noisy(1,j)=D_T_R1_xy_noisy(x_dir(j),y_dir(j));
    end
    
    
    [a b]=size(P_floor1);
    
    sparse1_true=zeros(a,b);
    
    
    for i=1:length(traj1)
        %for j=1:length(Ytri)
        sparse1_true(traj1(i),traj2(i))=Newer1_true(i);
        %end
    end
    
    
    
    obs_sparse1_true=sparse1_true(:);
    true_power=P_floor1(:);
    %
    
    %% Kalman estimates results
    
    % Measurement 1
    for j=1:length(kal_est1)
        %      New1(y_dir(j),x_dir(j))=abs(P_floor1(y_dir(j),x_dir(j)));
        Newer1(j)=P_floor(kal_est1(j),kal_est2(j));
        % Newer1_noisy(1,j)=D_T_R1_xy_noisy(x_dir(j),y_dir(j));
    end
    
    
    
    
    [a b]=size(P_floor);
    
    sparse1=zeros(a,b);
    
    
    for i=1:length(kal_est1)
        %for j=1:length(Ytri)
        sparse1(kal_est1(i),kal_est2(i))=Newer1(i);
        %end
    end
    
    
    %     if print_figures
    %
    %         %         figure
    %         %         surf(P_floor1)
    %         %         title('True Power')
    %         figure(2)
    %         drawnow
    %         subplot(121)
    %         surf(sparse1_true)
    %         title(['Trajectory #= ' num2str(k) ', Sparse 1 True '])
    %         subplot(122)
    %         surf(sparse1)
    %         title(['Trajectory #= ' num2str(k) ', Sparse 1'])
    %         %         frame=getframe(figure(2));
    %         %         writeVideo(aviobj,frame);
    %
    %
    %     end
    
    obs_sparse1=sparse1(:);
    true_power=P_floor1(:);
    
    sparse_kf_est(k,:,:)=sparse1;
end
end
%close(aviobj)
%save sparse_kf sparse_kf_est P_floor vlc_noise P_floor_all



