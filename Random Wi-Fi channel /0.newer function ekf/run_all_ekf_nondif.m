clc;clear all;close all;
%%
tic
%%
load nondiffusing.mat P_floor1
%%
load curvy_decimeters.mat X X1
%%
n_loss=4; % path loss coefficient
sigma=6; % mW std sigma
NUMBER_OF_OBS=100; %
NUMBER_OF_TRAJ=1000; %
seed_traj=10; % repeat the experiments
seed_kal=0; %
vlc_noise=2.50e-15; % VLC channel noise
qx=0.1; % process noise
%%
ekf_rmse1=[];
ekf_rmse_round1=[];
%%
%%
for i=1:seed_traj
    [P_log_dB_map,d_est,d]=log_normal_path_loss(n_loss,sigma,i); %% log normal shadowing

[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1, P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse1,CI1]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);
end
%%
vlc_noise=7.50e-16; % VLC channel noise
ekf_rmse2=[];
ekf_rmse_round2=[];

for i=1:seed_traj
        [P_log_dB_map,d_est,d]=log_normal_path_loss(n_loss,sigma,i); %% log normal shadowing

[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse2,CI2]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);

end
%%
vlc_noise=2.50e-16; % VLC channel noise
ekf_rmse3=[];
ekf_rmse_round3=[];

for i=1:seed_traj
        [P_log_dB_map,d_est,d]=log_normal_path_loss(n_loss,sigma,i); %% log normal shadowing

[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse3,CI3]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);
end

%%
vlc_noise=7.50e-17; % VLC channel noise
ekf_rmse4=[];
ekf_rmse_round4=[];

for i=1:seed_traj
        [P_log_dB_map,d_est,d]=log_normal_path_loss(n_loss,sigma,i); %% log normal shadowing

[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse4,CI4]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);

end
%%

vlc_noise=2.50e-17; % VLC channel noise
ekf_rmse5=[];
ekf_rmse_round5=[];

for i=1:seed_traj
        [P_log_dB_map,d_est,d]=log_normal_path_loss(n_loss,sigma,i); %% log normal shadowing

[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse5,CI5]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);

end

%%
vlc_noise=7.50e-18; % VLC channel noise
ekf_rmse6=[];
ekf_rmse_round6=[];

for i=1:seed_traj
        [P_log_dB_map,d_est,d]=log_normal_path_loss(n_loss,sigma,i); %% log normal shadowing

[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse6,CI6]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);
end

%%
vlc_noise=2.50e-18; % VLC channel noise
ekf_rmse7=[];
ekf_rmse_round7=[];

for i=1:seed_traj
        [P_log_dB_map,d_est,d]=log_normal_path_loss(n_loss,sigma,i); %% log normal shadowing

[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse7,CI7]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);

end
%%
vlc_noise=7.50e-19; % VLC channel noise
ekf_rmse8=[];
ekf_rmse_round8=[];

for i=1:seed_traj
[P_log_dB_map,d_est,d]=log_normal_path_loss(n_loss,sigma,i); %% log normal shadowing

[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse8,CI8]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
random_EKF_RMSE_VLC=[mean_rmse1 mean_rmse2 mean_rmse3 mean_rmse4...   %
            mean_rmse5 mean_rmse6 mean_rmse7 mean_rmse8];      %
random_confidence_ekf_nondif=[CI1; CI2; CI3; CI4; CI5; CI6; CI7; CI8];   %
SNR=[25 30 35 40 45 50 55 60];                                           
figure
plot(SNR,random_EKF_RMSE_VLC)%
save randoom_wifi_vlc_nondif_ekf random_EKF_RMSE_VLC random_confidence_ekf_nondif       %
toc                                                            %