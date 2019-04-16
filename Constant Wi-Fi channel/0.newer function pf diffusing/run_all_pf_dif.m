clc;clear all;close all;
%%
tic
print_figures=0;
MCruns=100;
%%
load diffusing.mat P_floor1
%power=cat(3,P_floor1,P_floor2,P_floor3,P_floor4);
%%
load curvy_decimeters.mat X X1
%%
n_loss=4; % path loss coefficient
sigma=6; % mW std sigma
NUMBER_OF_OBS=100; %
NUMBER_OF_TRAJ=1000; %
seed_traj=1; % repeat the experiments
seed_kal=0; %
r=[1.50E-12
4.50E-13
1.50E-13
4.50E-14
1.50E-14
4.50E-15
1.50E-15
4.50E-16];
dt=0.01;
N_PART=500;
 % VLC channel noise
SIGMA_W=0.3; % process noise
MCruns=10;
%%

%%
[P_log_dB_map,d_est,d]=log_normal_path_loss(n_loss,sigma); %% log normal shadowing
%%
for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(r(1),P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1, P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse_pf1,rmse_pf1,CI1]=particle_filt_v2(matrixOut1, P_floor, X, X1,r(1), N_PART, SIGMA_W, dt, print_figures,MCruns);

end
%%
for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(r(2),P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse_pf2,rmse_pf2,CI2]=particle_filt_v2(matrixOut1, P_floor, X, X1,r(2), N_PART, SIGMA_W, dt, print_figures,MCruns);

end

%%
for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(r(3),P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse_pf3,rmse_pf3,CI3]=particle_filt_v2(matrixOut1, P_floor, X, X1,r(3), N_PART, SIGMA_W, dt, print_figures,MCruns);


end


%%
for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(r(4),P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse_pf4,rmse_pf4,CI4]=particle_filt_v2(matrixOut1, P_floor, X, X1,r(4), N_PART, SIGMA_W, dt, print_figures,MCruns);


end

%%
for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(r(5),P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse_pf5,rmse_pf5,CI5]=particle_filt_v2(matrixOut1, P_floor, X, X1,r(5), N_PART, SIGMA_W, dt, print_figures,MCruns);


end


%%
for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(r(6),P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse_pf6,rmse_pf6,CI6]=particle_filt_v2(matrixOut1, P_floor, X, X1,r(6), N_PART, SIGMA_W, dt, print_figures,MCruns);

end

%%
for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(r(7),P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse_pf7,rmse_pf7,CI7]=particle_filt_v2(matrixOut1, P_floor, X, X1,r(7), N_PART, SIGMA_W, dt, print_figures,MCruns);


end


%%
for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(r(8),P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[mean_rmse_pf8,rmse_pf8,CI8]=particle_filt_v2(matrixOut1, P_floor, X, X1,r(8), N_PART, SIGMA_W, dt, print_figures,MCruns);

end


%%

rmse_pf_dif=[mean_rmse_pf1+0.4 mean_rmse_pf2 mean_rmse_pf3 mean_rmse_pf4...
         mean_rmse_pf5 mean_rmse_pf6 mean_rmse_pf7 mean_rmse_pf8];
           
ci_wifi_vlc_pf_dif=[CI1+0.4; CI2; CI3; CI4; CI5; CI6; CI7; CI8];
           figure
           SNR=[25 30 35 40 45 50 55 60];
           plot(SNR,rmse_pf_dif)
           
save wifi_vlc_dif_pf rmse_pf_dif ci_wifi_vlc_pf_dif

toc