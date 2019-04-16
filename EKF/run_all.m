clc;clear all;close all;
%%
tic
%%
load power_map_decimeters P_floor1
%%
load curvy_36x36_decimeters.mat X X1
%%
n_loss=4; % path loss coefficient
sigma=6; % mW std sigma
NUMBER_OF_OBS=100; %
NUMBER_OF_TRAJ=1000; %
seed_traj=100; % repeat the experiments
seed_kal=0; %
vlc_noise=2.50e-15; % VLC channel noise
qx=0.3; % process noise
%%
ekf_rmse1=[];
ekf_rmse_round1=[];
%%
[P_log_dB_map,d_est,d]=log_normal_path_loss(n_loss,sigma); %% log normal shadowing
%%
for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1, P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[XE,EKF_RMSE1,EKF_RMSE_ROUND1]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);
%%
ekf_rmse1=[ekf_rmse1 EKF_RMSE1];
ekf_rmse_round1=[ekf_rmse_round1 EKF_RMSE_ROUND1];

end
%%
fprintf('RMSE 25 dB=%.3f\n',mean(ekf_rmse1))
fprintf('RMSE round 25 dB=%.3f\n',mean(ekf_rmse_round1))


confi_int_25_db=ekf_rmse_round1;
SEM = std(confi_int_25_db)/sqrt(length(confi_int_25_db));               % Standard Error
ts = tinv([0.05  0.95],length(confi_int_25_db)-1);      % T-Score
CI_25 = mean(confi_int_25_db) + ts*SEM;

fprintf('CI round 25 dB=[%.3f %.3f]\n',CI_25)


%fprintf('RMSE true & RMSE round=[%.3f %.3f]\n',[mean(ekf_rmse) mean(ekf_rmse_round)])
%%
vlc_noise=7.50e-16; % VLC channel noise
ekf_rmse2=[];
ekf_rmse_round2=[];

for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[XE,EKF_RMSE2,EKF_RMSE_ROUND2]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);
%%
ekf_rmse2=[ekf_rmse2 EKF_RMSE2];
ekf_rmse_round2=[ekf_rmse_round2 EKF_RMSE_ROUND2];
end
fprintf('RMSE 30 dB=%.3f\n',mean(ekf_rmse2))
fprintf('RMSE round 30 dB=%.3f\n',mean(ekf_rmse_round2))
confi_int_30_db=ekf_rmse_round2;
SEM = std(confi_int_30_db)/sqrt(length(confi_int_30_db));               % Standard Error
ts = tinv([0.05  0.95],length(confi_int_30_db)-1);      % T-Score
CI_30 = mean(confi_int_30_db) + ts*SEM;

fprintf('CI round 30 dB=[%.3f %.3f]\n',CI_30)
%%
vlc_noise=2.50e-16; % VLC channel noise
ekf_rmse3=[];
ekf_rmse_round3=[];

for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[XE,EKF_RMSE3,EKF_RMSE_ROUND3]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);
%%
ekf_rmse3=[ekf_rmse3 EKF_RMSE3];
ekf_rmse_round3=[ekf_rmse_round3 EKF_RMSE_ROUND3];
end
fprintf('RMSE 35 dB=%.3f\n',mean(ekf_rmse3))
fprintf('RMSE round 35 dB=%.3f\n',mean(ekf_rmse_round3))

confi_int_35_db=ekf_rmse_round3;
SEM = std(confi_int_35_db)/sqrt(length(confi_int_35_db));               % Standard Error
ts = tinv([0.05  0.95],length(confi_int_35_db)-1);      % T-Score
CI_35 = mean(confi_int_35_db) + ts*SEM;

fprintf('CI round 35 dB=[%.3f %.3f]\n',CI_35)

%%
vlc_noise=7.50e-17; % VLC channel noise
ekf_rmse4=[];
ekf_rmse_round4=[];

for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[XE,EKF_RMSE4,EKF_RMSE_ROUND4]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);
%%
ekf_rmse4=[ekf_rmse4 EKF_RMSE4];
ekf_rmse_round4=[ekf_rmse_round4 EKF_RMSE_ROUND4];
end
fprintf('RMSE 40 dB=%.3f\n',mean(ekf_rmse4))
fprintf('RMSE round 40 dB=%.3f\n',mean(ekf_rmse_round4))

confi_int_40_db=ekf_rmse_round4;
SEM = std(confi_int_40_db)/sqrt(length(confi_int_40_db));               % Standard Error
ts = tinv([0.05  0.95],length(confi_int_40_db)-1);      % T-Score
CI_40 = mean(confi_int_40_db) + ts*SEM;

fprintf('CI round 40 dB=[%.3f %.3f]\n',CI_40)
%%

vlc_noise=2.50e-17; % VLC channel noise
ekf_rmse5=[];
ekf_rmse_round5=[];

for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[XE,EKF_RMSE5,EKF_RMSE_ROUND5]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);
%%
ekf_rmse5=[ekf_rmse5 EKF_RMSE5];
ekf_rmse_round5=[ekf_rmse_round5 EKF_RMSE_ROUND5];
end
fprintf('RMSE 45 dB=%.3f\n',mean(ekf_rmse5))
fprintf('RMSE round 45 dB=%.3f\n',mean(ekf_rmse_round5))

confi_int_45_db=ekf_rmse_round5;
SEM = std(confi_int_45_db)/sqrt(length(confi_int_45_db));               % Standard Error
ts = tinv([0.05  0.95],length(confi_int_45_db)-1);      % T-Score
CI_45 = mean(confi_int_45_db) + ts*SEM;

fprintf('CI round 45 dB=[%.3f %.3f]\n',CI_45)
%%


vlc_noise=7.50e-18; % VLC channel noise
ekf_rmse6=[];
ekf_rmse_round6=[];

for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[XE,EKF_RMSE6,EKF_RMSE_ROUND6]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);
%%
ekf_rmse6=[ekf_rmse6 EKF_RMSE6];
ekf_rmse_round6=[ekf_rmse_round6 EKF_RMSE_ROUND6];
end
fprintf('RMSE 50 dB=%.3f\n',mean(ekf_rmse6))
fprintf('RMSE round 50 dB=%.3f\n',mean(ekf_rmse_round6))

confi_int_50_db=ekf_rmse_round6;
SEM = std(confi_int_50_db)/sqrt(length(confi_int_50_db));               % Standard Error
ts = tinv([0.05  0.95],length(confi_int_50_db)-1);      % T-Score
CI_50 = mean(confi_int_50_db) + ts*SEM;

fprintf('CI round 50 dB=[%.3f %.3f]\n',CI_50)
%%



vlc_noise=2.50e-18; % VLC channel noise
ekf_rmse7=[];
ekf_rmse_round7=[];

for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[XE,EKF_RMSE7,EKF_RMSE_ROUND7]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);
%%
ekf_rmse7=[ekf_rmse7 EKF_RMSE7];
ekf_rmse_round7=[ekf_rmse_round7 EKF_RMSE_ROUND7];
end
fprintf('RMSE 55 dB=%.3f\n',mean(ekf_rmse7))
fprintf('RMSE round 55 dB=%.3f\n',mean(ekf_rmse_round7))

confi_int_55_db=ekf_rmse_round7;
SEM = std(confi_int_55_db)/sqrt(length(confi_int_55_db));               % Standard Error
ts = tinv([0.05  0.95],length(confi_int_55_db)-1);      % T-Score
CI_55 = mean(confi_int_55_db) + ts*SEM;

fprintf('CI round 55 dB=[%.3f %.3f]\n',CI_55)
%%
vlc_noise=7.50e-19; % VLC channel noise
ekf_rmse8=[];
ekf_rmse_round8=[];

for i=1:seed_traj
    
[tra]=rand_traj(NUMBER_OF_OBS,NUMBER_OF_TRAJ,i); % create n-random trajectories

[tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra); % Wifi trilateration

[kal_res,mean_rmse_kf]=lin_kal_filt(tri_res,tra,seed_kal); % Kalman filter trilateration results

[sparse_kf_est,P_floor_all,P_floor]=fit_process(vlc_noise,P_floor1,kal_res,tra); % Fit process

[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor); % Signal 

[matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res);

[XE,EKF_RMSE8,EKF_RMSE_ROUND8]=ekf_true_and_database(X, X1,vlc_noise,qx,matrixOut1,P_floor);
%%
ekf_rmse8=[ekf_rmse8 EKF_RMSE8];
ekf_rmse_round8=[ekf_rmse_round8 EKF_RMSE_ROUND8];
end
fprintf('RMSE 60dB=%.3f\n',mean(ekf_rmse8))
fprintf('RMSE round 60 dB=%.3f\n',mean(ekf_rmse_round8))

confi_int_60_db=ekf_rmse_round8;
SEM = std(confi_int_60_db)/sqrt(length(confi_int_60_db));               % Standard Error
ts = tinv([0.05  0.95],length(confi_int_60_db)-1);      % T-Score
CI_60 = mean(confi_int_60_db) + ts*SEM;

fprintf('CI round 60 dB=[%.3f %.3f]\n',CI_60)
%%
save wifi_vlc_nondif_ekf EKF_RMSE_ROUND1 EKF_RMSE_ROUND2 EKF_RMSE_ROUND3 EKF_RMSE_ROUND4...
    EKF_RMSE_ROUND5 EKF_RMSE_ROUND6 EKF_RMSE_ROUND7 EKF_RMSE_ROUND8
toc