clc;clear all;close all;
set(0,'DefaultFigureWindowStyle','normal')

%%
x_sigma=[1 2 4 6 8 10];
trila=[4.925501 6.080094 9.08358 13.28294 17.930062 22.436396];
trila_kal=[3.984 4.686 6.663 9.633 13.188 16.975];
ekf_30_dB=[4.773 5.056 5.103 5.203 5.33 5.506];
ekf_40_dB=[0.749 0.731 0.74 0.726 0.733 0.735];
ekf_50_dB=[0.459 0.458 0.461 0.457 0.459  0.463];
ekf_60_dB=[0.373 0.375 0.373 0.373 0.373  0.378];
figure
set(gca,'fontsize',12)
hold on
plot(x_sigma,trila,'-o',x_sigma,trila_kal,'-*',x_sigma,ekf_30_dB,'-hex',...
    x_sigma,ekf_40_dB,'-diamond',x_sigma,ekf_50_dB,'-penta',x_sigma,ekf_60_dB,'-square','linewidth',2,'markersize',12)
legend('Trilateration only (Wi-Fi)','Kalman filtered trilateration results (Wi-Fi)','Wi-fi aided VLC (Extended Kalman filter, Peak SNR 30 dB)','Wi-fi aided VLC (Extended Kalman filter, Peak SNR 40 dB)','Wi-fi aided VLC (Extended Kalman filter, Peak SNR 50 dB)','Wi-fi aided VLC (Extended Kalman filter, Peak SNR 30 dB)',...
    'location','NorthWest')
xlabel('Standard deviation of Wi-Fi Path Loss, X_{\sigma} (dB)')
ylabel('RMSE (decimeters)')
hold off
%%
PSNR=[25 30 35 40 45 50 55 60];
ekf_w_var=[8.235 4.773 1.145 0.749 0.569 0.459 0.394 0.373];
ekf_w_const=[4.845 1.419 0.948 0.681 0.545 0.444 0.387 0.368];
figure
set(gca,'fontsize',14)
hold on
plot(PSNR,ekf_w_var,'-*',PSNR,ekf_w_const,'-square','linewidth',2,'markersize',12)
legend('EKF, varying weights','EKF, constant weight')
xlabel('PSNR (dB)')
ylabel('RMSE (dm)')
hold off

%%
PSNR=[25 30 35 40 45 50 55 60];

final_fit_error_sigma_1=[9.66E-09 5.22E-09 2.95E-09 1.54E-09 8.21E-10 3.82E-10 1.51E-10 2.86E-11];
final_fit_error_sigma_10=[5.08E-08 2.96E-08 2.05E-08 1.67E-08 1.57E-08 1.56E-08 1.56E-08 1.56E-08];
figure
set(gca,'fontsize',14)
hold on
semilogy(PSNR,final_fit_error_sigma_1,'-hex',PSNR,final_fit_error_sigma_10,'-penta','linewidth',2,'markersize',12)
legend('Fit error, X_{\sigma}=1 dB','Fit error, X_{\sigma}=10 dB')
xlabel('PSNR (dB)')
ylabel('RMSE (W)')
hold off
%%
number_of_traj=[100 200 300 500 700 1000];
snr_30_db_w_vary=[16.51 8.52 8.411 6.698  6.698  4.351 ];
snr_40_db_w_vary=[2.896 0.955 0.723 0.747 0.762 0.749];
snr_50_db_w_vary=[1.259 0.52 0.461 0.458 0.461 0.459];
snr_60_db_w_vary=[1.219 0.466 0.373 0.372 0.372 0.373];
figure
set(gca,'fontsize',14)
hold on
plot(number_of_traj,snr_40_db_w_vary,'-hex',number_of_traj,snr_50_db_w_vary,'-^',...
    number_of_traj,snr_60_db_w_vary,'-o','linewidth',2,'markersize',12)
legend('EKF, PSNR 40 dB','EKF, PSNR 50 dB','EKF, PSNR 60 dB')
xlabel('Number of trajectories')
ylabel('RMSE (dm)')
hold off
%%
number_of_traj=[100 200 300 500 700 1000];
snr_30_db_w_constant=[20.188 5.304 1.405 1.389 1.469 1.419];
snr_40_db_w_constant=[10.055 1.442 0.75 0.671 0.655 0.651];
snr_50_db_w_constant=[3.794 0.989 0.513 0.445 0.445 0.444];
snr_60_db_w_constant=[3.137 0.957 0.481 0.376 0.372 0.368];
figure
set(gca,'fontsize',14)
hold on
plot(number_of_traj,snr_40_db_w_constant,'-hex',number_of_traj,snr_50_db_w_constant,'-^',...
    number_of_traj,snr_60_db_w_constant,'-o','linewidth',2,'markersize',12)
legend('EKF, PSNR 40 dB','EKF, PSNR 50 dB','EKF, PSNR 60 dB')
xlabel('Number of trajectories')
ylabel('RMSE (dm)')
hold off