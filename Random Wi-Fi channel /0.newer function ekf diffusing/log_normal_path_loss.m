function [P_log_dB_map,d_est,d]=log_normal_path_loss(n_loss,sigma,i)
% set(0,'DefaultFigureWindowStyle','normal')

randn('seed',9999+i);
Pt_mw=40; %Transmitted power in mW
Pt_dBm=10*log10(40);
Gt=1;%Gain of the Transmitted antenna in dBi 
Gr=1;%Gain of the Receiver antenna in dBi
f=10^9;%Transmitted signal frequency in Hertz 
lambda=3*10^8/f; %Wavelength in meters
%% Room Model
xr=linspace(1,50,50); %Create the power grid
yr=linspace(1,50,50);
xt=0;yt=0;
for i=1:length(xr)
    for j=1:length(yr)
d(i,j) = sqrt((xt-xr(i)).^2+(yt-yr(j)).^2); %Array of TX and RX distances in meters 
    end
end
[m n]=size(d);
%%
distance=reshape(d,[1 m*n]);
L=2; %Other System Losses, No Loss case L=1
d0=d(1,1); % Reference distance

%% From Friis model
Pr_d0=(Pt_mw*Gt*Gr*lambda^2)/((4*pi*d0)^2*L); % Reference power at reference distance d0 
Pr_d0_dBm=10*log10(Pr_d0); %Convert to dBm

%%
n_loss=4;%path loss exponent
%sigma=(1.25); %mW
sigma_dbm=10*log10(sigma*randn);  %dBm
nu=10/log(10); %constant
noise_sigma=(1); % dB
%% No noise log normal loss

P_log_dBm=Pr_d0_dBm-10*n_loss*log10(distance/d0); % Eq. 6
P_signal=10.^(P_log_dBm./10); %mW
%% Log normal shadowing loss
% X_path_dBm= (0 + sigma_dbm*randn(1,length(P_log_dBm))); %dBm
X_path=(0+sigma*randn(1,length(P_log_dBm))); %mW
%SNR=min(10.*log10(((10.^(P_log_dBm)./10)./(10.^(X./10)))));
% fprintf('SNR=%f',SNR)
% %% Randomness noise
%  X_mW=1;
%  X_dBm=10*log10(X_mW);
%  X_noise_dBm=X_dBm*randn(1,length(P_log_dBm)); % dB
% X_noise_mW=10.^((X_noise_dBm)./10); % mW
%% Path loss+shadowing
P_log_dBm=P_log_dBm+X_path;%% bug dbm mw FIXED
%%+X_noise_dBm; % Eq. 6
P_log=10.^((P_log_dBm)./10); %mW
%% Distance estimation
%d_est=d0*((P_log(1)/Pr_d0)^(-1/n_loss))*exp(-(sigma^2)/(2*n_loss^2*nu^2))
  for i=1:length(P_log)
  d_est(i)=d0*((P_log(i)/Pr_d0)^(-1/n_loss))*exp(-(sigma^2)/(2*n_loss^2*nu^2)); %% bug sigma dbm to sigma FIXED
  end
  d_est=reshape(d_est,[m n]);

%% SNR 
% SNR=10*log10(mean(P_log.^2)/mean(X_noise_mW.^2));
% fprintf('Signal-to-Noise Ratio=%.3f\n',mean(SNR))

%%

P_log_dB_map=reshape(P_log_dBm,[m n]);
end
% figure
% surf(-P_log_dB_map)
% title('Log Normal Path Loss')
% zlabel('dB')
% figure
% surf(d)
% title('True Euclidean Distance')
% zlabel('Decimeter')
% figure
% surf(d_est)
% title('Estimated Euclidean Distance')
% zlabel('Decimeter')
% 
% 
% % save log_normal_distance d d_est
% save log_normal_distance_noise d d_est
