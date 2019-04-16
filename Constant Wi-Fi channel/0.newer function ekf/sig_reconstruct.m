function[I0_lowpass_nearest]=sig_reconstruct(sparse_kf_est,P_floor_all,P_floor)
% close all;clc;clear all;
% set(0,'DefaultFigureWindowStyle','docked')
print_figures=0;
%load sparse_kf sparse_kf_est P_floor P_floor_all
[m n v]=size(sparse_kf_est);
sparse_for_rec=sparse_kf_est(1,:,:);
    
sparse1=reshape(sparse_for_rec, [n, v]);

 if print_figures
% figure
% surf(sparse1)
% 
figure
surf(P_floor)

end
% pause
%% 
[m n]=size(sparse1);
x=sparse1(:)';
x(x==0)=NaN;
%%
F1=abs(fillmissing(x,'spline'));
f_1=reshape(F1,[m n]);
P_floor=reshape(P_floor_all(1,:,:),[n v]);
rmse_spline=sqrt(immse(P_floor,f_1));
if print_figures
fprintf('MSE of Spline Method=%.3e\n',rmse_spline)
fprintf('-----------------------------------------------\n')
figure
% surf(P_floor1)
%hold on;
surf(f_1)
title('Spline Method')
end
%%
F2=abs(fillmissing(x,'pchip'));
f_2=reshape(F2,[m n]);
mse_pchip=immse(P_floor,f_2);
if print_figures
fprintf('MSE of Pchip Method=%.3e\n',mse_pchip)
fprintf('-----------------------------------------------\n')
figure
%surf(P_floor1)
%hold on;
surf(f_2)
title('Pchip Method')
end
%%
F3=abs(fillmissing(x,'linear'));
f_3=reshape(F3,[m n]);
mse_linear=immse(P_floor,f_3);
if print_figures
fprintf('MSE of Linear Method=%.3e\n',mse_linear)
fprintf('-----------------------------------------------\n')
figure
%surf(P_floor1)
%hold on;
surf(f_3)
title('Linear Method')
end
%%
F4=abs(fillmissing(x,'nearest'));
f_4=reshape(F4,[m n]);
rmse_nearest=sqrt(immse(P_floor,f_4));
% fprintf('MSE of Nearest Method=%.3e\n',rmse_nearest)
if print_figures
fprintf('MSE of Nearest Method=%.3e\n',rmse_nearest)
figure
%surf(P_floor1)
%hold on;
surf(f_4)
title('Nearest Method')
end
% %% Frequency Domain Filter Spline
% I0 = f_1;
% ff = fft2(f_1); % Take Fourier Transform 2D
% FF = fftshift(ff);
% FF = abs(FF); % Get the magnitude
% FF = 100*log(FF+1); % Use log, for perceptual scaling, and +1 since log(0) is undefined
% %figure, surf(F)
% K0 = abs(min(ff(:)));% cut-off frequency
% %F = 20*log(abs(fftshift(ff))); % Shift center; get log magnitude
% lp = fir1(18,K0);  % Generate the lowpass filter (order, cut-off frequency)
% %lp = fir1(144,K0); %for high res.
% lp_2D = ftrans2(lp);  % Convert to 2-dimensions
% I0_lowpass_spline = imfilter(I0,lp_2D,'replicate');
% spline_low_pass_mse=immse(P_floor,I0_lowpass_spline);
% spline_low_pass_rmse=sqrt(spline_low_pass_mse);
% if print_figures
% fprintf('-----------------------------------------------\n')
% fprintf('Spline Low Pass MSE=%.3e\n',spline_low_pass_mse)
% fprintf('Spline Low Pass RMSE=%.3e\n',spline_low_pass_rmse)
% figure
% surf(I0_lowpass_spline)
% title('Spline Method')
% figure
% surf(P_floor1-I0_lowpass_spline)
% title('Original and spline difference')
% fprintf('Is LPF_spline better than mse_spline=%d\n',spline_low_pass_mse<mse_spline)
% end
%% Frequency Domain Filter Pchip
I0 = f_2;
ff = fft2(f_2); % Take Fourier Transform 2D
FF = fftshift(ff);
FF = abs(FF); % Get the magnitude
FF = 100*log(FF+1); % Use log, for perceptual scaling, and +1 since log(0) is undefined
%figure, surf(F)
K0 = abs(min(ff(:)));% cut-off frequency
%F = 20*log(abs(fftshift(ff))); % Shift center; get log magnitude
lp = fir1(18,K0);  % Generate the lowpass filter (order, cut-off frequency)
%lp = fir1(144,K0); %for high res.
lp_2D = ftrans2(lp);  % Convert to 2-dimensions
I0_lowpass_pchip = imfilter(I0,lp_2D,'replicate');
pchip_low_pass_mse=immse(P_floor,I0_lowpass_pchip);
%pchip_low_pass_rmse=sqrt(pchip_low_pass_mse);
if print_figures
fprintf('-----------------------------------------------\n')
fprintf('Pchip Low Pass MSE=%.3e\n',pchip_low_pass_mse)
%fprintf('Spline Low Pass RMSE=%.3e\n',pchip_low_pass_rmse)

figure
surf(I0_lowpass_pchip)
title('Pchip Method')
fprintf('Is LPF_spline better than mse_pchip=%d\n',pchip_low_pass_mse<mse_pchip)
end
%% Frequency Domain Filter Linear
I0 = f_3;
ff = fft2(f_3); % Take Fourier Transform 2D
FF = fftshift(ff);
FF = abs(FF); % Get the magnitude
FF = 100*log(FF+1); % Use log, for perceptual scaling, and +1 since log(0) is undefined
%figure, surf(F)
K0 = abs(min(ff(:)));% cut-off frequency
%F = 20*log(abs(fftshift(ff))); % Shift center; get log magnitude
lp = fir1(18,K0);  % Generate the lowpass filter (order, cut-off frequency)
%lp = fir1(144,K0); %for high res.
lp_2D = ftrans2(lp);  % Convert to 2-dimensions
I0_lowpass_linear = imfilter(I0,lp_2D);
%linear_low_pass_mse=immse(pchip_low_pass_mse,I0_lowpass_linear);
if print_figures
fprintf('-----------------------------------------------\n')
%fprintf('Linear Low Pass MSE=%.3e\n',linear_low_pass_mse)
figure
surf(I0_lowpass_linear)
title('Linear Method')
%fprintf('Is LPF_linear better than mse_linear=%d\n',linear_low_pass_mse<mse_pchip)
end
%% Frequency Domain Filter Nearest
I0 = f_4;
ff = fft2(f_4); % Take Fourier Transform 2D
FF = fftshift(ff);
FF = abs(FF); % Get the magnitude
FF = 100*log(FF+1); % Use log, for perceptual scaling, and +1 since log(0) is undefined
%figure, surf(F)
K0 = abs(min(ff(:)));% cut-off frequency
%F = 20*log(abs(fftshift(ff))); % Shift center; get log magnitude
lp = fir1(18,K0);  % Generate the lowpass filter (order, cut-off frequency)
%lp = fir1(144,K0); %for high res.
lp_2D = ftrans2(lp);  % Convert to 2-dimensions
I0_lowpass_nearest = imfilter(I0,lp_2D);
%nearest_low_pass_mse=immse(pchip_low_pass_mse,I0_lowpass_nearest);
%nearest_low_pass_rmse=sqrt(nearest_low_pass_mse);
if print_figures
fprintf('-----------------------------------------------\n')
%fprintf('Nearest Low Pass MSE=%.3e\n',nearest_low_pass_mse)
%fprintf('Nearest Low Pass MSE=%.3e\n',nearest_low_pass_rmse)

figure
surf(I0_lowpass_nearest)
title('Nearest Method')
%fprintf('Is LPF_linear better than mse_linear=%d\n',nearest_low_pass_mse<mse_pchip)
end
% figure
% surf((P_floor1-I0_lowpass_spline)/mean2(P_floor1))
% title('Percentage Error of SPLINE + LPF For 60 dB in VLC channel')


% I0_lowpass_nearest(1:end,48:50)=0;
% I0_lowpass_nearest(48:50,1:end)=0;
end
%save reconstruct_nearest_traj1 I0_lowpass_nearest P_floor r P_floor_all
%save reconstruct_spline_traj1 I0_lowpass_spline P_floor1
