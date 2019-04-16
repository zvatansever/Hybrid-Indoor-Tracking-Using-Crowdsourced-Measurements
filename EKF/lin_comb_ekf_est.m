function [matrixOut1 P_floor]=lin_comb_ekf_est(I0_lowpass_nearest,P_floor_all,P_floor,sparse_kf_est,kal_res)

%clc;clear all;close all;
set(0,'DefaultFigureWindowStyle','docked')
print_figures=0;
%  aviobj = VideoWriter('lin_comb_one_by_one_40_dB','Motion JPEG AVI');
%  open(aviobj);
%%
%load power_map_decimeters.mat
%load reconstruct_nearest_traj1 I0_lowpass_nearest P_floor r P_floor_all

%load reconstruct_spline_traj1 I0_lowpass_spline

%load sparse_kf sparse_kf_est
%load kal_filt_tri_result kal_res


% if print_figures
%     figure
%     surf(I0_lowpass_nearest)
%     title('Initial Map')
% end
% pause
% if print_figures
%     figure
%     surf(P_floor)
%     title('Initial True Map')
% end
% pause

[a b]=size(P_floor);
[m n v]=size(kal_res);
%w=1:m; % weights
w=ones(1,m);
%% First Trajectory
for k=2:m
    
    kal_power=kal_res(k,:,:);
    
    kal_power=reshape(kal_power, [n, v]);
    
    kal_est1=round(kal_power(1,:));
    
    kal_est2=round(kal_power(2,:));
    
    
    
    x_dir=round(kal_est1);
    x_dir(x_dir<=0)=1;
    x_dir(x_dir>=50)=50;
    y_dir=round(kal_est2);
    y_dir(y_dir<=0)=1;
    y_dir(y_dir>=50)=50;
    % Low pass trilateration map Measurement 1
    for j=1:length(y_dir)
        %      New1(y_dir(j),x_dir(j))=abs(P_floor(y_dir(j),x_dir(j)));
        Newer1(j)=I0_lowpass_nearest(x_dir(j),y_dir(j));
        %Newer1_noisy(1,j)=D_T_R1_xy_noisy(x_dir(j),y_dir(j));
    end
    
    % Real Time Measurement 1
    for j=1:length(y_dir)
        %      New1(y_dir(j),x_dir(j))=abs(P_floor(y_dir(j),x_dir(j)));
        Meas1(j)=P_floor_all(k,x_dir(j),y_dir(j));
        %Newer1_noisy(1,j)=D_T_R1_xy_noisy(x_dir(j),y_dir(j));
    end
    %
    % Linear combination of LPT and RT maps
    for j=1:length(Newer1)
        Pnew1(j)=0.5^w(k)*Newer1(j)+(1-0.5^w(k))*Meas1(j);
    end
    rmse1=sqrt(mean((Newer1-Pnew1).^2));
    
    
    % Replace values in LPT with lin comb vals.
    for i=1:length(Newer1)
        %for j=1:length(Ytri)
        I0_lowpass_nearest(x_dir(i),y_dir(i))=Pnew1(i);
        %end
    end
    
    
    h = fspecial('gaussian', [3 3],0.5);
    
    matrixOut1 = imfilter(I0_lowpass_nearest,h,'symmetric');
    
    %matrixOut1 = medfilt2(I0_lowpass_nearest,[1 1],'symmetric');
    P_floor=reshape(P_floor_all(k,:,:),[a b]);
    rmse2=sqrt(mean((P_floor(:)-matrixOut1(:)).^2));
    
    
    %     figure
    %     plot(1:5:length(Pnew1),Pnew1(1:5:end),'-hex',1:5:length(Newer1),Newer1(1:5:end),'-square',1:5:length(Meas1),Meas1(1:5:end),'-penta',...
    %         'linewidth',2,'markersize',8)
    %     legend('Estimated map meas.','Spline map meas.','True meas.')
    %     title('Linear Combination n=1, $\alpha=0.5$','Interpreter', 'Latex')
    %     figure
    %     plot((Newer1-Pnew1)/mean(Newer1))
    %     title('Percentage error between Spline meas and Estimated Map')
    %     figure
    %     surf(I0_lowpass)
    %     title('Values in Spline map replaced with Real time meas')
    %     figure
    %     surf((I0_lowpass-matrixOut1)/mean2(I0_lowpass))
    %     title('Percentage error between Estimated Map and Updated-Estimated Map')
    %     figure
    %     surf((P_floor-matrixOut1)/mean2(P_floor))
    %     title('Percentage error between Real Map and Updated-Estimated Map')
    %     if print_figures
    %         figure(1)
    %         drawnow
    %         subplot 121
    %         surf(matrixOut1)
    %         title(['Updated-Estimated Map after trajectory # ' num2str(k) ' '])
    %         subplot 122
    %         surf((P_floor-matrixOut1)/mean2(P_floor))
    %         title(['Percentage error between Real Map and Updated-Estimated Map after trajectory # ' num2str(k) ' '])
    %         %             frame=getframe(figure(1));
    %         %             writeVideo(aviobj,frame);
    %         fprintf(['RMSE btwn. Low pass map and Lin. Comb. on trajectory # ' num2str(k) ' =%.3e\n'],rmse1)
    %         fprintf('\n')
    %
    %
    %         fprintf(['RMSE btwn True Power Map and Updated-Estimated on trajectory # ' num2str(k) '  =%.3e\n'],rmse2)
    %         fprintf('\n')
    %     end
    
    
end
end
% close(aviobj)
%save lin_comb_res matrixOut1 P_floor r