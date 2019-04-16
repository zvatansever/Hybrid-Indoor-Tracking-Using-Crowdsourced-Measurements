function [tri_res,mean_rmse_tri]=trilate_wifi(d_est,tra)
%print_figures=0;
% aviobj = VideoWriter('trilater','MPEG-4');
% open(aviobj);
%% Load distances
%load log_normal_distance_noise d d_est
%load rand_paths tra
[m n v]=size(tra);
%%
for i=1:m
    trajec=tra(i,:,:);
    
    trajectory=reshape(trajec, [n, v]);
    
    traj1=round(trajectory(1,:));
    traj1(traj1<=0)=1;
    traj1(traj1>=50)=50;
    
    traj2=round(trajectory(2,:));
    traj2(traj2<=0)=1;
    traj2(traj2>=50)=50;
    
    % trajectory=[traj1;traj2];
    % d1=d;
    
    % d_est1=d;
    % d_est2=flip(d);
    % d_est3=fliplr(d);
    
    d_est1=d_est;
    d_est2=flip(d_est);
    d_est3=fliplr(d_est2);
    
    
    %% Measurement 1
    %for i=1:length(trajectory(1,:))
    for j=1:length(traj1)
        dist1(j)=abs(d_est1(traj1(j),traj2(j)));
    end
    %end
    % %d2=flip(d);
    %  d_est2=flip(d_est);
    for j=1:length(traj1)
        dist2(j)=abs(d_est2(traj1(j),traj2(j)));
    end
    % %d3=fliplr(d);
    %  d_est3=fliplr(d_est2);
    for j=1:length(traj1)
        dist3(j)=abs(d_est3(traj1(j),traj2(j)));
    end
    %%
    %% 1
    AP=[0 0;
        50 0;
        50 50];
    
    Xtri = zeros(1,length(traj1));
    Ytri = zeros(1,length(traj2));
    
    for j=1:length(traj1)
        [x, y] = two_tri(AP(1,1), AP(2, 1), AP(3,1), AP(1, 2), AP(2, 2), AP(3, 2), dist1(j), dist2(j), dist3(j));
        
        Xtri(j)=x;
        Ytri(j)=y;
        
        Xtri(Xtri<1)=1;
        Xtri(Xtri>50)=50;
        
        Ytri(Ytri<1)=1;
        Ytri(Ytri>50)=50;
        
        %         if print_figures
        %             drawnow
        %             plot(y,x,'*','linewidth',2)
        %             hold on
        %         end
    end
    %%
        error_dist1(i)=sqrt(mean((traj1-Xtri).^2+(traj2-Ytri).^2));

%     if print_figures
%         figure(1)
%         drawnow
%         plot(traj1,traj2,'--b',Xtri,Ytri,'-.r','linewidth',2)
%         legend('True','Trilateration')
%         title(['Trilateration results ' num2str(i) ' '])
%         
%            fprintf('RMSE=%f\n',error_dist1);
% 
%  
%     end
%     
    trilate_res=[Xtri;Ytri];
    tri_res(i,:,:)=trilate_res;
    
%     frame=getframe(figure(1));
%     writeVideo(aviobj,frame);
    
end
mean_rmse_tri=mean(error_dist1);
%fprintf('Mean RMSE=%f\n',mean_rmse_tri);
%close(aviobj)
%save tri_pos tri_res error_dist1 mean_rmse_tri
end
