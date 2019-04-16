[tra]=function(NUMBER_OF_OBS,NUMBER_OF_TRAJ,seeds)
set(0,'DefaultFigureWindowStyle','normal')

% aviobj = VideoWriter('random_trajectories','Motion JPEG AVI');
% open(aviobj); 

print_figures=0;
X1_MIN = 0 ; % the state space domain
X1_MAX =  50 ; %          "
X2_MIN = 0 ; %          "
X2_MAX =  50 ; %          "

%NUMBER_OF_OBS = 100; % number of observations

clf;
axis([X1_MIN X1_MAX X2_MIN X2_MAX ]);
axis square;
hold on;
% the trajectory ----------------------------------------------------------
disp('   ');
disp('     - trajectory of the mobile:')
disp('           [left button to choose a point')
disp('           [right button for the last point')

x = [];
y = [];
n = 10;
% NUMBER_OF_TRAJ=1000;
tra=zeros(NUMBER_OF_TRAJ,2,NUMBER_OF_OBS+1);
rng(seeds,'twister')
for i=1:NUMBER_OF_TRAJ
    
    x = randi([0 50],n,1);
    y = randi([0 50],n,1);
    % end
    
    % Interpolate with two splines and finer spacing.
    t       = 1:n;
    delta   = (n-1)/NUMBER_OF_OBS; % 200 is the number of points of the
    % trajectory, i.e. the number of observations
    ts      = 1: delta: n;
    X1_TRAJ = spline(t,x,ts);
    X2_TRAJ = spline(t,y,ts);
    N_TRAJ  = length(ts);
    X1_TRAJ(X1_TRAJ<=0)=1;
    X1_TRAJ(X1_TRAJ>=50)=50;
    
    X2_TRAJ(X2_TRAJ<=0)=1;
    X2_TRAJ(X2_TRAJ>=50)=50;
    
    trajectory=[X1_TRAJ;X2_TRAJ];
    
%     if print_figures
%         figure(1)
%         drawnow
%         
%         plot(trajectory(1,:),trajectory(2,:),'r-');
%         title(['Random mobile user movement in the room ' num2str(i) ' '])
%         xlabel('Length of the room (dm)')
%         ylabel('Width of the room (dm)')
%     end
    tra(i,:,:)=trajectory;
%     frame=getframe(figure(1));
%     writeVideo(aviobj,frame);
end
   %close(aviobj)
%save rand_paths tra