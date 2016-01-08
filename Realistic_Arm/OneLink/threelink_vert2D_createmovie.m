% Plots three link robotic arm and captures the plot as a movie frame,
% over simulation data of the arm. Number of plots/frames depends on 
% length of simulation and desired frame rate.
%
% Intended for use with 'threelink_vert2D_sim' function. Two of the
% inputs to this function are the outputs of that function and another
% input to this function is also an input to that function.
%
% Function inputs:
%
%   state trajectory (x): 6 x length(0:dt:t_f) matrix,
%                         rows 1-3 are joint angles (theta) in rad 
%                         rows 4-6 are joint vels. (theta_dot) in rad/sec
%                           from 'threelink_horz2D_sim' function
%
%   simulation time (t) : 1 x length(t) vector,
%                           from 'threelink_horz2D_sim' function
%
%   link lengths (len)  : 3 x 1 vector, lengths in m
%                           from 'threelink_horz2D_sim' function
%
%   desired frame rate  : scalar, 
%   (fps_des)             number of plots/frames to display per sec
%                           example: fps = 30;
%
% Function outputs:
%
%   movie struct (M)    : 1 x length(M) structure array,
%                         each "column" in the array is a frame in the
%                         movie and is a struct having two fields ('cdata'
%                         and 'colormap'), each with one value (frame data)

function M = threelink_vert2D_createmovie(x,t,fps_des,len)

% Determine the number of frames for the movie from the desired fps
% and the simulation time length
num_frames = ceil(fps_des*t(end));

% Determine time step spacing between frames - must be an integer
% (rounding up to the next integer guarantees that the frame index will
% not exceed the simulation time steps index)
frame_space = ceil(length(t)/num_frames);

% Determine index for selecting which time steps are frames
% (always includes frames for t = 0 and t = t_final)
frame_index = 1:frame_space:length(t);
if frame_index(end) ~= length(t)
    frame_index = [1:frame_space:length(t),length(t)];
end

% Set link lengths
l1 = len(1);
%l2 = len(2);
%l3 = len(3);

% Initialize figure and set axes for frame capture
% (prevent the frames from being displayed during capture)
figure;
max_rad = l1; %+ l2 + l3;
%axis([-1.2*max_rad,1.2*max_rad,-1.2*max_rad,1.2*max_rad]);
xlabel('Horizontal plane');
ylabel('Vertical (normal) axis');
title('Side view of vertical planar arm');

% Initialize movie structure array
M(length(frame_index)) = struct('cdata',[],'colormap',[]);

% Plotting and frame capture loop
for i = 1:length(frame_index)
    
    % Set theta angles at the current time step index
    th1 = x(1,frame_index(i));
    %th2 = x(2,frame_index(i));
    %th3 = x(3,frame_index(i));
    
    % Calculate (x,y) for link end points at current time step index
    % (first column is end point connected to base or previous link)
    x1 = [0,l1*cos(th1)];
    y1 = [0,l1*sin(th1)];

    %x2 = [l1*cos(th1),l1*cos(th1)+l2*cos(th1+th2)];
    %y2 = [l1*sin(th1),l1*sin(th1)+l2*sin(th1+th2)];

    %x3 = [l1*cos(th1)+l2*cos(th1+th2),...
     %      l1*cos(th1)+l2*cos(th1+th2)+l3*cos(th1+th2+th3)];
    %y3 = [l1*sin(th1)+l2*sin(th1+th2),...
      %     l1*sin(th1)+l2*sin(th1+th2)+l3*sin(th1+th2+th3)];
    
    % Plot the robotic arm at the current time step index
    plot(x1,y1,'-bo','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor','b')
    %hold on
    %plot(x2,y2,'-go','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','g','MarkerFaceColor','g')
    %plot(x3,y3,'-ro','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r')
    %plot((l1+l2+l3)*cos(x_target(1)),(l1+l2+l3)*sin(x_target(1)),'*');
    %plot(l1*cos(x_target(1)),l1*sin(x_target(1)),'*');
    %hold off
    
    % Reset the axes and labels
    axis([-1.2*max_rad,1.2*max_rad,-1.2*max_rad,1.2*max_rad]);
    xlabel('Horizontal plane');
    ylabel('Vertical (normal) axis');
    title('Side view of vertical planar arm');
    
    % Capture the plot as the ith frame
    M(i) = getframe;
end