% This function plots the movement of a 3-link arm
% It checks the dimension of the input states to accomodate the
% reaching (one arm) or signaling (two arms) task
% Can choose whether to close the figure after showing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
% l:      3 by 1 vector containing the length of the three links
% x:      states of the 3-link arm (polar coordinate)
% target: reaching target (polar)
% lambda: in case of signaling task, the weight of signal
% t:      approximate time of the movie shown

function [] = real_arm_movie_2_link( l, x, target, Ns,r_rest,d,ss)
    % t is the approximate time of the movie
    n = size(x, 1)/2;
    N = size(x, 2);
    L = sum(l);

    joint = zeros(2,N);
    tip = zeros(2,N);
    targetCartesian = tip_fun(target);
    for k = 1:N
        joint(:,k) = joint_fun(x(:,k));
        tip(:,k) = tip_fun(x(:,k));
    end
    figure, 
    for k = floor(linspace(1, N, Ns))
        clf; % clear current figure
        x1=[0,joint(1,k)];
        y1=[0,joint(2,k)];
        x2=[joint(1,k),tip(1,k)];
        y2=[joint(2,k),tip(2,k)];       
        plot(x1,y1,'-bo','LineWidth',2,'MarkerSize',6);
        hold on
        plot(x2,y2,'-go','LineWidth',2,'MarkerSize',6);
        hold on
        %plot([0,joint(1,k),tip(1,k)], [0,joint(2,k),tip(2,k)]), hold on, 
        plot(targetCartesian(1),targetCartesian(2),'*'), grid, % plot target
        title({'2-link robot arm reach';...
           ['restlength = ', mat2str(r_rest)];['ss = ', mat2str(ss)];['d = ', mat2str(d)]}),
        axis([-L,L,-L,L]),
        daspect([1,1,1]);
        pause(0.001);
    end    
end

