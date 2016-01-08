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

function [ ] = robot_arm_movie_3_link( l, x, target, t)
    % t is the approximate time of the movie
    n = size(x, 1)/2;
    N = size(x, 2);
    L = sum(l);
    pause_time = 0.02;
    switch n
        case 3 % 3-link arm reaching
            joints_pos = zeros(2, 4);
            targetCartesian = polar2Cartesian_fun(target);
            figure, 
            for k = floor(linspace(1, N, floor(t/pause_time)))
                joints_pos(:, 1) = [0;0];
                joints_pos(:, 2) = l(1)*[cos(x(1,k)); sin(x(1,k))];
                joints_pos(:, 3) = joints_pos(:, 2) + ...
                                   l(2)*[cos(x(1,k)+x(2,k));sin(x(1,k)+x(2,k))];
                joints_pos(:, 4) = joints_pos(:, 3) + ...
                                   l(3)*[cos(x(1,k)+x(2,k)+x(3,k));sin(x(1,k)+x(2,k)+x(3,k))];
                clf; % clear current figure
                plot(joints_pos(1,:),joints_pos(2,:)), hold on, % plot arms and hold
                plot(targetCartesian(1),targetCartesian(2),'*'), grid, % plot target
                title({'2-link robot arm reach';['target = ', mat2str(target)]}),
                axis([-L,L,-L,L]),
                daspect([1,1,1]);
                pause(pause_time);
            end
        
        case 6 % 3-link arm signaling
            targetCartesian_p = polar2Cartesian_fun(target(1:n));
            targetCartesian_n = polar2Cartesian_fun(target((n+1):(2*n)));
            index = floor(linspace(1, N, floor(t/pause_time))); % sample points
            H = length(index); % Number of sample points
            joints_pos = zeros(4, 4, H);
            joints_tra = zeros(4, H, 3);

            for i = 1:H
                k = index(i);
                joints_pos(:, 1, i) = [0;0;0;0];
                joints_pos(:, 2, i) = l(1)*[cos(x(1,k)); sin(x(1,k)); ...
                                            cos(x(7,k)); sin(x(7,k))];
                joints_pos(:, 3, i) = joints_pos(:, 2, i) + ...
                                      l(2)*[cos(x(1,k)+x(2,k));sin(x(1,k)+x(2,k)); ...
                                            cos(x(7,k)+x(8,k));sin(x(7,k)+x(8,k))];
                joints_pos(:, 4, i) = joints_pos(:, 3, i) + ...
                                      l(3)*[cos(x(1,k)+x(2,k)+x(3,k));sin(x(1,k)+x(2,k)+x(3,k)); ...
                                            cos(x(7,k)+x(8,k)+x(9,k));sin(x(7,k)+x(8,k)+x(9,k))];
                joints_tra(:, i, 1) = joints_pos(:,2,i);
                joints_tra(:, i, 2) = joints_pos(:,3,i);
                joints_tra(:, i, 3) = joints_pos(:,4,i);
            end
            figure,
            for i = 1:H
                clf; % clear current figure
                plot(joints_pos(1,:,i), joints_pos(2,:,i), 'k'), hold on, % plot arms and hold
                plot(joints_pos(3,:,i), joints_pos(4,:,i), 'k'), hold on,
                plot(joints_tra(1,:,1), joints_tra(2,:,1), 'b'), hold on, % plot trajectories
                plot(joints_tra(3,:,1), joints_tra(4,:,1), 'b'), hold on,
                plot(joints_tra(1,:,2), joints_tra(2,:,2), 'C'), hold on,
                plot(joints_tra(3,:,2), joints_tra(4,:,2), 'c'), hold on,
                plot(joints_tra(1,:,3), joints_tra(2,:,3), 'g'), hold on,
                plot(joints_tra(3,:,3), joints_tra(4,:,3), 'g'), hold on,
                plot(targetCartesian_n(1),targetCartesian_n(2),'*r'), hold on,
                plot(targetCartesian_p(1),targetCartesian_p(2),'*r'), grid, % plot target
                title({'2-link robot arm reach';['target = ', mat2str(target)]}),
                axis([-L,L,-L,L]),
                daspect([1,1,1]);
                pause(pause_time);
            end
    end    
end

