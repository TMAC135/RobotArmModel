% This function plots the movement of a 2-link arm
% It checks the dimension of the input states to accomodate the
% reaching (one arm) or signaling (two arms) task
% Can choose whether to close the figure after showing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
% l:      2 by 1 vector containing the length of the two links
% x:      states of the 2-link arm (polar coordinate)
% target: reaching target (polar coordinate)
% lambda: in case of signaling task, the weight of signal
% t:      approximate time of the movie shown

function [ ] = robot_arm_movie( l, x, target, lambda, t)
    % t is the approximate time of the movie
    nx = size(x, 1)/2;
    N = size(x, 2);
    switch nx
        case 2 % 2-link arm
        cartesianTarget = l(1) * [cos(target(1)); sin(target(1))];
        cartesianTarget = cartesianTarget + ...
            l(2) * [cos(target(1)+target(2)); sin(target(1)+target(2))];
        joints_pos = zeros(2, 3);
        figure, 
        for k = floor(linspace(1, N, floor(t/0.02)))
            joints_pos(:, 1) = [0;0];
            joints_pos(:, 2) = l(1)*[cos(x(1,k)); sin(x(1,k))];
            joints_pos(:, 3) = joints_pos(:, 2) + ...
                                l(2)*[cos(x(1,k)+x(2,k));sin(x(1,k)+x(2,k))];
            clf; % clear current figure
            plot(joints_pos(1,:),joints_pos(2,:)), hold on, % plot arms and hold
            plot(cartesianTarget(1),cartesianTarget(2),'*'), grid, % plot target
            title({'2-link robot arm reach';['target = ', mat2str(cartesianTarget)]}),
            axis([-l(1)-l(2),l(1)+l(2), -l(1)-l(2),l(1)+l(2)]),
            daspect([1,1,1]);
            pause(1e-2);
        end
                    
        case 4 % 2-link arm
        cartesianTarget_p = l(1) * [cos(target(1)); sin(target(1))];
        cartesianTarget_p = cartesianTarget_p + ...
            l(2) * [cos(target(1)+target(2)); sin(target(1)+target(2))];
        cartesianTarget_n = l(1) * [cos(target(5)); sin(target(5))];
        cartesianTarget_n = cartesianTarget_n + ...
            l(2) * [cos(target(5)+target(6)); sin(target(5)+target(6))];
        joints_pos = zeros(4, 3);
        figure, 
        for k = floor(linspace(1, N, floor(t/0.02)))
            clf; % clear current figure
            joints_pos(:, 1) = [0;0;0;0];
            joints_pos(:, 2) = l(1)*[cos(x(1,k)); sin(x(1,k)); cos(x(5,k)); sin(x(5,k))];
            joints_pos(:, 3) = joints_pos(:, 2) + ...
                                l(2)*[cos(x(1,k)+x(2,k));sin(x(1,k)+x(2,k));cos(x(5,k)+x(6,k));sin(x(5,k)+x(6,k))];
            if k == 1
                tip = joints_pos(:, 3);
            else
                tip = [tip, joints_pos(:, 3)];
            end
            plot(cartesianTarget_p(1), cartesianTarget_p(2), '*'), hold on,
            plot(cartesianTarget_n(1), cartesianTarget_n(2), '*'), hold on,
            plot(joints_pos(1,:), joints_pos(2,:), tip(1,:), tip(2,:)), hold on, % plot arms and hold
            plot(joints_pos(3,:), joints_pos(4,:), tip(3,:), tip(4,:)), hold on, % plot arms and hold
            title({'2-link robot arm reach'; ['lambda = ', num2str(lambda)]}),
            axis([-l(1)-l(2), l(1)+l(2), -l(1)-l(2), l(1)+l(2)]),
            daspect([1,1,1]);
            pause(1e-2);
        end
%         close; % close the figure after showing the movement
        
        case 6
            
    end    
end

