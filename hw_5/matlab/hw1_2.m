close all; clear; clc;

figure; axis([0 100 0 100]); hold on;

path = [];
cnt = 0;

while true
    tmp = ginput(1);

    if size(tmp, 1) == 0
        break;
    else
        cnt = cnt + 1;
        path(cnt, :) = tmp;
        scatter(tmp(1), tmp(2));
    end

end

% Not visible
% path = ginput();

n_order = 7; % order of poly
n_seg = size(path, 1) - 1; % segment number
n_poly_perseg = n_order + 1; % coef number of perseg
ts = zeros(n_seg, 1); % time distribution

% calculate time distribution in proportion to distance between 2 points
dist = zeros(n_seg, 1);
dist_sum = 0;
T = n_seg;
t_sum = 0;

for i = 1:n_seg
    dist(i) = sqrt((path(i + 1, 1) - path(i, 1)) ^ 2 + (path(i + 1, 2) - path(i, 2)) ^ 2);
    dist_sum = dist_sum + dist(i);
end

for i = 1:n_seg - 1
    ts(i) = dist(i) / dist_sum * T;
    t_sum = t_sum + ts(i);
end

ts(n_seg) = T - t_sum;

% or you can simply set all time distribution as 1
% for i = 1:n_seg
%     ts(i) = 1.0;
% end

poly_coef_x = MinimumSnapCloseformSolver(path(:, 1), ts, n_seg, n_order);
poly_coef_y = MinimumSnapCloseformSolver(path(:, 2), ts, n_seg, n_order);

% display the trajectory
tstep = 0.01;
pts = T / tstep + 1;
X_n = zeros(pts, 1);
Y_n = zeros(pts, 1);
k = 1;

for i = 0:n_seg - 1
    % ##########################################################################
    %   STEP 3: get the coefficients of i-th segment of both x-axis and y-axis
    % ##########################################################################
    Pxi = poly_coef_x((n_order + 1) * i + 1:(n_order + 1) * i + n_order + 1);
    Pyi = poly_coef_y((n_order + 1) * i + 1:(n_order + 1) * i + n_order + 1);

    for t = 0:tstep:ts(i + 1)
        X_n(k) = polyval(flipud(Pxi), t);
        Y_n(k) = polyval(flipud(Pyi), t);
        k = k + 1;
    end

end

plot(X_n, Y_n, 'Color', [0 1.0 0], 'LineWidth', 2);

function poly_coef = MinimumSnapCloseformSolver(waypoints, ts, n_seg, n_order)
    start_cond = [waypoints(1), 0, 0, 0];
    end_cond = [waypoints(end), 0, 0, 0];

    % ##################################################
    %   you have already finished this function in hw1
    % ################################################
    Q = getQ(n_seg, n_order, ts);

    % #####################
    %   STEP 1: compute M
    % #####################
    M = getM(n_seg, n_order, ts);

    % ######################
    %   STEP 2: compute Ct
    % ######################
    Ct = getCt(n_seg, n_order);
    C = Ct';
    R = C * inv(M)' * Q * inv(M) * Ct;

    R_cell = mat2cell(R, [n_seg + 7, 3 * (n_seg - 1)], [n_seg + 7, 3 * (n_seg - 1)]);
    R_PP = R_cell{2, 2};
    R_FP = R_cell{1, 2};

    % ######################
    %   STEP 3: compute dF
    % ######################
    dF = zeros(8 + (n_seg - 1), 1);

    % start condition
    dF(1:4) = start_cond;

    % waypoint position
    for i = 1:n_seg - 1
        dF(4 + i) = waypoints(i + 1);
    end

    % end condition
    dF(end - 3:end) = end_cond;

    % calculate dP*
    dP = -inv(R_PP) * R_FP' * dF;

    % get poly
    poly_coef = inv(M) * Ct * [dF; dP];
end
