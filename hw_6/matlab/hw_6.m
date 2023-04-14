close all; clear; clc;

v_max = 400;
a_max = 400;
color = ['r', 'b', 'm', 'g', 'k', 'c', 'c'];

% specify the center points of the flight corridor and the region of corridor
path = [50, 50;
        100, 120;
        180, 150;
        250, 80;
        280, 0];
x_length = 100;
y_length = 100;

n_order = 7;
n_poly_perseg = n_order + 1;
n_seg = size(path, 1);

corridor = zeros(4, n_seg);

for k = 1:n_seg
    corridor(:, k) = [path(k, 1), path(k, 2), x_length / 2, y_length / 2]';
end

% specify ts for each segment
% NOTE: currently time duration for each Bezier curve only supports [0, 1], ts is not used...
ts = zeros(n_seg, 1);

% for k = 1:n_seg
%     ts(k, 1) = 0.1;
% end

poly_coef_x = MinimumSnapCorridorBezierSolver(1, path(:, 1), corridor, ts, n_seg, n_order, v_max, a_max);
poly_coef_y = MinimumSnapCorridorBezierSolver(2, path(:, 2), corridor, ts, n_seg, n_order, v_max, a_max);

% display the trajectory and cooridor
figure; hold on;
plot(path(:, 1), path(:, 2), '*r');

for k = 1:n_seg
    plot_rect([corridor(1, k); corridor(2, k)], corridor(3, k), corridor(4, k));
end

t_sample = 0.01;
x_pos = zeros(1, n_seg * (1 / t_sample + 1));
y_pos = zeros(1, n_seg * (1 / t_sample + 1));
idx = 1;

% #############################
%   STEP 4: draw bezier curve
% #############################
for k = 0:n_seg - 1

    for t = 0:t_sample:1

        for i = 0:n_order
            basis_p = nchoosek(n_order, i) * t ^ i * (1 - t) ^ (n_order - i);
            x_pos(idx) = x_pos(idx) + poly_coef_x(k * n_poly_perseg + i + 1) * basis_p;
            y_pos(idx) = y_pos(idx) + poly_coef_y(k * n_poly_perseg + i + 1) * basis_p;
        end

        idx = idx + 1;
    end

end

for k = 0:n_seg - 1
    scatter(poly_coef_x(n_poly_perseg * k + 1:n_poly_perseg * (k + 1)), poly_coef_y(n_poly_perseg * k + 1:n_poly_perseg * (k + 1)), 100, color(k + 1));
end

plot(x_pos', y_pos', "Linewidth", 2.0, "color", "k");

function poly_coef = MinimumSnapCorridorBezierSolver(axis, waypoints, corridor, ts, n_seg, n_order, v_max, a_max)
    start_cond = [waypoints(1), 0, 0];
    end_cond = [waypoints(end), 0, 0];

    % #################################
    %   STEP 1: compute Q_0 of c'Q_0c
    % #################################
    [Q, M] = getQM(n_seg, n_order, ts);
    Q_0 = M' * Q * M;
    Q_0 = nearestSPD(Q_0);

    % ###########################
    %   STEP 2: get Aeq and beq
    % ###########################
    [Aeq, beq] = getAbeq(n_seg, n_order, ts, start_cond, end_cond);

    % #####################################################
    %   STEP 3: get corridor_range, Aieq and bieq
    %   STEP 3.1: get corridor_range of x-axis or y-axis,
    %   you can define corridor_range as [p1_min, p1_max;
    %                                   p2_min, p2_max;
    %                                   ...,
    %                                   pn_min, pn_max ];
    % #####################################################
    corridor_range = zeros(n_seg, 2);

    for k = 1:n_seg
        corridor_range(k, :) = [corridor(axis, k) - corridor(axis + 2, k), corridor(axis, k) + corridor(axis + 2, k)];
    end

    % ###############################
    %   STEP 3.2: get Aieq and bieq
    % ###############################
    [Aieq, bieq] = getAbieq(n_seg, n_order, corridor_range, ts, v_max, a_max);

    f = zeros(size(Q_0, 1), 1);
    poly_coef = quadprog(Q_0, f, Aieq, bieq, Aeq, beq);
end

function plot_rect(center, x_r, y_r)
    p1 = center + [-x_r; -y_r];
    p2 = center + [-x_r; y_r];
    p3 = center + [x_r; y_r];
    p4 = center + [x_r; -y_r];
    plot_line(p1, p2);
    plot_line(p2, p3);
    plot_line(p3, p4);
    plot_line(p4, p1);
end

function plot_line(p1, p2)
    a = [p1(:), p2(:)];
    plot(a(1, :), a(2, :), 'b');
end
