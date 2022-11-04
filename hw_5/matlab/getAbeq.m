function [Aeq, beq]= getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond)
    n_poly_perseg = n_order + 1;
    n_all_poly = n_seg * n_poly_perseg;
    % ################################
    %   p,v,a,j constraint in start,
    % ################################
    Aeq_start = zeros(4, n_all_poly);
    % beq_start = zeros(4, 1);
    % #########################################################
    %   STEP 2.1: write expression of Aeq_start and beq_start
    % #########################################################
    Aeq_start(1:4, 1:8) = getCoeffCons(0);
    beq_start = start_cond';

    % ################################
    %   p,v,a,j constraint in end,
    % ################################
    Aeq_end = zeros(4, n_all_poly);
    % beq_end = zeros(4, 1);
    % #####################################################
    %   STEP 2.2: write expression of Aeq_end and beq_end
    % #####################################################
    Aeq_end(1:4, end - 7:end) = getCoeffCons(ts(end));
    beq_end = end_cond';
    
    % ##############################################
    %   position constrain in all middle waypoints
    % ##############################################
    Aeq_wp = zeros(n_seg - 1, n_all_poly);
    beq_wp = zeros(n_seg - 1, 1);
    % ###################################################
    %   STEP 2.3: write expression of Aeq_wp and beq_wp
    % ###################################################
    for k = 0: n_seg - 2
        beq_wp(k + 1, 1) = waypoints(k + 2);
        coeff = getCoeffCons(ts(k + 1));
        Aeq_wp(k + 1, k * 8 + 1: k * 8 + 8) = coeff(1, :);
    end

    % #########################################################
    %   position continuity constrain between each 2 segments
    % #########################################################
    Aeq_con_p = zeros(n_seg - 1, n_all_poly);
    beq_con_p = zeros(n_seg - 1, 1);
    Aeq_con_v = zeros(n_seg - 1, n_all_poly);
    beq_con_v = zeros(n_seg - 1, 1);
    Aeq_con_a = zeros(n_seg - 1, n_all_poly);
    beq_con_a = zeros(n_seg - 1, 1);
    Aeq_con_j = zeros(n_seg - 1, n_all_poly);
    beq_con_j = zeros(n_seg - 1, 1);
    % #########################################################
    %   STEP 2.4: write expression of Aeq_con_p and beq_con_p
    % #########################################################

    % #########################################################
    %   velocity continuity constrain between each 2 segments
    % #########################################################
    % Aeq_con_v = zeros(n_seg - 1, n_all_poly);
    % beq_con_v = zeros(n_seg - 1, 1);
    % #########################################################
    %   STEP 2.5: write expression of Aeq_con_v and beq_con_v
    % #########################################################

    % #############################################################
    %   acceleration continuity constrain between each 2 segments
    % #############################################################
    % Aeq_con_a = zeros(n_seg - 1, n_all_poly);
    % beq_con_a = zeros(n_seg - 1, 1);
    % #########################################################
    %   STEP 2.6: write expression of Aeq_con_a and beq_con_a
    % #########################################################

    % #####################################################
    %   jerk continuity constrain between each 2 segments
    % #####################################################
    % Aeq_con_j = zeros(n_seg - 1, n_all_poly);
    % beq_con_j = zeros(n_seg - 1, 1);
    % #########################################################
    %   STEP 2.7: write expression of Aeq_con_j and beq_con_j
    % #########################################################

    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a; Aeq_con_j];
    beq_con = [beq_con_p; beq_con_v; beq_con_a; beq_con_j];

    for k = 0: n_seg - 2
        Aeq_con(4 * k + 1: 4 * k + 4, 8 * k + 1: 8 * k + 8) = getCoeffCons(ts(k + 1)); % k's end
        Aeq_con(4 * k + 1: 4 * k + 4, 8 * (k + 1) + 1: 8 * (k + 1) + 8) = -getCoeffCons(0); % k + 1's start
    end

    % ##############################################
    %   combine all components to form Aeq and beq
    % ##############################################
    Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
    beq = [beq_start; beq_end; beq_wp; beq_con];
end
