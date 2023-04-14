function [Q, M] = getQM(n_seg, n_order, ts)
    n_poly_perseg = n_order + 1;
    d_order = n_poly_perseg / 2;
    Q = [];
    M = [];
    M_k = getM(n_order);

    for k = 1:n_seg
        % ##############################################
        %   STEP 2.1 calculate Q_k of the k-th segment
        % ##############################################
        Q_k = zeros(n_poly_perseg, n_poly_perseg);
        t_k = 1;

        for i = 4:n_order

            for j = 4:n_order
                Q_k(i + 1, j + 1) = factorial(i) / factorial(i - d_order) * factorial(j) / factorial(j - d_order) / (i + j - n_order) * t_k ^ (i + k - n_order);
            end

        end

        Q = blkdiag(Q, Q_k);
        M = blkdiag(M, M_k);
    end

end
