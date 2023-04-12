function M = getM(n_seg, n_order, ts)
    M = [];
    n_poly_perseg = n_order + 1;

    for k = 1:n_seg
        M_k = zeros(8, n_poly_perseg);

        % ###############################################
        %   STEP 1.1: calculate M_k of the k-th segment
        % ###############################################
        t = ts(k);
        M_k(1:4, :) = getCoeffCons(0);
        M_k(5:8, :) = getCoeffCons(t);

        M = blkdiag(M, M_k);
    end

end
