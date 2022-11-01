function M = getM(n_seg, n_order, ts)
    coeff = getCoeffCons(1);
    n_poly_perseg = n_order + 1;
    M = [];
    for k = 1: n_seg
        M_k = zeros(8, n_poly_perseg); % (start + end) x n_poly_perseg, i.e R^8 -> R^4(start) + R^4(end)
        % ###############################################
        %   STEP 1.1: calculate M_k of the k-th segment 
        % ###############################################
        t = ts(k);
        for i = 0: 3
            M_k(i + 1, i + 1) = coeff(i + 1, i + 1);
        end

        for i = 0: 3
            for j = 1: n_order
                if i == j
                    M_k(i + 4 + 1, j + 1) = coeff(i + 1, j + 1);
                else
                    M_k(i + 4 + 1, j + 1) = coeff(i + 1, j + 1) * t^(j - 1);
                end
            end
        end

        M = blkdiag(M, M_k);
    end
end