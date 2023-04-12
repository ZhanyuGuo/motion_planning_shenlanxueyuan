function Ct = getCt(n_seg, n_order)

    % #########################################
    %   STEP 2.1: finish the expression of Ct
    % #########################################
    d_order = (n_order + 1) / 2;
    m = 2 * d_order * n_seg;
    n = 2 * d_order + (n_seg - 1) * d_order;
    Ct = zeros(m, n);

    % hash map
    d_vector = zeros(m, 1);
    i = 1;

    for k = 1:n_seg
        % k-th piece
        for t = 0:1
            % start and end
            for d = 0:d_order -1
                % derivative
                d_vector(i) = getHash(k, t, d);
                i = i + 1;
            end

        end

    end

    % fix start condition
    col = 1;
    k = 1;
    t = 0;

    for d = 0:d_order - 1
        val = getHash(k, t, d);
        [row, ~] = find(d_vector == val);
        Ct(row, col) = 1;
        col = col + 1;
    end

    % fix waypoint position
    t = 1;
    d = 0;

    for k = 1:n_seg - 1
        val = getHash(k, t, d);
        [row, ~] = find(d_vector == val);
        Ct(row, col) = 1;

        val = getHash((k + 1), (t - 1), d);
        [row, ~] = find(d_vector == val);
        Ct(row, col) = 1;

        col = col + 1;
    end

    % fix end condition
    k = n_seg;
    t = 1;

    for d = 0:d_order - 1
        val = getHash(k, t, d);
        [row, ~] = find(d_vector == val);
        Ct(row, col) = 1;
        col = col + 1;
    end

    % other free variables
    t = 1;

    for k = 1:n_seg - 1

        for d = 1:d_order - 1
            val = getHash(k, t, d);
            [row, ~] = find(d_vector == val);
            Ct(row, col) = 1;

            val = getHash((k + 1), (t - 1), d);
            [row, ~] = find(d_vector == val);
            Ct(row, col) = 1;

            col = col + 1;
        end

    end

end

function val = getHash(k, t, d)
    val = k * 100 + t * 10 + d;
end
