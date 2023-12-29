classdef TP1
    properties
        name;
        ul_bu;
        ul_bl;
        ll_bu;
        ll_bl;
        ul_nobj = 2;
        ll_nobj = 2;
        ul_con = 1;
        ll_con = 1;
    end

    methods
        function obj = TP1()
            obj.name = 'TP1';
            obj.ul_bl = [0];
            obj.ul_bu = [1];
            obj.ll_bl = [-1, -1];
            obj.ll_bu = [1, 1];
        end

        function [f, c] = evaluate_u(obj, y, x)
            f(:, 1) = x(:, 1) - y;
            f(:, 2) = x(:, 2);

            c = (1 + x(:, 1) + x(:, 2)) * (-1);

        end


        function [f, c] = evaluate_l(obj, y, x)
            f(:, 1) = x(:, 1);
            f(:, 2) = x(:, 2);

            c = ( y.^2 - x(:, 1) .^2 - x(:, 2) .^2) * (-1);
        end

        function pf = PF_UL(obj, n2)
            n = n2 * 100;
            y = linspace(1/sqrt(2), 1, n);

            x2 = -0.5 - 0.25 * sqrt(8 * y .^2 - 4);
            x2(1) = -0.5;
            x1 = -1 - x2;

            f1 = x1 - y;
            f2 = x2;
            k1 = [f1', f2'];

            x2 = -0.5 + 0.25 * sqrt(8 * y .^ 2 - 4);
            x2(1) = -0.5;
            x1 = -1 - x2;

            f1 = x1 - y;
            f2 = x2;
            k2 = [f1', f2'];

            pf = [k1; k2];
            % pf = mirror_ref_selection(pf, n2);
            pf = sparse_pf_selection(pf, n2);
        end

        function pf = PF_LL(obj, n, y)
            n2 = n * 2;
            angles = linspace(pi, pi * 3 / 2, n2);
            pf(:, 1) = cos(angles) * y;
            pf(:, 2) = sin(angles) * y;
            pf = sparse_pf_selection(pf, n);
        end

        function ps = PS_LL(obj, n, y)
            angles = linspace(pi, pi * 3 / 2, n);
            x1 = cos(angles) * y;
            x2 = sin(angles) * y;
            ps = [x1', x2'];
        end

        function ul_FElimit = return_ulFE(obj)
            ul_FElimit = 9334;
        end

        function ll_FElimit = return_llFE(obj)
            ll_FElimit = 447632;
        end

        function literature_hvmean = return_meanHV(obj)
            literature_hvmean = 0.3632;
        end

        function literature_hvstd = return_stdHV(obj)
            literature_hvstd = 0.00369;
        end

        function thv = target_hv(obj)
            thv =  0.3632;
        end

        function tigd = target_igd(obj)
            tigd = 0.05226;
        end




    end
end