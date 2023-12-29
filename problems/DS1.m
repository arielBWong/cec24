classdef DS1
    properties
        K;
        name;
        ul_bu;
        ul_bl;
        ll_bu;
        ll_bl;
        r = 0.1;
        tau = 1;
        gamma = 1;
        alpha = 1;
        ul_nobj = 2;
        ll_nobj = 2;
        ul_con = 0;
        ll_con = 0;

    end
    methods
        function obj = DS1(K, tau)
            
            
            obj.tau = tau;
            if obj.tau > 0
                obj.name = 'DS1';
            else
                obj.name = 'DS1D';
            end
            
            obj.K = K;
            obj.ul_bu = [4, ones(1, K-1) .* K];
            obj.ul_bl = [1, ones(1, K-1) * (-K)];
            obj.ll_bu = ones(1, K) * K;
            obj.ll_bl = ones(1, K) * (-K);
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            J = ((2 : obj.K) - 1)/2;
            f(:,1) = 1 + obj.r - cos(obj.alpha * pi * xu(:, 1)) + sum((xu(:, 2:obj.K)- J) .^ 2, 2)...
                + obj.tau * sum((xl(:, 2:obj.K) - xu(:, 2:obj.K)).^2, 2) - obj.r * cos(obj.gamma * pi * xl(:, 1) ./ xu(:,1) / 2);
            
            f(:,2) = 1 + obj.r - sin(obj.alpha * pi * xu(:,1)) + sum((xu(:,2:obj.K) - J).^2,2)...
                + obj.tau * sum((xl(:, 2:obj.K) - xu(:, 2:obj.K)).^2, 2) - obj.r * sin(obj.gamma * pi * xl(:, 1) ./ xu(:,1) / 2);
            
            c = [];
        end
        
        function [f, c] = evaluate_l(obj, xu, xl)
            f(:, 1) = xl(:, 1).^2 + sum((xl(:, 2:obj.K) - xu(:,2:obj.K)).^2, 2) + sum(10*(1 - cos(pi * (xl(:,2:obj.K) - xu(:, 2:obj.K))/obj.K)), 2);
            f(:, 2) = sum((xl(:, 1:obj.K) - xu(:, 1:obj.K)).^2, 2) + sum(10 * (abs(sin(pi * (xl(:, 2:obj.K) - xu(:, 2:obj.K))/obj.K))), 2);
            c = [];
        end
        
        function pf = PF_LL(obj, n, xu)
            if size(xu, 1) > 1
                error('LL PF only uses one row vector ');
            end
            
            if xu(1) > obj.ll_bu(1)
                sample_bu = obj.ll_bu(1);
            else
                sample_bu = xu(1);
            end
            n2 = n * 2;
            s = linspace(0, sample_bu, n2);
            f1 = s .^ 2;
            f2 = (s - xu(1)) .^ 2;
            pf = [f1', f2'];  
            pf = sparse_pf_selection(pf, n);
        end

        function ps = PS_LL(obj, n, xu)
            if size(xu, 1) > 1
                error('LL PF only uses one row vector ');
            end
            
            if xu(1) > obj.ll_bu(1)
                sample_bu = obj.ll_bu(1);
            else
                sample_bu = xu(1);
            end

             xl_1 = linspace(0, sample_bu, n);
             xl_rest = repmat(xu(2:end), n, 1);
             ps = [xl_1', xl_rest];
        end
        
        function pf = PF_UL(obj, n)
            n2 = n * 2;
            angles = linspace(pi, pi * 3/2, n2);
            f1 = cos(angles) * (1 + obj.r) + (1 + obj.r);
            f2 = sin(angles) * (1 + obj.r) + (1 + obj.r);
            pf = [f1', f2'];
            pf = sparse_pf_selection(pf, n);
        end
        
        function ul_FElimit = return_ulFE(obj)
            ul_FElimit = 37950;
        end
        
        function ll_FElimit = return_llFE(obj)
            ll_FElimit = 1023122;
        end
        
        function literature_hvmean = return_meanHV(obj)
            literature_hvmean = 1.1520;
        end
        
        function literature_hvstd = return_stdHV(obj)
            literature_hvstd = 0.00267;
        end

        function thv = target_hv(obj)
            thv = 1.152;
        end

        function tigd = target_igd(obj)
             tigd = 0.07858;
        end
        
        
    end
end