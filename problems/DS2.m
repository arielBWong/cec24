classdef DS2
    properties
        K;
        name;
        ul_bu;
        ul_bl;
        ll_bu;
        ll_bl;
        r = 0.25;
        tau = 1;
        gamma = 4;
        ul_nobj = 2;
        ll_nobj = 2;
        ul_con = 0;
        ll_con = 0;
    end
    methods
        function obj = DS2(K, tau)
            
          
            obj.K = K;
            obj.tau = tau;
            
            if obj.tau > 0
                obj.name = 'DS2';
            else
                obj.name = 'DS2D';
            end
            
            obj.ul_bu = ones(1, K) .* K;
            obj.ul_bl = [0.001, ones(1, K-1) * (-K)];
            obj.ll_bu = ones(1, K) * K;
            obj.ll_bl = ones(1, K) * (-K);
            
        end
        
        function [v1, v2] = trade_off(obj, y)
            if size(y, 1) > 1
                error('can only process row vector');
            end
            
            if  y(1) <= 1 && y(1) > 0
                v1 = cos(0.2 * pi) * y(1)+ sin(0.2 * pi) * sqrt(abs(0.02 * sin(5 * pi * y(1))));
                v2 = -sin(0.2 * pi) * y(1) + cos(0.2 * pi) * sqrt(abs(0.02 * sin(5 * pi * y(1))));
            else                
                v1 = y(1) - (1 - cos(0.2* pi));
                v2 = 0.1 * (y(1) - 1) - sin(0.2 * pi);
            end
        end
        
        function [f, c] = evaluate_u(obj, y, x)
            n = size(y, 1);
            f = [];
            for i = 1: n
                
                middle_term = sum((y(i, 2:end).^2 + 10*(1-cos(pi/obj.K * y(i, 2:end)))), 2) + obj.tau * sum((x(i, 2:end) - y(i, 2:end)).^2);
                [v1, v2] = obj.trade_off(y(i, :));
                
                f1 = v1 + middle_term - obj.r * cos(obj.gamma * pi * x(i, 1)/(2 *y(i, 1)));
                f2 = v2 + middle_term - obj.r * sin(obj.gamma * pi * x(i, 1)/(2 * y(i, 1)));
                
                f = [f; [f1, f2]];
            end
            c = [];
        end
        
        function [f, c] = evaluate_l(obj, y, x)
            I = 1:obj.K;
            
            f1 = x(:, 1).^2 + sum((x(:, 2:end) - y(:, 2:end)).^2, 2);
            f2 = sum(I .* ((x-y).^2), 2);
            
            f = [f1, f2];
            c = [];            
        end
        
        function pf = PF_LL(obj, n, y)
            n2 = n*2;
            if size(y, 1) > 1
                error('LL PF only uses one row vector ');
            end
            s = linspace(0, y(1), n2);
            f1 = s .^ 2;
            f2 = (s - y(1)).^2;
            pf = [f1', f2'];             
            pf = sparse_pf_selection(pf, n);
        end

        function ps = PS_LL(obj, n, xu)
            xl_1 = linspace(0, xu(1), n);
            xl_rest = repmat(xu(2:end), n, 1);

            ps = [xl_1', xl_rest];            
        end
        
        
        function pf = PF_UL(obj, n2)
           
            y1 = [0.001, 0.2, 0.4, 0.6, 0.8, 1];
            F1 = [];
            F2 = [];
            n = n2*100;
            for i = 1:length(y1)
                x = linspace(0, y1(i), n);
                [v1, v2] = obj.trade_off([y1(i)]);
                f1 = v1 - obj.r * cos(obj.gamma * pi * x/(2 *y1(1)));
                f2 = v2 - obj.r * sin(obj.gamma * pi * x/(2 *y1(1)));
 
                F1 = [F1, f1];
                F2 = [F2, f2];
            end
            
            pf_candidate = [F1', F2'];

            
            [fronts, ~, ~] = nd_sort(pf_candidate, (1:size(pf_candidate, 1))');
            nd_ids = fronts(1).f;
            pf = pf_candidate(nd_ids, :);   
            pf = sparse_pf_selection(pf, n2);
        end
        
        function ul_FElimit = return_ulFE(obj)
            ul_FElimit = 51483;
        end
        
        function ll_FElimit = return_llFE(obj)
            ll_FElimit = 1668389;
        end
        
        function literature_hvmean = return_meanHV(obj)
            literature_hvmean = 0.5176;
        end
        
        function literature_hvstd = return_stdHV(obj)
            literature_hvstd = 0.0204;
        end

        function thv = target_hv(obj)
            thv = 0.5176;
        end

        function tigd = target_igd(obj)
             tigd = 0.1087;
        end
         
        
    end   
end
