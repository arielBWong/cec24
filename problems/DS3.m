classdef DS3
    properties
        K;
        name;
        ul_bu;
        ul_bl;
        ll_bu;
        ll_bl;
        r = 0.2;
        tau = 1;
        ul_nobj = 2;
        ll_nobj = 2;
        ul_con = 1;
        ll_con = 1;
    end
    methods
        function obj = DS3(K, tau)
            

            obj.K = K;
            obj.tau = tau;
            if obj.tau > 0
                obj.name = 'DS3';
            else
                obj.name = 'DS3D';
            end
            obj.ul_bu = ones(1, K) .* K;
            obj.ul_bl = zeros(1, K);
            obj.ll_bu = ones(1, K) * K;
            obj.ll_bl = ones(1, K) * (-K);            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            %xu1 should discrete, a multiple of 0.1
            xu(:,1)=round(xu(:,1),1);
            
            R = 0.1 + 0.15 * abs(sin(2 * pi * (xu(:,1) - 0.1))); 
            J = (3 : obj.K) / 2 ;
            
            f(:, 1) = xu(:, 1)+ sum((xu(:,3: obj.K)- J).^ 2, 2) + obj.tau * sum((xl(:, 3: obj.K) - xu(:, 3:obj.K)).^2, 2) ...
                - R .* cos(4 * atan((xu(:, 2) - xl(:, 2)) ./ (xu(:, 1) - xl(:, 1))));
            f(:, 2) = xu(:, 2) + sum((xu(:, 3:obj.K) - J).^2, 2) + obj.tau * sum((xl(:, 3: obj.K) - xu(:, 3:obj.K)).^2, 2)...
                - R .* sin(4 * atan((xu(:, 2) - xl(:, 2)) ./ (xu(:, 1) - xl(:, 1))));
            c = 1 - xu(:, 1).^2 - xu(:, 2); 
        end
        
        
        function [f, c] = evaluate_l(obj, xu, xl)
            
            f(:,1) = xl(:, 1) + sum((xl(:, 3:obj.K) - xu(:, 3:obj.K)).^2, 2);
            f(:,2) = xl(:, 2) + sum((xl(:, 3:obj.K) - xu(:, 3:obj.K)).^2, 2);
            c = (xl(:, 1) - xu(:, 1)) .^ 2 + (xl(:, 2) - xu(:, 2)) .^2 - obj.r^2;
        end
        
        function pf = PF_UL(obj, n2)
            
            n = n2 * 2;
            xu1 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4];
            xu2 = 1 - xu1 .^2;
            
            xu2_binary = xu2 < 0;
            xu2(xu2_binary) = 0;
            
            R = 0.1 + 0.15 * abs(sin(2 * pi * (xu1 - 0.1)));
            
            angles = linspace(pi, pi * 3/2, n);
            
            pf = [];
            for i = 1:length(xu1)
                x1 = obj.r * cos(angles) + xu1(i);
                x2 = obj.r * sin(angles) + xu2(i);
                
                f1 = xu1(i) - R(i) * cos(4 * atan((xu2(i) - x2) ./(xu1(i) - x1)));
                f2 = xu2(i) - R(i) * sin(4 * atan((xu2(i) - x2) ./(xu1(i) - x1)));
                
                pf = [pf; [f1', f2']];
            end
            
            [fronts, ~, ~] = nd_sort(pf, (1:size(pf, 1))');
            pf = pf(fronts(1).f, :);
            % pf = mirror_ref_selection(pf, n2);
            pf = sparse_pf_selection(pf, n2);
        end


        
        function pf = PF_LL(obj, n, xu)
            n2 = n * 2;
            angles = linspace(pi, pi * 3/2, n2);
            f1 = cos(angles) * obj.r + xu(1);
            f2 = sin(angles) * obj.r + xu(2);
            pf = [f1', f2'];
            
            pf = sparse_pf_selection(pf, n);
        end

        function ps = PS_LL(obj, n, xu)   
            try
                angles = linspace(pi, pi * 3/2, n);
                xl_1 = cos(angles) * obj.r + xu(1);
                xl_2 = sin(angles) * obj.r + xu(2);
                if size(xu, 2) > 2
                    xl_rest = repmat(xu(3: end), n, 1);
                else
                    xl_rest = [];
                end
                ps = [xl_1', xl_2', xl_rest];
            catch e
                fprintf(1,'The identifier was:\n%s',e.identifier);
                fprintf(1,'There was an error! The message was:\n%s',e.message);
                

            end
        end
        
        function ul_FElimit = return_ulFE(obj)
            ul_FElimit = 52149;
        end
        
        function ll_FElimit = return_llFE(obj)
            ll_FElimit = 1543789;
        end
        
        function literature_hvmean = return_meanHV(obj)
            literature_hvmean = 1.032;
        end
        
        function literature_hvstd = return_stdHV(obj)
            literature_hvstd = 0.0336;
        end

         function thv = target_hv(obj)
            thv = 1.032;
        end

        function tigd = target_igd(obj)
             tigd = 0.1667;
        end
        
    end
end