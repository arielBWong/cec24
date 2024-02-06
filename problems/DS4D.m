classdef DS4D
    properties
        K;
        L;
        name;
        ul_bu;
        ul_bl;
        ll_bu;
        ll_bl;
        nobj = 2;
        ul_nobj = 2;
        ll_nobj = 2;
        ul_con = 1;
        ll_con = 0;
    end
    methods
        function obj = DS4D(K, L)
            
            obj.name = 'DS4D';
            obj.K = K;
            obj.L = L;
            obj.ul_bu = [2];
            obj.ul_bl = [1];
            obj.ll_bu = [1, ones(1, K+L-1) * (K+L) ];
            obj.ll_bl = [0, ones(1, K+L-1) * (K+L) * (-1)];
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            f(:, 1) = (1 - xl(:, 1)) .* (1 + sum(xl(:, 2: obj.K).^2 , 2)) .* (1 - sum(xl(:, obj.K+1 : obj.L + obj.K).^2 , 2)) .* xu(:, 1);
            f(:, 2) =  xl(:, 1)  .* (1 + sum(xl(:, 2: obj.K).^2 , 2)) .* xu(:, 1);
            c = 1 -  (1 - xl(:, 1)) .* xu(:, 1) - 0.5 * xl(:, 1) .* xu(:, 1);
           
        end
        
        function [f, c] = evaluate_l(obj, xu, xl)
            f(:, 1) = (1 - xl(:, 1)) .* (1 + sum(xl(:, obj.K+1 : obj.L + obj.K).^2 , 2)) .* xu(:, 1);
            f(:, 2) =  xl(:, 1) .* (1 + sum(xl(:, obj.K+1 : obj.K + obj.L ).^2 , 2)) .* xu(:, 1);
            c = [];
        end
        
        function pf = PF_LL(obj, n, xu)
            n2 = n * 2;
            f1 = linspace(0, xu(1), n2);
            f2 = xu(1) -  f1;
            pf = [f1', f2'];
            pf = sparse_pf_selection(pf, n);
        end

       function ps = PS_LL(obj, n, xu)
           x1 = linspace(0, xu(1), n);
           xrest = zeros(n, length(obj.ll_bu));
           ps = [x1', xrest];
       end
       
        
        function pf = PF_UL(obj, n)
            n2 = n*2;
            s = linspace(1, 2, n2);
            f1 = 2 - s;
            f2 = 2 *(s-1);
            
            pf = [f1', f2'];
            pf = sparse_pf_selection(pf, n);
           
        end
        
        function ul_FElimit = return_ulFE(obj)
            ul_FElimit = 31887;
        end
        
        function ll_FElimit = return_llFE(obj)
            ll_FElimit = 751712;
        end
        
        
    end
end