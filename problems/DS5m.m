classdef DS5m
    properties
        K;
        L;
        name;
        ul_bu;
        ul_bl;
        ll_bu;
        ll_bl;
        nobj = 2;
    end
    methods
        function obj = DS5m(K, L)
            
            obj.name = 'DS5m';
            obj.K = K;
            obj.L = L;
            obj.ul_bu = [2, ones(1, K-1) * (K+L) ];
            obj.ul_bl = [1, ones(1, K-1) * (K+L) * (-1)];
            obj.ll_bu = [1, ones(1, L) * (K+L) ];
            obj.ll_bl = [0, ones(1, L) * (K+L) * (-1)];
            
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            f(:, 1) = (1 - xl(:, 1)) .* (1 + sum(xu(:, 2: obj.K).^2 , 2)) .* xu(:, 1);
            f(:, 2) =  xl(:, 1)  .* (1 + sum(xu(:, 2: obj.K).^2 , 2)) .* xu(:, 1);
            c = -1 * ((1 - xl(:, 1)) .* xu(:, 1) +  xl(:, 1) .* xu(:, 1) - 2 + 0.2 * floor(5 * ( 1 - xl(:, 1)) .* xu(:, 1) + 1));
            
        end
        
        function [f, c] = evaluate_l(obj, xu, xl)
            f(:, 1) = (1 - xl(:, 1)) .* (1 + sum(xl(:, 2: obj.L+1).^2 , 2)) .* xu(:, 1);
            f(:, 2) =  xl(:, 1) .* (1 + sum(xl(:, 2: obj.L +1 ).^2 , 2)) .* xu(:, 1);
            c = [];
        end
        
        function pf = PF_LL(obj, n, xu)
            n2 = n *2;
            f1 = linspace(0, xu(1), n2);
            f2 = xu(1) -  f1;
            pf = [f1', f2'];
            pf = sparse_pf_selection(pf, n);
        end
        
        function pf = PF_UL(obj, n)
            
            warning('This PF only generate 129 PF points');
            n2 = 1000;
            xu = [1, 1.2, 1.4, 1.6, 1.8];
            
            pf=[];
            for i = 1 : size(xu, 2)
                xl1_expand = [];
                xu_expand = [];
                xli = linspace(2*(1-1/xu(i)), 2*(1-0.9/xu(i)), n2);
                xl1_expand = [xl1_expand; xli'];
                xu_expand = [xu_expand; repmat(xu(i), n2, 1)];
                
                xl2_expand = zeros( size(xl1_expand, 1),(obj.K + obj.L-1));
                xl_expand = [xl1_expand, xl2_expand];
                pfi = obj.evaluate_u(xu_expand, xl_expand);
                pfi = sparse_pf_selection(pfi, 33);
                pf = [pf; pfi];
            end
        end
        
        
        function ul_FElimit = return_ulFE(obj)
            ul_FElimit = 35947;
        end
        
        function ll_FElimit = return_llFE(obj)
            ll_FElimit = 902194;
        end
        
        
    end
end