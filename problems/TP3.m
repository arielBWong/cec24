classdef TP3
    properties
        name;
        ul_bu;
        ul_bl;
        ll_bu;
        ll_bl;
        nobj = 2;
    end
    
    methods
        function obj = TP3()
            obj.ul_bu = [10];
            obj.ul_bl = [0];
            obj.ll_bu = [10, 10];
            obj.ll_bl = [0, 0];
            obj.name = 'TP3';
        end
        
        function [f, c] = evaluate_u(obj, xu, xl)
            f(:, 1) = xl(:, 1) + xl(:, 2) .^ 2 + xu(:, 1) + (sin(xu(:, 1) + xl(:, 1))) .^ 2;
            f(:, 2) = cos(xl(:, 2)) .* (0.1 + xu(:, 1)) .* exp(-(xl(:, 1)./(0.1 + xl(:, 2))));
            c(:, 1) = (xl(:, 1) - 0.5) .^ 2 + (xl(:, 2) - 5).^2 + (xu(:, 1) - 5).^2 - 16;
        end
        
        function [f, c] = evaluate_l(obj, xu, xl)
            f(:,1)= ((xl(:, 1) - 2).^2 + (xl(:, 2) - 1) .^ 2)/4 + ...
                (xl(:, 2) .* xu(:, 1) + (5 - xu(:, 1)).^ 2) / 16 + sin(xl(:, 2)/10);
            
            f(:,2)= (xl(:, 1).^2 + (xl(:, 2) - 6) .^ 4 - 2 * xl(:, 1) .* xu(:, 1) - (5 - xu(:, 1)) .^ 2) / 80;
            
            c(:,1) = xl(:, 1).^2 - xl(:, 2);
            c(:,2) = xl(:, 2) + 5 * xl(:, 1) .^ 2 - 10;
            c(:,3) = xl(:, 2) + xu(:, 1)/6 - 5;
            c(:,4) = -xl(:, 1);
        end
        
        function pf = PF_UL(obj, n)
            fprintf('no theoretical PF \n');
            pf = [];
        end
        
        function pf = PF_LL(obj, n, y)
            fprintf('no theoretical PF \n');
            pf = [];
        end
        
        function ul_FElimit = return_ulFE(obj)
            ul_FElimit = 15550;
        end
        
        function ll_FElimit = return_llFE(obj)
            ll_FElimit = 625862;
        end
        
        
    end
end