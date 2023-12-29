classdef TP2
    properties
        name;
        K;
        ul_bu;
        ul_bl;
        ll_bu;
        ll_bl;
        nobj = 2;
        ul_nobj = 2;
        ll_nobj = 2;
        ul_con = 0;
        ll_con = 0;
    end
    
    methods
        function obj = TP2(K)
            obj.name = 'TP2';
            obj.K = K;
            obj.ul_bl = [-1];
            obj.ul_bu = [2];
            obj.ll_bl = ones(1, obj.K) * (-1);
            obj.ll_bu = ones(1, obj.K) * 2; 
        end
        
        function [f, c] = evaluate_u(obj, y, x)
            f(:, 1) = (x(:, 1) - 1) .^2 + sum(x(:, 2:end) .^2, 2) + y.^2;
            f(:, 2) = (x(:, 1) - 1) .^2 + sum(x(:, 2:end) .^2, 2) + (y - 1).^2;
            c = [];
        end
 
        
        function [f, c] = evaluate_l(obj, y, x)
            f(:,1) = sum(x.^2, 2);
            f(:,2) = (x(:,1) - y).^2 + sum(x(:, 2:end) .^ 2, 2);
            c=[];
        end
        
        function pf = PF_UL(obj, n)
            n2 = n * 100;
            y = linspace(0.5, 1, n2);
            
            x1 = y;
            x_rest = zeros(n2, obj.K-1);
            x = [x1', x_rest];
            
            pf = obj.evaluate_u(y', x);
            
            figure(1)
            plot(pf(:, 1), pf(:, 2), 'k.'); hold on;
            pf = sparse_pf_selection(pf, n);
            plot(pf(:, 1), pf(:, 2), 'ro');
            close();
            
        end
        
        function pf = PF_LL(obj, n, y)
            n2 = n * 2;
            x1 = linspace(0, y, n2);
            x_rest = zeros(n2, obj.K-1);
            x = [x1', x_rest];
            y = repmat(y, n2, 1);
            
            pf = obj.evaluate_l(y, x);
            pf = sparse_pf_selection(pf, n);
        end  

        function ps = PS_LL(obj, n, y)
            x1 = linspace(0, y, n);
            x_rest = zeros(n, obj.K-1);
            ps = [x1', x_rest];
        end
        
        
        function ul_FElimit = return_ulFE(obj)
            ul_FElimit = 11867;
        end
        
        function ll_FElimit = return_llFE(obj)
            ll_FElimit = 152414;
        end
        
        function literature_hvmean = return_meanHV(obj)
            literature_hvmean = 0.2264;
        end
        
        function literature_hvstd = return_stdHV(obj)
            literature_hvstd = 0.00368;
        end

        function thv = target_hv(obj)
            thv = 0.2264 ;
        end

        function tigd = target_igd(obj)
             tigd = 0.02494;
        end
    end
end