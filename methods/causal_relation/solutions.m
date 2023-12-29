classdef solutions < handle
    properties
        xu = [];
        XL = [];
        XU = [];
        XUE = [];
        FU = [];
        FL = [];
        FC = [];
        FLC = [];
    end

    methods
        function obj = solutions()
            % this class has dependent on pop_sort method
            obj.XU = [];
            obj.XL = [];
            obj.xu = [];
            obj.XUE = [];
            obj.FU = [];
            obj.FL = [];
            obj.FC = [];
            obj.FLC = [];
        end

        function add(obj, xu, XL, FU, FL, FC, FLC)
            % this forms list of cell array
            obj.xu = [obj.xu, {xu}];
            obj.XL = [obj.XL, {XL}];

            if isempty(XL)  % single level problem
                obj.XU = obj.xu;
                obj.XUE = [];
            else
                obj.XU = [obj.XU, {repmat(xu, size(XL, 1), 1)}];
                obj.XUE = [obj.XUE, {[repmat(xu, size(XL, 1), 1), linspace(0, 1, size(XL, 1))']}];
            end
            obj.FU = [obj.FU, {FU}];
            obj.FC = [obj.FC, {FC}];
            obj.FL = [obj.FL, {FL}];
            obj.FLC = [obj.FLC, {FLC}];
        end

        function merge(obj,  new_solutions)

            obj.xu = [obj.xu, new_solutions.xu];
            obj.XL = [obj.XL, new_solutions.XL];
            obj.XU = [obj.XU, new_solutions.XU];
            obj.FU = [obj.FU, new_solutions.FU];
            obj.FL = [obj.FL, new_solutions.FL];
            obj.FC = [obj.FC, new_solutions.FC];
            obj.FLC = [obj.FLC, new_solutions.FLC];
            obj.XUE = [obj.XUE, new_solutions.XUE];
        end

        function copy(obj, another_solutions)
            if ~isempty(obj.xu)
                obj.XU = [];
                obj.XL = [];
                obj.xu = [];
                obj.XUE = [];
                obj.FU = [];
                obj.FL = [];
                obj.FC = [];
                obj.FLC = [];
            end
            obj.merge(another_solutions);
        end

        function nd_sort(obj, varargin)
            % convert cell array to matrix
            tmpFU = obj.FUs;
            tmpFL = obj.FLs;
            tmpFC = obj.FCs;
            tmpFLC = obj.FLCs;
            tmpXU = obj.XUs;
            tmpXL = obj.XLs;

            % ND sort wrapper
            pop.X = tmpXU;
            pop.F = tmpFU;
            pop.C = tmpFC;
            pop.dummy1 = tmpXL;
            pop.dummy2 = tmpFLC;
            pop.dummy3 = tmpFL;

            [pop, front_idx] = pop_sort(pop);
            nd_front_id = front_idx == 1;

            % keep ND
            ndXU = pop.X(nd_front_id, :);
            if isempty(varargin) 
                ndFU = pop.F(nd_front_id, :);
                if ~isempty(pop.C)
                    ndFC = pop.C(nd_front_id, :);
                else
                    ndFC = [];
                end

                if ~isempty(pop.dummy1)
                    ndXL = pop.dummy1(nd_front_id, :);
                else
                    ndXL = [];
                end

                if ~isempty(pop.dummy3)
                    ndFL = pop.dummy3(nd_front_id, :);
                else
                    ndFL = [];
                end

                if ~isempty(pop.dummy2)
                    ndFLC = pop.dummy2(nd_front_id, :);
                else
                    ndFLC = [];
                end

                % distribute back to object variables
                unique_ndXU = unique(ndXU, 'rows', 'stable');
                num_ndxu = size(unique_ndXU, 1);
                obj.clear_data;
                for ii = 1:num_ndxu
                    tmp_xu = unique_ndXU(ii, :);
                    ia = ismember(ndXU, tmp_xu, 'rows');
                    %xu, XL, FU, FL, FC, FLC
                    one_xu = tmp_xu;

                    if ~isempty(ndXL)
                        one_XL = ndXL(ia, :);
                    else
                        one_XL = [];
                    end

                    one_FU = ndFU(ia, :);

                    if ~isempty(ndFL)
                        one_FL = ndFL(ia,:);
                        
                    else
                        one_FL = [];
                    end

                    if ~isempty(ndFC)
                        one_FC = ndFC(ia,:);
                    else
                        one_FC = [];
                    end

                    if ~isempty(ndFLC)
                        one_FLC = ndFLC(ia, :);
                    else
                        one_FLC = [];
                    end

                    obj.add(one_xu, one_XL, one_FU, one_FL, one_FC, one_FLC);
                end
            else
                unique_ndXU = unique(ndXU, 'rows', 'stable');
                num_ndxu = size(unique_ndXU, 1);
                obj.clear_data;

                for ii = 1:num_ndxu
                    tmp_xu = unique_ndXU(ii, :);
                    ia = ismember(tmpXU, tmp_xu, 'rows');

                    one_xu = tmp_xu;
                     if ~isempty(tmpXL)
                        one_XL = tmpXL(ia, :);
                    else
                        one_XL = [];
                    end

                    one_FU = tmpFU(ia, :);

                    if ~isempty(tmpFL)
                        one_FL = tmpFL(ia,:);
                    else
                        one_FL = [];
                    end

                    if ~isempty(tmpFC)
                        one_FC = tmpFC(ia,:);
                    else
                        one_FC = [];
                    end

                    if ~isempty(tmpFLC)
                        one_FLC = tmpFLC(ia, :);
                    else
                        one_FLC = [];
                    end
                    obj.add(one_xu, one_XL, one_FU, one_FL, one_FC, one_FLC);
                end
            end
        end



        function FU =  FUs(obj)
            FU = cat(1, obj.FU{:});
        end

        function FL = FLs(obj)
            FL = cat(1, obj.FL{:});
        end

        function FC = FCs(obj)
            FC = cat(1, obj.FC{:});
        end

        function FLC = FLCs(obj)
            FLC = cat(1, obj.FLC{:});
        end

        function XU = XUs(obj)
            XU = cat(1, obj.XU{:});
        end

        function xu = xus(obj)
            xu = cat(1, obj.xu{:});
        end

        function XL = XLs(obj)
            XL = cat(1, obj.XL{:});
        end


        function clear_data(obj)
            obj.XU = [];
            obj.XL = [];
            obj.xu = [];
            obj.XUE = [];
            obj.FU = [];
            obj.FL = [];
            obj.FC = [];
            obj.FLC = [];
        end

    end
end