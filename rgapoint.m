classdef rgapoint < rga
    %RGAPOINT RGA Point
    %   Subclass of RGA
    methods
        function obj = rgapoint(x,y,z,w)
            %RGAPOINT Create RGA Point object
            switch nargin
                case 0
                    obj.m = [0 randn(1,4) zeros(1,11)];
                case 1
                    obj.m = [0 x(:)' zeros(1,11)];
                case 3
                    obj.m = [0 x y z 1 zeros(1,11)];
                case 4
                    obj.m = [0 x y z w zeros(1,11)];
                otherwise
                    error('Inputs not compatible')
            end
        end

        function obj = unitize(obj)
            %UNITIZE Unitize the point
            obj.m(5) = 1;
        end
    end
end