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
                    if isa(x,'rga')
                        v = vector(x);
                        obj.m = v.m;
                    else
                        obj.m = [0 x(:)' zeros(1,11)];
                    end
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

        function h = plot(obj)
            %PLOT Plot the point at w=1
            obj = unitize(obj); % not really needed
            h = plot3(obj.m(2),obj.m(3),obj.m(4),'*');
        end
    end

    methods(Access=protected)
        function obj = dotAssign(obj,indexOp,c)
            switch indexOp.Name
                case {'e0','eps1234'}
                    obj.m(1) = 0; warning('input ignored')
                case {'e1','eps423'}
                    obj.m(2) = c;
                case {'e2','eps431'}
                    obj.m(3) = c;
                case {'e3','eps412'}
                    obj.m(4) = c;
                case {'e4','eps321'}
                    obj.m(5) = 0; warning('input ignored')
                case {'e23','eps41'}
                    obj.m(6) = 0; warning('input ignored')
                case {'e31','eps42'}
                    obj.m(7) = 0; warning('input ignored')
                case {'e12','eps43'}
                    obj.m(8) = 0; warning('input ignored')
                case {'e43','eps12'}
                    obj.m(9) = 0; warning('input ignored')
                case {'e42','eps31'}
                    obj.m(10) = 0; warning('input ignored')
                case {'e41','eps23'}
                    obj.m(11) = 0; warning('input ignored')
                case {'e321','eps4'}
                    obj.m(12) = 0; warning('input ignored')
                case {'e412','eps3'}
                    obj.m(13) = 0; warning('input ignored')
                case {'e431','eps2'}
                    obj.m(14) = 0; warning('input ignored')
                case {'e423','eps1'}
                    obj.m(15) = 0; warning('input ignored')
                case {'e1234','eps0'}
                    obj.m(16) = 0; warning('input ignored')
                otherwise
                    error('basis not recognized')
            end
        end
    end
end