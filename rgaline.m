classdef rgaline < rga
    %RGALINE RGA Line
    % Subclass of RGA
    methods
        function obj = rgaline(varargin)
            %RGALINE Line specified by direction & moment vectors
            if nargin == 1 && isa(varargin{1},'rga')
                b = bivector(varargin{1});
                obj.m = b.m;
                obj.anti = varargin{1}.anti;
                obj.dispeps = varargin{1}.dispeps;
            else
                switch nargin
                    case 0
                        v = randn(3,1);
                        m = randn(3,1);
                        v = v - (v'*m)/(m'*m)*m;
                    case 1
                        v = varargin{1}([11 10 9]);
                        m = varargin{1}(6:8);
                    case 2
                        v = varargin{1};
                        m = varargin{2};
                    case 6
                        v = [varargin{1:3}];
                        m = [varargin{4:6}];
                    otherwise
                        error('Inputs not compatible')
                end
                v = v(:); m = m(:);
                if isnumeric(v) && isnumeric(m) && abs(v'*m) > 10*eps
                    error('Direction & Moment must be perpendicular')
                else
                    obj.m = [zeros(5,1); m; v([3 2 1]); zeros(5,1)];
                end
            end
        end

        function obj = unitize(obj)
            %UNITIZE Unitize the line
            obj.m(9:11) = obj.m(9:11)/norm(obj.m(9:11));
            obj.m(6:8) = obj.m(6:8) ...
                - dot(obj.m(6:8),obj.m([11 10 9])) ...
                / dot(obj.m([11 10 9]),obj.m([11 10 9])) ...
                * obj.m([11 10 9]);
        end

        function obj = direction(obj)
            %DIRECTION Return direction part
            obj.m([1:8 12:end]) = 0;
        end

        function obj = moment(obj)
            %MOMENT Return moment part
            obj.m([1:5 9:end]) = 0;
        end

        function h = plot(obj)
            %PLOT Plot the line
            % Project origin onto line, then use quiver
            o = rgapoint(0,0,0);
            poL = antiwedge(weightlc(obj)^o,obj);
            h = quiver3(poL.m(2)/poL.m(5),poL.m(3)/poL.m(5),poL.m(4)/poL.m(5),...
                obj.m(11),obj.m(10),obj.m(9));
        end
    end

    methods(Access=protected)
        function obj = dotAssign(obj,indexOp,c)
            switch indexOp.Name
                case {'e0','eps1234'}
                    obj.m(1) = 0; warning('input ignored')
                case {'e1','eps423'}
                    obj.m(2) = 0; warning('input ignored')
                case {'e2','eps431'}
                    obj.m(3) = 0; warning('input ignored')
                case {'e3','eps412'}
                    obj.m(4) = 0; warning('input ignored')
                case {'e4','eps321'}
                    obj.m(5) = 0; warning('input ignored')
                case {'e23','eps41'}
                    obj.m(6) = c;
                case {'e31','eps42'}
                    obj.m(7) = c;
                case {'e12','eps43'}
                    obj.m(8) = c;
                case {'e43','eps12'}
                    obj.m(9) = c;
                case {'e42','eps31'}
                    obj.m(10) = c;
                case {'e41','eps23'}
                    obj.m(11) = c;
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