classdef rgaline < rga
    %RGALINE RGA Line
    % Subclass of RGA
    methods
        function obj = rgaline(varargin)
            %RGALINE Line specified by direction & moment vectors
            if nargin == 1 && isa(varargin{1},'rga')
                b = bivector(varargin{1});
                obj.m = b.m;
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
end