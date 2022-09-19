classdef rgaline < rga
    %RGALINE RGA Line
    % Subclass of RGA
    methods
        function obj = rgaline(varargin)
            %RGALINE Line specified by direction & moment vectors
            switch nargin
                case 0
                    v = randn(3,1);
                    m = randn(3,1);
                    v = v - (v'*m)/(m'*m)*m;
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
    end
end