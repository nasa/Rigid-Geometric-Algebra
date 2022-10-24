classdef rgamotor < rga
    %RGAMOTOR RGA Motor: Rigid motion operator for rotation & translation
    % Subclass of RGA
    methods
        function obj = rgamotor(varargin)
            %RGAMOTOR Construct rigid motion operator for rotation & translation
            switch nargin
                case 0
                    r = randn(4,1);
                    u = randn(3,1);
                    u(4) = -1/r(4)*(dot(r(1:3),u(1:3)));
                case 1 % input is a multivector
                    obj.m = varargin{1}.m;
                case 2 % r and u
                    r = varargin{1};
                    u = varargin{2};
                case 3 % phi, d, and L
                    phi = varargin{1};
                    d = varargin{2}*rga("e0");
                    L = varargin{3};
                case 4 % phi, d, v, m
                    phi = varargin{1};
                    d = varargin{2}*rga("e0");
                    v = varargin{3};
                    m = varargin{4};
                    L = rgaline(v,m);
                case 8 % phi d v1 v2 v3 m1 m2 m3
                    phi = varargin{1};
                    d = varargin{2}*rga("e0");
                    v = [varargin{3:5}];
                    m = [varargin{6:8}];
                    L = rgaline(v,m);
            end
            if exist("L","var")
                Q = L*sin(phi) + rga("e1234")*cos(phi) ...
                    + antiwedgedot(d,L)*cos(phi) - d*sin(phi);
                obj.m = Q.m;
            elseif exist("r","var") && exist("u","var")
                r = r(:); u = u(:);
                if isnumeric(r) && isnumeric(u) && abs(r'*u) > 10*eps
                    error('Bulk & Weight must be perpendicular')
                else
                    obj.m = [u(4); zeros(4,1); u(1:3); r([3 2 1]); zeros(4,1); r(4)];
                end
            end
        end

        function obj = unitize(obj)
            %UNITIZE Unitize the motor
            obj.m([9:11 16]) = obj.m([9:11 16])/norm(obj.m([9:11 16]));
            %obj.m(1) = -1/obj.m(16)*dot(obj.m([11 10 9]),obj.m(6:8));
            %obj = obj - proj(bulk(obj),weight(obj));
            obj = rgamotor(obj);
        end
    end

    methods(Access=protected)
        function obj = dotAssign(obj,indexOp,c)
            switch indexOp.Name
                case {'e0','eps1234'}
                    obj.m(1) = c;
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
                    obj.m(16) = c;
                otherwise
                    error('basis not recognized')
            end
        end
    end

end