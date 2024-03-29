classdef rgamotor < rga
    %RGAMOTOR RGA Motor: Rigid motion operator for rotation & translation
    % Subclass of RGA
    % In RGA a motor Q is an object such that Q*x*~Q, where * indicates
    % antiwedgedot, and ~ indicates antireverse, rotates & translates the
    % RGA object x though angle 2\phi about axis L and distance 2d along L.
    % L can be any unitized line, i.e. it does not have to pass through the
    % origin.  A motor has 6 bivector elements along with a scalar and a 
    % pseudoscalar element.  By analogy with dual quaternion, four of the 
    % elements correspond to the real part, and four to the dual part.
    % The real and dual parts must be perpendicular, and the real part must
    % have unit norm for a well-posed rotation & translation operation.
    % Unitization ensures these contraints are met.
    methods
        function obj = rgamotor(varargin)
            %RGAMOTOR Construct rigid motion operator for rotation & translation
            %  Q = rgamotor creates a random RGA motor object
            %  Q = rgamotor(M) subclasses the RGA object M into a motor
            %  Q = rgamotor(r,u) creates a motor from dual quaternion r + eps*u
            %  Q = rgamotor(phi,d,L) uses the angle phi, distance d, and axis L
            %  Q = rgamotor(phi,d,v,m) uses the angle phi, distance d, and
            %      axis direction & moment, v & m
            %  Q = rgamotor(phi,d,vx,vy,vz,mx,my,mz) uses the angle phi,
            %      distance d, and axis direction & moment components
            % Defaults to anti basis.
            switch nargin
                case 0
                    r = randn(4,1);
                    u = randn(4,1);
                    u = u - (u'*r)/(r'*r)*r;
                case 1 % input is a multivector
                    %obj = varargin{1}; % doesn't preserve class
                    obj.m = varargin{1}.m;
                    obj.anti = varargin{1}.anti;
                    obj.dispeps = varargin{1}.dispeps;
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
                obj.anti = true;
            elseif exist("r","var") && exist("u","var")
                r = r(:); u = u(:);
                if isnumeric(r) && isnumeric(u) && abs(r'*u) > 100*eps
                    error('Bulk & Weight must be perpendicular')
                else
                    obj.m = [u(4); zeros(4,1); u(1:3); r([3 2 1]); zeros(4,1); r(4)];
                    obj.anti = true;
                end
            end
        end

        function obj = unitize(obj)
            %UNITIZE Unitize the motor
            obj.m([9:11 16]) = obj.m([9:11 16])/norm(obj.m([9:11 16]));
            if obj.m(16) < 0
                %obj.m([9:11 16]) = -obj.m([9:11 16]);
                obj.m = -obj.m;
                %warning('motor orientation sign swapped')
            end
            obj.m([6:8 1]) = obj.m([6:8 1]) ...
                - dot(obj.m([6:8 1]),obj.m([11 10 9 16])) ...
                / dot(obj.m([11 10 9 16]),obj.m([11 10 9 16])) ...
                * obj.m([11 10 9 16]);
            obj = rgamotor(obj);
        end

        function varargout = extract(obj)
            %EXTRACT Extract motor parameters from motor
            % [phi,d,v,m] = extract(obj)
            % [phi,d,L] = extract(obj)
            % where the motor rotates through 2*phi and translates through
            % 2*d along line L defined by direction v and moment m.
            if abs(obj.e1234) == 1
                phi = 0;
                t = [obj.e23 obj.e31 obj.e12];
                d = norm(t);
                if d == 0
                    v = zeros(1,3);
                else
                    v = t/d;
                end
                m = zeros(1,3);
            else
                phi = acos(obj.e1234);
                sph = sin(phi);
                d = -(obj.e0/sph);
                v = [obj.e41 obj.e42 obj.e43]/sph;
                m = ([obj.e23 obj.e31 obj.e12] - d*v*cos(phi))/sph;
            end
            varargout{1} = phi;
            varargout{2} = d;
            if nargout == 3
                L = rgaline(v,m);
                varargout{3} = L;
            else
                varargout{3} = v;
                varargout{4} = m;
            end
        end

        function h = plot(obj)
            %PLOT Plot the motor as a coordinate frame
            o = screw(obj,rgapoint(0,0,0,1));
            xax = screw(obj,rgapoint(1,0,0,1))-o;
            yax = screw(obj,rgapoint(0,1,0,1))-o;
            zax = screw(obj,rgapoint(0,0,1,1))-o;
            hold on
            ho = plot3(o.m(2),o.m(3),o.m(4),'k*');
            hx = quiver3(o.m(2),o.m(3),o.m(4),xax.m(2),xax.m(3),xax.m(4),'r');
            hy = quiver3(o.m(2),o.m(3),o.m(4),yax.m(2),yax.m(3),yax.m(4),'g');
            hz = quiver3(o.m(2),o.m(3),o.m(4),zax.m(2),zax.m(3),zax.m(4),'b');
            hy.Color = [0 .5 0];
            h = [ho;hx;hy;hz];
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