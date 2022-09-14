function Q = rgamotor(varargin)
%RGAMOTOR Rigid motion operator for rotation & translation
switch nargin
    case 0
        r = randn(4,1);
        u = randn(4,1);
        r = r - (r'*u)/(u'*u)*u;
    case 2 % r and u
        r = varargin{1};
        u = varargin{2};
    case 3 % phi, d, and L
        phi = varargin{1};
        d = varargin{2}*rga("e0");
        L = varargin{3};
    case 4 % phi, d, v, m
        phi = varargin{1};
        d = varargin{2};
        v = varargin{3};
        m = varargin{4};
        L = rgaline(v,m);
    case 8 % could be r1 r2 r3 r4 u1 u2 u3 u4 or phi d v1 v2 v3 m1 m2 m3
        phi = varargin{1};
        d = varargin{2};
        v = [varargin{3:5}];
        m = [varargin{6:8}];
        L = rgaline(v,m);
end
if exist("L","var")
    Q = L*sin(phi) + rga("e1234")*cos(phi) ...
            + antiwedgedot(d,L)*cos(phi) - d*sin(phi);
else
    if isnumeric(r) && isnumeric(u) && abs(r'*u) > 10*eps
        error('Bulk & Weight must be perpendicular')
    else
        Q = rga([r(1); zeros(4,1); r(1:3); u(1:3); zeros(4,1); u(4)]);
    end
end