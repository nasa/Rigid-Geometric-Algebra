function L = rgaline(varargin)
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
    L = rga([zeros(5,1); m; v([3 2 1]); zeros(5,1)]);
end