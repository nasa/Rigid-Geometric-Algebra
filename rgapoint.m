function p = rgapoint(x,y,z,w)
switch nargin
    case 0
        p = rga([0 randn(1,4) zeros(1,11)]);
    case 1
        p = rga([0 x(:)' zeros(1,11)]);
    case 3
        p = rga([0 x y z 1 zeros(1,11)]);
    case 4
        p = rga([0 x y z w zeros(1,11)]);
    otherwise
        error('Inputs not compatible')
end