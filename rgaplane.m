function f = rgaplane(x,y,z,w)
switch nargin
    case 0
        f = rga(fliplr([0 randn(1,4) zeros(1,11)]));
    case 1
        f = rga(fliplr([0 x(:)' zeros(1,11)]));
    case 3
        f = rga(fliplr([0 x y z 1 zeros(1,11)]));
    case 4
        f = rga(fliplr([0 x y z w zeros(1,11)]));
    otherwise
        error('Inputs not compatible')
end
f.anti = true;
end