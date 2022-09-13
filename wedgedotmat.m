function M = wedgedotmat(obj)
% WEDGEDOTMAT Matrix for wedgedot of two multivectors
I = eye(3); J = fliplr(I);
weks = @(x) fliplr(skew(x));
s = obj.m(1);
v3 = obj.m(2:4); v3f = flipud(v3(:)); v4 = obj.m(5);
b3 = obj.m(6:8); b3f = flipud(b3(:));
b4 = obj.m(9:11); b4f = flipud(b4(:));
t = obj.m(12:15);
t3 = obj.m(12);
t4 = obj.m(13:15); t4f = flipud(t4(:));
p = obj.m(16);
v3 = v3(:); b3 = b3(:); b4 = b4(:); t4 = t4(:); t = t(:);
O33 = zeros(3,3); O31 = zeros(3,1); O13 = O31';
M = [s,            v3',    0,          -b3',            O13, -t3,             O13,   0;
    v3,   s*I-skew(b3),  O31, t3*I-skew(v3),            O33,  b3,             O33, O31;
    v4,           b4f',    s,         -t4f',          -v3f',  -p,           -b3f',  t3;
    b3, -t3*I+skew(v3),  O31,  s*I-skew(b3),            O33, -v3,             O33, O31;
    b4,  v4*J+weks(t4), -v3f,  p*J+weks(b4),  s*I-skew(b3f)', t4, -t3*I-skew(v3f), b3f;
    t3,           -b3',    0,          -v3',            O13,   s,             O13,   0;
    t4,  -p*J-weks(b4),  b3f, v4*J+weks(t4), t3*I-skew(v3f)',-b4,   s*I+skew(b3f), v3f;
    p,           -t4f',  -t3,         -b4f',          -b3f',  v4,            v3f',   s];
end