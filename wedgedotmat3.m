function M = wedgedotmat3(obj)
arguments
    obj ga3
end
I = eye(3);
s = obj.m(1); v = obj.m(2:4); b = obj.m(5:7); p = obj.m(8);
v = v(:); b = b(:);
M = [s v' -b' -p;
    v,  s*I-skew(b), p*I-skew(v),  b;
    b, -p*I+skew(v), s*I-skew(b), -v;
    p -b' -v' s];
end