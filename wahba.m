function Q = wahba(M,N)
%WAHBA Q = wahba(M,N) finds motor Q such that N = Q*M*~Q
% M and N are corresonding arrays of multivectors containing any
% combination of points, lines, planes, and/or motors.
% If all objects in M & N are the same multivector type, e.g. they are all
% planes, all points, etc. then M & N can be ordinary arrays.  If there are
% mixed types within M & N, then M & N must be cell arrays.
% This function requires the rga class.

if nargin < 1
    wahba_test
    return
end
lenM = length(M);
C = 0; Cw = 0;
for i = 1:lenM
    B = Psi(N{i}.m) - Xi(M{i}.m);
    Bw = Psi(weight(N{i}).m) - Xi(weight(M{i}).m);
    C = C + B'*B;
    Cw = Cw + Bw'*Bw;
end

%% Solve for rotational part of motor
Hw(16,4) = 1; Hw(11:-1:9,1:3) = eye(3);
Aw = Hw'*Cw*Hw;
[~,Sw,Vw] = svd(Aw,'vector');
tol = max(size(Aw))*eps(norm(Aw));
k = find(Sw<tol);
if length(k) > 1
     warning('Rotational solution may be inaccurate')
end
qw = Hw*Vw(:,end); % min since svd returns in decreasing order
if qw(end)<0
    qw = -qw; % conventional to make scalar part positive
end

%% Solve for remaining part of motor
Hb(16,4) = 0; Hb([6:8 1],:) = eye(4); % s = Hb'*qb
%Cb = C*Hb;
qwhat = Hw'*qw;
A = [Hb'*C*Hb qwhat; qwhat' 0];
d = [-Hb'*C*Hw*qwhat; 0];
%A = [Cb'*Cb qwhat; qwhat' 0];
%d = [-Hb'*C'*C*Hw*qwhat; 0];
x = A\d;
qb = Hb*x(1:4);

%% Create output motor
Q = unitize(rgamotor(qw([11 10 9 16]),qb([6:8 1])));
end

%% Helper functions

function A = Psi(m)
% Psi(Q.m) is the productmat such that Psi(Q.m)*x.m' = Q*x, where * is the
% antiwedgedot product.
arguments
    m (1,16)
end
I = eye(3); J = fliplr(I);
skew = @(x) [0,-x(3),x(2);x(3),0,-x(1);-x(2),x(1),0];
weks = @(x) fliplr(skew(x));
s = m(1);
v3 = m(2:4); v3f = flipud(v3(:)); v4 = m(5);
b3 = m(6:8); b3f = flipud(b3(:));
b4 = m(9:11); b4f = flipud(b4(:));
t3 = m(12);
t4 = m(13:15); t4f = flipud(t4(:));
p = m(16);
v3 = v3(:); b3 = b3(:); b4 = b4(:); t4 = t4(:);
O33 = zeros(3,3); O31 = zeros(3,1); O13 = O31';
A = [p,   -t4f',          -t3, -b4f',           -b3f',           v4,   v3f',          s;
    -t4f, p*I+skew(b4f),   b3, -v4*I-skew(t4f),  t3*J+weks(v3),  b4f,  s*J+weks(b3),  v3;
    0,    O13,             p,   O13,            -t4',            0,   -b4',           v4;
    b4f,  v4*I+skew(t4f), -v3,  p*I+skew(b4f),   s*J+weks(b3),   t4f, -t3*J-weks(v3), b3;
    O31,  O33,            -t4,  O33,             p*I-skew(b4),   O31, -v4*I+skew(t4), b4;
    -v4,  -b4f',           s,   t4f',           -v3f',           p,   -b3f',          t3;
    O31,  O33,             b4,  O33,             v4*I-skew(t4),  O31,  p*I-skew(b4),  t4;
    0,    O13,            -v4,  O13,            -b4',            0,    t4',           p];
end

function B = Xi(m,tilde)
% Xi(Q.m) is commutor productmat s/t Xi(Q.m)*x.m' is equal to x*Q where *
% is the antiwedgedot product.
% Xi(Q.m,true) applies an antireverse to the input so that
% Xi(Q.m,true)*x.m' = Xi(antirev(Q).m)*x.m' = x*~Q;
% Note that Xi(Q.m)' = rga.productmat(lcomp(Q),'wedgedot')
arguments
    m (1,16)
    tilde logical = false
end
if ~tilde
    % This logic seems b/ward b/c B below was originally def'd for the case
    % when a single call to rga.productmat returned both Psi & Xi and the
    % use case was sandwich products Q*x*~Q = Psi(Q.m)*Xi(Q.m)*x.m, ie
    % there was an implicit assumption that input had already had
    % antireverse applied to it.  So we need to take out here:
    m = diag([1 -ones(1,10), ones(1,5)])*m';
end
I = eye(3); J = fliplr(I);
skew = @(x) [0,-x(3),x(2);x(3),0,-x(1);-x(2),x(1),0];
weks = @(x) fliplr(skew(x));
s = m(1);
v3 = m(2:4); v3f = flipud(v3(:)); v4 = m(5);
b3 = m(6:8); b3f = flipud(b3(:));
b4 = m(9:11); b4f = flipud(b4(:));
t3 = m(12);
t4 = m(13:15); t4f = flipud(t4(:));
p = m(16);
v3 = v3(:); b3 = b3(:); b4 = b4(:); t4 = t4(:);
O33 = zeros(3,3); O31 = zeros(3,1); O13 = O31';
B = [p,  t4f',           t3,  b4f',            b3f',           v4,   v3f',          s;
    t4f,  p*I+skew(b4f),  b3,  -v4*I-skew(t4f), t3*J+weks(v3),  -b4f, -s*J-weks(b3), -v3;
    0,    O13,            p,   O13,             -t4',           0,    b4',           -v4;
    -b4f, v4*I+skew(t4f), -v3, p*I+skew(b4f),   s*J+weks(b3),   -t4f, t3*J+weks(v3), -b3;
    O31,  O33,            -t4, O33,             p*I-skew(b4),   O31,  v4*I-skew(t4), -b4;
    -v4,  b4f',           -s,  -t4f',           v3f',           p,    -b3f',          t3;
    O31,  O33,            -b4, O33,             -v4*I+skew(t4), O31,  p*I-skew(b4),   t4;
    0,    O13,            v4,  O13,             b4',            0,    t4',            p];
end