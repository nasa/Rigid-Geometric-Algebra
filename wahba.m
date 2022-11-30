function Q = wahba(M,N)
%WAHBA Q = wahba(M,N) finds motor Q such that N = Q*M*~Q
% M and N are corresonding arrays of multivectors containing any
% combination of points, lines, planes, and/or motors.
% If all objects in M & N are the same multivector type, e.g. they are all
% planes, all points, etc. then M & N can be ordinary arrays.  If there are
% mixed types within M & N, then M & N must be cell arrays.

if nargin < 1
    test
    return
end

%% Transcribe multivectors to column vectors
lenM = length(M);
for i = lenM:-1:1
    if iscell(M)
        %m(:,i) = M{i}.m'; n(:,i) = N{i}.m';
        mw(:,i) = weight(M{i}).m'; nw(:,i) = weight(N{i}).m';
        mb(:,i) = bulk(M{i}).m'; nb(:,i) = bulk(N{i}).m';
    else
        %m(:,i) = M(i).m'; n(:,i) = N(i).m';
        mw(:,i) = weight(M(i)).m'; nw(:,i) = weight(N(i)).m';
        mb(:,i) = bulk(M(i)).m'; nb(:,i) = bulk(N(i)).m';
    end
end

%% Solve for rotational part of motor
C = 0;
for i = 1:lenM
    C = C + Psi(nw(:,i))'*Xi(mw(:,i));
end
C = (C+C')/2;
Ew(16,4) = 1; Ew(11:-1:9,1:3) = eye(3);
[V,e] = eig(Ew'*C*Ew,'vector');
V = real(V); e = real(e);
V(abs(V)<1e4*eps) = 0;
[~,k] = sort(e,'descend');
qw = Ew*V(:,k(1));

%% Solve for translation part of motor
U = eye(16); U(2:11,2:11)=-eye(10); % anti-reverse
C = 0;
d = 0;
for i = 1:lenM
    C = C + (Xi(qw)'*Xi(mw(:,i)) + Psi(qw)*Psi(mw(:,i))*U);
    d = d + nb(:,i) - Xi(qw)'*Psi(qw)*mb(:,i);
end
Eb(16,4) = 0; Eb([1 6:8],:) = eye(4);
qb = Eb*((Eb'*(C'*C)*Eb)\(Eb'*C'*d));

%% Enforce constraints and create output motor
qw = qw/norm(qw);
qb([1 6:8]) = reject(qb([1 6:8]),qw([16 11 10 9]));
Q = rgamotor(qw([11 10 9 16]),qb([6:8 1]));
end

%% Helper functions

function A = Psi(m)
% Psi(Q.m) is the productmat such that Psi(Q.m)*x.m' = Q*x, where * is the
% antiwedgedot product.
arguments
    m (1,16)
end
I = eye(3); J = fliplr(I);
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

function a2 = reject(a,b)
a2 = a - dot(a,b)/dot(b,b)*b;
end

function test
%% Test Wahba
Ltru = unitize(rgaline);
phi = pi*rand - pi;
d = 10*randn;
R = unitize(rgamotor(phi,0,Ltru)); % pure rotation about L
%R0 = unitize(rgamotor(phi,0,direction(Ltru))); % pure rotation about origin
S = rgamotor(0,d,Ltru); % pure translation along direction of L
%S0 = rgamotor(0,d,direction(Ltru)); % shuold be same as S
R.anti = true; S.anti = true; R0.anti = true; S0.anti = true;
%Q0 = rgamotor(R0*S0);
Qtru = rgamotor(R*S);
%Qtru = unitize(rgamotor); Q.anti = true;
disp('Qtru = ')
disp(Qtru)
[phi,d,v,m] = extract(Qtru)
n = 10;
for i = n:-1:1
    phi = 0.1*pi/180*rand;
    d = 0.1*randn;
    dR = unitize(rgamotor(phi,0,Ltru)); % rotation error
    dS = rgamotor(0,d,Ltru); % translation error
    dR.anti = true; dS.anti = true;
    dQ = rgamotor(dR*dS); % motor error
    Q = dQ*Qtru*~dQ;
    a(i) = unitize(rgapoint); a(i).anti = true;
    b(i) = rgapoint(Q*a(i)*~Q); b(i).anti = true;
    f(i) = unitize(rgaplane); % anti by default
    g(i) = rgaplane(Q*f(i)*~Q); 
    K(i) = unitize(rgaline); K(i).anti = true;
    L(i) = rgaline(Q*K(i)*~Q); L(i).anti = true;
    R(i) = unitize(rgamotor); R(i).anti = true;
    S(i) = rgamotor(Q*R(i)*~Q); S(i).anti = true;
    M{i,1} = a(i); N{i,1} = b(i);
    M{i,2} = f(i); N{i,2} = g(i);
    M{i,3} = K(i); N{i,3} = L(i);
    M{i,4} = R(i); N{i,4} = S(i);
end
M = M(:);
N = N(:);

Qhat = wahba(M,N);
Qhat.anti = true;
disp('Qhat =')
disp(Qhat)
[phi,d,v,m] = extract(Qhat)
disp('Qhat - Q =')
disp(Qhat - Q)
[phi,d,v,m] = extract(rgamotor(Qhat - Q))
end