function Q = wahba(M,N)
%WAHBA Q = wahba(M,N) finds motor Q such that N = Q*M*~Q
% M and N are corresonding arrays of multivectors containing any
% combination of points, lines, planes, and/or motors.
% If all objects in M & N are the same multivector type, e.g. they are all
% planes, all points, etc. then M & N can be ordinary arrays.  If there are
% mixed types within M & N, then M & N must be cell arrays.

if nargin < 1
    test
    Q = @test;
    return
end

%% Transcribe multivectors to column vectors
lenM = length(M);
for i = lenM:-1:1
    if iscell(M)
        Mi = M{i}; Ni = N{i};
    else
        Mi = M(i); Ni = N(i);
    end
    mw(:,i) = weight(Mi).m'; nw(:,i) = weight(Ni).m';
    mb(:,i) = bulk(Mi).m'; nb(:,i) = bulk(Ni).m';
end
U = eye(16); U(2:11,2:11)=-eye(10); % anti-reverse

%% Solve for rotational part of motor
C = 0;
for i = 1:lenM
    C = C + Psi(nw(:,i))'*Xi(mw(:,i));
    %C = C + Psi(nw(:,i)+nb(:,i))'*Xi(mw(:,i)+mb(:,i)); % works for unit directions
    %C = C + Psi(nw(:,i)+nb(:,i))'*Xi(mw(:,i)) + U*Xi(nw(:,i)+nb(:,i))'*Psi(mw(:,i))*U;
end
C = (C+C');
Ew(16,4) = 1; Ew(11:-1:9,1:3) = eye(3);
[V,e] = eig(Ew'*C*Ew,'vector');
V = real(V); e = real(e);
V(abs(V)<1e4*eps) = 0;
if all(all(abs(V-eye(4)) < 1e4*eps))
    qw(16,1) = 1;
    warning('Identity rotation found; possibly due to poor observability')
else
    [~,k] = sort(e,'descend');
    qw = Ew*V(:,k(1));
end

%% Solve for remaining part of motor
C = 0;
d = 0;
for i = 1:lenM
    C = C + (Xi(qw)'*Xi(mw(:,i)) + Psi(qw)*Psi(mw(:,i))*U);
    d = d + nb(:,i) - Xi(qw)'*Psi(qw)*mb(:,i);
end
Eb(16,4) = 0; Eb([6:8 1],:) = eye(4); % s = Eb'*qb
%qb = Eb*((Eb'*(C'*C)*Eb)\(Eb'*C'*d));
qb = Eb*lsqminnorm((Eb'*(C'*C)*Eb),(Eb'*C'*d));
%qb = Eb*pinv((Eb'*(C'*C)*Eb))*(Eb'*C'*d);
%qb = Eb*Eb'*lsqminnorm(C,d);
rnk = rank(Eb'*(C'*C)*Eb);
if rnk < 4
    warning(['Translational solution may be inaccurate; rank = ',...
        num2str(rnk), ' vs. 4 needed'])
end

%% Enforce constraints and create output motor
qw = qw/norm(qw);
qb([6:8 1]) = reject(qb([6:8 1]),qw([11 10 9 16]));
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

function test(opts)
%% Test Wahba
arguments
    opts.points (1,1) double = 1
    opts.point_sigma (1,1) double = 0
    opts.planes (1,1) double = 1
    opts.plane_sigma (1,1) = 0
    opts.directions (1,1) double = 1
    opts.direction_sigma (1,1) double = 0
    opts.lines (1,1) = 1
    opts.line_sigma (1,1) double = 0
    opts.motors (1,1) double = 1
    opts.motor_sigma (1,1) double = 0
    opts.rotation_sigma (1,1) double = 0
    opts.translation_sigma (1,1) double = 0
    opts.axis_sigma (1,1) double = 0;
end
%Ltru = unitize(rgaline);
%phi = pi*rand - pi;
%d = 10*randn;
%R = unitize(rgamotor(phi,0,Ltru)); % pure rotation about L
%R0 = unitize(rgamotor(phi,0,direction(Ltru))); % pure rotation about origin
%S = rgamotor(0,d,Ltru); % pure translation along direction of L
%S0 = rgamotor(0,d,direction(Ltru)); % should be same as S
%R.anti = true; S.anti = true; 
%R0.anti = true; S0.anti = true;
%Q0 = rgamotor(R0*S0);
%Qtru = rgamotor(R*S);
Qtru = unitize(rgamotor); Qtru.anti = true;
[~,~,vtru,mtru] = extract(Qtru);
n = opts.points + opts.planes + opts.directions + opts.lines + opts.motors;
for i = n:-1:1
    pherr = opts.rotation_sigma*randn;
    derr = opts.translation_sigma*randn;
    vper = vtru + opts.axis_sigma*randn(1,3);
    mper = mtru + opts.axis_sigma*randn(1,3);
    mper = reject(mper,vper);
    Lper = unitize(rgaline(vper,mper)); % perturbed axis (dir & mom)
    dR = unitize(rgamotor(pherr,0,Lper)); % rotation error
    dS = rgamotor(0,derr,Lper); % translation error
    dR.anti = true; dS.anti = true;
    dQ = rgamotor(dR*dS); % motor error
    Q = dQ*Qtru*~dQ;
    if opts.points
        M{i} = unitize(rgapoint); M{i}.anti = true;
        N{i} = rgapoint(Q*M{i}*~Q); N{i}.anti = true;
        if opts.point_sigma
            N{i}.m = N{i}.m + opts.point_sigma*randn(1,16).*N{i}.m;
            N{i} = unitize(N{i});
        end
        opts.points = opts.points - 1;
    elseif opts.planes
        M{i} = unitize(rgaplane); % anti by default
        N{i} = rgaplane(Q*M{i}*~Q);
        if opts.plane_sigma
            N{i}.m = N{i}.m + opts.plane_sigma*randn(1,16).*N{i}.m;
            N{i} = unitize(N{i});
        end
        opts.planes = opts.planes - 1;
    elseif opts.directions
        M{i} = rgaline(unit(randn(1,3)),zeros(1,3)); M{i}.anti = true;
        N{i} = rgaline(Q*M{i}*~Q); N{i}.anti = true;
        if opts.direction_sigma
            N{i}.m = N{i}.m + opts.direction_sigma*randn(1,16).*N{i}.m;
            N{i} = unitize(N{i});
        end
        opts.directions = opts.directions - 1;
    elseif opts.lines
        M{i} = unitize(rgaline); M{i}.anti = true;
        N{i} = rgaline(Q*M{i}*~Q); N{i}.anti = true;
        if opts.line_sigma
            N{i}.m = N{i}.m + opts.line_sigma*randn(1,16).*N{i}.m;
            N{i} = unitize(N{i});
        end
        opts.lines = opts.lines - 1;
    elseif opts.motors
        M{i} = unitize(rgamotor); M{i}.anti = true;
        N{i} = rgamotor(Q*M{i}*~Q); N{i}.anti = true;
        if opts.motor_sigma
            N{i}.m = N{i}.m + opts.motor_sigma*randn(1,16).*N{i}.m;
            N{i} = unitize(N{i});
        end
        opts.motors = opts.motors - 1;
    end
    %M{i,1} = a(i); N{i,1} = b(i); % points
    %M{i,2} = f(i); N{i,2} = g(i); % planes
    %M{i,3} = u(i); N{i,3} = v(i); % directions
    %M{i,4} = K(i); N{i,4} = L(i); % lines
    %M{i,5} = R(i); N{i,5} = S(i); % motors
end
%M = M(:);
%N = N(:);
%M = reshape(M(:,[2 4]),[],1);
%N = reshape(N(:,[2 4]),[],1);

Qest = wahba(M,N);
Qest.anti = true;
disp('Qtru = ')
disp(Qtru)
disp('Qest =')
disp(Qest)
disp('Qest - Q =')
disp(Qest - Q)
[phit,dt,vt,mt] = extract(Qtru);
[phie,de,ve,me] = extract(Qest);
disp('Rotation error [deg] (phie - phit) =')
disp(180/pi*(phie-phit))
disp('Translation error (de - dt) =')
disp(de - dt)%, if abs(de-dt)>1, keyboard, end
disp('True Direction & Moment of Screw Axis = ')
disp([vt mt])
disp('Est Direction & Moment of Screw Axis = ')
disp([ve me])

end