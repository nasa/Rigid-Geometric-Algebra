function Q = wahba(M,N,algorithm)
%WAHBA Q = wahba(M,N) finds motor Q such that N = Q*M*~Q
% M and N are corresonding arrays of multivectors containing any
% combination of points, lines, planes, and/or motors.
% If all objects in M & N are the same multivector type, e.g. they are all
% planes, all points, etc. then M & N can be ordinary arrays.  If there are
% mixed types within M & N, then M & N must be cell arrays.
% This function requires the rga class.
% Q = wahba(M,N,algorithm) allows a choice of solution method.
% If algorithm = 'one-step' wahba uses the approach described by
% Perwass, C., Geometric Algebra with Applications in Engineering,
% Springer, 2010.  This is the default.
% Otherwise, it uses a two-step approach that first solves for the
% rotational part of the motor, then "completes the screw" to solve for the
% translational part.  This was loosely inspired by
% Walker, Shao, and Volz, "Estimating 3-D Location Parameters Using Dual
% Number Quaternions," CVGIP: Image Understanding, 54(3) 358-367, 1991.

if nargin < 1
    wahba_test
    return
end
lenM = length(M);
if nargin < 3
    algorithm = 'one-step';
end

switch algorithm
    case 'one-step'
        %% Null-space approach
        % Solve for motor in one step, as in Perwass p.266
        H([1 6:11 16],:) = eye(8);
        C = 0;
        for i = 1:lenM
            %C = C + Xi(M{i}.m) - Psi(N{i}.m);
            B = Xi(M{i}.m) - Psi(N{i}.m);
            C = C + B'*B;
        end
        A = C*H;
        [~,S,V] = svd(A,'vector');
        tol = max(size(A))*eps(norm(A));
        k = find(S<tol);
        if isempty(k)
            if any(abs(diff(S))<tol) % check for repeats
                k = 7:8; % min pair since svd returns in decreasing order
            else
                k = 8; % no repeats so just one min that will be last
            end
        end
        lenk = length(k);
        if lenk > 1
            warning('Null space of Xi(M)-Psi(N) spans more than one motor')
        end
        %[~,k] = min(S); % not needed; svd returns them in decreasing order
        for j = lenk:-1:1
            Qsvd(j) = rgamotor(rga(H*V(:,k(j)),true));
            wtnormQ(j) = norm(Qsvd(j),'weight');
        end
        %Qsvd = unitize(rgamotor(rga(H*V(:,k),true))) if we unitize it, it's no
        %longer a null vector, so need to scale entire thing!
        keep = (wtnormQ>tol);
        Q = unitize(1/wtnormQ(keep)*Qsvd(keep));

    otherwise
        %% Transcribe multivectors to column vectors
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
        end
        C = (C+C')/2;
        Hw(16,4) = 1; Hw(11:-1:9,1:3) = eye(3);
        [V,e] = eig(Hw'*C*Hw,'vector');
        V = real(V); e = real(e);
        V(abs(V)<1e4*eps) = 0;
        if all(all(abs(V-eye(4)) < 1e4*eps))
            qw(16,1) = 1;
            warning('Identity rotation found; possibly due to poor observability')
        else
            [~,k] = sort(e,'descend');
            if abs(e(k(1))-e(k(2))) < 10*eps(e(k(1)))
                warning('Repeated eigenvalues in rotational solution')
            end
            qw = Hw*V(:,k(1));
        end

        %% Solve for remaining part of motor
        C = 0;
        d = 0;
        for i = 1:lenM
            C = C + (Xi(qw)'*Xi(mw(:,i)) + Psi(qw)*Psi(mw(:,i))*U);
            d = d + nb(:,i) - Xi(qw)'*Psi(qw)*mb(:,i);
        end
        Hb(16,4) = 0; Hb([6:8 1],:) = eye(4); % s = Hb'*qb

        Chat = C*Hb;
        qwhat = Hw'*qw;
        A = [Chat'*Chat; qwhat'];
        rnk = rank(A);
        if rnk < 4
            warning(['Translational solution may be inaccurate; rank = ',...
                num2str(rnk), ' vs. 4 needed'])
        end
        qb = Hb*(A\[Chat'*d; 0]);

        %% Create output motor
        Q = unitize(rgamotor(qw([11 10 9 16]),qb([6:8 1])));
end
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