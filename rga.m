classdef (InferiorClasses = {?sym}) rga < matlab.mixin.indexing.RedefinesDot
    % G(3,0,1) Geometric Algebra for homogeneous transformations of 3D
    % rigid bodies in 4D projective space, primarily based on Lengyal's
    % conventions:
    %    https://projectivegeometricalgebra.org/
    %    https://rigidgeometricalgebra.org/

    properties
        m (1,16) % coefficients of [e0 e1 e2 e3 e4 e23 e31 e12 e43 e42 e41 e321 e412 e431 e423 e1234]
        anti = false % set to true to use anti-basis elements; will "stick" to any obj it touches
        dispeps = 1e-15 % threshold for displaying small terms; set to 0 to show all
    end

    methods
        function obj = rga(v,anti)
            %RGA Construct rga multivector object
            %   obj = rga([e0 e1 e2 e3 e4 e23 e31 e12 e43 e42 e41 e321 e412 e431 e423 e1234])
            %   Symbolic inputs ok using above syntax.
            %   obj = rga('e0') or rga("e0") etc. also recognized for single elements
            %   Dual basis input e.g. 'eps0' etc. also works and creates e1234, etc.
            %   To create objects like 5e0 - 3e12 etc. first create a set
            %   of basis elements using RGA.BASES.
            if nargin == 0
                obj.m = randn(1,16);
            elseif ischar(v) || isstring(v)
                m = eye(1,16);
                switch v
                    case {'e0','eps1234'}
                        obj.m = m;
                    case {'e1','eps423'}
                        obj.m = circshift(m,1);
                    case {'e2','eps431'}
                        obj.m = circshift(m,2);
                    case {'e3','eps412'}
                        obj.m = circshift(m,3);
                    case {'e4','eps321'}
                        obj.m = circshift(m,4);
                    case {'e23','eps41'}
                        obj.m = circshift(m,5);
                    case {'e31','eps42'}
                        obj.m = circshift(m,6);
                    case {'e12','eps43'}
                        obj.m = circshift(m,7);
                    case {'e43','eps12'}
                        obj.m = circshift(m,8);
                    case {'e42','eps31'}'
                        obj.m = circshift(m,9);
                    case {'e41','eps23'}
                        obj.m = circshift(m,10);
                    case {'e321','eps4'}
                        obj.m = circshift(m,11);
                    case {'e412','eps3'}
                        obj.m = circshift(m,12);
                    case {'e431','eps2'}
                        obj.m = circshift(m,13);
                    case {'e423','eps1'}
                        obj.m = circshift(m,14);
                    case {'e1234','eps0'}
                        obj.m = circshift(m,15);
                    otherwise
                        error('string or char input not recognized')
                end
                if contains(v,'eps') && nargin == 1
                    obj.anti = true;
                end
            else
                if isa(v,"rga")
                    obj = v;
                else
                    obj.m = v;
                end
            end
            if nargin == 2
                obj.anti = anti;
            end
        end

        function iseq = eq(a,b)
            %EQ Overload == for rga objects
            if a.m == b.m
                iseq = true;
            else
                iseq = false;
            end
        end

        function obj = scalar(obj)
            %SCALAR Return scalar part
            obj.m(2:end) = 0;
        end

        function obj = vector(obj)
            %VECTOR Return vector part
            obj.m([1 6:end]) = 0;
        end

        function obj = bivector(obj)
            %BIVECTOR Return bivector part
            obj.m([1:5 12:end]) = 0;
        end

        function obj = trivector(obj)
            %TRIVECTOR Return trivector part
            obj.m([1:11 16]) = 0;
        end

        function obj = antivector(obj)
            %ANTIVECTOR Return antivector part
            obj = trivector(obj);
            obj.anti = true;
        end

        function obj = pseudoscalar(obj)
            %PSEUDOSCALAR Return pseudoscalar part
            obj.m(1:15) = 0;
        end

        function obj = antiscalar(obj)
            %ANTISCALAR Return antiscalar part
            obj = pseudoscalar(obj);
            obj.anti = true;
        end

        function obj = bulk(obj)
            %BULK Bulk - portion not containing factors of e4
            if isempty(obj)
                obj = diag([1 1 1 1 0 1 1 1 0 0 0 1 0 0 0 0]);
            else
                obj.m([5 9:11 13:end]) = 0;
            end
        end

        function obj = weight(obj)
            %WEIGHT Weight - portion containing factors of e4
            if isempty(obj)
                obj = diag([0 0 0 0 1 0 0 0 1 1 1 0 1 1 1 1]);
            else
                obj.m([1:4 6:8 12]) = 0;
            end
        end

        function g = gr(obj)
            %GR Grade - Number of basis vectors present in factorization
            nz = find(obj.m ~= 0);
            g = 0;
            if any(ismember(nz,[2 3 4 5])), g = 1; end
            if any(ismember(nz,[6 7 8 9 10 11])), g = 2; end
            if any(ismember(nz,[12 13 14 15])), g = 3; end
            if any(ismember(nz,16)), g = 4; end
        end

        function ag = antigr(obj)
            %ANTIGR Antigrade - Number of basis vectors absent in factorization
            ag = 4 - gr(obj);
        end

        function obj = rev(obj)
            %REV Reverse
            if isempty(obj)
                obj = diag([1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 1]);
            else
                obj.m(6:15) = -obj.m(6:15);
            end
        end

        function obj = antirev(obj)
            %ANTIREV Anti-reverse
            if isempty(obj)
                obj = diag([1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 1 1 1 1 1]);
            else
                obj.m(2:11) = -obj.m(2:11);
            end
        end

        function obj = not(obj)
            % Overload ~ for reverse & antireverse - DEPRECATED
            warning('~ is being deprecated in favor of '' and .''')
            if obj.anti
                obj = antirev(obj);
                obj.anti = true;
            else
                obj = rev(obj);
            end
        end

        function obj = ctranspose(obj)
            % Overload ' for reverse
            obj = rev(obj);
        end

        function obj = transpose(obj)
            % Overload .' for anti-reverse
            obj = antirev(obj);
        end

        function obj = rcomp(obj)
            %RCOMP Right complement
            obj = rev(obj); % if isempty(input obj) then LHS is a matrix
            if isa(obj,"rga")
                obj.m = flipud(obj.m(:));
            else 
                obj = flipud(obj); % rcomp(rga.empty)*a.m(:) <=> rcomp(a)
            end
        end

        function obj = lcomp(obj)
            %LCOMP Left complement
            obj = antirev(obj);
            if isa(obj,"rga")
                obj.m = flipud(obj.m(:));
            else
                obj = flipud(obj);
            end
        end

        function obj = plus(a,b)
            %PLUS RGA addition
            % If both inputs are the same type of object, output will
            % be the same type.
            % If at least one input uses the anti basis, output will too.
            % If only one input is a multivector, then this adds the other
            % input to the scalar component of the multivector, and the
            % output will be a generic rga object.
            obj = rga;
            arga = isa(a,"rga");
            brga = isa(b,"rga");
            if arga && brga
                obj.m = a.m + b.m;
                cla = class(a); clb = class(b);
                if matches(cla,clb) % Enforce common subclass
                    obj = feval(cla,obj);
                end
                aanti = a.anti;
                banti = b.anti;
            elseif ~arga && length(a)==1
                obj.m = b.m;
                obj.m(1) = a + obj.m(1);
                aanti = false;
                banti = b.anti;
            elseif ~brga && length(b)==1
                obj.m = a.m;
                obj.m(1) = obj.m(1) + b;
                aanti = a.anti;
                banti = false;
            else
                error('Inputs incompatible')
            end
            if aanti || banti % Enforce stickiness of anti property
                obj.anti = true;
            end
        end

        function obj = uminus(obj)
            %UMINUS RGA unary minus
            obj.m = -obj.m;
        end

        function obj = minus(a,b)
            %MINUS RGA subtraction
            obj = plus(a,-b);
        end

        function [A,B] = productmat(obj,type)
            % PRODUCTMAT Matrices for computing various multivector products
            % Common ones hardcoded here; less common ones
            % computed on the fly via createprodmat
            I = eye(3); J = fliplr(I);
            skew = @(x) [0,-x(3),x(2);x(3),0,-x(1);-x(2),x(1),0];
            weks = @(x) fliplr(skew(x));
            s = obj.m(1);
            v3 = obj.m(2:4); v3f = flipud(v3(:)); v4 = obj.m(5);
            b3 = obj.m(6:8); b3f = flipud(b3(:));
            b4 = obj.m(9:11); b4f = flipud(b4(:));
            t3 = obj.m(12);
            t4 = obj.m(13:15); t4f = flipud(t4(:));
            p = obj.m(16);
            v3 = v3(:); b3 = b3(:); b4 = b4(:); t4 = t4(:);
            O33 = zeros(3,3); O31 = zeros(3,1); O13 = O31';
            switch type
                case "wedgedot"
                    A = [s,            v3',    0,          -b3',            O13, -t3,             O13,   0;
                        v3,   s*I-skew(b3),  O31, t3*I-skew(v3),            O33,  b3,             O33, O31;
                        v4,           b4f',    s,         -t4f',          -v3f',  -p,           -b3f',  t3;
                        b3, -t3*I+skew(v3),  O31,  s*I-skew(b3),            O33, -v3,             O33, O31;
                        b4,  v4*J+weks(t4), -v3f,  p*J+weks(b4), s*I-skew(b3f)',  t4, -t3*I-skew(v3f), b3f;
                        t3,           -b3',    0,          -v3',            O13,   s,             O13,   0;
                        t4,  -p*J-weks(b4),  b3f, v4*J+weks(t4), t3*I-skew(v3f)',-b4,   s*I+skew(b3f), v3f;
                        p,           -t4f',  -t3,         -b4f',          -b3f',  v4,            v3f',   s];
                    if nargout == 2
                        B = createprodmat(obj,@wedgedot,true);
                    end
                case "wedge"
                    A = [s,            O13,    0,           O13,            O13,   0,             O13,   0;
                        v3,            s*I,  O31,           O33,            O33, O31,             O33, O31;
                        v4,            O13,    s,           O13,            O13,   0,             O13,   0;
                        b3,       skew(v3),  O31,           s*I,            O33, O31,             O33, O31;
                        b4,           v4*J, -v3f,           O33,            s*I, O31,             O33, O31;
                        t3,           -b3',    0,          -v3',            O13,   s,             O13,   0;
                        t4,      -weks(b4),  b3f,          v4*J,     -skew(v3f)',O31,             s*I, O31;
                        p,           -t4f',  -t3,         -b4f',          -b3f',  v4,            v3f',   s];
                    if nargout == 2
                        B = createprodmat(obj,@wedge,true);
                    end
                case "motor"
                    % Unitize attitude part first
                    r = [b4f;p];
                    r = r/sqrt(sum(r.^2));
                    b4f = r(1:3); p = r(4);
                    A = [2*I*(1-b4f'*b4f) + 2*(b4f*b4f') - I, 2*cross(b4f,b3);
                        O13, 1];
                    B = [2*skew(b4f)*p, 2*(p*b3 - s*b4f);
                        O13, 0];
                case "antiwedgedot"
                    A = [p,   -t4f',          -t3, -b4f',           -b3f',           v4,   v3f',          s;
                        -t4f, p*I+skew(b4f),   b3, -v4*I-skew(t4f),  t3*J+weks(v3),  b4f,  s*J+weks(b3),  v3;
                        0,    O13,             p,   O13,            -t4',            0,   -b4',           v4;
                        b4f,  v4*I+skew(t4f), -v3,  p*I+skew(b4f),   s*J+weks(b3),   t4f, -t3*J-weks(v3), b3;
                        O31,  O33,            -t4,  O33,             p*I-skew(b4),   O31, -v4*I+skew(t4), b4;
                        -v4,  -b4f',           s,   t4f',           -v3f',           p,   -b3f',          t3;
                        O31,  O33,             b4,  O33,             v4*I-skew(t4),  O31,  p*I-skew(b4),  t4;
                        0,    O13,            -v4,  O13,            -b4',            0,    t4',           p]; 
                    % B is s/t A(a)*B(a)*x is linear algebra equivalent to a.*x.*a.' 
                    B = [p,   t4f',           t3,  b4f',            b3f',           v4,   v3f',          s;
                        t4f,  p*I+skew(b4f),  b3,  -v4*I-skew(t4f), t3*J+weks(v3),  -b4f, -s*J-weks(b3), -v3;
                        0,    O13,            p,   O13,             -t4',           0,    b4',           -v4;
                        -b4f, v4*I+skew(t4f), -v3, p*I+skew(b4f),   s*J+weks(b3),   -t4f, t3*J+weks(v3), -b3;
                        O31,  O33,            -t4, O33,             p*I-skew(b4),   O31,  v4*I-skew(t4), -b4;
                        -v4,  b4f',           -s,  -t4f',           v3f',           p,    -b3f',          t3;
                        O31,  O33,            -b4, O33,             -v4*I+skew(t4), O31,  p*I-skew(b4),   t4;
                        0,    O13,            v4,  O13,             b4',            0,    t4',            p];

                case "antiwedge" 
                    A = [p, -t4f', -t3, -b4f',     -b3f', v4,  v3f',     s;
                        O31, p*I, O31,  -skew(t4f), t3*J,  b4f, weks(b3), v3;
                        0,   O13, p,    O13,       -t4',  0,   -b4',      v4;
                        O31, O33, O31,  p*I,       O33,   t4f, -t3*J,    b3;
                        O31, O33, O31,  O33,       p*I,   O31, skew(t4), b4;
                        %0,   O13, 0,  O13,       O13,   p,   O13,      t3;
                         0,   O13, 0,  O13,       O13,   p,   -b3',      t3;
                        O31, O33, O31,  O33,       O33,   O31, p*I,      t4;
                        0,   O13, 0,    O13,       O13,   0,   O13,      p]; 
                    if nargout == 2
                        B = createprodmat(obj,@antiwedge,true);
                    end
                otherwise
                    A = createprodmat(obj,eval(strcat("@",type)),false);
                    if nargout == 2
                        B = createprodmat(obj,eval(strcat("@",type)),true);
                    end
                    %error('product type not recognized')
            end
        end

        function obj = dot(a,b)
            %DOT RGA dot product
            % If at least one input uses the anti basis, output will too.
            M = diag([1 1 1 1 0 1 1 1 0 0 0 1 0 0 0 0]);
            arga = isa(a,"rga");
            if nargin == 1 && arga
                obj = a.m(:)'*M;
                return
            end
            brga = isa(b,"rga");
            if arga && brga
                s = a.m(:)'*M*b.m(:);
                obj = rga(s*eye(16,1));
            else
                error('Inputs incompatible')
            end
            if a.anti || b.anti % Enforce stickiness
                obj.anti = true;
            end
        end

        function obj = antidot(a,b)
            %ANTIDOT RGA anti-dot product
            % If at least one input uses the anti basis, output will too.
            M = diag(fliplr([1 1 1 1 0 1 1 1 0 0 0 1 0 0 0 0]));
            arga = isa(a,"rga");
            if nargin == 1 && arga
                obj = a.m(:)'*M;
                return
            end
            brga = isa(b,"rga");
            if arga && brga
                s = a.m(:)'*M*b.m(:);
                obj = rga(s*flipud(eye(16,1)));
            else
                error('Inputs incompatible')
            end
            if a.anti || b.anti % Enforce stickiness
                obj.anti = true;
            end
        end

        function obj = wedgedot(a,b)
            %WEDGEDOT RGA wedgedot product
            % If at least one input uses the anti basis, output will too.
            % If only one input is a multivector, the other input must be
            % an ordinary scalar, and the output type will be the type of
            % the input multivector.
            arga = isa(a,"rga");
            if nargin == 1 && arga
                obj = productmat(a,"wedgedot");
                return
            end
            brga = isa(b,"rga");
            obj = rga;
            if arga && brga
                A = productmat(a,"wedgedot");
                obj.m = A*b.m(:);
                aanti = a.anti;
                banti = b.anti;
            elseif ~arga && length(a)==1
                obj.m = a * b.m;
                obj = feval(class(b),obj);
                aanti = false;
                banti = b.anti;
            elseif ~brga && length(b)==1
                obj.m = b * a.m;
                obj = feval(class(a),obj);
                aanti = a.anti;
                banti = false;
            else
                error('Inputs incompatible')
            end
            if aanti || banti % Enforce stickiness of anti property
                obj.anti = true;
            end
        end

        function obj = mtimes(a,b)
            % Overload * as wedgedot product
            obj = wedgedot(a,b);
        end

        function obj = times(a,b)
            % Overload .* as antiwedgedot product
            obj = antiwedgedot(a,b);
        end

        function obj = wedge(a,b)
            %WEDGE RGA wedge product
            % If at least one input uses the anti basis, output will too.
            % If only one input is a multivector, the other input must be
            % an ordinary scalar, and the output type will be the type of
            % the input multivector.
            arga = isa(a,"rga");
            if nargin == 1 && arga
                obj = productmat(a,"wedge");
                return
            end
            brga = isa(b,"rga");
            obj = rga;
            if arga && brga
                A = productmat(a,"wedge");
                obj.m = A*b.m(:);
                aanti = a.anti;
                banti = b.anti;
            elseif ~arga && length(a)==1
                obj.m = a * b.m;
                obj = feval(class(b),obj);
                aanti = false;
                banti = b.anti;
            elseif ~brga && length(b)==1
                obj.m = b * a.m;
                obj = feval(class(a),obj);
                aanti = a.anti;
                banti = false;
            else
                error('Inputs incompatible')
            end
            if aanti || banti % Enforce stickiness
                obj.anti = true;
            end
        end

        function obj = mpower(a,b)
            % Overload ^ as wedge product
            obj = wedge(a,b);
        end

        function obj = power(a,b)
            % Overload .^ as wedge product
            obj = antiwedge(a,b);
        end

        function obj = antiwedge(a,b)
            %ANTIWEDGE Anti-wedge product
            if nargin == 1 && isa(a,"rga")
                obj = productmat(a,"antiwedge");
            else
                obj = lcomp(wedge(rcomp(a),rcomp(b)));
                %obj = rga;
                %obj.m = productmat(a,"antiwedge")*b.m(:);
            end
        end

        function obj = antiwedgedot(a,b)
            %ANTIWEDGEDOT Anti-wedgedot product
            if nargin == 1 && isa(a,"rga")
                % A(a)*b = antiwedgedot(a,b).m
                obj = productmat(a,"antiwedgedot");
            elseif nargin == 2 && isa(b,"rga") && isempty(a)
                % B*(antirev(b))*a = A(a)*b = antiwedgedot(a,b).m
                [~,obj] = productmat(b,"antiwedgedot");
            else
                obj = lcomp(wedgedot(rcomp(a),rcomp(b)));
            end
        end

        function obj = rint(a,b)
            %RINT Right interior product
            if nargin == 1 && isa(a,"rga")
                obj = antiwedge(a)*rcomp(rga.empty);
            else
                obj = antiwedge(a,rcomp(b));
            end
        end

        function obj = lint(a,b)
            %LINT Left interior product
            if nargin == 1 && isa(a,"rga")
                obj = antiwedge(lcomp(a));
            else
                obj = antiwedge(lcomp(a),b);
            end
        end

        function obj = antirint(a,b)
            %ANTIRINT Right interior antiproduct
            if nargin == 1 && isa(a,"rga")
                obj = wedge(a)*rcomp(rga.empty);
            else
                obj = wedge(a,rcomp(b));
            end
        end

        function obj = antilint(a,b)
            %ANTILINT Left interior antiproduct
            if nargin == 1 && isa(a,"rga")
                obj = wedge(lcomp(a));
            else
                obj = wedge(lcomp(a),b);
            end
        end

        function obj = proj(a,b)
            %PROJ Project a onto b
            obj = lint(rint(weight(b),a),b);
        end

        function obj = antiproj(a,b)
            %ANTIPROJ Antiproject a onto b
            obj = antilint(antirint(weight(b),a),b);
        end

        function obj = commutate(a,b,type)
            %COMMUTATE Four types:
            % commutate(a,b,"-^") = [a,b]-^ (default)
            % commutate(a,b,"+^") = [a,b]+^
            % commutate(a,b,"-v") = [a,b]-v
            % commutate(a,b,"+v") = [a,b]+v
            arguments
                a rga
                b rga
                type string = "-^"
            end
            switch type
                % GDC'21 talk shows reverses on 2nd argument, but poster
                % and wiki do not.  To get Euclidean distances shown on
                % poster, must not use reverses.  If they are included, get
                % Inf and NaN distances.
                case "-^"
                    obj = wedgedot(a,b) - wedgedot(b,a);
                    %obj = wedgedot(a,rev(b)) - wedgedot(b,rev(a));
                case "+^"
                    obj = wedgedot(a,b) + wedgedot(b,a);
                    %obj = wedgedot(a,rev(b)) + wedgedot(b,rev(a));
                case "-v"
                    obj = antiwedgedot(a,b) - antiwedgedot(b,a);
                    %obj = antiwedgedot(a,rev(b)) - antiwedgedot(b,rev(a));
                case "+v"
                    obj = antiwedgedot(a,b) + antiwedgedot(b,a);
                    %obj = antiwedgedot(a,rev(b)) + antiwedgedot(b,rev(a));
                otherwise
                    error('type not recognized')
            end
            obj.m = obj.m/2;
        end

        function n = norm(obj,type)
            %NORM Norm - four types:
            % bulk norm
            % weight norm
            % geometric norm
            % projected geometric norm
            arguments
                obj rga
                type string = "geom"
            end
            switch type
                case "bulk"
                    n2 = dot(obj,obj);
                    n = sqrt(n2.m(1));
                case "weight"
                    n2 = antidot(obj,obj);
                    n = sqrt(n2.m(end));
                case "geom"
                    n = norm(obj,"bulk") + norm(obj,"weight");
                case "projgeom"
                    n = norm(obj,"bulk")/norm(obj,"weight");
                otherwise
                    error('type not recognized')
            end
        end

        function d = dist(a,b)
            %DIST Euclidean distance between a & b, if defn'd.
            if (isa(a,'rgapoint') && isa(b,'rgapoint')) || ...
                    (isa(a,'rgapoint') && isa(b,'rgaplane')) || ...
                    (isa(a,'rgaplane') && isa(b,'rgapoint'))
                d = norm(commutate(a,b,"-^"),'weight') ...
                    /norm(commutate(a,b,"+v"),'weight');
            elseif (isa(a,'rgapoint') && isa(b,'rgaline')) || ...
                    (isa(a,'rgaline') && isa(b,'rgapoint'))
                d = norm(commutate(a,b,"+^"),'weight') ...
                    /norm(commutate(a,b,"+v"),'weight');
            elseif isa(a,'rgaline') && isa(b,'rgaline')
                d = norm(commutate(a,b,"+^"),'weight') ...
                    /norm(commutate(a,b,"-v"),'weight');
            else
                error('Euclidean distance not defined for these inputs')
            end
        end

        function obj = bulkrc(obj)
            %BULKRC Bulk right complement
            if isempty(obj)
                [~,B] = productmat(rga('e1234'),'wedgedot');
                obj = B*rev(obj);
            else
                obj = wedgedot(rev(obj),rga('e1234'));
            end
        end

        function obj = bulklc(obj)
            %BULKLC Bulk left complement
            if isempty(obj)
                [~,B] = productmat(rga('e1234'),'wedgedot');
                obj = B*antirev(obj);
            else
                obj = wedgedot(antirev(obj),rga('e1234'));
            end
        end

        function obj = weightrc(obj)
            %WEIGHTRC Weight right complement
            if isempty(obj)
                obj = antiwedgedot(rga('e0'))*rev(obj);
            else
                obj = antiwedgedot(rga('e0'),rev(obj));
            end
        end

        function obj = weightlc(obj)
            %WEIGHTLC Weight left complement
            if isempty(obj)
                obj = antiwedgedot(rga('e0'))*antirev(obj);
            else
                obj = antiwedgedot(rga('e0'),antirev(obj));
            end
        end

        function X = screw(M,X)
            %SCREW Translate and rotate the object X using the motor M
            arguments
                M (1,1) rgamotor
                X (1,1) rga
            end
            clX = class(X);
            anX = X.anti;
            M = unitize(M);
            if isa(X,'rgapoint')
                [A,B] = productmat(M,'motor');
                x = X.m(2:5);
                X = rga([0; (A+B)*x(:); zeros(11,1)]);
            else
                X = antiwedgedot(antiwedgedot(M,X),antirev(M));
            end
            X = feval(clX,X);
            if anX
                X.anti = true;
            end
        end

        function X = unscrew(M,X)
            %UNSCREW Untranslate and unrotate the object X using the motor M
            arguments
                M (1,1) rgamotor
                X (1,1) rga
            end
            clX = class(X);
            anX = X.anti;
            M = unitize(M);
            if isa(X,'rgapoint')
                [A,B] = productmat(M,'motor');
                x = X.m(2:5);
                X = rga([0; (A-B)*x(:); zeros(11,1)]);
            else
                X = antiwedgedot(antiwedgedot(antrev(M),X),M);
            end
            X = feval(clX,X);
            if anX
                X.anti = true;
            end
        end

        function obj = zapeps(obj,tol)
            %ZAPEPS Remove (zap) multivector components less than eps or tol
            arguments
                obj
                tol = 1e-16;
            end
            obj.m(abs(obj.m)<tol) = 0;
        end

        function dstr = char(obj)
            % CHAR String representation of object
            if obj.anti
                es = char(949);
            else
                es = 'e';
            end
            e0 = char([es,8320]);
            e1 = char([es,8321]);
            e2 = char([es,8322]);
            e3 = char([es,8323]);
            e4 = char([es,8324]);
            e23 = char([es,8322,8323]);
            e31 = char([es,8323,8321]);
            e12 = char([es,8321,8322]);
            e43 = char([es,8324,8323]);
            e42 = char([es,8324,8322]);
            e41 = char([es,8324,8321]);
            e321 = char([es,8323,8322,8321]);
            e412 = char([es,8324,8321,8322]);
            e431 = char([es,8324,8323,8321]);
            e423 = char([es,8324,8322,8323]);
            e1234 = char([es,8321,8322,8323,8324]);
            if obj.anti
                obj.m = flipud(obj.m(:));
            end
            e = char(e0,e1,e2,e3,e4,e23,e31,e12,e43,e42,e41,e321,e412,e431,e423,e1234);
            c_ = @(mi) char(strtrim(formattedDisplayText(mi)));
            if isnumeric(obj.m)
                c = c_;
                s = sign(obj.m);
                z = 0;
            else
                c = @(mi) char(['(',c_(mi),')']);
                s = ones(16,1);
                z = sym(0);
            end
            dstr = '';
            for i = 1:16
                mi = obj.m(i);
                if isequal(mi,z)
                    continue
                end
                if isnumeric(obj.m(i)) && abs(obj.m(i)) <= obj.dispeps
                    continue
                end
                if i>1 && s(i)>0
                    pm = '+ ';
                else
                    pm = '';
                end
                cmi = c(mi);
                dstr = [dstr,pm,cmi,e(i,:)]; %#ok<*AGROW>
                if i == 1 || i== 5 || i== 11 || i==15
                    dstr = [dstr,'\n'];
                end
                if length(cmi)>40
                    dstr = [dstr,'\n'];
                end
            end
        end

        function disp(obj)
            % DISP Display object on command line
            if numel(obj) == 1
                dstr = char(obj);
                if isempty(dstr)
                    fprintf(formattedDisplayText(0))
                else
                    fprintf([dstr,'\n'])
                end
            else
                fprintf(blanks(4))
                fprintf([repmat('%dx',1,ndims(obj)-1),'%d '],size(obj))
                fprintf('%s',class(obj))
                fprintf('\n')
            end
        end

    end

    methods (Hidden)
        function A = createprodmat(a,fnhandle,commute)
            %CREATEPRODMAT Create product matrices c = A(a)*b or c = B(b)*a
            e = cell(1,16);
            [e{:}] = rga.bases;
            for i = 16:-1:1
                if commute % c = B(b)*a
                    A(:,i) = fnhandle(e{i},a).m(:);
                else % c = A(a)*b
                    A(:,i) = fnhandle(a,e{i}).m(:);
                end
            end
        end
    end

    methods (Static)

        function [e0,e1,e2,e3,e4,e23,e31,e12,e43,e42,e41,e321,e412,e431,e423,e1234] = bases
            %BASES Create full set of basis elements for RGA
            e0 = rga('e0');
            e1 = rga('e1'); e2 = rga('e2'); e3 = rga('e3'); e4 = rga('e4');
            e23 = rga('e23'); e31 = rga('e31'); e12 = rga('e12');
            e43 = rga('e43'); e42 = rga('e42'); e41 = rga('e41');
            e321 = rga('e321'); e412 = rga('e412'); e431 = rga('e431'); e423 = rga('e423');
            e1234 = rga('e1234');
        end

        function [eps0,eps1,eps2,eps3,eps4,eps23,eps31,eps12,eps43,eps42,eps41,eps321,eps412,eps431,eps423,eps1234] = antibases
            %ANTIBASES Create full set of anti-basis elements for RGA
            [eps1234,eps423,eps431,eps412,eps321,eps41,eps42,eps43,...
                eps12,eps31,eps23,eps4,eps3,eps2,eps1,eps0] = rga.bases;
            eps1234.anti = true;
            eps423.anti = true; eps431.anti = true; eps412.anti = true;
            eps321.anti = true; 
            eps41.anti = true; eps42.anti = true; eps43.anti = true; 
            eps12.anti = true; eps31.anti = true; eps23.anti = true; 
            eps4.anti = true; eps3.anti = true; eps2.anti = true; 
            eps1.anti = true; eps0.anti = true;
        end

        function producttab(prod)
            % PRODUCTTAB Display basis element product table
            % Default is wedgedot but inputs can be "wedge", "antiwedge",
            % "antiwedgedot" as well.
            arguments
                prod string = ""
            end
            a = eye(1,16);
            b = eye(1,16);
            mtab = cell(16,16);
            for i = 1:16
                for j = 1:16
                    switch prod
                        case 'rint'
                            mtab{i,j} = strtrim(char(rint(rga(a),rga(b))));
                        case 'lint'
                            mtab{i,j} = strtrim(char(lint(rga(a),rga(b))));
                        case 'wedge'
                            mtab{i,j} = strtrim(char(wedge(rga(a),rga(b))));
                        case 'antiwedge'
                            mtab{i,j} = strtrim(char(antiwedge(rga(a),rga(b))));
                        case 'antiwedgedot'
                            mtab{i,j} = strtrim(char(antiwedgedot(rga(a),rga(b))));
                        otherwise
                            mtab{i,j} = strtrim(char(wedgedot(rga(a),rga(b))));
                    end
                    mtab{i,j}(isspace(mtab{i,j})) = [];
                    if isempty(mtab{i,j})
                        mtab{i,j} = ' 0';
                    end
                    if strcmp(mtab{i,j}(1),'1')
                        mtab{i,j} = ['+',mtab{i,j}];
                    end
                    mtab{i,j}(2) = [];
                    b = circshift(b,1);
                end
                a = circshift(a,1);
            end
            s = replace(string(mtab),"\n","");
            % have to transpose b/c fprintf runs down the columns:
            fprintf('%6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s\n',s')
        end

        function h = plotaxes
            %PLOTAXES Create orthogonal axes in plotspace
            h = axes('Projection','Perspective');%,...
            %'XLim',[-1 1],'YLim',[-1 1],'ZLim',[-1 1]);
            axis equal
            hold on
            view(37,30)
            grid on
            quiver3(-1,0,0,1,0,0,2,'color','#A2142F')
            quiver3(0,-1,0,0,1,0,2,'color','#77AC30')
            quiver3(0,0,-1,0,0,1,2,'color','#0072BD')
        end
    end

    methods (Access = protected)

        function c = dotReference(obj,indexOp)
            switch indexOp.Name
                case {'e0','eps1234'}
                    c = obj.m(1);
                case {'e1','eps423'}
                    c = obj.m(2);
                case {'e2','eps431'}
                    c = obj.m(3);
                case {'e3','eps412'}
                    c = obj.m(4);
                case {'e4','eps321'}
                    c = obj.m(5);
                case {'e23','eps41'}
                    c = obj.m(6);
                case {'e31','eps42'}
                    c = obj.m(7);
                case {'e12','eps43'}
                    c = obj.m(8);
                case {'e43','eps12'}
                    c = obj.m(9);
                case {'e42','eps31'}
                    c = obj.m(10);
                case {'e41','eps23'}
                    c = obj.m(11);
                case {'e321','eps4'}
                    c = obj.m(12);
                case {'e412','eps3'}
                    c = obj.m(13);
                case {'e431','eps2'}
                    c = obj.m(14);
                case {'e423','eps1'}
                    c = obj.m(15);
                case {'e1234','eps0'}
                    c = obj.m(16);
                otherwise
                    error('basis not recognized')
            end
        end

        function obj = dotAssign(obj,indexOp,c)
            switch indexOp.Name
                case {'e0','eps1234'}
                    obj.m(1) = c;
                case {'e1','eps423'}
                    obj.m(2) = c;
                case {'e2','eps431'}
                    obj.m(3) = c;
                case {'e3','eps412'}
                    obj.m(4) = c;
                case {'e4','eps321'}
                    obj.m(5) = c;
                case {'e23','eps41'}
                    obj.m(6) = c;
                case {'e31','eps42'}
                    obj.m(7) = c;
                case {'e12','eps43'}
                    obj.m(8) = c;
                case {'e43','eps12'}
                    obj.m(9) = c;
                case {'e42','eps31'}
                    obj.m(10) = c;
                case {'e41','eps23'}
                    obj.m(11) = c;
                case {'e321','eps4'}
                    obj.m(12) = c;
                case {'e412','eps3'}
                    obj.m(13) = c;
                case {'e431','eps2'}
                    obj.m(14) = c;
                case {'e423','eps1'}
                    obj.m(15) = c;
                case {'e1234','eps0'}
                    obj.m(16) = c;
                otherwise
                    error('basis not recognized')
            end
        end

        function n = dotListLength(obj,indexOp,indexContext) %#ok<INUSD> 
            n = 1;
        end

    end

end