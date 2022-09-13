classdef (InferiorClasses = {?sym}) ga3
    % G(3,0,0) 
    %   Geometric Algebra for ordinary 3D space w/Lengyal's basis elements.

    properties
        m % multivector represented as [e0 e1 e2 e3 e23 e31 e12 e321]
    end

    methods
        function obj = ga3(v)
            %GA3 Construct a ga3 object
            %   obj = ga3([e0 e1 e2 e3 e23 e31 e12 e321])
            %   Symbolic inputs ok using above syntax.
            %   obj = 'e0' or "e0" etc. also recognized for single elements
            %   To create objects like 5e0 - 3e12 etc. first create a set
            %   of basis elements using GA3BASES (but not for symbolic
            %   coefficients - use symbolic array syntax as above).
            if nargin == 0
                obj.m = zeros(1,8);
            elseif ischar(v) || isstring(v)
                m = eye(1,8);
                switch v
                    case 'e0'
                        obj.m = m;
                    case 'e1'
                        obj.m = circshift(m,1);
                    case 'e2'
                        obj.m = circshift(m,2);
                    case 'e3'
                        obj.m = circshift(m,3);
                    case 'e23'
                        obj.m = circshift(m,4);
                    case 'e31'
                        obj.m = circshift(m,5);
                    case 'e12'
                        obj.m = circshift(m,6);
                    case 'e321'
                        obj.m = circshift(m,7);
                    otherwise
                        error('string or char input not recognized')
                end
            else
                obj.m = v;
            end
        end

        function obj = scalar(obj)
            %SCALAR Return scalar part
            obj.m(2:end) = 0;
        end

        function obj = vector(obj)
            %SCALAR Return vector part
            obj.m([1 5:end]) = 0;
        end

        function obj = bivector(obj)
            %SCALAR Return bivector part
            obj.m([1:4 8]) = 0;
        end

        function obj = antivector(obj)
            %ANTIVECTOR Return antivector part
            obj = bivector(obj);
        end

        function obj = pseudoscalar(obj)
            %SCALAR Return pseudoscalar part
            obj.m(1:7) = 0;
        end

        function obj = antiscalar(obj)
            %ANTISCALAR Return antiscalar part
            obj = pseudoscalar(obj);
        end

        function obj = plus(a,b)
            %PLUS GA3 addition
            % If only one input is a multivector, then this adds the other
            % input to the scalar component of the multivector.
            obj = ga3;
            aga3 = isa(a,"ga3");
            bga3 = isa(b,"ga3");
            if aga3 && bga3
                obj.m = a.m + b.m;
            elseif ~aga3
                obj.m = b.m;
                obj.m(1) = a + obj.m(1);
            elseif ~bga3
                obj.m = a.m;
                obj.m(1) = obj.m(1) + b;
            else
                error('Inputs incompatible')
            end
        end

        function obj = uminus(obj)
            %UMINS GA3 unary minus
            obj.m = -obj.m;
        end

        function obj = minus(a,b)
            %MINUS GA3 subtraction
            obj = plus(a,-b);
        end

        function obj = mpower(a,b)
            %MPOWER GA3 wedge product
            aga3 = isa(a,"ga3");
            bga3 = isa(b,"ga3");
            if aga3 && bga3
                obj = a*b - a.*b;
            else
                error('Inputs incompatible')
            end
        end

        function obj = times(a,b)
            %TIMES GA3 dot product
            aga3 = isa(a,"ga3");
            bga3 = isa(b,"ga3");
            if aga3 && bga3
                M = blkdiag(eye(4),-eye(4));
                s = a.m(:)'*M*b.m(:);
                obj = ga3(s*eye(8,1));
            else
                error('Inputs incompatible')
            end
        end

        function obj = mtimes(a,b)
            %MTIMES GA3 wedgedot product
            obj = ga3;
            aga3 = isa(a,"ga3");
            bga3 = isa(b,"ga3");
            if aga3 && bga3
                A = wedgedotmat(a);
                obj.m = A*b.m(:);
            elseif ~aga3
                obj.m = a .* b.m;
            elseif ~bga3
                obj.m = b .* a.m;
            else
                error('Inputs incompatible')
            end
        end

        function obj = mldivide(a,b)
            %MLDIVIDE c = a\b for GA3
            obj = ga3;
            aga3 = isa(a,"ga3");
            bga3 = isa(b,"ga3");
            if aga3 && bga3
                A = wedgedotmat(a);
                obj.m = A\b.m(:);
            elseif ~aga3
                obj.m = a .\ b.m;
            elseif ~bga3
                obj.m = b .\ a.m;
            else
                error('Inputs incompatible')
            end            
        end

        function obj = inv(obj)
            %INV Inverse of the input
            obj = obj\ga3(eye(8,1));
        end

        function obj = mrdivide(a,b)
            %MRDIVIDE c = a/b for GA3
            obj = ga3;
            aga3 = isa(a,"ga3");
            bga3 = isa(b,"ga3");
            if aga3 && bga3
                obj= a*inv(b); %#ok<MINV> 
            elseif ~aga3
                obj.m = a ./ b.m;
            elseif ~bga3
                obj = b*inv(a); %#ok<MINV> 
            else
                error('Inputs incompatible')
            end            
        end

        function g = grade(obj)
            %GRADE Number of basis vectors present in factorization
            nz = find(obj.m ~= 0);
            g = 0;
            if any(ismember(nz,[2 3 4])), g = 1; end
            if any(ismember(nz,[5 6 7])), g = 2; end
            if any(ismember(nz,8)), g = 3; end
        end

        function ag = antigrade(obj)
            %ANTIGRADE Number of basis vectors absent in factorization
            ag = 3 - grade(obj);
        end

        function obj = reverse(obj)
            %REVERSE
            obj.m(5:8) = -obj.m(5:8);
        end

        function dstr = char(obj)
            % CHAR String representation of object
            e0 = char(['e',8320]);
            e1 = char(['e',8321]);
            e2 = char(['e',8322]);
            e3 = char(['e',8323]);
            e23 = char(['e',8322,8323]);
            e31 = char(['e',8323,8321]);
            e12 = char(['e',8321,8322]);
            e321 = char(['e',8323,8322,8321]);
            e = char(e0,e1,e2,e3,e23,e31,e12,e321);
            c_ = @(mi) char(strtrim(formattedDisplayText(mi)));
            if isnumeric(obj.m)
                c = c_;
                s = sign(obj.m);
                z = 0;
            else
                c = @(mi) char(['(',c_(mi),')']);
                s = ones(8,1);
                z = sym(0);
            end
            dstr = '';
            for i = 1:8
                mi = obj.m(i);
                if isequal(mi,z)
                    continue
                end
                if i>1 && s(i)>0
                    pm = '+ ';
                else
                    pm = '';
                end
                cmi = c(mi);
                dstr = [dstr,pm,cmi,e(i,:)]; %#ok<*AGROW> 
                if length(cmi)>40
                    dstr = [dstr,'\n'];
                end
            end
        end

        function disp(obj)
            % DISP Display object on command line
            dstr = char(obj);
            if isempty(dstr)
                fprintf(formattedDisplayText(0))
            else
                fprintf([dstr,'\n'])
            end
        end

    end

    methods (Static)
        
        function wedgedottab
            % WEDGEDOTTAB Display basis element wedgedot table
            a = eye(1,8);
            b = eye(1,8);
            mtab = cell(8,8);
            for i = 1:8
                for j = 1:8
                    mtab{i,j} = strtrim(char(ga3(a)*ga3(b)));
                    mtab{i,j}(isspace(mtab{i,j})) = [];
                    if strcmp(mtab{i,j}(1),'1')
                        mtab{i,j} = ['+',mtab{i,j}];
                    end
                    mtab{i,j}(2) = [];
                    b = circshift(b,1);
                end
                a = circshift(a,1);
            end
            % have to transpose b/c fprintf runs down the columns:
            fprintf('%6s %6s %6s %6s %6s %6s %6s %6s\n',string(mtab)')
        end

        function wedgetab
            % WEDGETAB Display basis element wedge table
            a = eye(1,8);
            b = eye(1,8);
            mtab = cell(8,8);
            for i = 1:8
                for j = 1:8
                    mtab{i,j} = strtrim(char(ga3(a)^ga3(b)));
                    mtab{i,j}(isspace(mtab{i,j})) = [];
                    if isempty(mtab{i,j})
                        mtab{i,j} = '0';
                    else
                        if strcmp(mtab{i,j}(1),'1')
                            mtab{i,j} = ['+',mtab{i,j}];
                        end
                        mtab{i,j}(2) = [];
                    end
                    b = circshift(b,1);
                end
                a = circshift(a,1);
            end
            % have to transpose b/c fprintf runs down the columns:
            fprintf('%6s %6s %6s %6s %6s %6s %6s %6s\n',string(mtab)')
        end

        function revtable
            % REVTABLE Display basis element reverse table
            a = eye(1,8);
            revtable = cell(1,8);
            for i = 1:8
                ga3a = ga3(a);
                revtable{1,i} = strtrim(char(reverse(ga3a)));
                revtable{1,i}(isspace(revtable{1,i})) = [];
                if strcmp(revtable{1,i}(1),'1')
                    revtable{1,i} = ['+',revtable{1,i}];
                end
                revtable{1,i}(2) = [];
                a = circshift(a,1);
            end
            % have to transpose b/c fprintf runs down the columns:
            fprintf('%6s %6s %6s %6s %6s %6s %6s %6s\n',string(revtable)')
        end

    end

    methods (Hidden)

        function M = wedgedotmat(obj)
            % WEDGEDOTMAT Matrix for wedgedot of two multivectors
            I = eye(3);
            s = obj.m(1); v = obj.m(2:4); b = obj.m(5:7); p = obj.m(8);
            v = v(:); b = b(:);
            M = [s v' -b' -p;
                v,  s*I-skew(b), p*I-skew(v),  b;
                b, -p*I+skew(v), s*I-skew(b), -v;
                p -b' -v' s];
        end

    end
end