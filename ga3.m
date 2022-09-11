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

        function obj = plus(a,b)
            %PLUS GA3 addition
            obj = ga3;
            aga3 = isa(a,"ga3");
            bga3 = isa(b,"ga3");
            if aga3 && bga3
                obj.m = a.m + b.m;
            elseif ~aga3
                obj.m = a + b.m;
            elseif ~bga3
                obj.m = a.m + b;
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

        function obj = mtimes(a,b)
            %MTIMES GA3 wedgedot product
            obj = ga3;
            aga3 = isa(a,"ga3");
            bga3 = isa(b,"ga3");
            if aga3 && bga3
                I = eye(3);
                as = a.m(1); av = a.m(2:4); ab = a.m(5:7); ap = a.m(8);
                av = av(:); ab = ab(:);
                X = [as av' -ab' -ap;
                    av,  as*I-skew(ab), ap*I-skew(av),  ab;
                    ab, -ap*I+skew(av), as*I-skew(ab), -av;
                    ap -ab' -av' as];
                obj.m = X*b.m(:);
            elseif ~aga3
                obj.m = a .* b.m;
            elseif ~bga3
                obj.m = b .* a.m;
            else
                error('Inputs incompatible')
            end
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
        
        function mtable
            % MTABLE Display basis element multiplication table
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

    end
end