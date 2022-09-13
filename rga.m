classdef (InferiorClasses = {?sym}) rga
    % G(3,0,1) 
    %   Geometric Algebra for homogeneous 4D space w/Lengyal's basis elements.

    properties
        m % multivector [e0 e1 e2 e3 e4 e23 e31 e12 e43 e42 e41 e321 e412 e431 e423 e1234]
    end

    methods
        function obj = rga(v)
            %RGA Construct a rga object
            %   obj = rga([e0 e1 e2 e3 e4 e23 e31 e12 e43 e42 e41 e321 e412 e431 e423 e1234])
            %   Symbolic inputs ok using above syntax.
            %   obj = 'e0' or "e0" etc. also recognized for single elements
            %   To create objects like 5e0 - 3e12 etc. first create a set
            %   of basis elements using RGABASES (but not for symbolic
            %   coefficients - use symbolic array syntax as above).
            if nargin == 0
                obj.m = zeros(1,16);
            elseif ischar(v) || isstring(v)
                m = eye(1,16);
                switch v
                    case 'e0'
                        obj.m = m;
                    case 'e1'
                        obj.m = circshift(m,1);
                    case 'e2'
                        obj.m = circshift(m,2);
                    case 'e3'
                        obj.m = circshift(m,3);
                    case 'e4'
                        obj.m = circshift(m,4);
                    case 'e23'
                        obj.m = circshift(m,5);
                    case 'e31'
                        obj.m = circshift(m,6);
                    case 'e12'
                        obj.m = circshift(m,7);
                    case 'e43'
                        obj.m = circshift(m,8);
                    case 'e42'
                        obj.m = circshift(m,9);
                    case 'e41'
                        obj.m = circshift(m,10);
                    case 'e321'
                        obj.m = circshift(m,11);
                    case 'e412'
                        obj.m = circshift(m,12);
                    case 'e431'
                        obj.m = circshift(m,13);
                    case 'e423'
                        obj.m = circshift(m,14);
                    case 'e1234'
                        obj.m = circshift(m,15);
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
            obj.m([1 6:end]) = 0;
        end

        function obj = bivector(obj)
            %SCALAR Return bivector part
            obj.m([1:5 12:end]) = 0;
        end

        function obj = trivector(obj)
            %TRIVECTOR Return trivector part
            obj.m([1:11 16]) = 0;
        end

        function obj = antivector(obj)
            %ANTIVECTOR Return antivector part
            obj = trivector(obj);
        end

        function obj = pseudoscalar(obj)
            %SCALAR Return pseudoscalar part
            obj.m(1:15) = 0;
        end

        function obj = antiscalar(obj)
            %ANTISCALAR Return antiscalar part
            obj = pseudoscalar(obj);
        end

        function obj = plus(a,b)
            %PLUS RGA addition
            % If only one input is a multivector, then this adds the other
            % input to the scalar component of the multivector.
            obj = ga3;
            arga = isa(a,"rga");
            brga = isa(b,"rga");
            if arga && brga
                obj.m = a.m + b.m;
            elseif ~arga
                obj.m = b.m;
                obj.m(1) = a + obj.m(1);
            elseif ~brga
                obj.m = a.m;
                obj.m(1) = obj.m(1) + b;
            else
                error('Inputs incompatible')
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

        function obj = mpower(a,b)
            %MPOWER RGA wedge product
            arga = isa(a,"rga");
            brga = isa(b,"rga");
            if arga && brga
                obj = a*b - a.*b;
            else
                error('Inputs incompatible')
            end
        end

        function obj = times(a,b)
            %TIMES RGA dot product
            arga = isa(a,"rga");
            brga = isa(b,"rga");
            if arga && brga
                M = blkdiag(eye(8),-eye(8));
                s = a.m(:)'*M*b.m(:);
                obj = rga(s*eye(16,1));
            else
                error('Inputs incompatible')
            end
        end

        function obj = mtimes(a,b)
            %MTIMES RGA wedgedot product
            obj = rga;
            arga = isa(a,"rga");
            brga = isa(b,"rga");
            if arga && brga
                A = wedgedotmat(a);
                obj.m = A*b.m(:);
            elseif ~arga
                obj.m = a .* b.m;
            elseif ~brga
                obj.m = b .* a.m;
            else
                error('Inputs incompatible')
            end
        end

        function g = grade(obj)
            %GRADE Number of basis vectors present in factorization
            nz = find(obj.m ~= 0);
            g = 0;
            if any(ismember(nz,[2 3 4 5])), g = 1; end
            if any(ismember(nz,[6 7 8 9 10 11])), g = 2; end
            if any(ismember(nz,[12 13 14 15])), g = 3; end
            if any(ismember(nz,16)), g = 4; end
        end

        function ag = antigrade(obj)
            %ANTIGRADE Number of basis vectors absent in factorization
            ag = 4 - grade(obj);
        end

        function obj = reverse(obj)
            %REVERSE
            obj.m(6:15) = -obj.m(6:16);
        end

        function obj = antireverse(obj)
            %ANTIREVERSE
            obj.m(2:11) = -obj.m(1:11);
        end

        function dstr = char(obj)
            % CHAR String representation of object
            e0 = char(['e',8320]);
            e1 = char(['e',8321]);
            e2 = char(['e',8322]);
            e3 = char(['e',8323]);
            e4 = char(['e',8324]);
            e23 = char(['e',8322,8323]);
            e31 = char(['e',8323,8321]);
            e12 = char(['e',8321,8322]);
            e43 = char(['e',8324,8323]);
            e42 = char(['e',8324,8322]);
            e41 = char(['e',8324,8321]);
            e321 = char(['e',8323,8322,8321]);
            e412 = char(['e',8324,8321,8322]);
            e431 = char(['e',8324,8323,8321]);
            e423 = char(['e',8324,8322,8323]);
            e1234 = char(['e',8321,8322,8323,8324]);
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
            a = eye(1,16);
            b = eye(1,16);
            mtab = cell(16,16);
            for i = 1:16
                for j = 1:16
                    mtab{i,j} = strtrim(char(rga(a)*rga(b)));
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

        function wedgetab
            % WEDGETAB Display basis element wedge table
            a = eye(1,16);
            b = eye(1,16);
            mtab = cell(16,16);
            for i = 1:16
                for j = 1:16
                    mtab{i,j} = strtrim(char(rga(a)^rga(b)));
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
            fprintf('%6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s\n',string(mtab)')
        end

    end

    methods (Hidden)

    end
end