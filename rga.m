classdef (InferiorClasses = {?sym}) rga
    % G(3,0,1) 
    %   Geometric Algebra for homogeneous 4D space w/Lengyal's basis elements.

    properties
        m (1,16) % multivector [e0 e1 e2 e3 e4 e23 e31 e12 e43 e42 e41 e321 e412 e431 e423 e1234]
        anti = false % set to true for anti
    end

    methods
        function obj = rga(v,anti)
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
            if nargin == 2
                obj.anti = anti;
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
            obj.anti = true;
        end

        function obj = pseudoscalar(obj)
            %SCALAR Return pseudoscalar part
            obj.m(1:15) = 0;
        end

        function obj = antiscalar(obj)
            %ANTISCALAR Return antiscalar part
            obj = pseudoscalar(obj);
            obj.anti = true;
        end

        function obj = bulk(obj)
            %BULK Bulk - portion not containing factors of e4
            obj.m([5 9:11 13:end]) = 0;
        end

        function obj = weight(obj)
            %WEIGHT Weight - portion containing factors of e4
            obj.m([1:4 6:8 12]) = 0;
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
            obj.m(6:15) = -obj.m(6:15);
        end

        function obj = antirev(obj)
            %ANTIREV Anti-reverse
            obj.m(2:11) = -obj.m(2:11);
        end

        function obj = not(obj)
            % Overload ~ for reverse & antireverse
            if obj.anti
                obj = antirev(obj);
            else
                obj = rev(obj);
            end
        end

        function obj = rcomp(obj)
            %RCOMP Right complement
            obj = rev(obj);
            obj.m = flipud(obj.m(:));
        end

        function obj = lcomp(obj)
            %LCOMP Left complement
            obj = antirev(obj);
            obj.m = flipud(obj.m(:));
        end

        function obj = plus(a,b)
            %PLUS RGA addition
            % If only one input is a multivector, then this adds the other
            % input to the scalar component of the multivector.
            obj = rga;
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

        function obj = dot(a,b)
            %DOT RGA dot product
            arga = isa(a,"rga");
            brga = isa(b,"rga");
            if arga && brga
                M = diag([1 1 1 1 0 -1 -1 -1 0 0 0 -1 0 0 0 0]);
                s = a.m(:)'*M*b.m(:);
                obj = rga(s*eye(16,1));
            else
                error('Inputs incompatible')
            end
        end

        function obj = antidot(a,b)
            %ANTIDOT RGA anti-dot product
            arga = isa(a,"rga");
            brga = isa(b,"rga");
            if arga && brga
                M = diag(fliplr([1 1 1 1 0 -1 -1 -1 0 0 0 -1 0 0 0 0]));
                s = a.m(:)'*M*b.m(:);
                obj = rga(s*flipud(eye(16,1)));
            else
                error('Inputs incompatible')
            end
        end

        function obj = times(a,b)
            % Overload .* as dot product
            if isa(a,"rga") && a.anti && isa(b,"rga") && b.anti
                obj = antidot(a,b);
            else
                obj = dot(a,b);
            end
        end

        function obj = wedgedot(a,b)
            %WEDGEDOT RGA wedgedot product
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

        function obj = mtimes(a,b)
            % Overload * as wedgedot product
            if isa(a,"rga") && a.anti && isa(b,"rga") && b.anti
                obj = antiwedgedot(a,b);
            else
                obj = wedgedot(a,b);
            end
        end

        function obj = wedge(a,b)
            %WEDGE RGA wedge product
            obj = rga;
            arga = isa(a,"rga");
            brga = isa(b,"rga");
            if arga && brga
                [~,A] = wedgedotmat(a);
                obj.m = A*b.m(:);
            elseif ~arga
                obj.m = a .* b.m;
            elseif ~brga
                obj.m = b .* a.m;
            else
                error('Inputs incompatible')
            end
        end

        function obj = mpower(a,b)
            % Overload ^ as wedge product
            if isa(a,"rga") && a.anti && isa(b,"rga") && b.anti
                obj = antiwedge(a,b);
            else
                obj = wedge(a,b);
            end
        end

        function obj = antiwedge(a,b)
            %ANTIWEDGE Anti-wedge product
            obj = lcomp(wedge(rcomp(a),rcomp(b)));
        end

        function obj = antiwedgedot(a,b)
            %ANTIWEDGEDOT Anti-wedgedot product
            obj = lcomp(wedgedot(rcomp(a),rcomp(b)));
        end

        function obj = rint(a,b)
            %RINT Right interior product
            obj = antiwedge(a,rcomp(b));
        end

        function obj = lint(a,b)
            %LINT Left interior product
            obj = antiwedge(lcomp(a),b);
        end

        function obj = antirint(a,b)
            %ANTIRINT Right interior antiproduct
            obj = wedge(a,rcomp(b));
        end

        function obj = antilint(a,b)
            %ANTILINT Left interior antiproduct
            obj = wedge(lcomp(a),b);
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
                case "-^"
                    obj = (wedge(a,rev(b)) - wedge(b,rev(a)));
                case "+^"
                    obj = (wedge(a,rev(b)) + wedge(b,rev(a)));
                case "-v"
                    obj = (antiwedge(a,rev(b)) - antiwedge(b,rev(a)));
                case "+v"
                    obj = (antiwedge(a,rev(b)) + antiwedge(b,rev(a)));
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
                    n2 = dot(obj,rev(obj));
                    n = sqrt(n2.m(1));
                case "weight"
                    n2 = antidot(obj,antirev(obj));
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
            %DIST Euclidean distance between a & b
            d = norm(commutate(a,b,"-^"),'weight') ...
                /norm(commutate(a,b,"+v"),'weight');
        end

        function obj = bulkrc(obj)
            %BULKRC Bulk right complement
            obj = wedgedot(rev(obj),rga('e1234'));
        end

        function obj = bulklc(obj)
            %BULKLC Bulk left complement
            obj = wedgedot(antirev(obj),rga('e1234'));
        end

        function obj = weightrc(obj)
            %WEIGHTRC Weight right complement
            obj = antiwedgedot(rga('e0'),rev(obj));
        end

        function obj = weightlc(obj)
            %WEIGHTLC Weight left complement
            obj = antiwedgedot(rga('e0'),antirev(obj));
        end

        function P = screw(M,P)
            %SCREW Translate & rotate the point P using the motor M
            [A,B] = productmat(M,'motor');
            x = P.m(2:5);
            P = rga([0; (A+B)*x(:); zeros(11,1)]);
        end

        function P = unscrew(M,P)
            %UNSCREW Untranslate & unrotate the point P using the motor M
            [A,B] = productmat(M,'motor');
            x = P.m(2:5);
            P = rga([0; (A-B)*x(:); zeros(11,1)]);
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

    end

end