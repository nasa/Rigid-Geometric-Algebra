classdef rgaplane < rga
    %RGAPLANE RGA Plane
    %   Subclass of RGA
    %   In RGA, a plane in 3D space is a four-element trivector whose
    %   elements specify the plane's normal and its distance from the
    %   origin.  It is the projection of a 4D volume spanned by the
    %   trivector bases onto the onto the plane with projective coordinate
    %   w=1.  A plane is dual to a point in that when expressed in the dual
    %   or anti basis elements, it looks like a point.  For this reason,
    %   the constructor defaults to anti basis.  A unitized point has unit
    %   normal vector.
    methods
        function obj = rgaplane(x,y,z,w)
            %RGAPLANE Construct RGA plane from its normal
            % F = rgaplane creates a random RGA plane object
            % F = rgaplane(M) subclasses the RGA object M into a plane
            % F = rgaplane(nx,ny,nz) creates a plane with normal n at w=1
            % F = rgaplane(nx,ny,nz,w) creates a plane at (nx,ny,nz,w)
            % Defaults to anti basis.
            switch nargin
                case 0
                    obj.m = fliplr([0 randn(1,4) zeros(1,11)]);
                case 1
                    if isa(x,"rga")
                        t = trivector(x);
                        obj.m = t.m;
                        obj.anti = x.anti;
                        obj.dispeps = x.dispeps;
                    else
                        obj.m = fliplr([0 x(:)' zeros(1,11)]);
                    end
                case 3
                    obj.m = fliplr([0 x y z 1 zeros(1,11)]);
                case 4
                    obj.m = fliplr([0 x y z w zeros(1,11)]);
                otherwise
                    error('Inputs not compatible')
            end
            obj.anti = true;
        end

        function obj = unitize(obj)
            %UNITIZE Unitize the plane
            obj.m(13:15) = obj.m(13:15)/norm(obj.m(13:15));
        end

        function h = plot(obj)
            %PLOT Plot the plane
            % use ax + by + cz + d = 0
            a = obj.m(15); b = obj.m(14); c = obj.m(13); d = obj.m(12);
            if c==0 && b==0 % just plot yz plane
                h = patch(-d/a*[1 1 1 1],[1 1 -1 -1],[1 -1 -1 1],'c','FaceAlpha',0.5);
            elseif c==0 && a==0 % just plot xz plane
                h = patch([1 -1 -1 1],-d/b*[1 1 1 1],[1 1 -1 -1],'c','FaceAlpha',0.5);
            elseif b==0 && a==0 % just plot xy plane
                h = patch([1 -1 -1 1],[1 1 -1 -1],-d/c*[1 1 1 1],'c','FaceAlpha',0.5);
            else
                x = [1 -1 -1 1];
                y = [1 1 -1 -1];
                z = -1/c*(a*x + b*y + d);
                h = patch(x,y,z,'c','FaceAlpha',0.5);
            end
            v = h.Vertices;
            c = h.EdgeColor;
            hold on
            for i = 1:4
                j = mod(i,4)+1;
                hq(i) = quiver3(v(i,1),v(i,2),v(i,3),v(j,1)-v(i,1),...
                    v(j,2)-v(i,2),v(j,3)-v(i,3),'color',c); %#ok<AGROW>
            end
            h.UserData = hq;
        end
    end

    methods(Access=protected)
        function obj = dotAssign(obj,indexOp,c)
            switch indexOp.Name
                case {'e0','eps1234'}
                    obj.m(1) = 0; warning('input ignored')
                case {'e1','eps423'}
                    obj.m(2) = 0; warning('input ignored')
                case {'e2','eps431'}
                    obj.m(3) = 0; warning('input ignored')
                case {'e3','eps412'}
                    obj.m(4) = 0; warning('input ignored')
                case {'e4','eps321'}
                    obj.m(5) = 0; warning('input ignored')
                case {'e23','eps41'}
                    obj.m(6) = 0; warning('input ignored')
                case {'e31','eps42'}
                    obj.m(7) = 0; warning('input ignored')
                case {'e12','eps43'}
                    obj.m(8) = 0; warning('input ignored')
                case {'e43','eps12'}
                    obj.m(9) = 0; warning('input ignored')
                case {'e42','eps31'}
                    obj.m(10) = 0; warning('input ignored')
                case {'e41','eps23'}
                    obj.m(11) = 0; warning('input ignored')
                case {'e321','eps4'}
                    obj.m(12) = c;
                case {'e412','eps3'}
                    obj.m(13) = c;
                case {'e431','eps2'}
                    obj.m(14) = c;
                case {'e423','eps1'}
                    obj.m(15) = c;
                case {'e1234','eps0'}
                    obj.m(16) = 0; warning('input ignored')
                otherwise
                    error('basis not recognized')
            end
        end
    end
end