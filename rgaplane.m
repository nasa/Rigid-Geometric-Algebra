classdef rgaplane < rga
    %RGAPLANE RGA Plane
    % Subclass of RGA
    methods
        function obj = rgaplane(x,y,z,w)
            %RGAPLANE Construct RGA plane from its normal
            switch nargin
                case 0
                    obj.m = fliplr([0 randn(1,4) zeros(1,11)]);
                case 1
                    if isa(x,"rga")
                        t = trivector(x);
                        obj.m = t.m;
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
    end
end