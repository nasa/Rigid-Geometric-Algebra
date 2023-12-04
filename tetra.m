classdef tetra
    %TETRA Tetrahedron
    %   Detailed explanation goes here

    properties
        vertices
        edges
        faces
    end

    methods
        function obj = tetra(specification,originOffset)
            %TETRA Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                specification = 2; % default to sideLength
                originOffset (1,3) = [0 0 0];
            end
            switch class(specification)
                case 'double' % sideLength: must be a single scalar
                    if length(specification) == 1
                        sideLength = specification;
                    else
                        error('sideLength must be a single scalar value')
                    end
                    v(4,:) = [0,0,1];
                    v(3,:) = [-sqrt(2/9),-sqrt(2/3),-1/3];
                    v(2,:) = [-sqrt(2/9),sqrt(2/3),-1/3];
                    v(1,:) = [sqrt(8/9),0,-1/3];
                    v = sideLength/2*v + originOffset(:)';
                    for i = 4:-1:1
                        P(i) = rgapoint(v(i,1),v(i,2),v(i,3));
                    end
                case 'rgapoint' % must have 4 elements
                    if length(specification) == 4
                        P = specification;
                    else
                        error('Must specify 4 vertices')
                    end
                case 'rgaline' % must have 6 elements
                    if length(specification) == 6
                        E = specification;
                    else
                        error('Must specify 6 edges')
                    end
                case 'rgaplane' % must have 4 elements
                    if length(specification) == 4
                        F = specification;
                    else
                        error('Must specify 4 faces')
                    end
                otherwise
            end
            C42 = nchoosek(1:4,2);
            C43 = nchoosek(1:4,3);
            if exist('P','var')
                for i = 6:-1:1
                    E(i) = unitize(rgaline(P(C42(i,1))^P(C42(i,2))));
                end
                for i = 4:-1:1
                    F(i) = unitize(rgaplane(P(C43(i,1))^P(C43(i,2))^P(C43(i,3))));
                end
            elseif exist('E','var')
                C62 = nchoosek(1:6,2); 
                C62 = C62([1 2 6 13],:); % don't need them all
                for i = 4:-1:1
                    F(i) = tetra.planefrom2lines(E(C62(i,1)),E(C62(i,2)));
                end
                for i = 4:-1:1
                    P(i) = unitize(rgapoint(F(C43(i,1)).^F(C43(i,2)).^F(C43(i,3))));
                    P(i).anti = false;
                end
            elseif exist('F','var')
                for i = 4:-1:1
                    P(i) = unitize(rgapoint(F(C43(i,1)).^F(C43(i,2)).^F(C43(i,3))));
                    P(i).anti = false;
                end
                for i = 6:-1:1
                    E(i) = unitize(rgaline(P(C42(i,1))^P(C42(i,2))));
                end
            else
                error('Something went wrong')
            end
            obj.vertices = P;
            obj.edges = E;
            obj.faces = F;
            if nargin == 0
                plot(obj);
            end
        end

        function h = plot(obj)
            %PLOT Plot the tetrahedron object
            P = obj.vertices;
            h(4) = patch('faces',[1 2 3],'vertices',[...
                P(4).e1,P(4).e2,P(4).e3;...
                P(1).e1,P(1).e2,P(1).e3;...
                P(2).e1,P(2).e2,P(2).e3],'facecolor','c');
            h(3) = patch('faces',[1 2 3],'vertices',[...
                P(4).e1,P(4).e2,P(4).e3;...
                P(2).e1,P(2).e2,P(2).e3;...
                P(3).e1,P(3).e2,P(3).e3],'facecolor','m');
            h(2) = patch('faces',[1 2 3],'vertices',[...
                P(4).e1,P(4).e2,P(4).e3;...
                P(3).e1,P(3).e2,P(3).e3;...
                P(1).e1,P(1).e2,P(1).e3],'facecolor','y');
            h(1) = patch('faces',[1 2 3],'vertices',[...
                P(1).e1,P(1).e2,P(1).e3;...
                P(2).e1,P(2).e2,P(2).e3;...
                P(3).e1,P(3).e2,P(3).e3],'facecolor','g');
        end
    end

    methods (Hidden, Static)
        function F = planefrom2lines(L1,L2)
            % First project origin onto one line to find a point on it,
            % then wedge this with the other line to find a plane
            o = rgapoint(0,0,0);
            q = rgapoint(L2.^(o^rcomp(weight(L2))));
            F = unitize(rgaplane(L1^q));
        end

    end
end