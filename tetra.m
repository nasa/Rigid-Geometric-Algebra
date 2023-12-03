classdef tetra
    %TETRA Tetrahedron
    %   Detailed explanation goes here

    properties
        vertices
        edges
        faces
    end

    methods
        function obj = tetra(sideLength,originOffset)
            %TETRA Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                sideLength (1,1) = 2;
                originOffset (1,3) = [0 0 0];
            end
            v(4,:) = [0,0,1];
            v(3,:) = [-sqrt(2/9),-sqrt(2/3),-1/3];
            v(2,:) = [-sqrt(2/9),sqrt(2/3),-1/3];
            v(1,:) = [sqrt(8/9),0,-1/3];
            v = sideLength/2*v + originOffset(:)';
            for i = 4:-1:1
                P{i} = rgapoint(v(i,1),v(i,2),v(i,3));
            end
            Ce = nchoosek(1:4,2);
            for i = 6:-1:1
                E{i} = rgaline(P{Ce(i,1)}^P{Ce(i,2)});
            end
            Cf = nchoosek(1:4,3);
            for i = 4:-1:1
                F{i} = rgaplane(P{Cf(i,1)}^P{Cf(i,2)}^P{Cf(i,3)});
            end
            obj.vertices = P;
            obj.edges = E;
            obj.faces = F;
            if nargin == 0
                plot(obj)
            end
        end

        function h = plot(obj)
            %PLOT Plot the tetrahedron object
            P = obj.vertices;
            for i = 1:4
                plot(P{i});
            end
            patch('faces',[1 2 3],'vertices',[...
                P{4}.e1,P{4}.e2,P{4}.e3;...
                P{1}.e1,P{1}.e2,P{1}.e3;...
                P{2}.e1,P{2}.e2,P{2}.e3],'facecolor','c');
            patch('faces',[1 2 3],'vertices',[...
                P{4}.e1,P{4}.e2,P{4}.e3;...
                P{2}.e1,P{2}.e2,P{2}.e3;...
                P{3}.e1,P{3}.e2,P{3}.e3],'facecolor','m');
            patch('faces',[1 2 3],'vertices',[...
                P{4}.e1,P{4}.e2,P{4}.e3;...
                P{3}.e1,P{3}.e2,P{3}.e3;...
                P{1}.e1,P{1}.e2,P{1}.e3],'facecolor','y');
            patch('faces',[1 2 3],'vertices',[...
                P{1}.e1,P{1}.e2,P{1}.e3;...
                P{2}.e1,P{2}.e2,P{2}.e3;...
                P{3}.e1,P{3}.e2,P{3}.e3],'facecolor','g');
        end
    end
end