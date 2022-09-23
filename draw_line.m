%% Line Drawing
clf
hold off
hold on

%% Define the line, and points lying on it, "traditionally"
px = 1; py = 0; pz = 1; p = [px;py;pz];
qx = 0; qy = 1; qz = 1; q = [qx;qy;qz];
v = (q - p)'
m = cross(p,q)'

%% Start with Matlab native drawing
plot3(px,py,pz,'*'), text(px,py,pz,'p')
plot3(qx,qy,qz,'*'), text(qx,qy,qz,'q')
line([px;qx],[py;qy],[pz;qz])
axis equal
view(37,30)
grid on
quiver3(-1,0,0,1,0,0,2,'k')
quiver3(0,-1,0,0,1,0,2,'k')
quiver3(0,0,-1,0,0,1,2,'k')

%% Next create GA versions
p = rgapoint(px,py,pz,1);
q = rgapoint(qx,qy,qz,1);
L = p^q
L = rgaline(L.m);
v = direction(L); v.anti = true
m = moment(L)


% % Create lines with zero moment representing the image plane axes
% I = rgaline([1 0 0],[0 0 0]);
% J = rgaline([0 1 0],[0 0 0]);
% % Next create a line, and then lines K_I and K_J that contain a 
% % point on each of the lines L and I, and L and J respectively.
% %L = rgaline; 
%p = rgapoint(1,0,1); q = rgapoint(0,1,1);
%L = p^q; L = rgaline(L.m);
% %L = unitize(L); % not really needed
% v = direction(L) % need to downcast as a vector...?
% K_I = commutate(L,I,"-v");
% K_J = commutate(L,J,"-v");
% % Finally extract the points on the I and J axes, respectively, & plot.
% poni = antiwedge(wedge(K_I,v),I); % shoudl this be e1 not e41?
% ponj = antiwedge(wedge(K_J,v),J); % should this be e2 not e42?
% line([nonzeros(poni.m);0],[0;nonzeros(ponj.m)])

% 1. Find perpendicular distance of line from origin
o = rgapoint(0,0,0);
doL = dist(o,L);
% 2. Create bounding box around origin big enough to contain the line
a = rgapoint(1,1,1,1); b = rgapoint(-1,1,1,1); c = rgapoint(-1,-1,1,1); d = rgapoint(1,-1,1,1);
e = rgapoint(1,-1,-1,1); f = rgapoint(1,1,-1,1); g = rgapoint(-1,1,-1,1); h = rgapoint(-1,-1,-1,1);
fpx = rgaplane(doL,0,0,1); fmx = rgaplane(-doL,0,0,1);
fpy = rgaplane(0,doL,0,1); fmy = rgaplane(0,-doL,0,1);
fpz = rgaplane(0,0,doL,1); fmz = rgaplane(0,0,-doL,1);
% 3. Search all 6 faces for intersections with line

% 4. Plot line using intersection points

% For planes: Project origin and at least two other points onto plane, then
% create patch from those points.
