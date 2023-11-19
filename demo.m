%% *Demonstrate capabilities of RGA class*

%% Introduction
% There are lots of methods implemented in RGA; only a few demo'd here.
% Click the documentation hyperlinks for more info.  For best enjoyment of
% this demo, try using the Matlab editor publish function, and/or save it
% and run it as a Live Script.
clear
clc
help rga

%% Creation of RGA objects
% There are lots of ways to create RGA objects.  Here are a few.  Note that
% $e_0$ is the scalar element, and $e_{1234}$ is the pseudoscalar.  Also
% note that some of the objects use $\epsilon$ rather than $e$ for basis
% elements.  More on that below.
P1 = rgapoint % create a random RGA point
L1 = rgaline % create a random RGA line
F1 = rgaplane % create random RGA plane
Q1 = rgamotor % create random RGA motor
M1 = rga % create random general RGA multivector
P2 = rgapoint([53.767, 183.39, -225.88, 1]) % create specified RGA point
F2 = rgaplane([0.52416, 1.3645, -1.2105, 0.070602]) % create RGA plane
try
    L2 = rgaline(randn(3,1),randn(3,1)) % not all inputs are possible
catch ME
    warning(ME.message)
    v = randn(3,1);
    m = randn(3,1);
    v = v - (v'*m)/(m'*m)*m;
    L2 = rgaline(v,m)
end

%% Using the anti property
% It is sometimes convenient to use an alternative set of dual basis
% vectors, or anti-basis elements.  Note that a slightly different symbol
% is used to indicate the dual basis.  The anti property is sticky or
% dominant in that if an anti multivector interacts with a non-anti
% multivector, the result will be anti.  This is a convention adopted for
% convenience within this implementation.
L2.anti = true
M2 = rga(randn(16,1),true) % set anti property via constructor (only for RGA, not subclasses!)

%% Using basis elements
% One may also create multivectors from RGA basis elements, in various
% ways:
[e0,e1,e2,e3,e4,e23,e31,e12,e43,e42,e41,e321,e412,e431,e423,e1234] = rga.bases;
M3 = e0 + e1 + e2 + e3 + 2*e23 + 2*e31 + 2*e12 + e321
M4 = 1*rga('e0') + 5*rga('e12') - 3*e431 % can create elements "on the fly"
M5 = -2*rga('eps2') + 13*rga('eps23') - 4*e1234 % can mix in anti; note "stickiness"

%% Unitizing
% It is often helpful to unitize RGA objects, which means something a
% little different for each type of object.  So it is not meaningful to
% unitize a generic RGA multivector object.
P3 = unitize(P1)
L3 = unitize(L1)
F3 = unitize(F1)
Q3 = unitize(Q1)

%% Symbolic RGA objects
% RGA objects support most symbolic operations
syms a0 a1 a2 a3 a23 a31 a12 a321 b0 b1 b2 b3 b23 b31 b12 b321 real
syms a4 a43 a42 a41 a412 a431 a423 a1234 b4 b43 b42 b41 b412 b431 b423 b1234 real
M6 = rga([a0 a1 a2 a3 a4 a23 a31 a12 a43 a42 a41 a321 a412 a431 a423 a1234])
M7 = rga([b0 b1 b2 b3 b4 b23 b31 b12 b43 b42 b41 b321 b412 b431 b423 b1234])
P4 = rgapoint(a1,a2,a3,a4)
P5 = rgapoint(b1,b2,b3,b4)

%% Using products
% Products include dot, antidot, wedge, antiwedge, wedgedot, and
% antiwedgedot.
B1 = dot(P1,P2) + wedge(P1,P2)
B2 = wedgedot(P1,P2) % B2 = B1 *only* for vectors
X1 = dot(L1,L2) + wedge(L1,L2)
X2 = wedgedot(L1,L2)
A1 = antidot(F1,F2) + antiwedge(F1,F2)
A2 = antiwedgedot(F1,F2) % A2 = A1 *only* for planes (anti-vectors)
Y1 = antidot(L1,L2) + antiwedge(L1,L2)
Y2 = antiwedgedot(L1,L2)
Z1 = wedge(P4,P5)
Z2 = antiwedgedot(P4,P5)

%% Using overloaded product operators
% The mpower method ('^') is overloaded to the exterior (wedge) product,
% and the mtimes method ('*') is overloaded to the geometric (wedgedot)
% product. The element-wise versions of these (power, '.^' and times, .*)
% are overloaded to the anti-wedge and anti-wedgedot products.
B3 = P1^P2
B4 = F1.^F2
B5 = P1*P2
B6 = F1.*F2

%% Using reverses and antireverses
% RGA motors perform screw operations using a sandwich antiproduct that
% involves a generalization of the conjugate called the antireverse.  The
% transpose operators (' and .') are overloaded to apply the reverse or
% antireverse, respectively
R1 = rev(Q3)
R2 = antirev(Q3)
N1 = antiwedgedot(Q3,antiwedgedot(L1,antirev(Q3)))
N2 = Q3.*L1.*Q3.'
[ph3,d3,v3,m3] = extract(Q3) % Extract angle, distance, direction & moment

%% Arrays of RGA objects
% Can use ordinary Matlab array indexing for objects of same type, but for
% heterogenous combinations, must use cell arrays.  At present, one cannot
% perform linear algebra on these arrays... this would be an interesting
% thing to pursue eventually.  There is a "Clifford Multivector Toolbox,"
% which partially inspired the RGA class, that does permit such operations,
% and it also supports arbitrary GA signatures.
M = [M1 M2; M3 M4]
M(1)
H = {P1 L1 F1 Q1}
H{4}

%% Now for some plotting and geometry
clf
hold off
hold on

%% Define a line, and points lying on it
px = 1; py = 0; pz = 0; p = [px;py;pz];
qx = 0; qy = 1; qz = 0; q = [qx;qy;qz];
v = (q - p)' % direction
m = cross(p,q)' % moment

%% Plot with Matlab native drawing
plot3(px,py,pz,'*'), text(px,py,pz,'p')
plot3(qx,qy,qz,'*'), text(qx,qy,qz,'q')
line([px;qx],[py;qy],[pz;qz])
axis equal
view(37,30)
grid on
quiver3(-1,0,0,1,0,0,2,'r') % x axis
quiver3(0,-1,0,0,1,0,2,'g') % y axis
quiver3(0,0,-1,0,0,1,2,'b') % z axis

%% Next create & plot GA versions
p = rgapoint(px,py,pz,1);
q = rgapoint(qx,qy,qz,1);
L = p^q % rga object from p to q
L = rgaline(L); % subclass into line object
v = direction(L); v.anti = true
m = moment(L)
plot(p);
plot(q);
plot(L);

%% Find perpendicular distance of line from origin
o = rgapoint(0,0,0);
doL = dist(o,L);

%% Create bounding box around origin big enough to contain the line
% Choose here to have all plane normals pointing "in" toward origin.
% (From eq. of a plane, distance from origin is oppositely directed from
% the plane's normal vector.)
fpx = rgaplane(-1,0,0,doL); hfpx = plot(fpx); hfpx.FaceColor = 'r';
fmx = rgaplane(1,0,0,doL);  hfmx = plot(fmx); hfmx.FaceColor = 'm';
fpy = rgaplane(0,-1,0,doL);  hfpy = plot(fpy); hfpy.FaceColor = 'g';
fmy = rgaplane(0,1,0,doL); hfmy = plot(fmy); hfmy.FaceColor = 'y';
fpz = rgaplane(0,0,-1,doL); hfpz = plot(fpz); hfpz.FaceColor = 'b';
fmz = rgaplane(0,0,1,doL); hfmz = plot(fmz); hfmz.FaceColor = 'c';

%% Search all 6 faces for intersections with line (at w=1)
ipx = unitize(rgapoint(commutate(L,fpx,'+v')));
imx = unitize(rgapoint(commutate(L,fmx,'+v')));
ipy = unitize(rgapoint(commutate(L,fpy,'+v')));
imy = unitize(rgapoint(commutate(L,fmy,'+v')));
ipz = unitize(rgapoint(commutate(L,fpz,'+v')));
imz = unitize(rgapoint(commutate(L,fmz,'+v')));

%% Plot line using intersection points
if ipx.m(5)==1, plot3(ipx.m(2),ipx.m(3),ipx.m(4),'r+'), end
if imx.m(5)==1, plot3(imx.m(2),imx.m(3),imx.m(4),'ro'), end
if ipy.m(5)==1, plot3(ipy.m(2),ipy.m(3),ipy.m(4),'k+'), end
if imy.m(5)==1, plot3(imy.m(2),imy.m(3),imy.m(4),'ko'), end
if ipz.m(5)==1, plot3(ipz.m(2),ipz.m(3),ipz.m(4),'b+'), end
if imz.m(5)==1, plot3(imz.m(2),imz.m(3),imz.m(4),'bo'), end

%% Solving a Generalized Wahba Problem in RGA
% Given a set of corresponding RGA objects, find the motor that relates
% them. (Here, directions are unit lines through the origin)
wahba_test % default test case with each type noise-free observation
wahba_test(points=1,planes=2,directions=0,lines=0,motors=0) % still noise-free
wahba_test(points=1,planes=2,directions=0,lines=0,motors=0,...
    point_sigma=0.1,plane_sigma=0.1) % add some noise w/specified sigmas
wahba_test(points=0,planes=0,directions=2,lines=0,motors=0,fullpose=false) % original Wahba
wahba_test(points=0,planes=0,directions=2,lines=0,motors=0,fullpose=false,...
    direction_sigma=0.1) % original Wahba w/noise

%% Display basis element product tables
% Each cell in the interior of a table shows the application of the
% specified product to the pair of elements specified in the corresponding
% 1st row and 1st column, in that order.  Zero results are blank.
rga.producttab('wedge')
rga.producttab('wedgedot')
rga.producttab('antiwedge')
rga.producttab('antiwedgedot')
rga.producttab('rint')
rga.producttab('lint')

%% Using linear algebra representation
% It is sometimes convenient work with a linear algebra representation of
% RGA objects and operations. Just as complex numbers and quaternions can
% be represented as 2- or 4-element linear algebra column or row vectors,
% respectively, RGA objects can be cast as 16-element linear algebra
% vectors.  These can be extracted from the "m" property of an RGA object
% contains the coefficients of the basis elements.  And just as products
% involving complex numbers, quaternions, and cross-products of ordinary
% vectors can be represented as a matrix-vector product, so can RGA
% products.  Keep in mind that the linear algebra vectors will always
% contain the coefficients in the same order as the regular basis vectors,
% not the reversed order of the anti-basis elements.  Many other RGA
% operations can be represented as linear algebra matrix-vector products,
% e.g. reverses and anti-reverses.
m1 = M1.m % Extract multivect coefficients into linear algebra vector
m6 = M6.m
productmat(M1,"wedge") % 16x16 matrix that can be used to multiply M1 by something
(productmat(M1,"wedge")*M3.m(:))' % = linear algebra version of M1^M3
wedge(M1,M3)
U = eye(16); U(2:11,2:11)=-eye(10); % anti-reverse matrix
(U*M1.m(:))'
antirev(M1)