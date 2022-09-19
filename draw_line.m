%% Line Drawing
clf
% First create lines with zero moment representing the image plane axes
I = rgaline([1 0 0],[0 0 0]);
J = rgaline([0 1 0],[0 0 0]);
% Next create a line, and then lines K_I and K_J that contain a 
% point on each of the lines L and I, and L and J respectively.
%L = rgaline; 
p = rgapoint(1,0,1); q = rgapoint(0,1,1);
L = p^q; L = rgaline(L.m);
%L = unitize(L); % not really needed
v = direction(L) % need to downcast as a vector...?
K_I = commutate(L,I,"-v");
K_J = commutate(L,J,"-v");
% Finally extract the points on the I and J axes, respectively, & plot.
poni = antiwedge(wedge(K_I,v),I); % shoudl this be e1 not e41?
ponj = antiwedge(wedge(K_J,v),J); % should this be e2 not e42?
line([nonzeros(poni.m);0],[0;nonzeros(ponj.m)])