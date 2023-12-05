%% Tutorial

clc
clear
%% Motivation
% It can be useful to extend the real numbers.  A common example is to add an 
% imaginary number, $i$ that squares to $-1$, to get the complex numbers.  Adding 
% two more of these, $j$ and $k$, along with a multiplication rule such as $ijk=-1$ 
% (Hamilton) or $kji = -1$(Shuster), gets the quaternions. Another example is 
% to add a number, $\epsilon$, that squares to $0$, to get the dual numbers, and 
% this can be generalized for example to the dual quaternions, which will involve 
% products of $\epsilon$ and $i, j,$ and $k$.
% 
% In any of these cases, one ends up with a kind of compound scalar quantity, 
% such as $a + ib$ (complex number), $a+\epsilon b$ (dual number), or $a + ib 
% + jc + kd$ (quaternion). It is often convenient to view these amalgams as vectors, 
% e.g. complex numbers as two-dimensional vectors, quaternions as four-dimenional 
% vectors, etc.  In this view, symbols such as $i$, $j$, $k$, and $\epsilon$ represent 
% basis vectors, and some convention is sometimes adopted to label a basis vector 
% for the "ordinary" real line, such as *1*. While often useful, this perspective 
% can obscure other important properties, such as the ability to support division, 
% which vectors do not have. Another important distinction that such a view misses 
% is that of covariance vs. contravariance, which has to do with whether or not 
% the coefficients of the basis vectors scale proportionally or inversely with 
% changes in the basis vectors.  There can also be confusion between associating 
% a vector with a point, which is fixed relative to an origin, and an arrow that 
% can move around as needed, for example to perform vector addition.
% 
% Geometric algebra is both an approach to generalizing algebra over the reals 
% to algebras over arbitrary combinations of real, imaginary, and dual numbers, 
% and an approach to definitively associating such numbers with geometrical objects, 
% such as points, lines, planes, and volumes, and their transformations.  The 
% former has its roots in the work of Clifford, and the latter in the works of 
% Hamilton and Grassman.  Geometric algebra may also be viewed as a form of compact 
% notation for multilinear operations that are often performed using tensor algebras.
%% Conventions
% As so often occurs, e.g. with the quaternions, a number of differing conventions 
% have evolved that concern representations and interpretations of various Geometric 
% Algebras.  This tutorial and the Matlab class it uses rely primarily on Lengyal's 
% conventions:
% 
% <https://projectivegeometricalgebra.org/ https://projectivegeometricalgebra.org/>
% 
% <https://rigidgeometricalgebra.org/ https://rigidgeometricalgebra.org/>
% 
% Let there be $p$ elements of the algebra that square to $+1$, $q$ elements 
% that square to $-1$, and $r$ elements that square to $0$.  Denote this algebra 
% as $G(p,q,r)$.  The set $(p,q,r)$ defines the "signature" of the algebra.  Let 
% these elements behave as a vector space (addition, scalar multiplication, etc.). 
% There are $n=p+q+r$ of these vectors, along with one scalar element.  It is 
% common to label the $n$ elements sequentially as $e_i$, where the first $p$ 
% of the $e_i$ are those that square to $+1$, etc.  It will sometimes also be 
% convenient to assign $e_0$ as a label for a sort of basis associated with the 
% scalar element, although this is a bit of an abuse of the convention. 
%% Wedge Product
% Grasssman's wedge (exterior) product is fundamental in all Geometric Algebras, 
% and it is from this that the basics of geometry come in.   The wedge product 
% of an element with itself is always zero (unless it is the scalar element).  
% The wedge product of $e_i$ and $e_j$, where $i\neq j$, is _not_ another vector, 
% but rather a new object called a bivector: $e_{ij} = e_i \wedge e_j$. The geometric 
% picture is that $e_{ij}$ is the plane spanned by $e_i$ and $e_j$.  And $e_{ijk} 
% = e_i \wedge e_j \wedge e_k$, a trivector, is the volume spanned by $e_i$, $e_j$, 
% and $e_k$.  This contruction can continue into higher geometric dimensions until 
% there are no more basis elements to consider.  When $n$ basis elements have 
% been wedged, the fullness of the space has been spanned. The "final" $n$-dimensional 
% hypervolume turns out to behave a lot like a scalar, so it is called the pseudoscalar.  
% The process of building up elements via wedge products will end up producing 
% $2^n$ basis elements, with 1 scalar, $n$ vectors, $n \choose 2$ bivectors, $n 
% \choose 3$ trivectors, etc., and finally the 1 pseudoscalar. A general linear 
% combination of basis elements is called a multivector.
% 
% The algebra treated by the RGA class is $G(3,0,1)$, which means there are 
% three elements that square to $+1$, zero elements that square to $-1$, and 1 
% element that squares to $0$.  Consequently, there are 16 total basis elements, 
% inclusive of the scalar $e_0$ and the pseudoscalar $e_{1234}$. The RGA class 
% has a static method to return these basis elements as variables, which can be 
% useful:

[e0,e1,e2,e3,e4,e23,e31,e12,e43,e42,e41,e321,e412,e431,e423,e1234] = rga.bases;
%% 
% The RGA class overloads the |mpower| method to implement the $\wedge$ symbol 
% as the wedge product.  Checking a few operations produces the expected results, 
% including for the operations on scalars:

e0^e0
%% 
% This is just scalar multiplication.  Similarly,

(6*e0)*(-3*e0)
%% 
% produces the expected result.  A scalar wedged with a vector just scales the 
% vector:

e0^e1
(.25*e0)^e2
%% 
% It's generally not necessary to include e0 though:

1 + e1
.25*e2
%% 
% A vector wedged with itself produces zero:

e1^e1
%% 
% A vector wedged with another vector produces a bivector:

e1^e2
%% 
% Note that bivectors anti-commute:

e2^e1
%% 
% Scalars also scale bivectors:

(1/3*e0)^e12
%% 
% A vector wedged with a trivector with which it does not share any common elements 
% produces the pseudoscalar:

e1^e423
%% 
% By contrast, here, $e_1$ is repeated, since $e_{321} = e_3 \wedge e_2 \wedge 
% e_1$:

e1^e321
%% 
% so the result is zero.
% 
% Note that the pseudoscalar wedged to itself is zero:

e1234^e1234
%% Geometric Product
% If the wedge product of a vector element with itself is zero, what does it 
% mean to say that vector elements "square to" $+1$, $-1$, or $0$?  Clearly the 
% expression "squares to" cannot be referencing the wedge product!  Stitching 
% these two threads together gives the geometric product.  Whenever there occurs 
% a product of a vector element with itself (even when it is embedded as part 
% of a bivector, trivector, etc.), the signature of the algebra gives the result.  
% Whenever there occurs a product of two different vector elements, the wedge 
% product gives the result.
% 
% The RGA class overloads the mtimes method to implement the $*$ symbol as the 
% geometric product.  A few operations can illustrate how the geometric product 
% works.
% 
% Scalar multiplication is the same as for the wedge product:

e0*e0
(13*e0)*e12
%% 
% Multiplying vectors with themselves gives the results expected from the signature:

e1*e1, e2*e2, e3*e3, e4*e4
%% 
% Compare these two operations to see what is going on:

e1*(e1^e2^e3), e1*e321
%% Duals
% Lengyal's approach highlights the concept of duality in Geometric Algebra 
% in a manner somewhat differently than most others.  Lengyal emphasizes that 
% for each basis element, one may view it in terms of the basis vectors that are 
% present, _or the basis vectors that are absent._ For example, in this view, 
% $e_1$ and $e_{423}$ are dual to each other, because in the triple wedge product 
% $e_{4}\wedge e_{2}\wedge e_{3}$ the basis vector $e_1$ is absent.  This idea 
% extends to the geometric objects in the algebra, i.e. a point and a plane are 
% dual to each other in this sense.  In Lengyal's approach, which he points out 
% was Grassman's as well, there exist duals to the products discussed above.  
% The dual of the wedge product is the anti-wedge product, indicated by $\vee$. 
% Rather than building up the basis elements that are present, the anti-wedge 
% builds up the basis vectors that are absent.  For example $e_{412} \vee e_{23}$ 
% combines the element missing from $e_{412}$, which is $e_3$, and the elements 
% missing from $e_{23}$, which are $e_4$ and $e_1$.  The result is the dual of 
% $e_3 \wedge e_4 \wedge e_1 = - e_4 \wedge e_3 \wedge e_1$, which is $-e_2$.  
% Similarly, the dual of the geometric product is the geometric anti-product.  
% Lengyal introduced the "wedge-dot" and "vee-dot" (or "anti-wedge-dot") notation, 
% $\wedge\!\!\!.$ and $\vee\!\!\!\!\dot$, respectively, to distinguish the geometric 
% product from the geometric anti-product. Since Matlab provides "dot" versions 
% of the "^" and "*" operators, RGA uses ".^" for the antiwedge and ".*" for the 
% anti-wedge-dot products.  A mnemonic to keep this straight is the "dot" preceding 
% the multiplication symbol stands for "down", which is the way the $\vee$ and 
% the $\vee\!\!\!\!\dot$ point.  Examples:

e412.^e23
%% 
% which gives the result just discussed.  Contrast the anti-wedge and anti-wedge-dot 
% products:

e412.^e12
%% 
% Here $e_3$ is missing on both sides of the anti-wedge, so the result is zero.

e412.*e12
%% 
% Here, $e_3$ is missing on the left, and $e_3$ and $e_4$ are missing on the 
% right. The geometric product of $e_3$ and $e_3\wedge e_4$ is $e_4$, and the 
% dual of $e_4$ is $e_{321}$.
% 
% It is sometimes easier to work with a dual set of bases elements, rather than 
% trying to remember which elements are missing.  All objects or the RGA class 
% have a property |anti| that one can set to true to express a multivector in 
% terms of the dual basis elements.  The symbol $\varepsilon$ indicates a dual 
% basis element. For example, here is a randomly-generated multivector, expressed 
% on both the standard and the dual basis:

X = rga
X.anti = true
%% 
% Similar with the basis elements, there is a static method to return the anti-basis 
% elements as variables:

[eps0,eps1,eps2,eps3,eps4,eps23,eps31,eps12,eps43,eps42,eps41,...
    eps321,eps412,eps431,eps423,eps1234] = rga.antibases
%% Geometric Operations
% Just as the wedge product builds up higher-dimensional objects by joining 
% lower-dimensional objects together, the anti-wedge product finds intersections 
% where higher-dimensional objects meet. To see some of this, it's helpful to 
% make use of subclasses the RGA class provides for points, lines, and planes.  
% First create two points.

p = rgapoint(1,1/2,1/3)
q = rgapoint(-1/3,-1/2,-1)
%% 
% This creates two points.  Notice they have an $e_4$ component as well, which 
% defaults to $+1$. In $G(3,0,1)$, it is customary to treat $e4$ as a basis vector 
% associated with $\infty$, so that this algebra becomes associated with projective 
% geometry.  In this view, the coordinate value $+1$ along the $e_4$ direction 
% is the hyperplane associated with physical space.  Each subclass corresponding 
% to a geometric object has a plot method.  Before making any plots, first create 
% a set of coordinate axes to give better context; the RGA class provides a conventient 
% static method for this.

h1 = figure(1);
rga.plotaxes;
%% 
% Now plot the points:

plot(p);
plot(q);
%% 
% Joining these two points with the wedge product creates a line:

L = p^q
%% 
% Although RGA doesn't "know" this is a line, so it remains a generic RGA object. 
% To gain access to methods that are unique to line objects, it's helpful to subclass 
% it into a line:

L = rgaline(L)
%% 
% The six bivector components of the line correspond to its Plucker coordinates, 
% which give the line's direction and moment.  The subclass provides methods to 
% extract this info:

v = direction(L)
m = moment(L)
%% 
% The |plot| method for the |rgaline| subclass allows for passing a point as 
% a second argument which will locate the end of the quiver it creates.

clf
rga.plotaxes;
plot(p); plot(q);
plot(L,p);
%% 
% Create another point:

r = rgapoint(-1,1,0)
%% 
% To create a plane, either join three unique points, or join a point and a 
% line that does not contain it:

f = p^q^r
g = L^r
%% 
% Reversing the order of the wedge operations gives the same result but with 
% all the signs changed.  This corresponds to a plane with the opposite "orientation 
% sense."

r^q^p
%% 
% Subclass the result into a plane, and notice that it defaults to the anti-basis 
% representation:

f = rgaplane(f)
%% 
% This makes clear that in the anti-basis, there is an antivector that corresponds 
% to what would be the plane's normal vector in a traditional setting. Plotting 
% it, notice that quivers around the boundary indicate the "orientation sense" 
% of the plane.

clf
rga.plotaxes;
plot(p); plot(q); plot(r); plot(L,p);
plot(f);
lighting none
%% 
% To see how intersections work, create a plane that corresponds to the $x,y$ 
% plane:

h = rgaplane(0,0,1,0)
%% 
% The anti-wedge product of f and |h| gives their intersection, which is the 
% line where they meet:

M = f.^h
%% 
% Notice that the |anti| property is sticky, even when subclassing |M|:

M = rgaline(M)
%% 
% To revert to the regular basis, set the anti property to false:

M.anti = false
direction(M)
moment(M)
%% 
% Plot everything to confirm the result makes sense:

clf
rga.plotaxes;
plot(f);
plot(h);
plot(M);
lighting none
%% 
% Note that in some other conventions, the $\wedge$ product is associated with 
% intersections ("meets") and the $\vee$ is associated with joins. This is a clue 
% that such conventions are making use of a dual basis.  For example, if one adopts 
% the convention that |f| and |h| are points, as they appear to be in the dual 
% basis representation, then it make sense to associate $\vee$ with the join operation, 
% since $f\vee h$ produces a line.
% 
% Note also that the anti-wedge $\vee$ of two lines does not actually give the 
% point where they intersect, because for example $(P_1 \wedge P_2) \vee (P_1 
% \wedge P_3)$ will contain $P_1 \vee P_1$ which is zero.  But if the two lines 
% do not intersect, the resulting scalar quantity gives information about the 
% way the two lines screw around one another.
%% Geometric Transformations
% While there are many other geometric manipulations that one may accomplish 
% with the wedge and anti-wedge products, it is the geometric product that provides 
% the mechanism for accomplishing transformations, such as reflections, rotations, 
% translations, etc. As with the quaternions and other similar entities, this 
% requires "sandwiching" the geometric object to be transformed, using a geometric 
% product or anti-product on each side, with another geometric object, for example
% 
% $$B = Q$$  $$\vee\!\!\!\!\dot$$ $$A$$ $$\vee\!\!\!\!\dot$$ $$\underset{\widetilde}Q$$
% 
% Notice another bit of notation: placing the tilde above a geometric algebra 
% object indicates the _reverse_ operation, and placing it under the object indicates 
% the anti-reverse.  In geometric algebra, the reverse is analagous to a conjugate 
% or tranpose operation.  What happens with the reverse is that any elements that 
% have been wedge together get their order reversed, so that for example $e_{12}$ 
% becomes $e_{21} = -e_{12}$.  The anti-reverse works the same way but using the 
% anti-wedge.
% 
% The ability of the geometric (anti-)product to perform transformations is 
% easiest to see by starting with reflections of a tetrahedron.  Create a starting 
% tetrahedron with the |tetra| class, which is affiliated with the RGA class.  
% For a starting tetrahedron, specify a sidelength of 1 and an offest from the 
% origin along the $x$ and $z$ axes:

clf
rga.plotaxes;
T1 = tetra(1,[2;0;1])
plot(T1);
%% 
% Extract the vertices of the tetrahedron:

U = T1.vertices;
%% 
% To reflect the vertices about the $x,y$ plane, use an anti-geometric product 
% sandwich.  First create a plane corresponding to the $x,y$ plane:

R1 = rgaplane(0,0,1,0); 
%% 
% Next sandwich the vertices using the anti-wedgedot and the anti-reverse.  
% The RGA class overloads the |ctranpose| and |tranpose| operators for the reverse 
% and anti-reverse, respectively.

for i = 4:-1:1
    V(i) = rgapoint(R1.*U(i).*R1.');
end
%% 
% Pass these new vertices into the tetra constructor to create a new tetrahedron:

T2 = tetra(V);
clf
rga.plotaxes;
plot(T1);
plot(R1);
plot(T2);
lighting none
%% 
% Now notice what happens with a second reflection about a second plane:

R2 = rgaplane(sqrt(2)/2,0,sqrt(2)/2,0);
for i = 4:-1:1
    W(i) = rgapoint(R2.*V(i).*R2.');
end
T3 = tetra(W);
clf
rga.plotaxes;
plot(T1);
plot(R1);
plot(T2);
plot(R2);
plot(T3);
lighting none
%% 
% Notice that the tetrahedron has been rotated through 90 degrees, which is 
% twice the angle separating the two planes, and that this rotation occurred around 
% the $y$ axis, which is the line where the two planes intersect.  This is easiet 
% to see from the green face, which started parallel to the $x,y$ plane, and ended 
% up parallel to the $y,z$ plane.  Combining the two sandwich products into one,
% 
% $$W_i = (R_2$$ $$\vee\!\!\!\!\dot$$ $$R_1)$$ $$\vee\!\!\!\!\dot$$ $$U_i$$ 
% $$\vee\!\!\!\!\dot$$ $$\underset{\widetilde}{(R_2 \vee\!\!\!\!\dot{}\ R_1})$$
% 
% Letting $Q = R_2 \vee\!\!\!\!\dot{}\ R_1$,

Q = R2.*R1; Q.anti = false
%% 
% it is clear that $Q$ looks like a quaternion, with a (pseudo-) scalar part 
% and a (bi-) vector part.  It's typical to think of a quaternion as rotating 
% a rigid body about its center of mass, which may seem different from the diagram 
% here.  If so, just imagine the tetrahedron is on one end of a pendulum that 
% swings about the origin.
% 
% By moving $R_2$ away from the origin, the line of intersection with $R_1$ 
% will no longer contain the origin, but all of the above still works, so already 
% there is something here that is more general than a quaternion:

R2.eps4 = 1
Q = R2.*R1; Q.anti = false
%% 
% Notice how the bivector part of $Q$ now contains both the direction and the 
% moment of the line of rotation.

for i = 4:-1:1
    W(i) = rgapoint(Q.*U(i).*Q.');
end
T4 = tetra(W);
clf
rga.plotaxes;
plot(T1);
plot(R1);
plot(T2);
plot(R2);
plot(T4);
lighting none
%% 
% 
% 
% Letting the line of intersection of the planes move to infinity results in 
% parallel planes:

R3 = rgaplane(0,0,1,2)
clf
rga.plotaxes;
plot(T1);
plot(R1);
plot(R3);
lighting none
Q = R3.*R1; Q.anti = false
%% 
% Now the bivector part of $Q$ contains only the moment, since a line at infinity 
% lies in all directions perpendicular to the moment.  Performing the double reflection,

for i = 4:-1:1
    W(i) = rgapoint(Q.*U(i).*Q.');
end
T5 = tetra(W);
clf
rga.plotaxes;
plot(T1);
plot(R1);
plot(T2);
plot(R3);
plot(T5);
lighting none
%% 
% One sees the double reflection results in a translation of twice the distance 
% between the planes.


%% 
%