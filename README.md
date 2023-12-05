# Rigid Geometric Algebra
A collection of Matlab classes implementing the Geometric Algebra G(3,0,1) for homogeneous transformations of 3D rigid bodies in 4D projective space, primarily based on Lengyal's convention.

# Motivation
There are many tremendous existing resources providing various implementations of Geometric Algebra.  One of the early ones, GABLE, was started in Matlab but moved on to C++ and Java. Another one that is still maintained is the Clifford Multivector Toolbox, which provides a really general set of capabilities.  There are many more.  What seems to be missing is a set of Matlab tools that implement Lengyal's conventions. Lengyal's approach highlights the concept of duality in Geometric Algebra in a manner somewhat differently than most others; some may find Lengyal's treatment to be more satisfying.  This set of classes is aimed for such users, who desire to work in a Matlab environment.

# Getting Started
There is a tutorial script and a demo script.  Both of these will be best utilized by saving the .m file as a Live Script, or using the Publish functionality of the editor to render a rich-text static version.  The tutorial is the place to start if Geometric Algebra, and/or Lengyal's approach to it are new to you.  If you already have some relevant background, the demo will show you how to do things you already understand using the objects and methods that this class provides.

# Applications
One application is included, in the script wahba.m.  This implements the pose solution of a generalized version of Wahba's Problem, using Lengyal's conventions and the null-space method of Perwass, as well as a two-step method inspired by a classic paper by Walker, Shao, and Volz that used dual quaternions.  The demo script illustrates this application at the end.

# References
If Geometric Algebra is new to you, there are some great resources online, and a number of textbooks. A few suggestions:

https://en.wikipedia.org/wiki/Geometric_algebra

https://projectivegeometricalgebra.org/

https://bivector.net/

https://enkimute.github.io/ganja.js/

Perwass, C., Geometric Algebra with Applications in Engineering, Springer, 2010, ISBN: 978-3540891680.

MacDonald, A., Geometric Algebra and Calculus, Vols. 1 & 2, CreateSpace Independent Publishing, 2021, ISBN: 978-1453854938, 978-1480132450.

Lengyal, E., Projective Geometric Algebra Illuminated (forthcoming as of Dec 2023), Terathon Software, 2023, ISBN: 979-8985358247.
