# Rigid Geometric Algebra
A collection of Matlab classes implementing the Geometric Algebra G(3,0,1) for homogeneous transformations of 3D rigid bodies in 4D projective space, primarily based on Lengyel's convention.

Copyright © 2023 United States Government as represented by the Administrator of the National Aeronautics and Space Administration. No copyright is claimed in the United States under Title 17, U.S. Code. All Other Rights Reserved.

Released under the NASA Open Source Agreement Version 3.1; refer to the included file NOSA GSC-19268-1.pdf for details.

# Motivation
There are many tremendous existing resources providing various implementations of Geometric Algebra.  One of the early ones, GABLE, was started in Matlab but moved on to C++ and Java. Another one that is still maintained is the Clifford Multivector Toolbox, which provides a really general set of capabilities.  There are many more.  What seems to be missing is a set of Matlab tools that implement Lengyel's conventions. Lengyel's approach highlights the concept of duality in Geometric Algebra in a manner somewhat differently than most others; some may find Lengyal's treatment to be more satisfying.  This set of classes is aimed for such users, who desire to work in a Matlab environment.

# Getting Started
There is a tutorial script and a demo script.  Both of these will be best utilized by saving the .m file as a Live Script, or using the Publish functionality of the editor to render a rich-text static version.  The tutorial is the place to start if Geometric Algebra, and/or Lengyel's approach to it are new to you.  If you already have some relevant background, the demo will show you how to do things you already understand using the objects and methods that this class provides.

# Applications
One application is included, in the script wahba.m.  This implements the pose solution of a generalized version of Wahba's Problem, using Lengyel's conventions and the null-space method of Perwass, as well as a two-step method inspired by a classic paper by Walker, Shao, and Volz that used dual quaternions.  The demo script illustrates this application at the end.

# State of This Project
This project was undertaken in the "spare time" of the original developer, motivated solely by academic interest.  As of the initial public release, it does not purport to be a complete and error-free implementation of Lengyel's methods.  As such it is provided "as is" without any warranty, express, or implied.  As a recipient, you acknowledge that you waive any claims against, and indemnify and holds harmless, NASA and its contractors and subcontractors.  The NASA civil service technical point of contact for this project (who is the original developer), is

J. Russell Carpenter

Russell "dot" Carpenter "at" nasa "dot" gov

NASA Goddard Space Flight Center Space Science Mission Operations

Also, in the Fall of 2023 Lengyel began making several changes to the online resources that formed the reference material for this implementation (https://rigidgeometricalgebra.org/), which for example changed the definition of the dot and anti-dot products relative to prior material, among other changes. Minimal effort has been made to ensure all of these changes have been properly incorporated.

# How to Contribute
This project is open to collaborations. If you find this repository useful, even a little, please fork or star the repository! It provides positive feedback to the developer as well as helpful analytics data to the platform maintainers. 

# References
If Geometric Algebra is new to you, there are some great resources online, and a number of textbooks. A few suggestions:

https://en.wikipedia.org/wiki/Geometric_algebra

https://projectivegeometricalgebra.org/

https://bivector.net/

https://enkimute.github.io/ganja.js/

Perwass, C., Geometric Algebra with Applications in Engineering, Springer, 2010, ISBN: 978-3540891680.

MacDonald, A., Geometric Algebra and Calculus, Vols. 1 & 2, CreateSpace Independent Publishing, 2021, ISBN: 978-1453854938, 978-1480132450.

Lengyel, E., Projective Geometric Algebra Illuminated (forthcoming as of Dec 2023), Terathon Software, 2023, ISBN: 979-8985358247.
