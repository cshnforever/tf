# tf
Transfer Function for Linear System

#Tfalg 
#-> Just Polynomials! -> You can put Any Roots or coefficients(But Cannot affect the roots)
#To find roots for arbitrary coefficients, use find_roots()
#It's just a OPTIMIZATION PROBLEM.
#If f(s)=a_ns^n+....+a_1s+a_0=h(s)+ik(s) -> Minimize |f(s)|^2!!
#1st, 2nd order is just calculation...
#But over 3rd order, it is Numerically calculated (Basic Iteration is 100000) -> So, Overlapped Roots cannot be detected easily!!

#(NOTE)
#If (Absolute value of Imaginary part) is <0.01, This code makes it zero! (to use only real parts!)
#So, Please modify this value to get a more accurate roots.

#Tf
#-> You can put Zeros, and Poles
#-> But, if you put coefficients directly, it doesn't calculates its zeros and poles....
#-> So, Use above Tfalg's find_roots()
#Also, It can make (Partial Fractions)!! -> Use tf_divide()
#And, it can make Response in Time!(From Laplace Transformation)
