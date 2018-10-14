# root-approximations
Functions for finding the roots of floating point numbers, using magic numbers similar to the famous inverse square root function

The magic numbers found are based upon the method described in Hackers Delight.
It should be noted, while they are approximations, they are not nessacarily faster.

std::sqrt seems to consistantly be faster than the approximation, however by contrast the cube root appears to trump both std::pow and std::cbrt:
http://quick-bench.com/ZgMLOIAzaCS1E2Y3EqP2e0l5u9Q

A lot of the magic numbers are really here more out of intrest, rather than having a practical application.