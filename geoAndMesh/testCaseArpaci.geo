Point(1) = {0, 0, 0, 0.5};
Point(2) = {0.5, 0, 0, 0.5};
Point(3) = {0.5, 0.5, 0, 0.5};
Point(4) = {0, 0.5, 0, 0.5};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};
Transfinite Surface {6} = {4, 3, 2, 1};
Transfinite Line {2, 4} = 3 Using Progression 1;
Transfinite Line {3, 1} = 3 Using Progression 1;
Recombine Surface {6};
Physical Line(10)={2,3};
Physical Line(11)={4,1};
Physical Surface(100) = {6};
