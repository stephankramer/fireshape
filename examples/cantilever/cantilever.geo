Point(1) = {0, -4, 0, 0.5};
Point(2) = {0, -2, 0, 0.5};
// Point(7) = {10, 0, 0, 0.5};
Point(3) = {0, 2, 0, 0.5};
Point(4) = {0, 4, 0, 0.5};
Point(5) = {20, 1, 0, 0.5};
Point(6) = {20, -1, 0, 0.5};

Line(1) = {4, 3};
Line(2) = {2, 1};
Line(3) = {5, 6};
Line(4) = {3, 2};
Line(5) = {1, 6};
Line(6) = {4, 5};
// Line(7) = {3, 7};

Curve Loop(1) = {5, -3, -6, 1, 4, 2};
Plane Surface(1) = {1};
Physical Surface(7) = {1};

Physical Curve(1) = {2, 1};
Physical Curve(2) = {3};
Physical Curve(3) = {5, 6, 4};
