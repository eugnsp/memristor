// At least GMSH 4.0 is required

lc1 = 1;
lc2 = 5;

Point(1) = {0, 0, 0, lc1};
Point(2) = {1, 0, 0, lc1};
Point(3) = {1, 5, 0, lc1};
Point(4) = {0, 5, 0, lc1};
Point(5) = {15, 0, 0, lc2};
Point(6) = {15, 40, 0, lc2};
Point(7) = {0, 40, 0, lc1};
Point(8) = {0, 23, 0, lc1};
Point(9) = {0, 17, 0, lc1};
Point(10) = {1.5, 20, 0, lc1};
Point(11) = {0, 20, 0, lc1};
Point(12) = {0, 60, 0, lc2};
Point(13) = {15, 60, 0, lc2};

Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 6};
Line(4) = {6, 7};
Line(5) = {7, 8};
Line(6) = {8, 9};
Line(7) = {9, 4};
Line(8) = {4, 1};
Line(9) = {2, 3};
Line(10) = {3, 4};
Ellipse(11) = {8, 11, 10, 10};
Ellipse(12) = {10, 11, 10, 9};
Line(13) = {7, 12};
Line(14) = {12, 13};
Line(15) = {13, 6};

Curve Loop(1) = {2, 3, 4, 5, 11, 12, 7, -10, -9};
Curve Loop(2) = {1, 9, 10, 8};
Curve Loop(3) = {6, -12, -11};
Curve Loop(4) = {-4, -15, -14, -13};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
