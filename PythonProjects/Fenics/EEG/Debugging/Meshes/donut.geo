dh = 0.1;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {-1, 0, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0, -1, 0, 1.0};
//+
Circle(1) = {3, 1, 3};
//+
Circle(2) = {4, 1, 2};
//+
Circle(3) = {2, 1, 5};
//+
Circle(4) = {5, 1, 3};
//+
Point(6) = {dh, 0, 0, 1.0};
//+
Point(7) = {-dh, 0, 0, 1.0};
//+
Point(8) = {0, dh, 0, 1.0};
//+
Point(9) = {0, -dh, 0, 1.0};
//+
Circle(5) = {7, 1, 8};
//+
Circle(6) = {8, 1, 6};
//+
Circle(7) = {6, 1, 9};
//+
Circle(8) = {9, 1, 7};//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {6, 7, 8, 5};
//+
Plane Surface(1) = {1, 2};
//+
Physical Surface("mainSurface", 10) = {1};
