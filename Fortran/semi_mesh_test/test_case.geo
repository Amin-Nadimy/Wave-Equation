//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {2, 0, 0, 1.0};
//+
Point(3) = {1, 1.732, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 1};
//+
Line Loop(1) = {1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Surface(1) = {1};
//+
Physical Line(2) = {3};
//+
Physical Line(3) = {1};
//+
Physical Line(4) = {2};
