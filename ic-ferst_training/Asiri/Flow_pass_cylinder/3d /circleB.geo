d=0.01;
Li=5*d; //12.5*d;
Lo=25*d;
W=3.145*d;
H=30*d;
e=d/30;
e2=2*d;

Point (1) = {Li, H/2, 0, e};

Point (2) = {Li+d/2, H/2, 0, e};
Point (3) = {Li, H/2+d/2, 0, e};
Point (4) = {Li-d/2, H/2, 0, e};
Point (5) = {Li, H/2-d/2, 0, e};

Point (6) = {0, 0, 0, e2};
Point (7) = {Li+Lo, 0, 0, e2};
Point (8) = {Li+Lo, H, 0, e2};
Point (9) = {0, H, 0, e2};


Circle (1) = {2, 1, 3};
Circle (2) = {3, 1, 4};
Circle (3) = {4, 1, 5};
Circle (4) = {5, 1, 2};


Line (5) = {9, 6};
Line (6) = {9, 8};
Line (7) = {6, 7};
Line (8) = {8, 7};


Line Loop (10) = {6, 8, -7, -5};
Line Loop (11) = {1,2,3,4};

Ruled Surface(102) = {10,11};

//+
Extrude {0, 0, W} {
	Surface{102}; 
}
//+
//Extrude {0, 0, W} {
//	Surface{102}; Layers{20}; 
//}

//+
Physical Surface(1) = {102}; //back (z)
Physical Surface(2) = {144}; //front (z)
Physical Surface(3) = {115}; //top
Physical Surface(4) = {123}; //floor
Physical Surface(5) = {127}; //inlet
Physical Surface(6) = {119}; //outlet
Physical Surface(7) = {131,143,135,139}; //cylinder

Physical Volume(1) = {1}; //VOL
