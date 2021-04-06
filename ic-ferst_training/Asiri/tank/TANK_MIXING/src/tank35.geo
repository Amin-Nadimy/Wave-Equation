R_t=25;
R_c=0.5;
R_c2=0.3;

in=Asin(R_c/R_t);
in2=Asin(R_c2/R_t);
in_mid=((R_t^2-R_c^2)^(0.5));
in_circ=1;


d1=(90*Pi)/180; //i1
d2=(130*Pi)/180; //o
d3=(45*Pi)/180;  //i2

e1=0.5;
e2=0.05;

Point(1)={0,0,0,e1};

Point(100)={R_t*Cos(d1-in),R_t*Sin(d1-in),0,e1};
Point(2)={R_t*Cos(d1+in),R_t*Sin(d1+in),0,e1};

Point(102)={R_t*Cos(d2-in),R_t*Sin(d2-in),0,e1};
Point(103)={R_t*Cos(d2+in),R_t*Sin(d2+in),0,e1};

Point(3)={-R_t,0,0,e1};
Point(4)={0,-R_t, 0,e1};
Point(5)={R_t,0,0,e1};

Point(104)={R_t*Cos(d3-in),R_t*Sin(d3-in),0,e1};
Point(105)={R_t*Cos(d3+in),R_t*Sin(d3+in),0,e1};

Line(1)={100,2};
Circle(2)={2, 1, 102};
Line(3)={102,103};
Circle(4)={103, 1, 3};
Circle(5)={3,1,4};
Circle(6)={4,1,5};
Circle(7)={5,1,104};
Line(8)={104,105};
Circle(9)={105,1,100};

Line Loop(10)={1,2,3,4,5,6,7,8,9};


Plane Surface(6) = {10};

Extrude {0,0,6} {
  Surface{6};
}




Point(200)={in_mid*Cos(d1), in_mid*Sin(d1), in_circ, e2};
Point(201)={in_mid*Cos(d1), in_mid*Sin(d1),in_circ-0.3, e2};
Point(202)={in_mid*Cos(d1), in_mid*Sin(d1), in_circ+0.3, e2};
Point(203)={R_t*Cos(d1-in2),R_t*Sin(d1-in2),in_circ,e2};
Point(204)={R_t*Cos(d1+in2),R_t*Sin(d1+in2),in_circ,e2};

Circle(205)={201,200,203};
Circle(206)={203, 200, 202};
Circle(207)={202, 200, 204};
Circle(208)={204, 200, 201};

//Line(205)={201,203};
//Line(206)={203, 202};
//Line(207)={202, 204};
//Line(208)={204,  201};

Point(300)={in_mid*Cos(d2), in_mid*Sin(d2), in_circ, e2};
Point(301)={in_mid*Cos(d2), in_mid*Sin(d2), in_circ-0.3, e2};
Point(302)={in_mid*Cos(d2), in_mid*Sin(d2), in_circ+0.3, e2};
Point(303)={R_t*Cos(d2-in2),R_t*Sin(d2-in2),in_circ,e2};
Point(304)={R_t*Cos(d2+in2),R_t*Sin(d2+in2),in_circ,e2};

Circle(305)={301,300,303};
Circle(306)={303, 300, 302};
Circle(307)={302, 300, 304};
Circle(308)={304, 300, 301};

//Line(305)={301,303};
//Line(306)={303,  302};
//Line(307)={302, 304};
//Line(308)={304,  301};

Point(400)={in_mid*Cos(d3), in_mid*Sin(d3), in_circ, e2};
Point(401)={in_mid*Cos(d3), in_mid*Sin(d3), in_circ-0.3, e2};
Point(402)={in_mid*Cos(d3), in_mid*Sin(d3), in_circ+0.3, e2};
Point(403)={R_t*Cos(d3-in2),R_t*Sin(d3-in2),in_circ,e2};
Point(404)={R_t*Cos(d3+in2),R_t*Sin(d3+in2),in_circ,e2};

Circle(405)={401,400,403};
Circle(406)={403, 400, 402};
Circle(407)={402, 400, 404};
Circle(408)={404, 400, 401};

//Line(405)={401, 403};
//Line(406)={403,  402};
//Line(407)={402,  404};
//Line(408)={404,  401};



Line Loop(409) = {405, 406, 407, 408}; //i2
Plane Surface(410) = {409};

Line Loop(411) = {205, 206, 207, 208}; //inlet 1
Plane Surface(412) = {411};

Line Loop(413) = {305, 306, 307, 308}; //o2
Plane Surface(414) = {413};


Line Loop(415) = {47, 19, -51, -8};
Plane Surface(416) = {409, 415};

Line Loop(417) = {22, 12, -23, -1};
Plane Surface(418) = {411, 417};

Line Loop(419) = {27, 14, -31, -3};
Plane Surface(420) = {413, 419};

//PHSYICAL SURFACES
Physical Surface(1) = {44, 40, 36, 48, 56, 28, 6, 416, 418, 420}; //walls
Physical Surface(2) = {57};//top
Physical Surface(3) = {412}; //i1
Physical Surface(4) = {410}; //i2
Physical Surface(5) = {414}; //o2


Delete {
  Volume{1};
}
Delete {
  Surface{52, 24, 32};
}

Surface Loop(429) = {40, 6, 418, 412, 28, 420, 414, 36, 57, 44, 48, 416, 410, 56};
Volume(430) = {429};
Physical Volume(431) = {430};


