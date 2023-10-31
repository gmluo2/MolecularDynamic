// functions in Z-matrix to define a new atom cordinate

// r1 -- adjacent atom; d -- shift vector
void ZMatrix(VECTOR3 &r1, VECTOR3 &d, VECTOR3 &res);

// r1, r2 are atoms' position, bl -- bound distance 
// axis -- normal to new-1-2 plane, R_1-2 X R_new-1
// ba -- bound angle between new-1-2
void ZMatrix(VECTOR3 &r1, VECTOR3 &r2, VECTOR3 &axis, float bl, float ba, VECTOR3 &res);

// bl -- bound distance; ba -- bound angle between new-1-2
// ta -- torsion angle between new-1-2 & 3-2-1 along axis 2-1
// r1, r2 & r3 -- atoms' position
void ZMatrix(VECTOR3 &r1, VECTOR3 &r2, VECTOR3 &r3, float bl, float ba, float ta, VECTOR3 &res);


// Torsion Angle -- to get r1 with r2-r3-r4 axis
// reverse procedure of ZMatrix(r1, r2, r3, bound-length, bound-angle)
double TorsionAngle(VECTOR3 &r1, VECTOR3 &r2, VECTOR3 &r3, VECTOR3 &r4); // arc degree [-PI, PI]
// bound angle between r21 and r23 -- 1-2-3
double BoundAngle(VECTOR3 &r1, VECTOR3 &r2, VECTOR3 &r3); // arc degree [0, PI]

