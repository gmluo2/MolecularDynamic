class DIHEDRAL_INTERACT {
public:
	VECTOR3 *r1, *r2, *r3, *r4;
	double phi; // dihedral angle
	double B; // B = cos(phi), see DL_POLY_3 manual

	SVECTOR<double, 3> r12, r23, r34;
	SVECTOR<double, 3> rp123, rp234;
	double Rp123, R2p123, Rp234, R2p234, inv_Rp;

	SVECTOR<double, 3> com_12_12, com_12_23, com_12_34;
	SVECTOR<double, 3> com_23_23, com_23_34, com_34_34;
	SVECTOR<double, 3> dB1, dB2, dB3, dB4;

	SVECTOR<double, 3> virial; // diagonal of virial matrix

	void init();
	void derive_B();
	void cal_virial();

	
};

#define DIH_V3_COMM_0(v1, v2) v1[1] * v2[1] + v1[2] * v2[2]
#define DIH_V3_COMM_1(v1, v2) v1[0] * v2[0] + v1[2] * v2[2]
#define DIH_V3_COMM_2(v1, v2) v1[0] * v2[0] + v1[1] * v2[1]

#define DELTA(i, j) (i == j ? 1 : 0)



class HINGE_CONSTRAINT_FORCE {
public:
	SVECTOR<double, 3> r00, r01, r02; // cordinates of the three atoms (or virtual positions) on CLUSTER #0, relative to position 10
	SVECTOR<double, 3> r11, r12; // cordinates of the two atoms (or virtual positions) on CLUSTER #1, relative to position 10
	SVECTOR<double, 3> f, torq; // force and torque from cluster #1 to cluster #0
	double f_10_00, f_10_01, f_10_02; // forces on atom 00, 01, 02 from 10
	double f_11_00, f_12_00; // forces on atom 00 from 11 and 12
	SVECTOR<double, 3> u_10_00, u_10_01, u_10_02; // unit vector from 10==>00, 10==>01 and 10==>02
	SVECTOR<double, 3> u_11_00, u_12_00; // unit vector from 11==>00, 12==>00

	SVECTOR<double, 3> u_torq, u_f; // the direction of f and torq
	double F, Torq;

	MATRIX<3> m, invm;
	SVECTOR<double, 3> b;

	// important: HINGE is between atoms 00 & 10, the parent hinge of cluster #0, from #0 to #1
	// the force and torque are from cluster #1 to cluster #0
	// So, r10 is the connecting point of cluster #0 to cluster #1,
	// r00, r01, r02 are coordinates in cluster #0, relative to position 10
	// r11 and r12 are the coordinates of atoms in cluster #1, relative to position 10 also

	//      11     12
	//       \    /
	//        \  /
	//         10
	//         ||
	//         ||  <== HINGE
	//         ||
	//         00
	//        /  \
	//       /    \
	//      01     02

	void init(VECTOR3 *r00, VECTOR3 *r01, VECTOR3 *r02, VECTOR3 *r10, VECTOR3 *r11, VECTOR3* r12) {
		// here r00, r01 and r02 are the inertial coordinates in cluster #0
		// here r10, r11 and r12 are the inertial (or global) coordinates in cluster #1
		// so, we need to convert the coordinates of 10, 11, 12 to those relative to 10
		double r;
		V32V3((*r00), this->r00) V32V3((*r01), this->r01) V32V3((*r02), this->r02)
		VECT3((*r10), (*r11), this->r11) VECT3((*r10), (*r12), this->r12)

		V32V3(this->r00, u_10_00) V3ABS2(u_10_00, r) r = 1.0 / sqrt(r); u_10_00.v[0] *= r; u_10_00.v[1] *= r; u_10_00.v[2] *= r;
		V32V3(this->r01, u_10_01) V3ABS2(u_10_01, r) r = 1.0 / sqrt(r); u_10_01.v[0] *= r; u_10_01.v[1] *= r; u_10_01.v[2] *= r;
		V32V3(this->r02, u_10_02) V3ABS2(u_10_02, r) r = 1.0 / sqrt(r); u_10_02.v[0] *= r; u_10_02.v[1] *= r; u_10_02.v[2] *= r;
		VECT3(this->r11, this->r00, u_11_00) V3ABS2(u_11_00, r) r = 1.0 / sqrt(r); u_11_00.v[0] *= r; u_11_00.v[1] *= r; u_11_00.v[2] *= r;
		VECT3(this->r12, this->r00, u_12_00) V3ABS2(u_12_00, r) r = 1.0 / sqrt(r); u_12_00.v[0] *= r; u_12_00.v[1] *= r; u_12_00.v[2] *= r;
	};
	void hinge_force(VECTOR<6> *fHinge) {
		double *v = fHinge->v;
		torq.v[0] = v[0]; torq.v[1] = v[1]; torq.v[2] = v[2]; f.v[0] = v[3]; f.v[1] = v[4]; f.v[2] = v[5];
		V3ABS2(f, F) 
		if (F > 0) {
			F = sqrt(F);  u_f.v[0] = f.v[0] / F; u_f.v[1] = f.v[1] / F; u_f.v[2] = f.v[2] / F;
		} 
		else {F = 0; u_f.v[0] = 1; u_f.v[1] = 0; u_f.v[2] = 0;}
		V3ABS2(torq, Torq); 
		if (Torq > 0) {
			Torq = sqrt(Torq); u_torq.v[0] = torq.v[0] / Torq; u_torq.v[1] = torq.v[1] / Torq; u_torq.v[2] = torq.v[2] / Torq;
		} 
		else {Torq = 0; u_torq.v[0] = 1; u_torq.v[1] = 0; u_torq.v[2] = 0;}
	};
	bool cal_constraint_force();
};

