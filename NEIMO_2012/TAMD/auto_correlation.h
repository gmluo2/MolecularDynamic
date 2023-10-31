namespace _velocity_correlation_ {

	using namespace _distribution_;

	class VATOM {
	public:
		int uid;
		//double v[3];
	};
	class VCELL {
	public:
		ARRAY<VATOM> a;
	};

	void init_vcell(MMOL_MD_CELL &mdcell, VCELL &vcell) {
	};
	void update_vcell(MMOL_MD_CELL &mdcell, VCELL &vcell);

	class V_AutoCorr {
	public:
		int mtype, mindx, cindx;
		DISTRIBT<float, float> dt;
	};

	bool read_vhist(char *ftitle, int istart, MMOL_MD_CELL &mdcell, int mtype, int mindx, int cindx, BUFFER< SVECTOR<float, 3> > &vhist);
	void auto_correlation_float3(BUFFER< SVECTOR<float, 3> > &vhist, DISTRIBT<int, double> &dt);

} // end of namespace _velocity_correlation_
