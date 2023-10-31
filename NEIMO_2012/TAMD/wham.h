namespace _WHAM_ {
	class DISTRIBT {
	public:
		_PROFILE_1D_<float> x, y, s; // s is the error bar
		_PROFILE_1D_<float> wb; // biased potential
		void release() {x.release(); y.release(); s.release(); wb.release();};
		void set(int npts) {
			if (x.N != npts) x.set(npts);
			if (y.N != npts) y.set(npts);
			if (s.N != npts) s.set(npts);
			if (wb.N != npts) wb.set(npts);
		};
	};
	class WHAM_DISTRIBT : public DISTRIBT {
	public:
		int Nc; // number of samples / configurations used in the distribution-calculation
		float F; // the related free energy shift / constant

		WHAM_DISTRIBT() {Nc = 0; F = 0;};
	};

	template <class T, class INDX> class CombPt {
	public:
		T v;
		ARRAY<INDX> indx;
	};

	class PtIndx {
	public:
		int iw, ipt;
	};

	class WHAM {
	public:
		ARRAY<WHAM_DISTRIBT> fi; // the distributions with different biased windows
		_PROFILE_EB_< float, CombPt<float, PtIndx>, float > f; // combined distribution
		_PROFILE_1D_<float> Fb; // backup for F of each window

		bool init_unbaised_prof();
		void construct_unbaised_prof();
		void calculate_F();
		void backup_F();
		float max_diff_F();
	};

	bool IterateWHAM(WHAM &wham, float dF_max, int max_loops);

	bool wham_read(char *parfname, WHAM &wham);
}
