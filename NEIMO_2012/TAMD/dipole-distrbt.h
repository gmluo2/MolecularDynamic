namespace _dipole_distribution_ {
using namespace _distribution_;
using namespace _MM_rdf_;

template <class MM, class CLUSTER> class DIPOLE_DISTRBT_THREAD_VAR : public NEIGHBOR_RDF_THREAD_VAR<MM, CLUSTER> {
public:
	NDISTRIBT<float, long> mu_dt, theta_dt;

	void *func; // the function to calculate the interaction 

	DIPOLE_DISTRBT_THREAD_VAR() {func = NULL;};
	~DIPOLE_DISTRBT_THREAD_VAR() {func = NULL;};
};

class DIPOLE_DISTRIBT_VAR : public NEIGHBOR_RDF_VAR {
public:
	NDISTRIBT<float, float> mu_dt, theta_dt;

	DIPOLE_DISTRIBT_VAR() {};
	~DIPOLE_DISTRIBT_VAR() {};
};

} // end of namespace _dipole_distribution_
