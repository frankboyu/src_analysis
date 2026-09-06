// C++ interface to the original Fortran edved_ implementation.
// Compile and link this file with get_edved_wkng_pol_callable.f.

#include <stdexcept>
#include <string>

namespace edved {

struct Parameters {
    double bgamman = 0.0;
    double sgamman = 0.0;
    double bphin = 0.0;
    double sphin = 0.0;
};

struct Result {
    double crs0 = 0.0;
    double crs = 0.0;
    double tcrs0 = 0.0;
    double tcrs = 0.0;
    double pd = 0.0;
    double thd = 0.0;
    double pvm = 0.0;
    double thvm = 0.0;
    double crs0_m1 = 0.0;
    double crs0_0 = 0.0;
    double crs0_1 = 0.0;
    double crs_m1 = 0.0;
    double crs_0 = 0.0;
    double crs_1 = 0.0;
};

// The Fortran source uses default REAL, which is normally 4-byte float.
extern "C" void edved_(
    int* in, int* ivm,
    float* ei, float* q2, float* q0, float* epsl, float* t,
    float* crs0, float* crs, float* tcrs0, float* tcrs,
    float* pd, float* thd, float* pvm, float* thvm,
    float* beamenergy, float* bgamman, float* sgamman,
    float* bphin, float* sphin);

class EdvedModel {
public:
    explicit EdvedModel(std::string) {}
    EdvedModel() = default;

    void initialize(int meson, const Parameters& parameters) {
        if (meson != 1 && meson != 3 && meson != 4)
            throw std::invalid_argument("meson must be 1 (rho), 3 (phi), or 4 (J/psi)");

        int in = 1;
        int ivm = meson;
        float zero = 0.0f;
        float beamenergy = 0.0f;
        float bgamman = static_cast<float>(parameters.bgamman);
        float sgamman = static_cast<float>(parameters.sgamman);
        float bphin = static_cast<float>(parameters.bphin);
        float sphin = static_cast<float>(parameters.sphin);
        callFortran(in, ivm, zero, zero, zero, zero, zero, zero, zero, zero, zero,
                    zero, zero, zero, zero, beamenergy, bgamman, sgamman, bphin, sphin);
        initialized_ = true;
        meson_ = meson;
    }

    Result calculate(double ei, double q2, double q0, double epsl, double t) const {
        if (!initialized_)
            throw std::logic_error("EdvedModel::initialize must be called first");

        int in = 0;
        int ivm = meson_;
        float fei = static_cast<float>(ei);
        float fq2 = static_cast<float>(q2);
        float fq0 = static_cast<float>(q0);
        float fepsl = static_cast<float>(epsl);
        float ft = static_cast<float>(t);
        float crs0 = 0.0f, crs = 0.0f, tcrs0 = 0.0f, tcrs = 0.0f;
        float pd = 0.0f, thd = 0.0f, pvm = 0.0f, thvm = 0.0f;
        float zero = 0.0f;
        callFortran(in, ivm, fei, fq2, fq0, fepsl, ft, crs0, crs, tcrs0, tcrs,
                    pd, thd, pvm, thvm, zero, zero, zero, zero, zero);

        Result result;
        result.crs0 = crs0;
        result.crs = crs;
        result.tcrs0 = tcrs0;
        result.tcrs = tcrs;
        result.pd = pd;
        result.thd = thd;
        result.pvm = pvm;
        result.thvm = thvm;
        return result;
    }

private:
    static void callFortran(
        int& in, int& ivm, float& ei, float& q2, float& q0, float& epsl, float& t,
        float& crs0, float& crs, float& tcrs0, float& tcrs,
        float& pd, float& thd, float& pvm, float& thvm,
        float& beamenergy, float& bgamman, float& sgamman, float& bphin, float& sphin) {
        edved_(&in, &ivm, &ei, &q2, &q0, &epsl, &t, &crs0, &crs, &tcrs0, &tcrs,
               &pd, &thd, &pvm, &thvm, &beamenergy, &bgamman, &sgamman, &bphin, &sphin);
    }

    bool initialized_ = false;
    int meson_ = 3;
};

// Evaluate phi photoproduction for one (photon energy, -t) point.
// Parameter order follows the requested interface: sgamman, bgamman,
// sphin, bphin.
Result getPhiResult(double photonEnergy, double minusT,
                   double sgamman, double bgamman,
                   double sphin, double bphin) {
    if (photonEnergy <= 0.0)
        throw std::invalid_argument("photonEnergy must be positive");
    if (minusT < 0.0)
        throw std::invalid_argument("minusT must be non-negative");

    EdvedModel model("input");
    model.initialize(3, {bgamman, sgamman, bphin, sphin});
    return model.calculate(0.0, 0.0, photonEnergy, 0.0, -minusT);
}

double getPhiCrossSection(double photonEnergy, double minusT,
                          double sgamman, double bgamman,
                          double sphin, double bphin) {
    return getPhiResult(photonEnergy, minusT, sgamman, bgamman, sphin, bphin).crs;
}

} // namespace edved
