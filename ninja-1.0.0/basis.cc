// Definition of the constructors of the class Basis.


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ninja/momentum.hh>
#include <ninja/spinors.hh>

#include <basis.hh>
using namespace ninja;

namespace {
  const Real BASIS_ESPSILON = 1.e-9;
}

namespace ninja {

  // returns x multiplied by the sign of the second argument
  inline Real sign(const Real & x, const Real & sign)
  {
    return sign < 0. ? -x : x;
  }

  Basis::Basis(const RealMomentum & k1, const RealMomentum & k2)
    : e1(k1), e2(k2), e3(), e4(), r1(0.), r2(0.), mpee12()
  {
    Real k1q = mp2(k1), k2q = mp2(k2), k1k2=mp(k1,k2);
    bool k1n = std::abs(k1q) < BASIS_ESPSILON;
    bool k2n = std::abs(k2q) < BASIS_ESPSILON;
    
    if (k1n && k2n) {
      mpee12 = k1k2;
    } else if (k1n) {
      r2 = 0.5*k2q/k1k2;
      e2 -= r2*k1;
      mpee12 = k1k2;
    } else if (k2n) {
      r1 = 0.5*k1q/k1k2;
      e1 -= r1*k2;
      mpee12 = k1k2;
    } else {
      //Real gamma = k1k2 * (ONE + std::sqrt(1 - (k1q/k1k2) * (k2q/k1k2)));
      Real gamma = k1k2 + sign(ONE,k1k2)*std::sqrt(k1k2*k1k2-k1q*k2q);
      r1 = k1q/gamma;
      r2 = k2q/gamma;
      Real den = ONE - r1*r2;
      e1 -= r1*k2;
      e1 /= den;
      e2 -= r2*k1;
      e2 /= den;
      mpee12 = k1k2/(ONE+r1*r2);
    }
    Spinor sp1(e1);
    Spinor sp2(e2);
    e3 = momentumFromSpinors(sp1,sp2);
    e4 = momentumFromSpinors(sp2,sp1);
  }

  Basis::Basis(const RealMomentum & k)
    : e1(k), e2(), e3(), e4(), r1(0.), r2(0.), mpee12()
  {
    RealMomentum k2(  sign(ONE,k[0])     , -sign(INVSQRT3,k[1]),
                     -sign(INVSQRT3,k[2]), -sign(INVSQRT3,k[3])   );
    Real kq = mp2(k), k1k2=mp(k,k2);
    bool kn = std::abs(kq) < BASIS_ESPSILON;
    if (kn) {
      e2 = k2;
      mpee12 = k1k2;
    } else {
      r1 = 0.5*kq/k1k2;
      e1 -= r1*k2;
      e2 = k2;
      //mpee12 = mp(e1,k2);
      mpee12 = k1k2;
    }
    Spinor sp1(e1);
    Spinor sp2(e2);
    e3 = momentumFromSpinors(sp1,sp2);
    e4 = momentumFromSpinors(sp2,sp1);
  }

} // namespace ninja
