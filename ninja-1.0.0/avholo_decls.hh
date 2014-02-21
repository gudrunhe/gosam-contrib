#ifndef AVHOLO_DECLS_HH
#define AVHOLO_DECLS_HH

#include <tmp_utils.hh>
#include <ninja/types.hh>

// real-type used by OneLoop
typedef ninja::Real AvhOloReal;

// complex-type used by OneLoop
typedef ninja::Complex AvhOloComplex;


// Declarations of functions in OneLoop

extern "C" {
  
  // set the threshold to distinguish between IR-divergent and
  // IR-finite cases
  void avh_olo_onshell_(const AvhOloReal & thrs);

  // set the renormalization scale
  void avh_olo_mu_set_(const AvhOloReal & mu);

  // 4-point MI
  void avh_olo_d0m_(AvhOloComplex rslt[3],
                    const AvhOloReal & p1, const AvhOloReal & p2,
                    const AvhOloReal & p3, const AvhOloReal & p4,
                    const AvhOloReal & p12, const AvhOloReal & p23,
                    const AvhOloReal & m1, const AvhOloReal & m2,
                    const AvhOloReal & m3, const AvhOloReal & m4);
  void avh_olo_d0c_(AvhOloComplex rslt[3],
                    const AvhOloComplex & p1, const AvhOloComplex & p2,
                    const AvhOloComplex & p3, const AvhOloComplex & p4,
                    const AvhOloComplex & p12, const AvhOloComplex & p23,
                    const AvhOloComplex & m1, const AvhOloComplex & m2,
                    const AvhOloComplex & m3, const AvhOloComplex & m4);

  // 3-point MI
  void avh_olo_c0m_(AvhOloComplex rslt[3],
                    const AvhOloReal & p1,const AvhOloReal & p2,
                    const AvhOloReal & p3,
                    const AvhOloReal & m1, const AvhOloReal & m2,
                    const AvhOloReal & m3);
  void avh_olo_c0c_(AvhOloComplex rslt[3],
                    const AvhOloComplex & p1,const AvhOloComplex & p2,
                    const AvhOloComplex & p3,
                    const AvhOloComplex & m1, const AvhOloComplex & m2,
                    const AvhOloComplex & m3);

  // 2-point MIs
  void avh_olo_b0m_(AvhOloComplex rslt[3], const AvhOloReal & p1,
                    const AvhOloReal & m1, const AvhOloReal & m2 );
  void avh_olo_b11m_(AvhOloComplex b11[3], AvhOloComplex b00[3],
                     AvhOloComplex b1[3], AvhOloComplex b0[3],
                     const AvhOloReal & p1, const AvhOloReal & m1,
                     const AvhOloReal & m2);
  void avh_olo_b0c_(AvhOloComplex rslt[3], const AvhOloComplex & p1,
                    const AvhOloComplex & m1, const AvhOloComplex & m2 );
  void avh_olo_b11c_(AvhOloComplex b11[3], AvhOloComplex b00p[3],
                     AvhOloComplex b1[3], AvhOloComplex b0[3],
                     const AvhOloComplex & p1, const AvhOloComplex & m1,
                     const AvhOloComplex & m2);
  
  // 1-point MI
  void avh_olo_a0m_(AvhOloComplex rslt[3], const AvhOloReal & m0);
  void avh_olo_a0c_(AvhOloComplex rslt[3], const AvhOloComplex & m0);

} // extern "C"

#endif // AVHOLO_DECLS_HH
