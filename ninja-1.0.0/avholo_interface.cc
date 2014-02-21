#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <ninja/avholo.hh>

namespace ninja {

  AvHOneLoop avh_olo;
  Real AvHOneLoop::mur2_ = 0;
  Real AvHOneLoop::ir_threshold_ = 1.0e-10;

}


#ifdef NINJA_USE_ONELOOP
#include <avholo_decls.hh>

#ifdef NINJA_USE_ONELOOP_WITH_CACHE
#include <ninja_hash_table.hh>
#endif

using namespace ninja;

// Common methods of cached and uncached OneLoop

namespace {

  inline void
  write_avh_result(const AvhOloComplex avh_res[3], Complex res[3])
  {
    res[0] = avh_res[0];
    res[1] = avh_res[1];
    res[2] = avh_res[2];
  }

} // namespace

namespace ninja {

  void AvHOneLoop::setInfraredThreshold(Real threshold)
  {
    ir_threshold_ = threshold;
    avh_olo_onshell_(ir_threshold_);
  }

  void AvHOneLoop::init(Real muRsq)
  {
    avh_olo_onshell_(ir_threshold_);
    avh_olo_mu_set_(std::sqrt(muRsq));
    if (mur2_ != muRsq) {
      clearIntegralCache();
      mur2_ = muRsq;
    }
  }

  void
  AvHOneLoop::getBubbleIntegralCM(Complex result[3],
                                  Real p1,
                                  const Complex & m1, const Complex & m2)
  {
    avh_olo_b0c_(result,p1,m1,m2);
  }

  void
  AvHOneLoop::getBubbleIntegralRM(Complex result[3],
                                  Real p1,
                                  Real m1, Real m2)
  {
    avh_olo_b0m_(result,p1,m1,m2);
  }

}  // namespace ninja

#ifdef NINJA_USE_ONELOOP_WITH_CACHE

// methods for AvHOlo with cache

namespace {

  struct  MIsResult {
    Complex data[3];
  };

  struct MIsRank2BubbleResult {
    Complex data11[3];
    Complex data1[3];
    Complex data0[3];
  };


  struct BoxArgsCM {
    Real arg1[6];
    Complex arg2[4];
  };
  inline bool operator== (const BoxArgsCM & x, const BoxArgsCM & y)
  {
    return (((x).arg1[0] == (y).arg1[0]) &&    
            ( (x).arg1[1] == (y).arg1[1] ) &&
            ( (x).arg1[2] == (y).arg1[2] ) &&
            ( (x).arg1[3] == (y).arg1[3] ) &&
            ( (x).arg1[4] == (y).arg1[4] ) &&
            ( (x).arg1[5] == (y).arg1[5] ) &&
            ( (x).arg2[0] == (y).arg2[0] ) &&
            ( (x).arg2[1] == (y).arg2[1] ) &&
            ( (x).arg2[2] == (y).arg2[2] ) &&
            ( (x).arg2[3] == (y).arg2[3] ));
  }

  struct BoxArgsRM {
    Real arg1[6];
    Real arg2[4];
  };
  inline bool operator== (const BoxArgsRM & x, const BoxArgsRM & y)
  {
    return (((x).arg1[0] == (y).arg1[0]) &&    
            ( (x).arg1[1] == (y).arg1[1] ) &&
            ( (x).arg1[2] == (y).arg1[2] ) &&
            ( (x).arg1[3] == (y).arg1[3] ) &&
            ( (x).arg1[4] == (y).arg1[4] ) &&
            ( (x).arg1[5] == (y).arg1[5] ) &&
            ( (x).arg2[0] == (y).arg2[0] ) &&
            ( (x).arg2[1] == (y).arg2[1] ) &&
            ( (x).arg2[2] == (y).arg2[2] ) &&
            ( (x).arg2[3] == (y).arg2[3] ));
  }

  struct BoxArgsNM {
    Real arg1[6];
  };
  inline bool operator== (const BoxArgsNM & x, const BoxArgsNM & y)
  {
    return (((x).arg1[0] == (y).arg1[0]) &&    
            ( (x).arg1[1] == (y).arg1[1] ) &&
            ( (x).arg1[2] == (y).arg1[2] ) &&
            ( (x).arg1[3] == (y).arg1[3] ) &&
            ( (x).arg1[4] == (y).arg1[4] ) &&
            ( (x).arg1[5] == (y).arg1[5] ));
  }

  struct  TriangleArgsCM {
    Real arg1[3];
    Complex arg2[3];
  };
  inline bool
  operator== (const TriangleArgsCM & x, const TriangleArgsCM & y)
  {
    return (((x).arg1[0] == (y).arg1[0]) &&
            ( (x).arg1[1] == (y).arg1[1] ) &&
            ( (x).arg1[2] == (y).arg1[2] ) &&
            ( (x).arg2[0] == (y).arg2[0] ) &&
            ( (x).arg2[1] == (y).arg2[1] ) &&
            ( (x).arg2[2] == (y).arg2[2] ));
  }

  struct TriangleArgsRM {
    Real arg1[3];
    Real arg2[3];
  };
  inline bool
  operator== (const TriangleArgsRM & x, const TriangleArgsRM & y)
  {
    return (((x).arg1[0] == (y).arg1[0]) &&
            ( (x).arg1[1] == (y).arg1[1] ) &&
            ( (x).arg1[2] == (y).arg1[2] ) &&
            ( (x).arg2[0] == (y).arg2[0] ) &&
            ( (x).arg2[1] == (y).arg2[1] ) &&
            ( (x).arg2[2] == (y).arg2[2] ));
  }

  struct  TriangleArgsNM {
    Real arg1[3];
  };
  inline bool
  operator== (const TriangleArgsNM & x, const TriangleArgsNM & y)
  {
    return (((x).arg1[0] == (y).arg1[0]) &&
            ( (x).arg1[1] == (y).arg1[1] ) &&
            ( (x).arg1[2] == (y).arg1[2] ));
  }

  struct BubbleArgsCM {
    Real arg1;
    Complex arg2[2];
  };
  inline bool operator== (const BubbleArgsCM & x, const BubbleArgsCM & y)
  {
    return (((x).arg1 == (y).arg1) &&
            ( (x).arg2[0] == (y).arg2[0] ) &&
            ( (x).arg2[1] == (y).arg2[1] ));
  }

  struct  BubbleArgsRM {
    Real arg1;
    Real arg2[2];
  };
  inline bool operator== (const BubbleArgsRM & x, const BubbleArgsRM & y)
  {
    return (((x).arg1 == (y).arg1) &&
            ( (x).arg2[0] == (y).arg2[0] ) &&
            ( (x).arg2[1] == (y).arg2[1] ));
  }

  struct BubbleArgsNM {
    Real arg1;
  };
  inline bool operator== (const BubbleArgsNM & x, const BubbleArgsNM & y)
  {
    return (((x).arg1 == (y).arg1));
  }

  struct TadpoleArgsCM {
    Complex arg2;
  };
  inline bool operator== (const TadpoleArgsCM & x, const TadpoleArgsCM & y)
  {
    return (( (x).arg2 == (y).arg2 ));
  }

  struct  TadpoleArgsRM {
    Real arg2;
  };
  inline bool operator== (const TadpoleArgsRM & x, const TadpoleArgsRM & y)
  {
    return (( (x).arg2 == (y).arg2 ));
  }


  // Hash Tables
  HashTable<BoxArgsCM,MIsResult> ht_4cm;
  HashTable<BoxArgsRM,MIsResult> ht_4rm;
  HashTable<TriangleArgsCM,MIsResult> ht_3cm;
  HashTable<TriangleArgsRM,MIsResult> ht_3rm;
  HashTable<BubbleArgsCM,MIsRank2BubbleResult> ht_2cm;
  HashTable<BubbleArgsRM,MIsRank2BubbleResult> ht_2rm;
  HashTable<TadpoleArgsCM,MIsResult> ht_1cm;
  HashTable<TadpoleArgsRM,MIsResult> ht_1rm;
#ifdef NINJA_MASSLESS
  HashTable<BoxArgsNM,MIsResult> ht_4nm;
  HashTable<TriangleArgsNM,MIsResult> ht_3nm;
  HashTable<BubbleArgsNM,MIsRank2BubbleResult> ht_2nm;
#endif

  // Initial number of buckets.  As long as the cache is cleared
  // (hence maintaining the buckets in memory) and not freed, changing
  // these numbers shouldn't make much difference.
  const std::size_t N_BOXES = 30; // 750
  const std::size_t N_TRIANGLES = 30; // 350
  const std::size_t N_BUBBLES = 30; // 75
  const std::size_t N_TADPOLES = 1; // 2

} // namespace


namespace ninja {

  void
  AvHOneLoop::getBoxIntegralCM(Complex result[3],
                               Real p1, Real p2,
                               Real p3, Real p4,
                               Real p12, Real p23,
                               const Complex & m1, const Complex & m2,
                               const Complex & m3, const Complex & m4)
  {
    BoxArgsCM args = {{p1,p2,p3,p12,p23},{m1,m2,m3,m4}};
    MIsResult *val;

    if (ht_4cm.empty())
      ht_4cm.resize(N_BOXES);
  
    if (ht_4cm.find(args,val)) {
      result[0] = val->data[0];
      result[1] = val->data[1];
      result[2] = val->data[2];
    } else {
      avh_olo_d0c_(result,p1,p2,p3,p4,p12,p23,m1,m2,m3,m4);
      val->data[0] = result[0];
      val->data[1] = result[1];
      val->data[2] = result[2];
    }
  }

  void
  AvHOneLoop::getBoxIntegralRM(Complex result[3],
                               Real p1, Real p2,
                               Real p3, Real p4,
                               Real p12, Real p23,
                               Real m1, Real m2,
                               Real m3, Real m4)
  {
    BoxArgsRM args = {{p1,p2,p3,p12,p23},{m1,m2,m3,m4}};
    MIsResult *val;

    if (ht_4rm.empty())
      ht_4rm.resize(N_BOXES);
  
    if (ht_4rm.find(args,val)) {
      result[0] = val->data[0];
      result[1] = val->data[1];
      result[2] = val->data[2];
    } else {
      avh_olo_d0m_(result,p1,p2,p3,p4,p12,p23,m1,m2,m3,m4);
      val->data[0] = result[0];
      val->data[1] = result[1];
      val->data[2] = result[2];
    }
  }

  void
  AvHOneLoop::getBoxIntegralNM(Complex result[3],
                               Real p1, Real p2, Real p3, Real p4,
                               Real p12, Real p23)
  {
    const AvhOloReal msq(0);
    BoxArgsNM args = {{p1,p2,p3,p12,p23}};
    MIsResult *val;

    if (ht_4nm.empty())
      ht_4nm.resize(N_BOXES);
 
    if (ht_4nm.find(args,val)) {
      result[0] = val->data[0];
      result[1] = val->data[1];
      result[2] = val->data[2];
    } else {
      avh_olo_d0m_(result,p1,p2,p3,p4,p12,p23,msq,msq,msq,msq);
      val->data[0] = result[0];
      val->data[1] = result[1];
      val->data[2] = result[2];
    }
  }


  void
  AvHOneLoop::getTriangleIntegralCM(Complex result[3],
                                    Real p1, Real p2,
                                    Real p3,
                                    const Complex & m1, const Complex & m2,
                                    const Complex & m3)
  {
    TriangleArgsCM args = {{p1,p2,p3},{m1,m2,m3}};
    MIsResult *val;

    if (ht_3cm.empty())
      ht_3cm.resize(N_TRIANGLES);
  
    if (ht_3cm.find(args,val)) {
      result[0] = val->data[0];
      result[1] = val->data[1];
      result[2] = val->data[2];
    } else {
      avh_olo_c0c_(result,p1,p2,p3,m1,m2,m3);
      val->data[0] = result[0];
      val->data[1] = result[1];
      val->data[2] = result[2];
    }
  }

  void
  AvHOneLoop::getTriangleIntegralRM(Complex result[3],
                                    Real p1, Real p2,
                                    Real p3,
                                    Real m1, Real m2,
                                    Real m3)
  {
    TriangleArgsRM args = {{p1,p2,p3},{m1,m2,m3}};
    MIsResult *val;

    if (ht_3rm.empty())
      ht_3rm.resize(N_TRIANGLES);
  
    if (ht_3rm.find(args,val)) {
      result[0] = val->data[0];
      result[1] = val->data[1];
      result[2] = val->data[2];
    } else {
      avh_olo_c0m_(result,p1,p2,p3,m1,m2,m3);
      val->data[0] = result[0];
      val->data[1] = result[1];
      val->data[2] = result[2];
    }
  }

  void
  AvHOneLoop::getTriangleIntegralNM(Complex result[3],
                                    Real p1, Real p2,
                                    Real p3)
  {
    const AvhOloReal msq(0);
    TriangleArgsNM args = {{p1,p2,p3}};
    MIsResult *val;

    if (ht_3nm.empty())
      ht_3nm.resize(N_TRIANGLES);
  
    if (ht_3nm.find(args,val)) {
      result[0] = val->data[0];
      result[1] = val->data[1];
      result[2] = val->data[2];
    } else {
      avh_olo_c0m_(result,p1,p2,p3,msq,msq,msq);
      val->data[0] = result[0];
      val->data[1] = result[1];
      val->data[2] = result[2];
    }
  }


  void
  AvHOneLoop::getRank2BubbleIntegralCM(Complex b11[3],
                                       Complex b1[3], Complex b0[3],
                                       Real p1,
                                       const Complex & m1, const Complex & m2)
  {
    BubbleArgsCM args = {p1,{m1,m2}};
    MIsRank2BubbleResult *val;

    if (ht_2cm.empty())
      ht_2cm.resize(N_BUBBLES);
  
    if (ht_2cm.find(args,val)) {
      b11[0] = val->data11[0];
      b11[1] = val->data11[1];
      b11[2] = val->data11[2];
      b1[0] = val->data1[0];
      b1[1] = val->data1[1];
      b1[2] = val->data1[2];
      b0[0] = val->data0[0];
      b0[1] = val->data0[1];
      b0[2] = val->data0[2];
    } else {
      AvhOloComplex avh_res00[3];
      avh_olo_b11c_(b11,avh_res00,b1,b0,p1,m1,m2);
      val->data11[0] = b11[0];
      val->data11[1] = b11[1];
      val->data11[2] = b11[2];
      val->data1[0] = b1[0];
      val->data1[1] = b1[1];
      val->data1[2] = b1[2];
      val->data0[0] = b0[0];
      val->data0[1] = b0[1];
      val->data0[2] = b0[2];
    }
  }

  void
  AvHOneLoop::getRank2BubbleIntegralRM(Complex b11[3],
                                       Complex b1[3], Complex b0[3],
                                       Real p1,
                                       Real m1, Real m2)
  {
    BubbleArgsRM args = {p1,{m1,m2}};
    MIsRank2BubbleResult *val;

    if (ht_2rm.empty())
      ht_2rm.resize(N_BUBBLES);
  
    if (ht_2rm.find(args,val)) {
      b11[0] = val->data11[0];
      b11[1] = val->data11[1];
      b11[2] = val->data11[2];
      b1[0] = val->data1[0];
      b1[1] = val->data1[1];
      b1[2] = val->data1[2];
      b0[0] = val->data0[0];
      b0[1] = val->data0[1];
      b0[2] = val->data0[2];
    } else {
      AvhOloComplex avh_res00[3];
      avh_olo_b11m_(b11,avh_res00,b1,b0,p1,m1,m2);
      val->data11[0] = b11[0];
      val->data11[1] = b11[1];
      val->data11[2] = b11[2];
      val->data1[0] = b1[0];
      val->data1[1] = b1[1];
      val->data1[2] = b1[2];
      val->data0[0] = b0[0];
      val->data0[1] = b0[1];
      val->data0[2] = b0[2];
    }
  }

  void
  AvHOneLoop::getRank2BubbleIntegralNM(Complex b11[3],
                                     Complex b1[3], Complex b0[3],
                                     Real p1)
  {
    const AvhOloReal msq(0);
    BubbleArgsNM args = {p1};
    MIsRank2BubbleResult *val;

    if (ht_2nm.empty())
      ht_2nm.resize(N_BUBBLES);
  
    if (ht_2nm.find(args,val)) {
      b11[0] = val->data11[0];
      b11[1] = val->data11[1];
      b11[2] = val->data11[2];
      b1[0] = val->data1[0];
      b1[1] = val->data1[1];
      b1[2] = val->data1[2];
      b0[0] = val->data0[0];
      b0[1] = val->data0[1];
      b0[2] = val->data0[2];
    } else {
      AvhOloComplex avh_res00[3];
      avh_olo_b11m_(b11,avh_res00,b1,b0,p1,msq,msq);
      val->data11[0] = b11[0];
      val->data11[1] = b11[1];
      val->data11[2] = b11[2];
      val->data1[0] = b1[0];
      val->data1[1] = b1[1];
      val->data1[2] = b1[2];
      val->data0[0] = b0[0];
      val->data0[1] = b0[1];
      val->data0[2] = b0[2];
    }
  }


  void
  AvHOneLoop::getTadpoleIntegralCM(Complex result[3], const Complex & m0)
  {
    TadpoleArgsCM args = {m0};
    MIsResult *val;

    if (ht_1cm.empty())
      ht_1cm.resize(N_TADPOLES);
  
    if (ht_1cm.find(args,val)) {
      result[0] = val->data[0];
      result[1] = val->data[1];
      result[2] = val->data[2];
    } else {
      avh_olo_a0c_(result,m0);
      val->data[0] = result[0];
      val->data[1] = result[1];
      val->data[2] = result[2];
    }
  }

  void
  AvHOneLoop::getTadpoleIntegralRM(Complex result[3], Real m0)
  {
    TadpoleArgsRM args = {m0};
    MIsResult *val;

    if (ht_1rm.empty())
      ht_1rm.resize(N_TADPOLES);
  
    if (ht_1rm.find(args,val)) {
      result[0] = val->data[0];
      result[1] = val->data[1];
      result[2] = val->data[2];
    } else {
      avh_olo_a0m_(result,m0);
      val->data[0] = result[0];
      val->data[1] = result[1];
      val->data[2] = result[2];
    }
  }


  void AvHOneLoop::clearIntegralCache()
  {
    if (!ht_4cm.empty())
      ht_4cm.clear();
    if (!ht_3cm.empty())
      ht_3cm.clear();
    if (!ht_2cm.empty())
      ht_2cm.clear();
#if 0 // clearing tadpole-MIs doesn't make much sense
    if (!ht_1cm.empty())
      ht_1cm.clear();
#endif

    if (!ht_4rm.empty())
      ht_4rm.clear();
    if (!ht_3rm.empty())
      ht_3rm.clear();
    if (!ht_2rm.empty())
      ht_2rm.clear();
#if 0 // clearing tadpole-MIs doesn't make much sense
    if (!ht_1rm.empty())
      ht_1rm.clear();
#endif

    if (!ht_4nm.empty())
      ht_4nm.clear();
    if (!ht_3nm.empty())
      ht_3nm.clear();
    if (!ht_2nm.empty())
      ht_2nm.clear();
  }


  void AvHOneLoop::freeIntegralCache()
  {
    if (!ht_4cm.empty())
      ht_4cm.free();
    if (!ht_3cm.empty())
      ht_3cm.free();
    if (!ht_2cm.empty())
      ht_2cm.free();
    if (!ht_1cm.empty())
      ht_1cm.free();

    if (!ht_4rm.empty())
      ht_4rm.free();
    if (!ht_3rm.empty())
      ht_3rm.free();
    if (!ht_2rm.empty())
      ht_2rm.free();
    if (!ht_1rm.empty())
      ht_1rm.free();

    if (!ht_4nm.empty())
      ht_4nm.free();
    if (!ht_3nm.empty())
      ht_3nm.free();
    if (!ht_2nm.empty())
      ht_2nm.free();
  }

} // namespace ninja

#else // NINJA_USE_ONELOOP_WITH_CACHE

// OneLoop without cache

namespace ninja {

  void
  AvHOneLoop::getBoxIntegralCM(Complex result[3],
                               Real p1, Real p2,
                               Real p3, Real p4,
                               Real p12, Real p23,
                               const Complex & m1, const Complex & m2,
                               const Complex & m3, const Complex & m4)
  {
    avh_olo_d0c_(result,p1,p2,p3,p4,p12,p23,m1,m2,m3,m4);
  }

  void
  AvHOneLoop::getBoxIntegralRM(Complex result[3],
                               Real p1, Real p2,
                               Real p3, Real p4,
                               Real p12, Real p23,
                               Real m1, Real m2,
                               Real m3, Real m4)
  {
    avh_olo_d0m_(result,p1,p2,p3,p4,p12,p23,m1,m2,m3,m4);
  }

  void
  AvHOneLoop::getTriangleIntegralCM(Complex result[3],
                                    Real p1, Real p2,
                                    Real p3,
                                    const Complex & m1, const Complex & m2,
                                    const Complex & m3)
  {
    avh_olo_c0c_(result,p1,p2,p3,m1,m2,m3);
  }

  void
  AvHOneLoop::getTriangleIntegralRM(Complex result[3],
                                    Real p1, Real p2,
                                    Real p3,
                                    Real m1, Real m2,
                                    Real m3)
  {
    avh_olo_c0m_(result,p1,p2,p3,m1,m2,m3);
  }

  void
  AvHOneLoop::getRank2BubbleIntegralCM(Complex b11[3],
                                       Complex b1[3], Complex b0[3],
                                       Real p1,
                                       const Complex & m1, const Complex & m2)
  {
    AvhOloComplex avh_res00[3];
    avh_olo_b11c_(b11,avh_res00,b1,b0,p1,m1,m2);
  }

  void
  AvHOneLoop::getRank2BubbleIntegralRM(Complex b11[3],
                                       Complex b1[3], Complex b0[3],
                                       Real p1,
                                       Real m1, Real m2)
  {
    AvhOloComplex avh_res00[3];
    avh_olo_b11m_(b11,avh_res00,b1,b0,p1,m1,m2);
  }

  void
  AvHOneLoop::getTadpoleIntegralCM(Complex result[3], const Complex & m0)
  {
    avh_olo_a0c_(result,m0);
  }

  void
  AvHOneLoop::getTadpoleIntegralRM(Complex result[3], Real m0)
  {
    avh_olo_a0m_(result,m0);
  }

  void AvHOneLoop::clearIntegralCache() {}
  void AvHOneLoop::freeIntegralCache() {}

} // namespace ninja

#endif // ! NINJA_USE_ONELOOP_WITH_CACHE
#else // NINJA_USE_ONELOOP

// No OneLoop library --> dummy implementation

#include <iostream>
#include <stdexcept>
using namespace std;

namespace {
  
  void avh_olo_disabled()
  {
    cerr << "ERROR IN NINJA: "
         << "No intergral library has been specified.\n"
         << "ERROR IN NINJA: Pleas enable OneLoop and/or LoopTools"
         << " during configuration.\n"
         << "ERROR IN NINJA: Alternatively, specify another custom library."
         << endl;
    NINJA_THROW (logic_error("No integral library."));
  }

}

namespace ninja {

  void AvHOneLoop::init(Real)
  {
    avh_olo_disabled();
  }

  void
  AvHOneLoop::getBubbleIntegralCM(Complex *,
                                  Real,
                                  const Complex &, const Complex &)
  {
    avh_olo_disabled();
  }

  void
  AvHOneLoop::getBubbleIntegralRM(Complex *,
                                  Real, Real, Real)
  {
    avh_olo_disabled();
  }

  void
  AvHOneLoop::getBoxIntegralCM(Complex *,
                               Real, Real, Real, Real, Real, Real,
                               const Complex &, const Complex &,
                               const Complex &, const Complex &)
  {
    avh_olo_disabled();
  }

  void
  AvHOneLoop::getBoxIntegralRM(Complex *,
                               Real, Real, Real, Real, Real, Real,
                               Real, Real, Real, Real)
  {
    avh_olo_disabled();
  }

  void
  AvHOneLoop::getTriangleIntegralCM(Complex *,
                                    Real, Real, Real,
                                    const Complex &, const Complex &,
                                    const Complex &)
  {
    avh_olo_disabled();
  }

  void
  AvHOneLoop::getTriangleIntegralRM(Complex *,
                                    Real, Real, Real,
                                    Real, Real, Real)
  {
    avh_olo_disabled();
  }

  void
  AvHOneLoop::getRank2BubbleIntegralCM(Complex *, Complex *,
                                       Complex *,
                                       Real,
                                       const Complex &, const Complex &)
  {
    avh_olo_disabled();
  }

  void
  AvHOneLoop::getRank2BubbleIntegralRM(Complex *, Complex *,
                                       Complex *,
                                       Real, Real, Real)
  {
    avh_olo_disabled();
  }

  void
  AvHOneLoop::getTadpoleIntegralCM(Complex *, const Complex &)
  {
    avh_olo_disabled();
  }

  void
  AvHOneLoop::getTadpoleIntegralRM(Complex *, Real)
  {
    avh_olo_disabled();
  }

  void AvHOneLoop::clearIntegralCache()
  {
    avh_olo_disabled();
  }

  void AvHOneLoop::freeIntegralCache()
  {
    avh_olo_disabled();
  }

} // namespace ninja

#endif // ! NINJA_USE_ONELOOP


#if !defined(NINJA_USE_ONELOOP_WITH_CACHE) || !defined(NINJA_USE_ONELOOP)

namespace ninja {

  void
  AvHOneLoop::getBoxIntegralNM(Complex rslt[3],
                               Real s21, Real s32, Real s43,
                               Real s14, Real s31, Real s42)
  {
    getBoxIntegralRM(rslt, s21, s32, s43, s14, s31, s42,
                     Real(), Real(), Real(), Real());
  }

  void
  AvHOneLoop::getTriangleIntegralNM(Complex rslt[3],
                                    Real s21, Real s32, Real s13)
  {
    getTriangleIntegralRM(rslt, s21, s32, s13, Real(), Real(), Real());
  }

  void
  AvHOneLoop::getRank2BubbleIntegralNM(Complex b11[3],
                                       Complex b1[3], Complex b0[3],
                                       Real s21)
  {
    getRank2BubbleIntegralRM(b11, b1, b0, s21, Real(), Real());
  }

} // namespace ninja

#endif // ! NINJA_USE_ONELOOP_WITH_CACHE
