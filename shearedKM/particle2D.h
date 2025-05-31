#pragma once
#include "vect.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


template<typename TRan, typename TPar>
void ini_particles(std::vector<TPar>& p_arr, TRan& myran, const std::string& ini_mode,
                   int n_par, const Vec_2<double>& gl_l, double sigma) {
  p_arr.reserve(n_par);
  if (ini_mode == "bimodal") {
    for (int i = 0; i < n_par; i++) {
      p_arr.emplace_back(myran, gl_l, Vec_2<double>());
      if (i < n_par / 2) {
        p_arr[i].omega = sigma;
      } else {
        p_arr[i].omega = -sigma;
      }
    }
  }
}

template<typename TRan, typename TPar, typename TSnap>
void ini_particles(std::vector<TPar>& p_arr, TRan& myran, const std::string& ini_mode,
                   int n_par, const Vec_2<double>& gl_l, double sigma, TSnap& snap) {
  p_arr.reserve(n_par);
  if (ini_mode == "rand" || ini_mode == "ordered") {
    for (int i = 0; i < n_par; i++) {
      p_arr.emplace_back(myran, gl_l, Vec_2<double>());
      if (i < n_par / 2) {
        p_arr[i].omega = sigma;
      } else {
        p_arr[i].omega = -sigma;
      }
    }
    if (ini_mode == "ordered") {
      for (auto& p : p_arr) {
        p.psi = 0;
      }
    }
  } else if (ini_mode == "resume") {
    snap.read_last_frame(p_arr);
  } else {
    std::cout << "IC must be one of 'rand', 'ordered', 'resume'." << std::endl;
    exit(1);
  }
}

class PassiveOscillator_2 {
public:
  PassiveOscillator_2() : pos(), psi() {}
  PassiveOscillator_2(const Vec_2<double>& pos0) : pos(pos0), psi() {}
  PassiveOscillator_2(const Vec_2<double>& pos0, double theta0): pos(pos0), psi() {}
  PassiveOscillator_2(const Vec_2<double>& pos0, double theta0, double psi0): pos(pos0), psi(psi0) {}
  template <typename TRan>
  PassiveOscillator_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o);

  double get_psi() const { return psi; }

  Vec_2<double> pos;
  double psi;         // phase for oscillation
  double omega = 0;   // intrinsic frequence
  double tau_psi = 0;
};


template<typename TRan>
PassiveOscillator_2::PassiveOscillator_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o) {
  pos.x = myran.doub() * l.x + o.x;
  pos.y = myran.doub() * l.y + o.y;
  psi = myran.doub() * M_PI * 2;
}


/*
class ActiveOscillator_2 {
public:
  ActiveOscillator_2() : pos(), theta(), psi() {}
  ActiveOscillator_2(const Vec_2<double>& pos0) : pos(pos0), theta(), psi() {}
  ActiveOscillator_2(const Vec_2<double>& pos0, double theta0)
    : pos(pos0), theta(theta0), psi() {}
  ActiveOscillator_2(const Vec_2<double>& pos0, double theta0, double psi0)
  : pos(pos0), theta(theta0), psi(psi0) {}

  double get_theta() const { return theta; }
  double get_psi() const { return psi; }


  template <typename TRan>
  ActiveOscillator_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o);

  template <typename TRan>
  ActiveOscillator_2(TRan& myran, double theta0, const Vec_2<double>& l, const Vec_2<double>& o);

  Vec_2<double> pos;
  double theta;     // moving direction
  double psi;       // phase for oscillation, decoupled from theta
  double omega = 0;
  double tau_psi = 0;
};


template<typename TRan>
ActiveOscillator_2::ActiveOscillator_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o) {
  pos.x = myran.doub() * l.x + o.x;
  pos.y = myran.doub() * l.y + o.y;
  theta = myran.doub() * M_PI * 2;
  psi = myran.doub() * M_PI * 2;
}

template<typename TRan>
ActiveOscillator_2::ActiveOscillator_2(TRan& myran, double theta0, const Vec_2<double>& l, const Vec_2<double>& o) {
  pos.x = myran.doub() * l.x + o.x;
  pos.y = myran.doub() * l.y + o.y;
  theta = theta0;
  psi = myran.doub() * M_PI * 2;
}
*/
