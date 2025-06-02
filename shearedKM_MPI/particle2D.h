#pragma once
#include "vect.h"
#include <iostream>
#include <cstring>
#define M_PI 3.14159265358979323846

template <typename T>
void tangle_1D(T& x, T L) {
  if (x < 0) {
    x += L;
  } else if (x >= L) {
    x -= L;
  }
}

class PassiveOscillator_2 {
public:
  PassiveOscillator_2() : pos(), psi() {}
  PassiveOscillator_2(const Vec_2<double>& pos0) : pos(pos0), psi() {}
  PassiveOscillator_2(const Vec_2<double>& pos0, double theta0) : pos(pos0), psi() {}
  PassiveOscillator_2(const Vec_2<double>& pos0, double theta0, double psi0) : pos(pos0), psi(psi0) {}
  template <typename TRan>
  PassiveOscillator_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o);

  double get_psi() const { return psi; }

  template <typename TInt>
  void copy_pos_psi_to(char* dest, TInt& idx) const;

  template <typename TInt>
  void copy_pos_psi_from(const char* src, TInt& idx);

  template <typename TInt>
  void copy_to(char* dest, TInt& idx) const;

  template <typename TInt>
  void copy_from(const char* dest, TInt& idx);

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


template<typename TInt>
void PassiveOscillator_2::copy_pos_psi_to(char* dest, TInt& idx) const {
  std::memcpy(dest + idx, &pos.x, 24);
  idx += 24;
}

template<typename TInt>
void PassiveOscillator_2::copy_pos_psi_from(const char* src, TInt& idx) {
  std::memcpy(&pos.x, src + idx, 24);
  idx += 24;
}

template<typename TInt>
void PassiveOscillator_2::copy_to(char* dest, TInt& idx) const {
  std::memcpy(dest + idx, &pos.x, 32);
  idx += 32;
}

template<typename TInt>
void PassiveOscillator_2::copy_from(const char* dest, TInt& idx) {
  std::memcpy(&pos.x, dest + idx, 32);
  idx += 32;
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

  template <typename TRan>
  ActiveOscillator_2(TRan& myran, const Vec_2<double>& l, const Vec_2<double>& o);

  template <typename TRan>
  ActiveOscillator_2(TRan& myran, double theta0, const Vec_2<double>& l, const Vec_2<double>& o);


  double get_theta() const { return theta; }
  double get_psi() const { return psi; }

  template <typename TInt>
  void copy_pos_psi_to(char* dest, TInt& idx) const;

  template <typename TInt>
  void copy_pos_psi_from(const char* src, TInt& idx);

  template <typename TInt>
  void copy_to(char* dest, TInt& idx) const;

  template <typename TInt>
  void copy_from(const char* dest, TInt& idx);

  Vec_2<double> pos;
  double psi;       // phase for oscillation, decoupled from theta
  double theta;     // moving direction
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

template<typename TInt>
void ActiveOscillator_2::copy_pos_psi_to(char* dest, TInt& idx) const {
  std::memcpy(dest + idx, &pos.x, 24);
  idx += 24;
}

template<typename TInt>
void ActiveOscillator_2::copy_pos_psi_from(const char* src, TInt& idx) {
  std::memcpy(&pos.x, src + idx, 24);
  idx += 24;
}

template<typename TInt>
void ActiveOscillator_2::copy_to(char* dest, TInt& idx) const {
  std::memcpy(dest + idx, &pos.x, 40);
  idx += 40;
}

template<typename TInt>
void ActiveOscillator_2::copy_from(const char* dest, TInt& idx) {
  std::memcpy(&pos.x, dest + idx, 40);
  idx += 40;
}
*/