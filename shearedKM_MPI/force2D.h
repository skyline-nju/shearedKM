#pragma once
#include "vect.h"
#include <cmath>
#include "communicator2D.h"

class AlignKernal {
public:
  AlignKernal(double r_cut) : r_cut_(r_cut), r_cut_square_(r_cut* r_cut) {}

  template <typename TPar>
  void align(TPar& p1, TPar& p2) const {
    double tau = sin(p2.psi - p1.psi);
    p1.tau_psi += tau;
    p2.tau_psi -= tau;
  }

  template <typename TPar>
  void operator ()(TPar& p1, TPar& p2) const;

  template <typename TPar, typename BdyCondi>
  void operator ()(TPar& p1, TPar& p2, const BdyCondi& bc) const;

private:
  double r_cut_;
  double r_cut_square_;
};


template<typename TPar>
void AlignKernal::operator()(TPar& p1, TPar& p2) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  if (r12.square() < r_cut_square_) {
    align(p1, p2);
  }
}

template<typename TPar, typename BdyCondi>
void AlignKernal::operator()(TPar& p1, TPar& p2, const BdyCondi& bc) const {
  Vec_2<double> r12 = p2.pos - p1.pos;
  bc.untangle(r12);
  if (r12.square() < r_cut_square_) {
    align(p1, p2);
  }
}


template <typename TNode, typename TFunc, typename TDomain>
void cal_force(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl,
               Communicator_2& comm, TFunc for_all_pair_force, const TDomain& dm) {
  int n_ghost = 0;
  comm.comm_before_cal_force(p_arr, cl, n_ghost, dm);
  for_all_pair_force();
  comm.clear_padded_particles(cl, p_arr, n_ghost);
}
