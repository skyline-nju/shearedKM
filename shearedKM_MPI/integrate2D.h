#pragma once
#include <cmath>
#include "particle2D.h"
#include "node.h"


class MotileOscillatorEM {
public:
  MotileOscillatorEM(double h, double D_theta, double D_psi, double v0)
    : h_(h), D_theta_(std::sqrt(24 * D_theta * h)), D_psi_(std::sqrt(24 * D_psi * h)), v0_(v0) {}

  template <class TPar, class TDomain, class TRan>
  void update(TPar& p, const TDomain& dm, TRan& myran) const;

  template <class TPar, class TDomain, class TRan, class TCellList>
  void update_par_cellList(TPar& p, const TDomain& dm, TRan& myran, TCellList& cl) const;

protected:
  double h_;
  double D_theta_;
  double D_psi_;
  double v0_;
};


template<class TPar, class TDomain, class TRan>
void MotileOscillatorEM::update(TPar& p, const TDomain& dm, TRan& myran) const {
  double d_theta = (myran.doub() - 0.5) * D_theta_;
  double d_psi = (p.omega + p.tau_psi) * h_ + (myran.doub() - 0.5) * D_psi_;
  p.theta += d_theta;
  p.psi += d_psi;

  double ux = cos(p.theta) * v0_;
  double uy = sin(p.theta) * v0_;
  p.pos.x += ux * h_;
  p.pos.y += uy * h_;
  dm.tangle(p.pos);

  p.tau_psi = 0;

  if (p.theta >= M_PI) {
    p.theta -= M_PI * 2;
  } else if (p.theta <= -M_PI) {
    p.theta += M_PI * 2;
  }

  if (p.psi >= M_PI) {
    p.psi -= M_PI * 2;
  } else if (p.psi < -M_PI) {
    p.psi += M_PI * 2;
  }
}


template<class TPar, class TDomain, class TRan, class TCellList>
void MotileOscillatorEM::update_par_cellList(TPar& p, const TDomain& dm, TRan& myran, TCellList& cl) const {
  auto ic_old = cl.get_ic(p);
  update(p, dm, myran);
  auto ic_new = cl.get_ic(p);
  if (ic_old != ic_new) {
    cl.update(p, ic_old, ic_new);
  }
}


// recreate cell lists when all particle have moved forward one step
template <typename TNode, typename UniFunc>
void integrate(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl,
               UniFunc f_move, Communicator_2& comm) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    f_move(*it);
  }
  cl.recreate(p_arr);
  comm.comm_after_integration(p_arr, cl);
}

// update cell list once one particle has moved from one cell to another cell
template <typename TNode, typename UniFunc>
void integrate2(std::vector<TNode>& p_arr, CellListNode_2<TNode>& cl,
                UniFunc f_move, Communicator_2& comm) {
  const auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    int ic_old = cl.get_ic(*it);
    f_move(*it);
    int ic_new = cl.get_ic(*it);
    if (ic_old != ic_new) {
      cl.update(*it, ic_old, ic_new);
    }
  }
  comm.comm_after_integration(p_arr, cl);
}
