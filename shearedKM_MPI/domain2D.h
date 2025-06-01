#pragma once
#include "config.h"
#include "vect.h"
#include "particle2D.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

Vec_2<int> get_proc_rank(const Vec_2<int>& proc_size, MPI_Comm group_comm);

template <typename T>
bool within(T x, T y, const Vec_2<T>& o, const Vec_2<T>& end) {
  return x >= o.x && x < end.x && y >= o.y && y < end.y;
}

/**
 * @brief The simulation domain is divided into grids.
 *
 * The mesh size of the grids should be equal to or slightly larger than the
 * cutoff radius of interactions between particles. The grids described by
 * this class are real grids, while the grids in celllist class may contain
 * ghost grids.
 *
 */
class Grid_2 {
public:
  Grid_2(const Vec_2<double>& gl_l, double r_cut,
         const Vec_2<int>& proc_size, MPI_Comm group_comm);

  const Vec_2<int>& n() const { return n_; }
  const Vec_2<int>& gl_n() const { return gl_n_; }
  const Vec_2<int>& origin() const { return origin_; }
  const Vec_2<double>& lc() const { return lc_; }
  const Vec_2<double>& inverse_lc() const { return inverse_lc_; }
private:
  Vec_2<int> n_;
  Vec_2<int> gl_n_;
  Vec_2<int> origin_;
  Vec_2<double> lc_; // length of each cell
  Vec_2<double> inverse_lc_;
};

/*
 * @brief Simulation domain in 2D
 *
 * The simulation domain is usually a rectangle box, with length l_(lx, ly) and
 * origin o_(o_x, o_y). If using MPI, the domain with size gl_l_ is divided into
 * subdomains with size l_ and origin o_.
 *
 * ! The subdomains initialized by the constructor function is of the same size.
 * Use function adjust() to specify the size and origin of each subdomain.
 */
class Domain_2 {
public:
  Domain_2(const Vec_2<double>& gl_l, const Grid_2& grid,
           const Vec_2<int>& proc_size, MPI_Comm group_comm);

  const Vec_2<int>& proc_size() const { return proc_size_; }
  const Vec_2<int>& proc_rank() const { return proc_rank_; }
  const Vec_2<double>& gl_l() const { return gl_l_; }
  const Vec_2<double>& l() const { return l_; }
  const Vec_2<double> origin() const { return o_; }
  MPI_Comm comm() const { return comm_; }

  void find_neighbor(int(*neighbor)[2]) const;


  bool contain_particle(double x, double y) const {
    return within(x, y, o_, o_ + l_);
  }

  template <typename TPar>
  bool contain_particle(const TPar& p) const {
    return contain_particle(p.pos.x, p.pos.y);
  }

protected:
#ifdef USE_MPI
  MPI_Comm comm_;
#endif
  Vec_2<int> proc_size_;
  Vec_2<int> proc_rank_;
  Vec_2<double> gl_l_;
  Vec_2<double> l_;
  Vec_2<double> o_;
};

/**
 * @brief Domain with periodic boundary condition in two dimension.
 *
 * Use flag_PBC_ to indicate whether the periodic bondary condition is applied
 * along x and y direction.
 */
class PeriodicDomain_2 : public Domain_2 {
public:
  PeriodicDomain_2(const Vec_2<double>& gl_l, const Grid_2& grid,
                   const Vec_2<int>& proc_size, MPI_Comm group_comm);

  void tangle(Vec_2<double>& pos) const;
  void untangle(Vec_2<double>& v) const;

  void set_PBC() {} //TODO reset flag_PBC_ when there are walls 
protected:
  Vec_2<double> gl_half_l_;
  Vec_2<bool> flag_PBC_;
};

inline void PeriodicDomain_2::tangle(Vec_2<double>& pos) const {
  if (flag_PBC_.x) {
    tangle_1D(pos.x, gl_l_.x);
  }
  if (flag_PBC_.y) {
    tangle_1D(pos.y, gl_l_.y);
  }
}

inline void PeriodicDomain_2::untangle(Vec_2<double>& r12_vec) const {
  if (flag_PBC_.x) {
    if (r12_vec.x < -gl_half_l_.x) {
      r12_vec.x += gl_l_.x;
    } else if (r12_vec.x > gl_half_l_.x) {
      r12_vec.x -= gl_l_.x;
    }
  }
  if (flag_PBC_.y) {
    if (r12_vec.y < -gl_half_l_.y) {
      r12_vec.y += gl_l_.y;
    } else if (r12_vec.y > gl_half_l_.y) {
      r12_vec.y -= gl_l_.y;
    }
  }
}