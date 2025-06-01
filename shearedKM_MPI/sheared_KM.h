#pragma once

#include <iostream>
#include "domain2D.h"
#include "cellList2D.h"
#include "particle2D.h"
#include "force2D.h"
#include "integrate2D.h"
#include "io2D.h"
#include "config.h"
#include "communicator2D.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

template <typename TDomain, typename TPar, typename TRan, typename TCellList>
void ini(std::vector<TPar>& p_arr, const TDomain& dm, TCellList& cl,
         double rho0, double sigma, const std::string& ini_mode,
         TRan& myran, io::Snap_GSD_2& gsd, double amplification) {
  int my_rank;
  MPI_Comm_rank(dm.comm(), &my_rank);
  double area_gl = dm.gl_l().x * dm.gl_l().y;
  double area = dm.l().x * dm.l().y;
  int n_gl = int(round(rho0 * area_gl));

  size_t n_max_per_core = size_t(rho0 * area * amplification);
  if (n_max_per_core > n_gl) {
    n_max_per_core = n_gl;
  }
  p_arr.reserve(n_max_per_core);

  double* x = new double[n_gl] {};
  double* y = new double[n_gl] {};
  double* theta = new double[n_gl] {};
  double* psi = new double[n_gl] {};
  double* omega = new double[n_gl] {};

  if (my_rank == 0) {
    if (ini_mode == "bimodal") {
      for (size_t i = 0; i < n_gl; i++) {
        x[i] = myran.doub() * dm.gl_l().x;
        y[i] = myran.doub() * dm.gl_l().y;
        theta[i] = M_PI * myran.doub() * 2;
        psi[i] = M_PI * myran.doub() * 2;
        if (i < n_gl / 2) {
          omega[i] = sigma;
        } else {
          omega[i] = -sigma;
        }
      }
      std::cout << "Create " << n_gl << " particles with random pos and ori" << std::endl;
    } else if (ini_mode == "resume") {
      gsd.read_last_frame(x, y, theta, psi, omega);
    } else {
      std::cout << "ini_mode need be one of rand or resume" << std::endl;
      exit(1);
    }

  }
  MPI_Barrier(dm.comm());
  MPI_Bcast(x, n_gl, MPI_DOUBLE, 0, dm.comm());
  MPI_Bcast(y, n_gl, MPI_DOUBLE, 0, dm.comm());
  MPI_Bcast(theta, n_gl, MPI_DOUBLE, 0, dm.comm());
  MPI_Bcast(psi, n_gl, MPI_DOUBLE, 0, dm.comm());
  MPI_Bcast(omega, n_gl, MPI_DOUBLE, 0, dm.comm());


  for (size_t i = 0; i < n_gl; i++) {
    if (dm.contain_particle(x[i], y[i])) {
      size_t j = p_arr.size();
      p_arr.push_back(TPar());
      p_arr[j].pos.x = x[i];
      p_arr[j].pos.y = y[i];
      p_arr[j].theta = theta[i];
      p_arr[j].psi = psi[i];
      p_arr[j].omega = omega[i];
    }
  }

  int my_np = p_arr.size();
  int tot_np;
  MPI_Reduce(&my_np, &tot_np, 1, MPI_INT, MPI_SUM, 0, dm.comm());
  if (my_rank == 0) {
    if (tot_np != n_gl) {
      std::cout << "Error when creating particles with tot_np=" << tot_np
        << " while n_gl=" << n_gl << std::endl;
      exit(1);
    }
  }
  cl.create(p_arr);

  MPI_Barrier(dm.comm());
  delete[] x;
  delete[] y;
  delete[] theta;
  delete[] psi;
  delete[] omega;
}




void run(const Vec_2<double>& gl_l,
         double rho0, double sigma,
         double D_theta, double D_psi,
         double h, double v0,
         int n_step, int snap_dt, double snap_log_sep,
         const std::string& ini_mode,
         unsigned long long seed,
         const Vec_2<int>& proc_size,
         MPI_Comm group_comm);
