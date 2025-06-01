#include "sheared_KM.h"

void run(const Vec_2<double>& gl_l,
         double rho0, double sigma,
         double D_theta, double D_psi,
         double h, double v0,
         int n_step, int snap_dt, double snap_log_sep,
         const std::string& ini_mode,
         unsigned long long seed,
         const Vec_2<int>& proc_size,
         MPI_Comm group_comm) {
  typedef BiNode<ActiveOscillator_2> node_t;
  int my_rank, tot_proc;
  MPI_Comm_rank(group_comm, &my_rank);
  MPI_Comm_size(group_comm, &tot_proc);
  std::vector<node_t> p_arr;
  const double r_cut = 1.0;

  Grid_2 grid(gl_l, r_cut, proc_size, group_comm);
  PeriodicDomain_2 dm(gl_l, grid, proc_size, group_comm);
  CellListNode_2<node_t> cl(dm, grid);
  Communicator_2 comm(dm, grid, rho0, 20.);

  // cal force
  AlignKernal kernal(r_cut);
  auto f1 = [&kernal](node_t* p1, node_t* p2) {
    kernal(*p1, *p2);
    };
  auto f2 = [&kernal, &dm](node_t* p1, node_t* p2) {
    kernal(*p1, *p2, dm);
    };

  //auto f3 = [&kernal](node_t* p1, node_t* p2, const Vec_2<double>& offset) {
  //  kernal(*p1, *p2, p2->pos - p1->pos + offset);
  //  };
  auto for_all_pair_force = [&cl, &f1, &f2]() {
    cl.for_each_pair(f1, f2);
    //cl.for_each_pair_fast(f1, f3);
    };



  // set output
  char basename[255];
  char log_file[255];
  char op_file[255];
  char gsd_file[255];
#ifdef _MSC_VER
  char folder[] = "data\\";
#else
  char folder[255];
  snprintf(folder, 255, "/home/ps/data/motile_KM/bimodal/L%g/", gl_l.x);
#endif

  char log_folder[255];
  char op_folder[255];
  snprintf(log_folder, 255, "%slog%s", folder, delimiter.c_str());
  snprintf(op_folder, 255, "%sop%s", folder, delimiter.c_str());
  if (my_rank == 0) {
    mkdir(log_folder);
    mkdir(op_folder);
  }

  snprintf(basename, 255, "L%g_%g_r%g_v%g_T%g_s%g_D%.4f_h%g_S%d",
           gl_l.x, gl_l.y, rho0, v0, D_psi, sigma, D_theta, h, seed);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);

  int log_dt = 10000;
  int gl_np = round(rho0 * gl_l.x * gl_l.y);
  int start = 0;
  io::Snap_GSD_2 gsd(gsd_file, n_step, snap_dt, start, h, snap_log_sep, gl_l, gl_np, ini_mode, group_comm);

  int log_interval = 10000;
  int op_interval = 100;
  snprintf(log_file, 255, "%s%s_t%d.dat", log_folder, basename, start);
  io::LogExporter log(log_file, n_step, log_dt, start, gl_np, group_comm);
  snprintf(op_file, 255, "%s%s_t%d.dat", op_folder, basename, start);
  io::OrderParaExporter op(op_file, start, n_step, op_interval, gl_np, group_comm);

  Ranq1 myran(seed + my_rank + start);


  // initialize particles and celllist
  ini(p_arr, dm, cl, rho0, sigma, ini_mode, myran, gsd, 20.);

  // ini integrator
  MotileOscillatorEM integrator(h, D_theta, D_psi, v0);

  auto one_par_move = [&integrator, &dm, &myran, h](node_t& p) {
    integrator.update(p, dm, myran);
  };

  // run
  for (int t = 1; t <= n_step; t++) {
    cal_force(p_arr, cl, comm, for_all_pair_force);

#ifdef POS_OMEGA
    for (auto& p : p_arr) {
      if (p.pos.x < gl_l.x * 0.5) {
        p.omega = sigma;
      } else {
        p.omega = -sigma;
      }
    }
#endif

    integrate(p_arr, cl, one_par_move, comm);
    gsd.dump(t, p_arr);
    log.record(t);
    op.dump(t, p_arr);
  }

  if (my_rank == 0) {
    std::cout << "Finish simulation!" << std::endl;
  }
  MPI_Barrier(group_comm);
}