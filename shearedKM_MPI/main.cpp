#include "sheared_KM.h"

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  double Lx = 512;
  double Ly = 256;

  double rho0 = atof(argv[1]);
  double T = atof(argv[2]);
  double v0 = 1;
  double D_theta = atof(argv[3]);
  double sigma = atof(argv[4]);

  double h = atof(argv[5]);
  int n_step = atoi(argv[6]);
  //int snap_interval = int(round(200 / h * 0.1));
  double snap_log_sep = atof(argv[7]);

  int seed = atoi(argv[8]);
  std::string ini_mode = argv[9];  // should be "bimodal" or "resume"

  int snap_interval = 1;
  if (snap_log_sep > 1) {
    snap_interval = int(snap_log_sep);
    snap_log_sep = -1;
  }
  int n_par = int(rho0 * Lx * Ly);


  Vec_2<int> proc_size = Vec_2<int>(2, 2);


  Vec_2<double> gl_l(Lx, Ly);

  run(gl_l, rho0, sigma,
      D_theta, T, h, v0, 
      n_step, snap_interval, snap_log_sep,
      ini_mode, seed, proc_size,
      MPI_COMM_WORLD);

  MPI_Finalize();

}
