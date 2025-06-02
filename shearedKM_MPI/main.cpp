#include "sheared_KM.h"

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  double Lx = 128;
  double Ly = 128;

  double rho0 = atof(argv[1]);
  double D_psi = atof(argv[2]);
  double D_t = atof(argv[3]);
  double gamma = atof(argv[4]);
  double sigma = atof(argv[5]);

  double h = atof(argv[6]);
  int n_step = atoi(argv[7]);
  //int snap_interval = int(round(200 / h * 0.1));
  double snap_log_sep = atof(argv[8]);

  int seed = atoi(argv[9]);
  std::string ini_mode = argv[10];  // should be "rand", "ordered" or "resume"

  int snap_interval = 1;
  if (snap_log_sep > 1) {
    snap_interval = int(snap_log_sep);
    snap_log_sep = -1;
  }
  int n_par = int(rho0 * Lx * Ly);


  Vec_2<int> proc_size = Vec_2<int>(1, 4);


  Vec_2<double> gl_l(Lx, Ly);

  run(gl_l, rho0, D_t, D_psi, gamma, sigma, h, n_step, snap_interval, snap_log_sep, ini_mode, seed, proc_size, MPI_COMM_WORLD);

  MPI_Finalize();
}
