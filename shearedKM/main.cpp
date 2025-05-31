#include <iostream>
#include "domain2D.h"
#include "cellList2D.h"
#include "particle2D.h"
#include "force2D.h"
#include "integrate2D.h"
#include "exporter2D.h"
#include "comn.h"

int main(int argc, char* argv[]) {
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
  std::string ini_mode = argv[10];  // should be "bimodal" or "resume"

  int snap_interval = 1;
  if (snap_log_sep > 1) {
    snap_interval = int(snap_log_sep);
    snap_log_sep = -1;
  }
  int n_par = int(rho0 * Lx * Ly);

  typedef BiNode<PassiveOscillator_2> node_t;
  Ranq2 myran(seed);
  Vec_2<double> gl_l(Lx, Ly);
  double r_cut = 1;
  Grid_2 grid(gl_l, r_cut);
  LeesEdwardsDomain_2 dm(gl_l);
  CellListNode_2<node_t> cl(dm, grid);
  std::vector<node_t> p_arr;

  // ini integrator
  ShearedOscillatorEM integrator(h, D_t, D_psi, gamma);

  // cal force
  AlignKernal kernal(r_cut);
  auto f1 = [&kernal](node_t* p1, node_t* p2) {
    kernal(*p1, *p2);
    };
  auto f2 = [&kernal, &dm](node_t* p1, node_t* p2) {
    kernal(*p1, *p2, dm);
    };

  // set output
  char basename[255];
  char log_file[255];
  char op_file[255];
  char gsd_file[255];
#ifdef _MSC_VER
  char folder[] = "data\\";
#else
  char folder[] = "/mnt/d/code/KM/data/";
#endif

  char log_folder[255];
  char op_folder[255];
  snprintf(log_folder, 255, "%slog%s", folder, delimiter.c_str());
  snprintf(op_folder, 255, "%sop%s", folder, delimiter.c_str());
  mkdir(log_folder);
  mkdir(op_folder);

  snprintf(basename, 255, "L%g_%g_r%g_Dt%g_T%g_g%g_s%g_h%g_S%d",
    Lx, Ly, rho0, D_t, D_psi, gamma, sigma, h, seed);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);

  int start = 0;

  exporter::Snap_GSD_2 gsd(gsd_file, n_step, snap_interval, start, h, snap_log_sep, gl_l, ini_mode);

  int log_interval = 10000;
  int op_interval = 100;
  snprintf(log_file, 255, "%s%s_t%d.dat", log_folder, basename, start);
  exporter::LogExporter log(log_file, start, n_step, log_interval, n_par);

  snprintf(op_file, 255, "%s%s_t%d.dat", op_folder, basename, start);
  exporter::OrderParaExporter op(op_file, start, n_step, op_interval);

  // ini particles
  ini_particles(p_arr, myran, ini_mode, n_par, gl_l, sigma, gsd);
  cl.create(p_arr);

  for (int t = 1; t <= n_step; t++) {
    dm.update_dx(t * h, gamma);
    cl.for_each_pair(f1, f2, dm.get_dx());
    for (int i = 0; i < n_par; i++) {
      integrator.update(p_arr[i], dm, myran);
      //integrator.update_par_cellList(p_arr[i], pdm, myran, cl);
    }
    cl.recreate(p_arr);
    gsd.dump(t, p_arr);
    op.dump(t, p_arr);
    log.record(t);
  }
}
