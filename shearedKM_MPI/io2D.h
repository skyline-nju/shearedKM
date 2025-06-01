#pragma once
#include <vector>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include "config.h"
#include "particle2D.h"
#include "gsd.h"
#include "comn.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

namespace io {

/**
 * @brief Basic class for exporting data.
 *
 * Define the timming to dump data.
 */
class ExporterBase {
public:
  ExporterBase(int n_step, int sep, int start, MPI_Comm group_comm);

  bool need_export(const int i_step) {
    if (log_sep_ < 0) {
      return i_step % sep_ == 0;
    } else {
      return need_export_log_scale(i_step);
    }
  }

  bool need_export_log_scale(const int i_step) {
    bool res = false;
    if (i_step + start_ == frames_[cur_frame_]) {
      cur_frame_++;
      res = true; 
    }
    return res;
  }


  void set_log_scale_frames(double h, double log_sep=0.1);

protected:
  int n_step_;    // total steps to run
  int sep_;
  int start_ = 0; // The first step 
  MPI_Comm comm_;
  int my_rank_ = 0;
  int tot_proc_ = 1;
  std::vector<int> frames_;
  int cur_frame_ = 0;
  double log_sep_ = -1;
};

/**
 * @brief Exporter to output log
 *
 * Output the parameters after the initialization.
 * Output the beginning and endding time of the simulation.
 * Record time every certain time steps.
 */
class LogExporter : public ExporterBase {
public:
  LogExporter(const std::string& outfile,
              int n_step, int sep, int start,
              int np_gl, MPI_Comm group_comm);

  ~LogExporter();

  void record(int i_step);

  std::ofstream fout;
private:
  std::chrono::time_point<std::chrono::system_clock> t_start_;
  int n_par_;
  int step_count_ = 0;
};


class Snap_GSD_2 : public ExporterBase {
public:
  Snap_GSD_2(const std::string& filename,
             int n_step, int sep, int &start,
             double h, double log_sep,
             const Vec_2<double>& gl_l,
             size_t n_par_gl,
             const std::string& open_flag,
             MPI_Comm group_comm);

  ~Snap_GSD_2();

  template <typename TPar>
  void get_data_from_par(const std::vector<TPar>& p_arr, float* pos, float* psi, float* omega);

  uint64_t get_time_step();

  int reset_start_time_step();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);

  template <typename TFloat>
  void read(int i_frame, TFloat* x, TFloat* y, TFloat* theta, TFloat* psi, TFloat* omega);

  template <typename TFloat>
  void read_last_frame(TFloat* x, TFloat* y, TFloat* theta, TFloat* psi, TFloat* omega);

  void ini(const std::string& filename, const std::string& open_flag, const Vec_2<double>& gl_l);

private:
  gsd_handle* handle_ = nullptr;
  Vec_2<double> gl_l_;
  size_t gl_np_;
};



class OrderParaExporter : public ExporterBase {
public:
  OrderParaExporter(const std::string& outfile, int start, int n_step, int sep, int gl_np, MPI_Comm group_comm);

  ~OrderParaExporter() {if (my_rank_ == 0) fout_.close();};

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& birds);

private:
  std::ofstream fout_;
  int gl_np_;
};

template <typename TPar>
void OrderParaExporter::dump(int i_step, const std::vector<TPar>& birds) {
  if (need_export(i_step)) {
    double gl_u[2] = {0};
    double my_u[2] = {0};
    int sum_n = 0;
    int my_n = 0;

    for (auto& p : birds) {
      double angle = p.get_psi();
      my_u[0] += cos(angle);
      my_u[1] += sin(angle);
      my_n++;
    }

    MPI_Reduce(my_u, gl_u, 2, MPI_DOUBLE, MPI_SUM, 0, comm_);
    MPI_Reduce(&my_n, &sum_n, 1, MPI_INT, MPI_SUM, 0, comm_);

    if (my_rank_ == 0) {
      if (sum_n != gl_np_) {
        std::cout << "Error, sum_n=" << sum_n << ", while gl_np_=" << gl_np_ << std::endl;
        exit(1);
      }
      double ux = gl_u[0] / gl_np_;
      double uy = gl_u[1] / gl_np_;
      double phi = sqrt(ux * ux + uy * uy);
      double theta = atan2(uy, ux);

      fout_ << i_step << "\t" << std::setprecision(8)
            << phi << "\t" << theta << "\n";
    }  
  }
}


template <typename TPar>
void Snap_GSD_2::get_data_from_par(const std::vector<TPar>& p_arr,
                                   float* pos,
                                   float* psi,
                                   float* omega) {
  int n_par = p_arr.size();
  int* n_par_arr = new int[tot_proc_]{};
  int* displs = new int[tot_proc_]{};
  MPI_Gather(&n_par, 1, MPI_INT, n_par_arr, 1, MPI_INT, 0, comm_);
  if (my_rank_ == 0) {
    size_t count = 0;
    for (int i = 0; i < tot_proc_; i++) {
      displs[i] += count;
      count += n_par_arr[i];
    }
    if (count != gl_np_) {
      std::cout << "Error, total particles = " << count << ", unequal to "
                << gl_np_ << std::endl;
      exit(1);
    }
  }

  float* my_pos = new float[n_par * 3];
  float* my_psi = new float[n_par];
  float* my_omega = new float[n_par];

  double half_Lx = gl_l_.x * 0.5;
  double half_Ly = gl_l_.y * 0.5;
  for (int j = 0; j < n_par; j++) {
    int j3 = j * 3;
    my_pos[j3    ] = p_arr[j].pos.x - half_Lx;
    my_pos[j3 + 1] = p_arr[j].pos.y - half_Ly;
    my_pos[j3 + 2] = p_arr[j].get_theta();
    my_psi[j] = p_arr[j].get_psi();
    my_omega[j] = p_arr[j].omega;
  }

  MPI_Gatherv(my_psi, n_par, MPI_FLOAT,
              psi, n_par_arr, displs, MPI_FLOAT,
              0, comm_);
  MPI_Gatherv(my_omega, n_par, MPI_FLOAT,
              omega, n_par_arr, displs, MPI_FLOAT,
              0, comm_);

  for (int i = 0; i < tot_proc_; i++) {
    n_par_arr[i] *= 3;
    displs[i] *= 3;
  }

  MPI_Gatherv(my_pos, n_par * 3, MPI_FLOAT,
              pos, n_par_arr, displs, MPI_FLOAT,
              0, comm_);

  delete[] n_par_arr;
  delete[] displs;
  delete[] my_pos;
  delete[] my_psi;
  delete[] my_omega;
}

template<typename TPar>
void Snap_GSD_2::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    float* pos = nullptr;
    float* psi = nullptr;
    float* omega = nullptr;
    
    if (my_rank_ == 0) {
      pos = new float[gl_np_ * 3];
      psi = new float[gl_np_];
      omega = new float[gl_np_];
    }

    get_data_from_par(p_arr, pos, psi, omega);
    uint32_t n_par = gl_np_;
    if (my_rank_ == 0) {
      //uint64_t step = get_time_step();
      uint64_t step = start_ + i_step;
      int nframes = gsd_get_nframes(handle_);

      std::cout << "dump frame " << nframes << " at time step " << step << std::endl;
      gsd_write_chunk(handle_, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &step);
      gsd_write_chunk(handle_, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &n_par);
      gsd_write_chunk(handle_, "particles/position", GSD_TYPE_FLOAT, n_par, 3, 0, pos);
      gsd_write_chunk(handle_, "particles/charge", GSD_TYPE_FLOAT, n_par, 1, 0, psi);
      gsd_write_chunk(handle_, "particles/mass", GSD_TYPE_FLOAT, n_par, 1, 0, omega);
      gsd_end_frame(handle_);
    }
    delete[] pos;
    delete[] psi;
    delete[] omega;
  }
}


template<typename TFloat>
void io::Snap_GSD_2::read(int i_frame,
                          TFloat* x, TFloat* y, TFloat* theta,
                          TFloat* psi, TFloat* omega) {
  if (my_rank_ == 0) {
    uint32_t n_par;
    const gsd_index_entry* chunk = gsd_find_chunk(handle_, i_frame, "particles/N");
    gsd_read_chunk(handle_, &n_par, chunk);
    if (n_par == gl_np_) {
      std::cout << "frame " << i_frame << ": find " << n_par << " particles" << std::endl;
    } else {
      std::cout << "Error for loading gsd file, wrong partiles number" << std::endl;
      exit(1);
    }
    float* pos = new float[size_t(n_par) * 3];
    chunk = gsd_find_chunk(handle_, i_frame, "particles/position");
    gsd_read_chunk(handle_, pos, chunk);

    float* psi_in = new float[n_par];
    chunk = gsd_find_chunk(handle_, i_frame, "particles/charge");
    gsd_read_chunk(handle_, psi_in, chunk);

    float* omega_in = new float[n_par];
    chunk = gsd_find_chunk(handle_, i_frame, "particles/mass");
    gsd_read_chunk(handle_, omega_in, chunk);

    double half_Lx = gl_l_.x * 0.5;
    double half_Ly = gl_l_.y * 0.5;
    for (int i = 0; i < n_par; i++) {
      x[i] = pos[i * 3 + 0] + half_Lx;
      y[i] = pos[i * 3 + 1] + half_Ly;
      theta[i] = pos[i * 3 + 2];
      tangle_1D(x[i], gl_l_.x);
      tangle_1D(y[i], gl_l_.y);

      psi[i] = psi_in[i];
      omega[i] = omega_in[i];
    }
    delete[] pos;
    delete[] psi_in;
    delete[] omega_in;
  }
}


template<typename TFloat>
void io::Snap_GSD_2::read_last_frame(TFloat* x, TFloat* y, 
                                     TFloat* theta, TFloat* psi,
                                     TFloat* omega) {
  if (my_rank_ == 0) {
    int nframes = gsd_get_nframes(handle_);
    if (nframes < 1) {
      std::cout << "Error, nframes=" << nframes << std::endl;
      exit(1);
    } else {
      read(nframes - 1, x, y, theta, psi, omega);
    }
  }
}


}

