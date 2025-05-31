#pragma once
#include <vector>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "config.h"
#include "particle2D.h"
#include "gsd.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

namespace exporter {

#ifdef _MSC_VER
const std::string delimiter("\\");
#else
const std::string delimiter("/");
#endif

/**
 * @brief Basic class for exporting data.
 *
 * Define the timming to dump data.
 */
class ExporterBase {
public:
  ExporterBase(int n_step, int sep, int start) : n_step_(n_step), sep_(sep), start_(start) {}

  bool need_export(const int i_step);

  void set_log_scale_frames(double h, double log_sep = 0.1);

protected:
  int n_step_;    // total steps to run
  int sep_;
  int start_ = 0; // The first step 
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
  LogExporter(const std::string& outfile, int start, int n_step, int sep, int np);

  ~LogExporter();

  void record(int i_step);

  std::ofstream fout;
private:
  std::chrono::time_point<std::chrono::system_clock> t_start_;
  int n_par_;
  int step_count_ = 0;
};

class OrderParaExporter : public ExporterBase {
public:
  OrderParaExporter(const std::string& outfile, int start, int n_step, int sep)
    : ExporterBase(n_step, sep, start), fout_(outfile) {}

  ~OrderParaExporter() { fout_.close(); }

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& birds);

private:
  std::ofstream fout_;
};

template <typename TPar>
void OrderParaExporter::dump(int i_step, const std::vector<TPar>& birds) {
  if (need_export(i_step)) {
    double phi, theta;
    
    double ux = 0;
    double uy = 0;

    for (auto& p : birds) {
      ux += cos(p.psi);
      uy += sin(p.psi);
    }

    ux /= birds.size();
    uy /= birds.size();

    phi = sqrt(ux * ux + uy * uy);
    theta = atan2(uy, ux);

    fout_ << i_step << "\t" << std::setprecision(8)
          << phi << "\t" << theta << "\n";
  }
}


class Snap_GSD_2 : public ExporterBase {
public:
  Snap_GSD_2(const std::string& filename,
             int n_step, int sep, int & start,
             double h, double log_sep,
             const Vec_2<double>& gl_l,
             const std::string& open_flag);

  ~Snap_GSD_2();

  template <typename TPar>
  void get_data_from_par(const std::vector<TPar>& p_arr, float* pos);

  template <typename TPar>
  void get_data_from_par(const std::vector<TPar>& p_arr, float* pos, float* omega);

  uint64_t get_time_step();

  int reset_start_time_step();

  template <typename TPar>
  void dump(int i_step, const std::vector<TPar>& p_arr);

  template <typename TPar>
  void read(int i_frame, std::vector<TPar>& p_arr);

  template <typename TPar>
  void read_last_frame(std::vector<TPar>& p_arr);

private:
  gsd_handle* handle_ = nullptr;
  Vec_2<double> half_l_;
};


template <typename TPar>
void Snap_GSD_2::get_data_from_par(const std::vector<TPar>& p_arr, float* pos) {
  size_t n_par = p_arr.size();
  for (size_t j = 0; j < n_par; j++) {
    size_t j3 = j * 3;
    pos[j3    ] = p_arr[j].pos.x - half_l_.x;
    pos[j3 + 1] = p_arr[j].pos.y - half_l_.y;
    pos[j3 + 2] = p_arr[j].get_psi();
  }
}


template<typename TPar>
void exporter::Snap_GSD_2::get_data_from_par(const std::vector<TPar>& p_arr, float* pos, float* omega) {
  size_t n_par = p_arr.size();
  for (size_t j = 0; j < n_par; j++) {
    size_t j3 = j * 3;
    pos[j3] = p_arr[j].pos.x - half_l_.x;
    pos[j3 + 1] = p_arr[j].pos.y - half_l_.y;
    pos[j3 + 2] = p_arr[j].get_psi();
    omega[j] = p_arr[j].omega;
  }
}


template <typename TPar>
void Snap_GSD_2::dump(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    uint32_t n_par = p_arr.size();
    float* pos = new float[n_par * 3];
    float* omega = new float[n_par];
    get_data_from_par(p_arr, pos, omega);
    //uint64_t step = get_time_step();
    uint64_t step = start_ + i_step;

    int nframes = gsd_get_nframes(handle_);

    std::cout << "dump frame " << nframes << " at time step " << step << std::endl;
    gsd_write_chunk(handle_, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &step);
    gsd_write_chunk(handle_, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &n_par);
    gsd_write_chunk(handle_, "particles/position", GSD_TYPE_FLOAT, n_par, 3, 0, pos);
    gsd_write_chunk(handle_, "particles/mass", GSD_TYPE_FLOAT, n_par, 1, 0, omega);

    gsd_end_frame(handle_);
    delete[] pos;
    delete[] omega;
  }
}

template <typename TPar>
void Snap_GSD_2::read(int i_frame, std::vector<TPar>& p_arr) {
  uint32_t n_par;
  const gsd_index_entry* chunk = gsd_find_chunk(handle_, i_frame, "particles/N");
  gsd_read_chunk(handle_, &n_par, chunk);
  std::cout << "frame " << i_frame  <<": find " << n_par << " particles" << std::endl;

  float* pos = new float[n_par * 3];
  chunk = gsd_find_chunk(handle_, i_frame, "particles/position");
  gsd_read_chunk(handle_, pos, chunk);

  float* omega = new float[n_par];
  chunk = gsd_find_chunk(handle_, 0, "particles/mass");
  gsd_read_chunk(handle_, omega, chunk);

  p_arr.reserve(n_par);

  double Lx = half_l_.x * 2;
  double Ly = half_l_.y * 2;

  for (int j = 0; j < n_par; j++) {
    TPar p;
    size_t j3 = j * 3;
    double x = pos[j3] + half_l_.x;
    double y = pos[j3 + 1] + half_l_.y;
    double psi = pos[j3 + 2];

    
    if (x < 0) {
      x += Lx;
    } else if (x >= Lx) {
      x -= Lx;
    }
    if (y < 0) {
      y += Ly;
    } else if (y >= Ly) {
      y -= Ly;
    }

    p.pos = Vec_2<double>(x, y);
    p.psi = psi;
    p.omega = omega[j];
    p_arr.push_back(p);
  }

  delete[] pos;
  delete[] omega;
}


template <typename TPar>
void Snap_GSD_2::read_last_frame(std::vector<TPar>& p_arr) {
  int nframes = gsd_get_nframes(handle_);
  if (nframes < 1) {
    std::cout << "Error, nframes=" << nframes << std::endl;
    exit(1);
  } else {
    read(nframes-1, p_arr);
  }
}

}

