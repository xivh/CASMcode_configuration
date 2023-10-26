#ifndef CASM_irreps_misc
#define CASM_irreps_misc

#include "casm/casm_io/Log.hh"
#include "casm/configuration/irreps/definitions.hh"

namespace CASM {
namespace irreps {

/// \brief Round entries that are within tol of being integer to that integer
/// value
inline Eigen::MatrixXcd prettyc(const Eigen::MatrixXcd &M, double tol = 1e-10) {
  Eigen::MatrixXcd Mp(M);
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < M.cols(); j++) {
      if (std::abs(std::round(M(i, j).real()) - M(i, j).real()) < tol) {
        Mp(i, j).real(std::round(M(i, j).real()));
      }
      if (std::abs(std::round(M(i, j).imag()) - M(i, j).imag()) < tol) {
        Mp(i, j).imag(std::round(M(i, j).imag()));
      }
    }
  }
  return Mp;
}

/// \brief Round entries that are within tol of being integer to that integer
/// value
inline Eigen::MatrixXd pretty(const Eigen::MatrixXd &M, double tol = 1e-10) {
  Eigen::MatrixXd Mp(M);
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < M.cols(); j++) {
      if (std::abs(std::round(M(i, j)) - M(i, j)) < tol) {
        Mp(i, j) = std::round(M(i, j));
      }
    }
  }
  return Mp;
}

inline void prettyp(std::ostream &sout, const Eigen::MatrixXcd &M,
                    double tol = 1e-10) {
  sout << "shape=(" << M.rows() << ", " << M.cols() << ")" << std::endl;
  sout << prettyc(M, tol) << std::endl;
}

inline void prettyp(std::ostream &sout, const Eigen::MatrixXd &M,
                    double tol = 1e-10) {
  sout << "shape=(" << M.rows() << ", " << M.cols() << ")" << std::endl;
  sout << pretty(M, tol) << std::endl;
}

template <int _required_verbosity = Log::standard>
void prettypc(Log &log, std::string what, Eigen::MatrixXcd const &M,
              double tol = 1e-10) {
  log.begin_section<_required_verbosity>();
  if (log.print()) {
    log.indent() << what << ": shape=(" << M.rows() << ", " << M.cols() << ")"
                 << std::endl;
    std::stringstream ss;
    ss << pretty(M, tol);
    log.verbatim(ss.str(), false);
    log.indent() << std::endl;
  }
  log.end_section();
}

template <int _required_verbosity = Log::standard>
void prettyp(Log &log, std::string what, Eigen::MatrixXd const &M,
             double tol = 1e-10) {
  log.begin_section<_required_verbosity>();
  if (log.print()) {
    log.indent() << what << ": shape=(" << M.rows() << ", " << M.cols() << ")"
                 << std::endl;
    std::stringstream ss;
    ss << pretty(M, tol);
    log.verbatim(ss.str(), false);
    log.indent() << std::endl;
  }
  log.end_section();
}

}  // namespace irreps
}  // namespace CASM

#endif
