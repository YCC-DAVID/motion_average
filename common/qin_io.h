#ifndef QIN_IO_H
#define QIN_IO_H

#include <Eigen/Core>
#include <iostream>
#include <string>

struct QinPose {
  std::string imgname;
  double f;
  double cx;
  double cy;
  int pxsz_x = 1;
  int pxsz_y = 1;
  int width;
  int height;
  double x, y, z;
  double omega, phi, kappa;

  Eigen::Matrix3d GetK() const {
    Eigen::Matrix3d K;
    K << f / pxsz_x, 0, double(width - 1) / 2.0 + cx / pxsz_x, 0, f / pxsz_y, double(height - 1) / 2.0 - cy / pxsz_y, 0, 0, 1;
    return K;
  }

  void SetK(const Eigen::Matrix3d &K) {
      // Qin file use x-> right y->up convension
    f = K(0, 0) * pxsz_x; // = K(1,1) * pxsz_y
    cx = (K(0, 2) - double(width-1) / 2.0) * pxsz_x;
    cy = (double(height-1) / 2.0 - K(1, 2)) * pxsz_y;
  }

  /**
   Implicitly convert Photogrammetry CS to ComputerVision CS
    [ cos(kappa)*cos(phi), cos(omega)*sin(kappa) + cos(kappa)*sin(omega)*sin(phi),   sin(kappa)*sin(omega) - cos(kappa)*cos(omega)*sin(phi)]
    [ cos(phi)*sin(kappa), sin(kappa)*sin(omega)*sin(phi) - cos(kappa)*cos(omega), - cos(kappa)*sin(omega) - cos(omega)*sin(kappa)*sin(phi)]
    [           -sin(phi),                                    cos(phi)*sin(omega),                                     -cos(omega)*cos(phi)]
   * @return
   */
  Eigen::Matrix3d GetR() const {
    Eigen::Matrix3d Rcv;
    Rcv << cos(kappa) * cos(phi), cos(omega) * sin(kappa) + cos(kappa) * sin(omega) * sin(phi), sin(kappa) * sin(omega) - cos(kappa) * cos(omega) * sin(phi), cos(phi) * sin(kappa),
        sin(kappa) * sin(omega) * sin(phi) - cos(kappa) * cos(omega), -cos(kappa) * sin(omega) - cos(omega) * sin(kappa) * sin(phi), -sin(phi), cos(phi) * sin(omega), -cos(omega) * cos(phi);
    return Rcv;
  }

  void SetR(const Eigen::Matrix3d &Rph) {
    double _phi, _omega, _kappa;
    _phi = std::asin(-Rph(2, 0));
    _omega = std::atan2(Rph(2, 1), -Rph(2, 2));
    _kappa = std::atan2(Rph(1, 0), Rph(0, 0));
    /// TODO: Singular handling while cos(phi) ~= 0
    phi = _phi;
    omega = _omega;
    kappa = _kappa;
  }

  Eigen::Vector3d GetC() const {
    Eigen::Vector3d C(x, y, z);
    return C;
  }

  void SetC(const Eigen::Vector3d &c) {
    x = c[0];
    y = c[1];
    z = c[2];
  }
};
std::istream &operator>>(std::istream &is, QinPose &qin) {
  is >> qin.imgname >> qin.f >> qin.cx >> qin.cy >> qin.pxsz_x >> qin.pxsz_y >> qin.width >> qin.height >> qin.x >> qin.y >> qin.z >> qin.omega >> qin.phi >> qin.kappa;
  return is;
}

std::ostream &operator<<(std::ostream &os, const QinPose &qin) {
  os << qin.imgname << " " << qin.f << " " << qin.cx << " " << qin.cy << " " << qin.pxsz_x << " " << qin.pxsz_y << " "
     << qin.width << " " << qin.height << " " << qin.x << " " << qin.y << " " << qin.z << " "
     << qin.omega<< " " << qin.phi << " " << qin.kappa << "\n";
  return os;
}

#endif  // QIN_IO_H