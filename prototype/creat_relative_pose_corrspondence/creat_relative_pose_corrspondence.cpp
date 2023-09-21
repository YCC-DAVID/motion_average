#include <spdlog/spdlog.h>
#include <Eigen/Dense>
#include <cxxopts.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <Eigen/Core>
#include<opencv2/opencv.hpp>
#include<opencv2/core/eigen.hpp>
#include "XYZGrid.h"
#include "gdal.h"
#include "mystructs.h"
#include "myutils.h"
#include "qin_io.h"
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>


template <typename T>
class StatVar {
  std::string _name;
  T _E;
  T _sqE;
  uint32_t _cnt;

 public:
  StatVar(std::string name) : _name(name) { _cnt = 0; }
  void measure(T _v) {
    if (_cnt == 0) {
      _E = _v;
      _sqE = _v.cwiseProduct(_v);
    } else {
      _E += (_v - _E) / (_cnt + 1);
      _sqE += (_v.cwiseProduct(_v) - _sqE) / (_cnt + 1);
    }
    _cnt += 1;
  }

  T getE() const { return _E; }
  T getVar() const { return _sqE - _E.cwiseProduct(_E); }
  uint32_t getCount() const { return _cnt; }

  void print(std::ostream& out) { out << "E(" << _name << ")=" << getE().transpose() << ", D(" << _name << ")=" << getVar().transpose() << " of " << getCount() << " measurements" << std::endl; }

  void print(spdlog::logger& logger) {
    std::stringstream ss;
    ss << "E(" << _name << ")=" << getE().transpose() << ", D(" << _name << ")=" << getVar().transpose() << " of " << getCount() << " measurements";
    logger.info(ss.str());
  }
};


std::vector<QinPose>readQinv2File(std::string input_qin) {
  std::ifstream ifs(input_qin);
  size_t num_pose;
  ifs >> num_pose;
  std::vector<QinPose> poses(num_pose);
  for (size_t i = 0; i < num_pose; ++i) ifs >> poses[i];
  return poses;
}

namespace fs = std::filesystem;
namespace logg = spdlog;
using namespace std;
using namespace motionavg::TranslateND;

bool read_correspondence_pts_pix_binary_num(const char* filename, unsigned int& numpts) {
  std::fstream fs(filename, std::ios::in | std::ios::binary);
  if (!fs.is_open()) {
    std::cout << "file opening error" << std::endl;
    return false;
  }
  fs.read((char*)&numpts, sizeof(unsigned int));
  fs.close();
}
//
bool read_correspondence_pts_pix_binary(const char* filename, float* datapointerout) {
  std::fstream fs(filename, std::ios::in | std::ios::binary);
  if (!fs.is_open()) {
    std::cout << "file opening error" << std::endl;
    return false;
  }

  unsigned int numpt;
  fs.read((char*)&numpt, sizeof(unsigned int));
  fs.read((char*)datapointerout, sizeof(float) * numpt * 7);
  fs.close();

  return true;
}

cxxopts::Options parseOptions(std::string exepath = "") {
  std::string exename = fs::path(exepath).filename().string();
  cxxopts::Options options(exename, "create_graph from folder");
  // clang-format off
  options.add_options()
    ("i,input_folder", "input folder", cxxopts::value<std::string>())
	("qinv2file", "qinv2 file", cxxopts::value<std::string>())
    ("gcpfile", "gcp json file", cxxopts::value<std::string>()->default_value(""))
	("o,output_name", "output filename", cxxopts::value<std::string>()->default_value("graph.txt"))
	("minrays", "minimum num of rays", cxxopts::value<int>()->default_value("3"))
	("minmatches", "minimum num of matches for valid pair", cxxopts::value<int>()->default_value("100"))
    ("sigma", "sigma of gaussian kernel for confidence [meter]", cxxopts::value<float>()->default_value("1"))
	("h,help", "Print help")
	("v,verbose", "verbose level (trace - 0, debug - 1, info - 2, warn - 3, error - 4, critical - 5, off - 6)", cxxopts::value<int>()->default_value("2"))
	;
  // clang-format on

  options.parse_positional({"input_folder", "qinv2file", "output_name"});
  options.positional_help("input_folder qinv2file output_name [options]");
  return options;
}
float INVALID_VALUE = std::numeric_limits<float>::infinity();

StatVar<Eigen::Vector2d> total_stat_reprojection_error("total reproj");


bool check_pose(cv::Mat R1_cand, cv::Mat R2_cand, cv::Mat t_cand, std::vector<cv::Point2f> inlier1, std::vector<cv::Point2f> inlier2, cv::Mat* R_res, cv::Mat* t_res) {
  cv::Mat R1(R1_cand), R2(R2_cand), t(t_cand), neg_t;  // 2 rotation and 1 translation get from essential matrix
  neg_t = -t;
  cv::Mat R1_t(R1.t()), R2_t(R2.t());
  cv::Mat P1 = cv::Mat::eye(3, 4, CV_64F);  // projection matrix of first camera
  cv::Mat P2(3, 4, CV_64F);                 // projection matrix of second camera
  cv::Mat CurrentR,CurrentT;

  int bestSolution = -1;
  int maxInFront = 0;
  int sample_rate=100;
  cv::Mat Mat_inlier1(2, inlier1.size(), CV_64F);
  cv::Mat Mat_inlier2(2, inlier2.size(), CV_64F);
  for (size_t i = 0; i < inlier1.size() / sample_rate; ++i) {
    Mat_inlier1.at<double>(0, i) = inlier1[i * sample_rate].x;
    Mat_inlier1.at<double>(1, i) = inlier1[i * sample_rate].y;

    Mat_inlier2.at<double>(0, i) = inlier2[i * sample_rate].x;
    Mat_inlier2.at<double>(1, i) = inlier2[i * sample_rate].y;
  }

  // check the possible pose
  for (int solution = 0; solution < 4; solution++) {
    if (solution == 0) {
      R1.copyTo(P2(cv::Rect(0, 0, 3, 3)));
      t.copyTo(P2.col(3));
      CurrentR = R1;
      CurrentT = t;
    } else if (solution == 1) {
      R1.copyTo(P2(cv::Rect(0, 0, 3, 3)));
      neg_t.copyTo(P2.col(3));
      CurrentR = R1;
      CurrentT = neg_t;
    } else if (solution == 2) {
      R2.copyTo(P2(cv::Rect(0, 0, 3, 3)));
      t.copyTo(P2.col(3));
      CurrentR = R2;
      CurrentT = t;
    /*} else if (solution == 3) {
      R1_t.copyTo(P2(cv::Rect(0, 0, 3, 3)));
      t.copyTo(P2.col(3));
      CurrentR = R1_t;
      CurrentT = t;
    } else if (solution == 4) {
      R2_t.copyTo(P2(cv::Rect(0, 0, 3, 3)));
      t.copyTo(P2.col(3));
      CurrentR = R2_t;
      CurrentT = t;*/
    } else {
      R2.copyTo(P2(cv::Rect(0, 0, 3, 3)));
      neg_t.copyTo(P2.col(3));
      CurrentR = R2;
      CurrentT = neg_t;
    }
    cv::Mat point3D(1, 4, CV_64F);
    int forward_num = 0;
    int point1_positive = 0;
    int point2_positive = 0;

    cv::triangulatePoints(P1, P2, Mat_inlier2, Mat_inlier1, point3D);

    // homography to 3d
    for (int j = 0; j < inlier1.size() / sample_rate; j++) {
      cv::Mat point3DNormalized = point3D.col(j) / point3D.col(j).at<double>(3);

      // identify if the point before the image
     
      cv::Mat point3DNormalized2 = CurrentR * point3DNormalized.rowRange(0, 3) + CurrentT;
      if (point3DNormalized.at<double>(2) > 0) point1_positive++;
      if (point3DNormalized2.at<double>(2) > 0) point2_positive++;
      if (point3DNormalized.at<double>(2) > 0 && point3DNormalized2.at<double>(2) > 0) forward_num++;
    
    }
    if (forward_num > maxInFront) {
        maxInFront = forward_num;
        bestSolution = solution;
        //break;
    }
  }

  if (bestSolution == -1) {
    std::cerr << "No valid solution found!" << std::endl;
    return false;
  } else {
    std::cout << "Valid solution is: " << bestSolution << std::endl;
  }
  if (bestSolution == 0) {
    R1.copyTo(*R_res);
    t.copyTo(*t_res);
  } else if (bestSolution == 1) {
    R1.copyTo(*R_res);
    neg_t.copyTo(*t_res);
  } else if (bestSolution == 2) {
    R2.copyTo(*R_res);
    t.copyTo(*t_res);
  //} else if (bestSolution == 3) {
  //  R1_t.copyTo(*R_res);
  //  t.copyTo(*t_res);
  //} else if (bestSolution == 4) {
  //  R2_t.copyTo(*R_res);
  //  t.copyTo(*t_res);
  } else {
    R2.copyTo(*R_res);
    neg_t.copyTo(*t_res);
  }
  return true;
}
bool check_pose_temp(cv::Mat R1_cand, cv::Mat R2_cand, cv::Mat t_cand, cv::Mat ori_R, cv::Mat ori_t,cv::Mat * R_final, cv::Mat* t_final) { 
  cv::Mat R1(R1_cand), R2(R2_cand), t(t_cand), neg_t(-t_cand),orig_R(ori_R),orig_t(ori_t);
  cv::Mat  Rot_final;
  cv::Scalar Rotdif;
  Rotdif = cv::sum(cv::abs(orig_R - R1));
  cout << "rot diff" << cv::sum(cv::abs(orig_R - R1)) << endl;
  //Rotdif2 = cv::sum(cv::sum(cv::abs(ori_R - R2_cand)));
  if (Rotdif.val[0] < 1) {
    R1.copyTo(*R_final);
  } else {
    R2.copyTo(*R_final);
  }
  Rot_final = *R_final;
 cv::Scalar transdiff;
  transdiff = cv::sum(cv::abs(orig_t + Rot_final.t() * t));
 cout << "trans diff" << transdiff << endl;
  if (transdiff.val[0] < 1) {
   t.copyTo(*t_final);
 } else {
   neg_t.copyTo(*t_final);
 }
 return true;
}

double pointToLineDistance(const cv::Mat& line, const cv::Point2f& point) {
  double uppat = std::abs(line.at<double>(0, 0) * point.x + line.at<double>(1, 0) * point.y + line.at<double>(2, 0));
  double lowpat = std::sqrt(line.at<double>(0, 0) * line.at<double>(0, 0) + line.at<double>(1,0) * line.at<double>(1, 0));
  double dis = uppat / lowpat;
  return dis;
        
}

bool create_edge_from_corrspondence_pts(const QinPose& srcPose, const QinPose& tgtPose, const float* matches, int numpts, float* offset_n_weight,
                                        float invsigmasq_2, int minimum_matches, double* xfm, double* cov, spdlog::logger* logger = nullptr) {
  // calcualting the reproject error
  auto srcK = srcPose.GetK();
  auto srcR = srcPose.GetR();
  auto srcC = srcPose.GetC();
  auto tgtK = tgtPose.GetK();
  auto tgtR = tgtPose.GetR();
  auto tgtC = tgtPose.GetC();  // from Qinpose, the pose in world coordinate

  Eigen::Vector3d init_bij = -srcC + tgtC;
  Eigen::Vector3d expected_bij{0, 0, 0};
  int num_measure = 0;
  float totalWeight = 0;

  StatVar<Eigen::Vector2d> stat_src2tgt_error("src2tgt_reproj");
  StatVar<Eigen::Vector2d> stat_tgt2tgt_error("tgt2tgt_reproj");
  StatVar<Eigen::Vector2d> stat_cortgt_error("cortgt_reproj");
  StatVar<Eigen::Vector3d> stat_bij("bij");

  // ceres::Problem problem;
  // ceres::Problem problem_o;
  double posi[3], posi_o[3];
  double Quan[4];
  double Rotation[9];
  double Rot[4];

  Eigen::Matrix3d src2tgtR = tgtR * srcR.transpose();  // rotation from src camera to tgt camera
  Eigen::Vector3d srcTrans = srcR * init_bij;          // translation in src coordinate

  std::vector<cv::Point3f> Matchpoints3d;
  std::vector<cv::Point2f> Matchpointsimg_beg;
  std::vector<cv::Point2f> Matchpointsimg_end;
  std::vector<cv::Point2f> cvReprj_px, ceReprj_px;

  for (int i = 0; i < 3; ++i) {
    posi[i] = srcTrans[i];
    posi_o[i] = posi[i];
  }
  
  Eigen::Quaterniond quan(src2tgtR);
  quan.normalize();
  

  cv::Mat vec_rota;
  cv::eigen2cv(src2tgtR, vec_rota);
  cv::Mat rvec;
  cv::Rodrigues(vec_rota, rvec);
  cv::Mat tvec(3, 1, CV_64F, posi_o);

  for (int k = 0; k < numpts; k++)  // k-th matched points
  {
    std::fill_n(offset_n_weight + 4 * k, 3, INVALID_VALUE);
    offset_n_weight[4 * k + 3] = 0;
    const float& srcX = matches[7 * k + 0];
    const float& srcY = matches[7 * k + 1];
    const float& tgtX = matches[7 * k + 2];
    const float& tgtY = matches[7 * k + 3];  // just use the pixel coordinate not care 3d coordinate

    Matchpointsimg_beg.push_back(cv::Point2f(srcX, srcY));
    Matchpointsimg_end.push_back(cv::Point2f(tgtX, tgtY));

    num_measure++;
  }

  if (num_measure < minimum_matches) return false;

  cv::Mat IntK(3, 3, CV_64F);
  cv::eigen2cv(tgtK, IntK);

  //std::cout << "camera matrix" << std::endl << IntK << std::endl;
  cv::Mat distCof(1, 5, CV_64F);
  distCof = cv::Mat(0, 0, 0, 0, 0);

  // calculate the fundamental matrix
  cv::Mat inlier;
  //cv::Mat F = cv::findFundamentalMat(Matchpointsimg_beg, Matchpointsimg_end, cv::FM_RANSAC, 0.001, 0.99, inlier);
  cv::Mat Ess = cv::findEssentialMat(Matchpointsimg_beg, Matchpointsimg_end, IntK,cv::RANSAC,0.99,0.001,inlier );
  //Ess = Ess.t();
  //std::cout << "essential Matrix:" << std::endl << Ess << std::endl;
  //std::cout << "Fundamental Matrix:" << std::endl << F << std::endl;

  // find the inlier
  std::vector<cv::Point2f> inliers1, inliers2;
  for (size_t i = 0; i < inlier.rows; ++i) {
    if (inlier.at<uchar>(i)) {
      inliers1.push_back(Matchpointsimg_beg[i]);
      inliers2.push_back(Matchpointsimg_end[i]);
    }
  }
  int inlier_num = inliers1.size();
  double inlier_rt = double((inlier_num+0.0) / num_measure);
  std::cout << "Number of inliers: " << inliers1.size() << " Ratio of inliers: " << inlier_rt << std::endl;

  // calculate the rotation and translation

  //cv::Mat E = IntK.t() * F * IntK;
  cv::Mat cvRot;
  cv::Rodrigues(rvec, cvRot);

  double trans_length = norm(tvec);
  cv::Mat R1, R2, t,t_l;
  cv::decomposeEssentialMat(Ess, R1, R2, t);
  cv::Mat t_uni = t/ norm(t);
  t_l = t_uni * trans_length;
  Eigen::Matrix3d R1_eig;
  cv::cv2eigen(R1,R1_eig);
  Eigen::Quaterniond quat_R1(R1_eig);

  Eigen::Matrix3d R2_eig;
  cv::cv2eigen(R2 ,R2_eig);
  Eigen::Quaterniond quat_R2(R2_eig);

  Eigen::Matrix3d R1_eig_t;
  cv::cv2eigen(R1.t(), R1_eig_t);
  Eigen::Quaterniond quat_R1_t(R1_eig_t);

  Eigen::Matrix3d R2_eig_t;
  cv::cv2eigen(R2.t(), R2_eig_t);
  Eigen::Quaterniond quat_R2_t(R2_eig_t);

  cout << "R1 : [" << quat_R1.w() << "," << quat_R1.x() << "," << quat_R1.y() << "," << quat_R1.z() << "]" << endl;
  cout << "R2 : [" << quat_R2.w() << "," << quat_R2.x() << "," << quat_R2.y() << "," << quat_R2.z() << "]" << endl;
  //cout << "R1 transpose : [" << quat_R1_t.w() << "," << quat_R1_t.x() << "," << quat_R1_t.y() << "," << quat_R1_t.z() << "]" << endl;
  //cout << "R2 transpose: [" << quat_R2_t.w() << "," << quat_R2_t.x() << "," << quat_R2_t.y() << "," << quat_R2_t.z() << "]" << endl;
  cv::Mat R_final, t_final_,t_final;

  // trangluate to test which pose is correct
  /*bool findFlag=check_pose(R1, R2, t_l, inliers1, inliers2, &R_final, &t_final_);
  if (!findFlag) return false;*/
  check_pose_temp(R1, R2, t_l, cvRot, tvec, &R_final, &t_final_);
  t_final = -R_final.t()*t_final_;
  

  double img2_err=0;
  
  for (size_t i = 0; i < inliers1.size(); ++i) {
    cv::Mat pt1 = (cv::Mat_<double>(3, 1) <<inliers1[i].x, inliers1[i].y, 1.0);
    cv::Mat line = IntK.t().inv() * Ess * IntK.inv() * pt1;
    //cv::Mat line = pt1.t() * F.t();
  
    double distance = pointToLineDistance(line, inliers2[i]);

    img2_err = img2_err + distance;
    //std::cout << "Distance from point " << i << " to its epipolar line: " << distance << std::endl;
  }

  Eigen::Matrix3d Final_Rota;
  cv::cv2eigen(R_final, Final_Rota);
  Eigen::Quaterniond quat_res(Final_Rota);
  quat_res.normalize();

 

  cv::Mat t_uni_vect = t_final / norm(t_final);
  cv::Mat rel_Trans = t_uni_vect * trans_length;
  cv::Mat Tvec;
  Tvec = -cvRot.t() * tvec;
  totalWeight = img2_err;

  cv::Mat Rotdif;
  Rotdif = cvRot - R_final;

  Eigen::Vector3d transdiff;
 
  double T_final[3],Rot_final[4];
  for (int i = 0; i < 3; i++) {
    transdiff[i] = t_final.at<double>(i) - posi[i];
    T_final[i] = t_final.at<double>(i);
  }
  Rot_final[0] = quat_res.w();
  Rot_final[1] = quat_res.x();
  Rot_final[2] = quat_res.y();
  Rot_final[3] = quat_res.z();
  

   Eigen::Vector3d Posi(T_final);
  /* totalWeight = summary.final_cost;*/
  cout << "Original trans : ["  << tvec.at<double>(0) << "," << tvec.at<double>(1) << "," << tvec.at<double>(2) << "]"<<endl;
  std::ostringstream oss;
  oss << "Corrspondence trans : [" << Posi.x() << "," << Posi.y() << "," << Posi.z() << "]";
  string trans_info = oss.str();
  logger->info(oss.str());
  oss.str("");
  oss.clear();
  cout << "Original trans : [" << quan.w() << "," << quan.x() << "," << quan.y() << "," << quan.z() << "]" << endl;
  oss << "Corrspondence Rot : [" << quat_res.w() << "," << quat_res.x() << "," << quat_res.y() << "," << quat_res.z() <<"]";
  logger->info(oss.str());
  oss.str("");
  oss.clear();
  cout << "Rot diff(corrspond-original): "
       << "\n"
       << Rotdif << "\n ";
  cout << "\n";
  cout << "sum of trans error:"
       << "\n"
       << transdiff.array().abs().sum() << "\n";
  cout << "sum of Rot error:"
       << "\n"
       << cv::sum(cv::abs(Rotdif)) << "\n";
  oss << "the average epipolar line error is [" << totalWeight / inlier_num << "], with [" << inlier_num << "] corrspondence points"<< '\n';
  logger->info(oss.str());
  oss.str("");
  oss.clear();
 
  bool valid = true;
  /*double diff_rate = (transdiff.norm()) / (Posi.norm());
  if (diff_rate > 0.01 && (transdiff.norm()) > 1) {
    valid = false;
    cout << "the result is not stable with the rate of:" << diff_rate << endl;
  }*/
 

  cout << "\n\n";
  std::vector<double> Pose_trans(T_final, T_final + sizeof(T_final) / sizeof(T_final[0]));
  std::vector<double> Pose_rota(Rot_final, Rot_final + sizeof(Rot_final) / sizeof(Rot_final[0]));
  std::vector<double> Pose;
  Pose.reserve(Pose_trans.size() + Pose_rota.size());
  Pose.insert(Pose.end(), Pose_trans.begin(), Pose_trans.end());
  Pose.insert(Pose.end(), Pose_rota.begin(), Pose_rota.end());

  // totalWeight = 1;
  std::fill_n(cov, 36, 0.);
  for (int _d = 0; _d < 7; ++_d) {
    xfm[_d] = Pose[_d];

    if (std::isnan(xfm[_d]) || std::isinf(xfm[_d])) valid = false;
  }
  for (int _d = 0; _d < 6; ++_d) {
    cov[7 * _d] = 1. / (totalWeight / inlier_num);
  }
  return valid;
}


using GraphT = TranslateGraph3WithGCP;

int main(int argc, char** argv) {
  cout << argv[0] << endl;
  GDALAllRegister();

  cxxopts::Options options = parseOptions(argv[0]);
  cxxopts::ParseResult args = options.parse(argc, argv);

  if (args.count("help") != 0) {
    cout << options.help() << endl;
    return 0;
  }
  if (args.count("input_folder") == 0 || args.count("qinv2file") == 0) {
    cout << options.help() << endl;
    return 0;
  }
  bool hasGCP = args.count("gcpfile") > 0;
  fs::path rootdirpath(args["input_folder"].as<string>());
  fs::path qinv2path(args["qinv2file"].as<string>());
  fs::path gcppath(args["gcpfile"].as<string>());
  fs::path outfilename(args["output_name"].as<string>()), outfilepath;
  int minimum_rays = args["minrays"].as<int>();
  bool preloadgrids = false;
  bool minimum_matches = args["minmatches"].as<int>();
  float sigma = args["sigma"].as<float>();
  float invsigmasq_2 = 2.f / std::pow(sigma, 2.f);
  if (outfilename.is_absolute())
    outfilepath = outfilename;
  else
    outfilepath = rootdirpath / outfilename;

  auto console_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
  console_sink->set_level(spdlog::level::debug);

  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(outfilepath.string() + ".create_graph.log");
  file_sink->set_level(spdlog::level::trace);

  spdlog::logger logger("LOG", {console_sink, file_sink});
  logger.set_level(static_cast<logg::level::level_enum>(args["verbose"].as<int>()));

  if (!fs::exists(qinv2path)) {
    logger.critical("QinPose Not found: {}", qinv2path.string());
    return 1;
  }

  if (hasGCP && !fs::exists(gcppath)) {
    logger.critical("GCP file Not Found: {}", gcppath.string());
    return 1;
  }

  fs::path lasdirpath = rootdirpath / "Point_clouds" / "point_clouds";
  fs::path tempfolderdirpath = rootdirpath / "Point_clouds" / "temp_folder_cluster";
  fs::path Neighborhood_file_path = lasdirpath / "Neighborhood_file.bin";
  if (!fs::exists(Neighborhood_file_path)) {
    logger.critical("File {} not found!", Neighborhood_file_path.string());
    return 1;
  }

  // Build graph;
  // using GraphT = TranslateGraph<3>;
  GraphT g;
  g.basepath = rootdirpath.string();
  g.baseID = -1;

  std::vector<QinPose> poses = readQinv2File(qinv2path.string());

  // Add Nodes
  for (size_t i = 0; i < poses.size(); ++i) {
    GraphT::Node n;
    n.name = poses[i].imgname;
    fs::path tmpfolder = tempfolderdirpath / fs::path(n.name).stem();
    if (!fs::exists(tmpfolder)) {
      logger.error("Tmp folder {} not found!", tmpfolder.string());
      return 1;
    }
    n.path = fs::relative(tmpfolder, rootdirpath).string();
    n.xfm[0] = poses[i].x;
    n.xfm[1] = poses[i].y;
    n.xfm[2] = poses[i].z;
    auto Rot = poses[i].GetR();
    /* double rot[9];
     for (int i = 0; i < 3; ++i) {
       for (int j = 0; j < 3; ++j) {
         rot[i * 3 + j] = Rot(i, j);
       }
     }
     cout << Rot << endl;
     double quat[4];*/
    Eigen::Quaterniond quan(Rot);
    // Eigen::Quaterniond quan(quat[0], quat[1], quat[2], quat[3]);
    quan.normalize();
    n.xfm[3] = quan.w();
    n.xfm[4] = quan.x();
    n.xfm[5] = quan.y();
    n.xfm[6] = quan.z();
    g.insertNode(n);
  }

  if (hasGCP) {
    ifstream ifs;
    // Add GCP
    using json = nlohmann::json;
    ifs.open(gcppath.string());
    json gcpfilejson = json::parse(ifs);
    json gcpjson = gcpfilejson["GCP"];
    for (auto& [gcpname, gcpval] : gcpjson.items()) {
      GCPGraph::GCPPoint gcppt;
      gcppt.name = gcpname;
      gcppt.x = gcpval["x"].get<double>();
      gcppt.y = gcpval["y"].get<double>();
      gcppt.z = gcpval["z"].get<double>();

      gcppt.ex = gcpval.value("ex", 999.f);
      gcppt.ey = gcpval.value("ey", 999.f);
      gcppt.ez = gcpval.value("ez", 999.f);

      int gcpid = g.gcps.size();
      g.gcps.push_back(gcppt);

      for (auto& [key, view] : gcpval["views"].items()) {
        // find viewid
        std::string viewname = key;
        auto nodeit = std::find_if(g.nodes.begin(), g.nodes.end(), [&viewname](const GraphT::Node& n) { return viewname == n.name; });
        int viewid = std::distance(g.nodes.begin(), nodeit);
        GCPGraph::GCPLink gcplink;
        gcplink.gcpid = gcpid;
        gcplink.viewid = viewid;
        gcplink.u = view["u"].get<double>();
        gcplink.v = view["v"].get<double>();
        g.gcplinks.push_back(gcplink);
      }
    }

    // Read GCP constraints
    for (auto& l : g.gcplinks) {
      auto& view = g.nodes[l.viewid];
      auto& gcppt = g.gcps[l.gcpid];

      string viewName = fs::path(view.name).stem().string();

      // clang-format off
      fs::path viewgridpath[3] = {
        tempfolderdirpath / viewName / fmt::format("{}_Xgrid.tif", viewName),
        tempfolderdirpath / viewName / fmt::format("{}_Ygrid.tif", viewName),
        tempfolderdirpath / viewName / fmt::format("{}_Zgrid.tif", viewName)};
      // clang-format on
      XYZGrid viewXYZ(false);
      viewXYZ.open(viewgridpath, fs::path());

      uint16_t dummy;
      viewXYZ.sample(l.u, l.v, l.dx, l.dy, l.dz, dummy);
      // TODO: The erro could be derived from grids and correspondences.
      l.ex = 0.01;
      l.ey = 0.01;
      l.ez = 0.01;
      double viewX = l.dx + view.xfm[0];
      double viewY = l.dy + view.xfm[1];
      double viewZ = l.dz + view.xfm[2];

      logger.trace("GCP[{}]-View[{}]: {} {} {}", gcppt.name, view.name, viewX - gcppt.x, viewY - gcppt.y, viewZ - gcppt.z);
    }
  }

  // Add Edges
  int numimg, maxpair;
  std::vector<int> neighborbin;
  {
    ifstream ifs;
    ifs.open(Neighborhood_file_path.string(), ios::binary);
    ifs.read(reinterpret_cast<char*>(&numimg), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&maxpair), sizeof(int));
    logger.debug("neighbor bin: {} {}", numimg, maxpair);
    neighborbin.resize(numimg * (maxpair + 1));
    ifs.read(reinterpret_cast<char*>(neighborbin.data()), numimg * (maxpair + 1) * sizeof(int));
    ifs.close();
  }

  std::vector<std::vector<int>> decoded_neighbors(numimg);
  for (int i = 0; i < numimg; ++i)  // i-th image
  {
    const int* ids = &neighborbin[i * (maxpair + 1)];
    for (const int* _pj = ids + 1; _pj < ids + maxpair + 1; ++_pj) {
      if (*_pj == -1) break;
      decoded_neighbors[i].push_back(*_pj);
    }
  }

  for (int i = 0; i < numimg; ++i)  // i-th image
  {
    const int* ids = &neighborbin[i * (maxpair + 1)];
    int srcId = i;
    string srcName = fs::path(g.nodes[srcId].name).stem().string();

    for (int j : decoded_neighbors[i]) {
      int tgtId = j, order = 0;
      // if (srcId != 3 || tgtId != 4) continue;
      string tgtName = fs::path(g.nodes[tgtId].name).stem().string();
      GraphT::Edge e;
      e.source = srcId;
      e.target = tgtId;

      fs::path correspondencepixbin = tempfolderdirpath / srcName / fmt::format("{}_{}_correspondence_pix.bin", srcName, tgtName);
      fs::path correspondenceoffsetbin = tempfolderdirpath / srcName / fmt::format("{}_{}_correspondence_offset.bin", srcName, tgtName);

      e.name = fs::relative(correspondencepixbin, rootdirpath).string();
      if (!fs::exists(correspondencepixbin)) {
        logger.error("Cannot find correspondence file: {}", correspondencepixbin.string());
        continue;
      }

      unsigned int numpts;
      read_correspondence_pts_pix_binary_num(correspondencepixbin.string().c_str(), numpts);
      logger.trace("numpts {}", numpts);
      if (numpts < minimum_matches) continue;
      std::unique_ptr<float> matches(new float[7 * numpts]);          //
      std::unique_ptr<float> offset_n_weight(new float[7 * numpts]);  // what this variable's use
      read_correspondence_pts_pix_binary(correspondencepixbin.string().c_str(), matches.get());

      // Now, make combination of src-neighbors and tgt-neighbors, except src-tgt
      std::vector<std::pair<int, int>> xyzgrid_pairs;
      /* for (int adj_srcId : decoded_neighbors[srcId])
         for (int adj_tgtId : decoded_neighbors[tgtId]) {
           if (adj_srcId != tgtId && adj_tgtId != srcId) xyzgrid_pairs.emplace_back(adj_srcId, adj_tgtId);
         }*/
      for (int adj_srcId : decoded_neighbors[srcId])
        if (adj_srcId != tgtId) xyzgrid_pairs.emplace_back(adj_srcId, tgtId);

      logger.info("Edge[{},{}] requires {} grid pairs", e.source, e.target, xyzgrid_pairs.size());

      Eigen::Matrix<double, 7, 1> all_xfm = Eigen::Matrix<double, 7, 1>::Zero();
      Eigen::Matrix<double, 6, 6> all_cov = Eigen::Matrix<double, 6, 6>::Identity();
      all_cov.fill(10000);  //

      //for (auto& p : xyzgrid_pairs) {
        //int adj_srcId, adj_tgtId;
        //std::tie(adj_srcId, adj_tgtId) = p;  // assign value to adj_srcId and ADJ_tgtId
        //std::string adj_srcName = fs::path(g.nodes[adj_srcId].name).stem().string();
        // std::string adj_tgtName = fs::path(g.nodes[adj_tgtId].name).stem().string();
        // clang-format off
    //#if 1
    //    
        logger.info("Edge[{},{}]... ", e.source, e.target);
    //    fs::path srctempgridpath[3] = {
    //      tempfolderdirpath / srcName / fmt::format("{}_{}_temgridX.tif", srcName, adj_srcName),
    //      tempfolderdirpath / srcName / fmt::format("{}_{}_temgridY.tif", srcName, adj_srcName),
    //      tempfolderdirpath / srcName / fmt::format("{}_{}_temgridZ.tif", srcName, adj_srcName)};
    //    fs::path tgttempgridpath[3] = {
    //      tempfolderdirpath / srcName / fmt::format("{}_{}_temgridX.tif", srcName, tgtName),
    //      tempfolderdirpath / srcName / fmt::format("{}_{}_temgridY.tif", srcName, tgtName),
    //      tempfolderdirpath / srcName / fmt::format("{}_{}_temgridZ.tif", srcName, tgtName)};
    //      /*tempfolderdirpath / tgtName / fmt::format("{}_{}_temgridX.tif", tgtName, srcName),
    //      tempfolderdirpath / tgtName / fmt::format("{}_{}_temgridY.tif", tgtName, srcName),
    //      tempfolderdirpath / tgtName / fmt::format("{}_{}_temgridZ.tif", tgtName, srcName)};*/
    //#else
    //    fs::path srctempgridpath[3] = {
    //      tempfolderdirpath / srcName / fmt::format("{}_{}_temgridX.tif", srcName, tgtName),
    //      tempfolderdirpath / srcName / fmt::format("{}_{}_temgridY.tif", srcName, tgtName),
    //      tempfolderdirpath / srcName / fmt::format("{}_{}_temgridZ.tif", srcName, tgtName)};
    //    fs::path tgttempgridpath[3] = {
    //      tempfolderdirpath / tgtName / fmt::format("{}_{}_temgridX.tif", tgtName, srcName),
    //      tempfolderdirpath / tgtName / fmt::format("{}_{}_temgridY.tif", tgtName, srcName),
    //      tempfolderdirpath / tgtName / fmt::format("{}_{}_temgridZ.tif", tgtName, srcName)};
    //#endif
        // clang-format on
        // for (int _d = 0; _d < 3; ++_d) {
        //  if (!fs::exists(srctempgridpath[_d])) {
        //    logger.error("File Not Found: {}", srctempgridpath[_d].string());
        //    return 1;
        //  }
        //  if (!fs::exists(tgttempgridpath[_d])) {
        //    logger.error("File Not Found: {}", tgttempgridpath[_d].string());
        //    return 1;
        //  }
        //}  // for (int _d = 0; _d < 3; ++_d)

        // XYZGrid srcTmpXYZ(preloadgrids);
        // srcTmpXYZ.open(srctempgridpath, fs::path());
        // XYZGrid tgtTmpXYZ(preloadgrids);
        // tgtTmpXYZ.open(tgttempgridpath, fs::path());  // what function of this part

        Eigen::Matrix<double, 7, 1> _pair_xfm;
        Eigen::Matrix<double, 6, 6> _pair_cov;

        bool valid = create_edge_from_corrspondence_pts( poses[srcId], poses[tgtId], matches.get(), numpts, offset_n_weight.get(), invsigmasq_2, minimum_matches,
                                                        _pair_xfm.data(), _pair_cov.data(), &logger);
        if (!valid) continue;
        double update_weight = 1. / all_cov(0, 0) + 1. / _pair_cov(0, 0);

        // use in translation average to make each one have same weight
        // all_xfm = (all_xfm * 1. / all_cov(0, 0) + _pair_xfm * 1. / _pair_cov(0, 0)) / update_weight;
        // all_cov = Eigen::Matrix<double, 7, 7>::Identity() * 1. / update_weight;
        // break;
        //}  // for (auto& p : xyzgrid_pairs)
        double sum = 0;
        for (int i = 0; i < 7; i++) sum += _pair_xfm.data()[i];
        if (abs(sum) <= 1) continue;
        memcpy(e.xfm, _pair_xfm.data(), sizeof(double) * 7);
        memcpy(e.cov, _pair_cov.data(), sizeof(double) * 36);
        e.order = order;
        double initXfm[3] = {g.nodes[tgtId].xfm[0] - g.nodes[srcId].xfm[0], g.nodes[tgtId].xfm[1] - g.nodes[srcId].xfm[1],
                             g.nodes[tgtId].xfm[2] - g.nodes[srcId].xfm[2]};  // according the xfm adjust relative translation

        double diffXfm[3] = {std::abs(initXfm[0] - e.xfm[0]), std::abs(initXfm[1] - e.xfm[1]), std::abs(initXfm[2] - e.xfm[2])};

        logger.debug("Init: {} {} {}", initXfm[0], initXfm[1], initXfm[2]);
        logger.debug("Edge: {} {} {}", e.xfm[0], e.xfm[1], e.xfm[2]);
        logger.debug("Diff: {} {} {} cm {}, {}", diffXfm[0] * 100, diffXfm[1] * 100, diffXfm[2] * 100, srcName, tgtName);
        //#pragma omp critical
        g.insertEdge(e);

        order++;
      
      // break;
    }

    // break;
  }

  cout << "prepare to save"
       << "\n";
  total_stat_reprojection_error.print(logger);
  ofstream ofs(outfilepath.string());
  ofs << &g;
  ofs.close();
  logger.info("Done");

  return 0;
}

bool IsTriangulatedPointInFrontOfCameras(const Eigen::Vector3d& pt1, const Eigen::Vector3d& pt2, const Eigen::Matrix3d& R, const Eigen::Vector3d& t) {
  const Eigen::Vector3d c = -R.transpose() * t;

  const Eigen::Vector3d dir1 = pt1;
  const Eigen::Vector3d dir2 = R.transpose() * pt2;

  const double dir1_sq = dir1.squaredNorm();
  const double dir2_sq = dir2.squaredNorm();
  const double dir1_dir2 = dir1.dot(dir2);
  const double dir1_pos = dir1.dot(c);
  const double dir2_pos = dir2.dot(c);

  return (dir2_sq * dir1_pos - dir1_dir2 * dir2_pos > 0 && dir1_dir2 * dir1_pos - dir1_sq * dir2_pos > 0);
}