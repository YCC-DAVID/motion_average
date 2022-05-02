#include "mystructs.h"

#include <iomanip>
#include <string>

#include "myutils.h"

using namespace std;

namespace motionavg {
namespace Affine2D {
void BBox::add(double x, double y) {
  ptMin[0] = std::min<double>(x, ptMin[0]);
  ptMin[1] = std::min<double>(y, ptMin[1]);
  ptMax[0] = std::max<double>(x, ptMax[0]);
  ptMax[1] = std::max<double>(y, ptMax[1]);
}

BBox BBox::transform(const double* xfm) const {
  BBox result;
  double tmpPt[2];
  XfmApply(xfm, ptMin[0], ptMin[1], tmpPt[0], tmpPt[1]);
  result.add(tmpPt[0], tmpPt[1]);
  XfmApply(xfm, ptMin[0], ptMax[1], tmpPt[0], tmpPt[1]);
  result.add(tmpPt[0], tmpPt[1]);
  XfmApply(xfm, ptMax[0], ptMin[1], tmpPt[0], tmpPt[1]);
  result.add(tmpPt[0], tmpPt[1]);
  XfmApply(xfm, ptMax[0], ptMax[1], tmpPt[0], tmpPt[1]);
  result.add(tmpPt[0], tmpPt[1]);
  return result;
}

BBox BBox::intersect(const BBox& other) const {
  BBox result;
  result.ptMin[0] = std::max<double>(ptMin[0], other.ptMin[0]);
  result.ptMin[1] = std::max<double>(ptMin[1], other.ptMin[1]);
  result.ptMax[0] = std::min<double>(ptMax[0], other.ptMax[0]);
  result.ptMax[1] = std::min<double>(ptMax[1], other.ptMax[1]);
  return result;
}

bool AffineModel::load(std::string p) {
  bool succ = true;
  filename = fs::path(p).filename().stem().string();

  // TODO: the filename encoding is not program friendly. Add a better deliminator?
  target_file = filename.substr(0, (filename.size() - 1) / 2);
  source_file = filename.substr(target_file.size() + 1, (filename.size() - 1) / 2);

  // Read content
  ifstream ifs(p);
  string line;
  std::getline(ifs, line);
  if (line.compare(0, 10, "affine parameter pixels", 10) == 0) {
    std::getline(ifs, line);
    auto tokens = split(line, ' ');
    if (tokens.size() >= 6) {
      for (int i = 0; i < 6; ++i) pixelParam[i] = std::stod(tokens[i]);
    } else {
      succ = false;
    }
  } else {
    succ = false;
  }

  std::getline(ifs, line);
  if (line.compare(0, 10, "affine parameter GEO", 10) == 0) {
    std::getline(ifs, line);
    auto tokens = split(line, ' ');
    if (tokens.size() >= 6) {
      for (int i = 0; i < 6; ++i) geoParam[i] = std::stod(tokens[i]);
    } else {
      succ = false;
    }
  } else {
    succ = false;
  }

  ifs.close();
  valid = succ;
  return succ;
}

bool PairwiseData::load(std::string p) {
  bool succ = affine.load(p);
  // then find feature pair
  if (succ) {
    fs::path parampath(p);
    fs::path rootdir = parampath.parent_path();
    string filename = parampath.filename().stem().string();
    fs::path ptspath = rootdir / (filename + "_feature_points.txt");

    tiepoints.clear();

    if (fs::exists(ptspath)) {
      ifstream ifs(ptspath.string());
      for (std::string line; std::getline(ifs, line);) {
        auto tokens = split(line, ' ');
        if (tokens.size() == 4) {
          TiePoint tp;
          tp.target[0] = std::stod(tokens[2]);
          tp.target[1] = std::stod(tokens[3]);
          tp.source[0] = std::stod(tokens[0]);
          tp.source[1] = std::stod(tokens[1]);
          tiepoints.push_back(tp);
        }
      }
      ifs.close();
    }
  }
  return succ;
}

size_t XfmGraph::insertNode(XfmGraph::Node& n) {
  nodes.push_back(n);
  return nodes.size() - 1;
}

void XfmGraph::rebase(int i) {
  if (i < 0) {
    double invBase[6];
    XfmInv(baseGeo, invBase);
    rebase(invBase);
    baseID = -1;
  } else if (i >= nodes.size()) {
    return;
  } else {
    rebase(nodes[i].poseXfm);
    baseID = i;
  }
}

void XfmGraph::rebase(const double const* newbase) {
  double invNewBase[6], rebaseMat[6], newbaseCopy[6];

  XfmComposite(baseGeo, newbase, newbaseCopy);

  XfmInv(newbaseCopy, invNewBase);
  XfmComposite(invNewBase, baseGeo, rebaseMat);

  double tmpGeo[6];
  for (size_t i = 0; i < nodes.size(); ++i) {
    XfmComposite(rebaseMat, nodes[i].poseXfm, tmpGeo);
    std::copy_n(tmpGeo, 6, nodes[i].poseXfm);
  }
  std::copy_n(newbaseCopy, 6, baseGeo);
}

size_t PoseGraph::insertEdge(PoseGraph::Edge& e) {
  edges.push_back(e);
  return edges.size() - 1;
}

size_t BundleGraph::insertEdge(BundleGraph::Edge& e) {
  edges.push_back(e);
  return edges.size() - 1;
}
}  // namespace Affine2D
}  // namespace motionavg

std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::AffineModel& m) {
  if (m.target_file.empty() || m.source_file.empty())
    ofs << "Name: " << m.filename << "\n";
  else
    ofs << "Target: " << m.target_file << "\nSource: " << m.source_file << "\n";

  ofs << " PixParam " << m.pixelParam[0] << " " << m.pixelParam[1] << " " << m.pixelParam[2] << " " << m.pixelParam[3] << " " << m.pixelParam[4] << " " << m.pixelParam[5] << "\n";
  ofs << " GeoParam " << m.geoParam[0] << " " << m.geoParam[1] << " " << m.geoParam[2] << " " << m.geoParam[3] << " " << m.geoParam[4] << " " << m.geoParam[5] << "\n";
  return ofs;
}

// -----------------------------------------

std::istream& operator>>(std::istream& ifs, motionavg::Affine2D::XfmGraph::Node& n) {
  char dummy;
  ifs >> n.name >> dummy >> n.path >> dummy >> n.left >> dummy >> n.top >> dummy >> n.right >> dummy >> n.bottom >> dummy;
  for (int i = 0; i < 6; ++i) ifs >> n.poseXfm[i] >> dummy;
  return ifs;
}

std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::XfmGraph::Node& n) {
  ofs << n.name << '\n' << n.path << '\n' << n.left << ' ' << n.top << ' ' << n.right << ' ' << n.bottom << '\n';

  for (int i = 0; i < 6; ++i) ofs << n.poseXfm[i] << ((i < 5) ? ',' : '\n');
  return ofs;
}

std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::XfmGraph& g) {
  ofs << std::fixed << std::setprecision(15);
  ofs << "XfmGraph\n";
  ofs << g.basepath << '\n';
  ofs << g.baseID << '\n';
  for (int i = 0; i < 6; ++i) {
    if (i > 0) ofs << ' ';
    ofs << g.baseGeo[i];
  }
  return ofs;
}

std::istream& operator>>(std::istream& ifs, motionavg::Affine2D::XfmGraph& g) {
  ifs.unsetf(std::ios_base::skipws);
  char dummy;
  string format_identifier;
  ifs >> format_identifier >> dummy;

  ifs >> g.basepath >> dummy;
  ifs >> g.baseID >> dummy;
  for (int i = 0; i < 6; ++i) ifs >> g.baseGeo[i] >> dummy;
  size_t num_nodes;
  ifs >> num_nodes >> dummy;
  g.nodes.resize(num_nodes);
  for (int i = 0; i < num_nodes; ++i) ifs >> g.nodes[i];
  return ifs;
}

// -----------------------------------------

std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::PoseGraph& g) {
  ofs << std::fixed << std::setprecision(15);
  ofs << "PoseGraph\n";
  ofs << g.basepath << '\n';
  ofs << g.baseID << '\n';
  for (int i = 0; i < 6; ++i) {
    if (i > 0) ofs << ' ';
    ofs << g.baseGeo[i];
  }
  ofs << '\n';
  ofs << g.nodes.size() << '\n';
  for (size_t i = 0; i < g.nodes.size(); ++i) ofs << g.nodes[i];
  ofs << g.edges.size() << '\n';
  for (size_t i = 0; i < g.edges.size(); ++i) ofs << g.edges[i];
  return ofs;
}

std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::PoseGraph::Edge& e) {
  ofs << e.name << '\n' << e.target << ',' << e.source << '\n';
  for (int i = 0; i < 6; ++i) ofs << e.regXfm[i] << ((i < 5) ? ',' : '\n');
  for (int i = 0; i < 36; ++i) ofs << e.covXfm[i] << ((i < 35) ? ',' : '\n');
  return ofs;
}

std::istream& operator>>(std::istream& ifs, motionavg::Affine2D::PoseGraph& g) {
  ifs.unsetf(std::ios_base::skipws);
  char dummy;
  string format_identifier;
  ifs >> format_identifier >> dummy;
  if (format_identifier != "PoseGraph") return ifs;
  ifs >> g.basepath >> dummy;
  ifs >> g.baseID >> dummy;
  for (int i = 0; i < 6; ++i) ifs >> g.baseGeo[i] >> dummy;
  size_t num_nodes, num_edges;
  ifs >> num_nodes >> dummy;
  g.nodes.resize(num_nodes);
  for (int i = 0; i < num_nodes; ++i) ifs >> g.nodes[i];
  ifs >> num_edges >> dummy;
  g.edges.resize(num_edges);
  for (int i = 0; i < num_edges; ++i) ifs >> g.edges[i];
  return ifs;
}

std::istream& operator>>(std::istream& ifs, motionavg::Affine2D::PoseGraph::Edge& e) {
  char dummy;
  ifs >> e.name >> dummy >> e.target >> dummy >> e.source >> dummy;
  for (int i = 0; i < 6; ++i) ifs >> e.regXfm[i] >> dummy;
  for (int i = 0; i < 36; ++i) ifs >> e.covXfm[i] >> dummy;
  return ifs;
}

// -----------------------------------------

std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::BundleGraph& g) {
  ofs << std::fixed << std::setprecision(15);
  ofs << "BundleGraph\n";
  ofs << g.basepath << '\n';
  ofs << g.baseID << '\n';
  for (int i = 0; i < 6; ++i) {
    if (i > 0) ofs << ' ';
    ofs << g.baseGeo[i];
  }
  ofs << '\n';
  ofs << g.nodes.size() << '\n';
  for (size_t i = 0; i < g.nodes.size(); ++i) ofs << g.nodes[i];
  ofs << g.edges.size() << '\n';
  for (size_t i = 0; i < g.edges.size(); ++i) ofs << g.edges[i];
  return ofs;
}

std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::BundleGraph::Edge& e) {
  ofs << e.name << '\n' << e.target << ',' << e.source << '\n';
  ofs << e.tiepoints.size() << '\n';
  for (size_t i = 0; i < e.tiepoints.size(); ++i) ofs << e.tiepoints[i].target[0] << ' ' << e.tiepoints[i].target[1] << ' ' << e.tiepoints[i].source[0] << ' ' << e.tiepoints[i].source[1] << '\n';
  return ofs;
}

std::istream& operator>>(std::istream& ifs, motionavg::Affine2D::BundleGraph& g) {
  ifs.unsetf(std::ios_base::skipws);
  char dummy;
  string format_identifier;
  ifs >> format_identifier >> dummy;
  if (format_identifier != "BundleGraph") return ifs;
  ifs >> g.basepath >> dummy;
  ifs >> g.baseID >> dummy;
  for (int i = 0; i < 6; ++i) ifs >> g.baseGeo[i] >> dummy;
  size_t num_nodes, num_edges;
  ifs >> num_nodes >> dummy;
  g.nodes.resize(num_nodes);
  for (int i = 0; i < num_nodes; ++i) ifs >> g.nodes[i];
  ifs >> num_edges >> dummy;
  g.edges.resize(num_edges);
  for (int i = 0; i < num_edges; ++i) ifs >> g.edges[i];
  return ifs;
}

std::istream& operator>>(std::istream& ifs, motionavg::Affine2D::BundleGraph::Edge& e) {
  char dummy;
  ifs >> e.name >> dummy >> e.target >> dummy >> e.source >> dummy;
  size_t num_points;
  ifs >> num_points >> dummy;
  e.tiepoints.resize(num_points);
  for (size_t i = 0; i < num_points; ++i) {
    ifs >> e.tiepoints[i].target[0] >> dummy >> e.tiepoints[i].target[1] >> dummy >> e.tiepoints[i].source[0] >> dummy >> e.tiepoints[i].source[1] >> dummy;
  }
  return ifs;
}