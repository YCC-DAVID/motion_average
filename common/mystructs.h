#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
namespace motionavg {
namespace TranslateND {
template <int N>
struct Node {
  std::string name;
  std::string path;
  double xfm[N];
  Node() { std::fill_n(xfm, N, 0.); }
};
template <int N>
struct Edge {
  std::string name;
  int source;
  int target;
  double xfm[N];
  double cov[N * N];
  Edge() {
    std::fill_n(xfm, N, 0.);
    std::fill_n(cov, N * N, 0.);
    for (int i = 0; i < N; ++i) cov[N * i + i] = 1.;
  }
};

template <int N>
struct TranslateGraph {
  static constexpr int DIM = N;
  using Node = Node<N>;
  using Edge = Edge<N>;

  std::string basepath;
  int baseID = -1;
  double baseGeo[N];
  std::vector<Node> nodes;
  std::vector<Edge> edges;

  TranslateGraph() { std::fill_n(baseGeo, N, 0.); }
  size_t insertNode(Node& n) {
    nodes.push_back(n);
    return nodes.size() - 1;
  }
  size_t insertEdge(Edge& e) {
    edges.push_back(e);
    return edges.size() - 1;
  }
  void rebase(int i) {
    if (i < 0) {
      double invBase[N];
      for (int i = 0; i < N; ++i) invBase[i] = -baseGeo[i];
      rebase(invBase);
      baseID = -1;
    } else if (i >= nodes.size()) {
      return;
    } else {
      rebase(nodes[i].xfm);
      baseID = i;
    }
  }
  void rebase(const double* newbase) {
    double newbaseCopy[N];
    for (int _d = 0; _d < N; ++_d) newbaseCopy[_d] = baseGeo[_d] + newbase[_d];

    for (size_t i = 0; i < nodes.size(); ++i) {
      for (int _d = 0; _d < N; ++_d) nodes[i].xfm[_d] = baseGeo[_d] + nodes[i].xfm[_d] - newbaseCopy[_d];
    }
    std::copy_n(newbaseCopy, N, baseGeo);
  }
};

// Store 3D point table and
struct GCPGraph {
  struct GCPPoint {
    GCPPoint() : x(0.f), y(0.f), z(0.f), ex(999.f), ey(999.f), ez(999.f){};
    std::string name;
    double x, y, z;
    double ex, ey, ez;  // std deviation
  };
  struct GCPLink {
    GCPLink() : gcpid(-1), viewid(-1), u(0.f), v(0.f), dx(0.f), dy(0.f), dz(0.f), ex(999.f), ey(999.f), ez(999.f){};
    int gcpid;
    int viewid;
    float u, v;         // u (col), v (row) on image (top left origin)
    float dx, dy, dz;   // 3d translation from camera center to GCP
    double ex, ey, ez;  // std deviation of 3d translation
  };

  GCPGraph() : baseGCP{0.f} {}

  std::vector<GCPPoint> gcps;
  std::vector<GCPLink> gcplinks;
  double baseGCP[3];
  void rebaseGCP(const double* newbase) {
    double newbaseCopy[3];
    std::copy_n(newbase, 3, newbaseCopy);
    for (auto& gcp : gcps) {
      gcp.x += baseGCP[0] - newbaseCopy[0];
      gcp.y += baseGCP[1] - newbaseCopy[1];
      gcp.z += baseGCP[2] - newbaseCopy[2];
    }
    std::copy_n(newbaseCopy, 3, baseGCP);
  }
};

struct TranslateGraph3WithGCP : public TranslateGraph<3>, public GCPGraph {
  
  void rebase(int i) {
    TranslateGraph<3>::rebase(i);
    GCPGraph::rebaseGCP(baseGeo);
  }
  void rebase(const double* newbase) {
    TranslateGraph<3>::rebase(newbase);
    GCPGraph::rebaseGCP(baseGeo);
    //assert((baseGeo[0] == baseGCP[0]) && (baseGeo[1] == baseGCP[1]) && (baseGeo[2] == baseGCP[2]));
  }
};
}  // namespace TranslateND

namespace Affine2D {
struct BBox {
  void add(double x, double y);
  BBox transform(const double* xfm) const;
  BBox intersect(const BBox& other) const;

  double ptMin[2] = {DBL_MAX, DBL_MAX};
  double ptMax[2] = {-DBL_MAX, -DBL_MAX};
};

struct AffineModel {
  bool load(std::string p);
  std::string filename;
  std::string target_file;
  std::string source_file;
  double pixelParam[6];
  double geoParam[6];
  bool valid = false;
};

struct TiePoint {
  double target[2];
  double source[2];
};

struct PairwiseData {
  bool load(std::string p);

  AffineModel affine;
  std::vector<TiePoint> tiepoints;

  bool valid = false;
};

struct Node {
  const static int DIM = 6;
  std::string name;
  std::string path;
  int left, right, top, bottom;
  double poseXfm[DIM];
};

struct XfmGraph {
  using Node = Node;
  const static int DIM = Node::DIM;
  std::string basepath;
  int baseID = -1;
  double baseGeo[DIM] = {0, 1, 0, 0, 0, 1};
  std::vector<Node> nodes;
  size_t insertNode(Node& n);

  void rebase(int i);
  void rebase(const double* newbase);
};

struct PoseGraph : public XfmGraph {
  struct Edge {
    std::string name;
    int target;
    int source;
    double regXfm[DIM] = {0, 1, 0, 0, 0, 1};
    double covXfm[DIM * DIM] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};

    Edge() : target(-1), source(-1) {}
  };

  std::vector<Edge> edges;
  size_t insertEdge(Edge& e);
};

struct BundleGraph : public XfmGraph {
  struct Edge {
    std::string name;
    int target;
    int source;
    std::vector<TiePoint> tiepoints;

    Edge() : target(-1), source(-1) {}
  };

  std::vector<Edge> edges;
  size_t insertEdge(Edge& e);
};

/*
 * Composition of Affine Transformations
 * lMat = [	l(1)	l(2)	l(0);	rMat = [	r(1)	r(2)	r(0);
 *			l(4)	l(5)	l(3);				r(4)	r(5)	r(3);
 *			0		0		1	]				0		0		1	]
 *
 * cMat = [	c(1)	c(2)	c(0);
 *			c(4)	c(5)	c(3);
 *			0		0		1	]
 *
 *	= [	l(1)*r(1)+l(2)*r(4)	l(1)*r(2)+l(2)*r(5)	l(1)*r(0)+l(2)*r(3)+l(0)	;
 *		l(4)*r(1)+l(5)*r(4)	l(4)*r(2)+l(5)*r(5)	l(4)*r(0)+l(5)*r(3)+l(3)	;
 *		0					0					1							]
 *
 */
template <typename LT, typename RT, typename CT>
void XfmComposite(const LT* const l, const RT* const r, CT* c) {
  c[0] = l[1] * r[0] + l[2] * r[3] + l[0];
  c[1] = l[1] * r[1] + l[2] * r[4];
  c[2] = l[1] * r[2] + l[2] * r[5];
  c[3] = l[4] * r[0] + l[5] * r[3] + l[3];
  c[4] = l[4] * r[1] + l[5] * r[4];
  c[5] = l[4] * r[2] + l[5] * r[5];
}

/*
 * Inverse of affine transformation
 * vMat =
 *
 * D = (v1*v5 - v2*v4)
 *
 * inv(vMat) = [	i(1)	i(2)	i(0);
 *				i(4)	i(5)	i(3);
 *				0		0		1	]
 *	=	[  v5/D		-v2/D	-(v0*v5 - v2*v3)/D	;
 *			-v4/D	v1/D	(v0*v4 - v1*v3)/D	;
 *			0		0		1					]
 */
template <typename IT, typename OT>
void XfmInv(const IT* v, OT* i) {
  IT D = (v[1] * v[5] - v[2] * v[4]);
  i[0] = -(v[0] * v[5] - v[2] * v[3]) / D;
  i[1] = v[5] / D;
  i[2] = -v[2] / D;
  i[3] = (v[0] * v[4] - v[1] * v[3]) / D;
  i[4] = -v[4] / D;
  i[5] = v[1] / D;
}

template <typename PT, typename IT, typename OT>
void XfmApply(const PT* param, const IT x, const IT y, OT& ox, OT& oy) {
  ox = param[0] + param[1] * x + param[2] * y;
  oy = param[3] + param[4] * x + param[5] * y;
}
}  // namespace Affine2D
}  // namespace motionavg

std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::AffineModel& m);

std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::XfmGraph& g);
std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::XfmGraph::Node& n);
std::istream& operator>>(std::istream& ifs, motionavg::Affine2D::XfmGraph& g);
std::istream& operator>>(std::istream& ifs, motionavg::Affine2D::XfmGraph::Node& n);

std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::PoseGraph& g);
std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::PoseGraph::Edge& e);
std::istream& operator>>(std::istream& ifs, motionavg::Affine2D::PoseGraph& g);
std::istream& operator>>(std::istream& ifs, motionavg::Affine2D::PoseGraph::Edge& e);

std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::BundleGraph& g);
std::ostream& operator<<(std::ostream& ofs, const motionavg::Affine2D::BundleGraph::Edge& e);
std::istream& operator>>(std::istream& ifs, motionavg::Affine2D::BundleGraph& g);
std::istream& operator>>(std::istream& ifs, motionavg::Affine2D::BundleGraph::Edge& e);

template <int N>
inline std::ostream& operator<<(std::ostream& ofs, const motionavg::TranslateND::Node<N>& n) {
  ofs << n.name << '\n' << n.path << '\n';
  for (int i = 0; i < N; ++i) ofs << n.xfm[i] << ((i < N - 1) ? ',' : '\n');
  return ofs;
}

template <int N>
inline std::ostream& operator<<(std::ostream& ofs, const motionavg::TranslateND::Edge<N>& e) {
  ofs << e.name << '\n' << e.target << ',' << e.source << '\n';
  for (int i = 0; i < N; ++i) ofs << e.xfm[i] << ((i < N - 1) ? ',' : '\n');
  for (int i = 0; i < N * N; ++i) ofs << e.cov[i] << ((i < N * N - 1) ? ',' : '\n');
  return ofs;
}

template <int N>
inline std::ostream& operator<<(std::ostream& ofs, const motionavg::TranslateND::TranslateGraph<N>* g) {
  ofs << std::fixed << std::setprecision(15);
  ofs << "TranslateGraph " << N << '\n';
  ofs << g->basepath << '\n';
  ofs << g->baseID << '\n';
  for (int i = 0; i < N; ++i) {
    if (i > 0) ofs << ' ';
    ofs << g->baseGeo[i];
  }
  ofs << '\n';
  ofs << g->nodes.size() << '\n';
  for (size_t i = 0; i < g->nodes.size(); ++i) ofs << g->nodes[i];
  ofs << g->edges.size() << '\n';
  for (size_t i = 0; i < g->edges.size(); ++i) ofs << g->edges[i];
  return ofs;
}

inline std::ostream& operator<<(std::ostream& ofs, const motionavg::TranslateND::GCPGraph::GCPPoint& p) {
  ofs << p.name << '\n' << p.x << ',' << p.y << ',' << p.z << '\n';
  ofs << p.ex << ',' << p.ey << ',' << p.ez << '\n';
  return ofs;
}

inline std::ostream& operator<<(std::ostream& ofs, const motionavg::TranslateND::GCPGraph::GCPLink& l) {
  ofs << l.gcpid << ',' << l.viewid << ',' << l.u << ',' << l.v << '\n';
  ofs << l.dx << ',' << l.dy << ',' << l.dz << '\n';
  ofs << l.ex << ',' << l.ey << ',' << l.ez << '\n';
  return ofs;
}

inline std::ostream& operator<<(std::ostream& ofs, const motionavg::TranslateND::GCPGraph* g) {
  ofs << std::fixed << std::setprecision(15);
  ofs << "GCPGraph\n";
  ofs << g->gcps.size() << '\n';
  for (size_t i = 0; i < g->gcps.size(); ++i) ofs << g->gcps[i];
  ofs << g->gcplinks.size() << '\n';
  for (size_t i = 0; i < g->gcplinks.size(); ++i) ofs << g->gcplinks[i];
  return ofs;
}

inline std::ostream& operator<<(std::ostream& ofs, const motionavg::TranslateND::TranslateGraph3WithGCP* g) {
  ofs << "Translate3GraphWithGCP\n";
  ofs << static_cast<const motionavg::TranslateND::TranslateGraph<3>*>(g);
  ofs << static_cast<const motionavg::TranslateND::GCPGraph*>(g);
  return ofs;
}

////

template <int N>
inline std::istream& operator>>(std::istream& ifs, motionavg::TranslateND::Node<N>& n) {
  char dummy;
  ifs >> n.name >> dummy >> n.path >> dummy;
  for (int i = 0; i < N; ++i) ifs >> n.xfm[i] >> dummy;
  return ifs;
}

template <int N>
inline std::istream& operator>>(std::istream& ifs, motionavg::TranslateND::Edge<N>& e) {
  char dummy;
  ifs >> e.name >> dummy >> e.target >> dummy >> e.source >> dummy;
  for (int i = 0; i < N; ++i) ifs >> e.xfm[i] >> dummy;
  for (int i = 0; i < N * N; ++i) ifs >> e.cov[i] >> dummy;
  return ifs;
}

template <int N>
inline std::istream& operator>>(std::istream& ifs, motionavg::TranslateND::TranslateGraph<N>* g) {
  ifs.unsetf(std::ios_base::skipws);
  char dummy;
  std::string format_identifier;
  int dim;
  ifs >> format_identifier >> dummy >> dim >> dummy;
  if (format_identifier != "TranslateGraph") return ifs;
  if (dim != N) return ifs;
  ifs >> g->basepath >> dummy;
  ifs >> g->baseID >> dummy;
  for (int i = 0; i < N; ++i) ifs >> g->baseGeo[i] >> dummy;
  size_t num_nodes, num_edges;
  ifs >> num_nodes >> dummy;
  g->nodes.resize(num_nodes);
  for (int i = 0; i < num_nodes; ++i) ifs >> g->nodes[i];
  ifs >> num_edges >> dummy;
  g->edges.resize(num_edges);
  for (int i = 0; i < num_edges; ++i) ifs >> g->edges[i];
  return ifs;
}

inline std::istream& operator>>(std::istream& ifs, motionavg::TranslateND::GCPGraph::GCPPoint& p) {
  char dummy;
  ifs >> p.name >> dummy >> p.x >> dummy >> p.y >> dummy >> p.z >> dummy;
  ifs >> p.ex >> dummy >> p.ey >> dummy >> p.ez >> dummy;
  return ifs;
}

inline std::istream& operator>>(std::istream& ifs, motionavg::TranslateND::GCPGraph::GCPLink& l) {
  char dummy;
  ifs >> l.gcpid >> dummy >> l.viewid >> dummy >> l.u >> dummy >> l.v >> dummy;
  ifs >> l.dx >> dummy >> l.dy >> dummy >> l.dz >> dummy;
  ifs >> l.ex >> dummy >> l.ey >> dummy >> l.ez >> dummy;
  return ifs;
}

inline std::istream& operator>>(std::istream& ifs, motionavg::TranslateND::GCPGraph* g) {
  ifs.unsetf(std::ios_base::skipws);
  char dummy;
  std::string format_identifier;
  ifs >> format_identifier >> dummy;
  if (format_identifier != "GCPGraph") return ifs;
  size_t num_gcps, num_links;
  ifs >> num_gcps >> dummy;
  g->gcps.resize(num_gcps);
  for (int i = 0; i < num_gcps; ++i) ifs >> g->gcps[i];
  ifs >> num_links >> dummy;
  g->gcplinks.resize(num_links);
  for (int i = 0; i < num_links; ++i) ifs >> g->gcplinks[i];
  return ifs;
}

inline std::istream& operator>>(std::istream& ifs, motionavg::TranslateND::TranslateGraph3WithGCP* g) {
  ifs.unsetf(std::ios_base::skipws);
  char dummy;
  std::string format_identifier;
  ifs >> format_identifier >> dummy;
  if (format_identifier != "Translate3GraphWithGCP") return ifs;
  ifs >> static_cast<motionavg::TranslateND::TranslateGraph<3>*>(g);
  ifs >> static_cast<motionavg::TranslateND::GCPGraph*>(g);
  return ifs;
}