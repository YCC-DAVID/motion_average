#pragma once

#include <iostream>
#include <fstream>
#include <vector>
namespace motionavg {

	struct BBox {
		void add(double x, double y);
		BBox transform(const double* xfm) const;
		BBox intersect(const BBox& other) const;

		double ptMin[2] = { DBL_MAX, DBL_MAX };
		double ptMax[2] = { -DBL_MAX, -DBL_MAX };
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
		std::string name;
		std::string path;
		int left, right, top, bottom;
		double poseXfm[6];
	};

	struct XfmGraph {
		using Node = Node;
		std::string basepath;
		int baseID = -1;
		double baseGeo[6] = { 0, 1, 0, 0, 0, 1 };
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
			double regXfm[6] = { 0,1,0,0,0,1 };
			double covXfm[36] = { 1,0,0,0,0,0,
								0,1,0,0,0,0,
								0,0,1,0,0,0,
								0,0,0,1,0,0,
								0,0,0,0,1,0,
								0,0,0,0,0,1 };

			Edge() :target(-1), source(-1) {}
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

			Edge() :target(-1), source(-1) {}
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
	void XfmComposite(const LT* const l, const RT* const r, CT* c)
	{
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
	void XfmInv(const IT* v, OT* i)
	{
		IT D = (v[1] * v[5] - v[2] * v[4]);
		i[0] = -(v[0] * v[5] - v[2] * v[3]) / D;
		i[1] = v[5] / D;
		i[2] = -v[2] / D;
		i[3] = (v[0] * v[4] - v[1] * v[3]) / D;
		i[4] = -v[4] / D;
		i[5] = v[1] / D;
	}

	template<typename PT, typename IT, typename OT>
	void XfmApply(const PT* param, const IT x, const IT y, OT& ox, OT& oy) {
		ox = param[0] + param[1] * x + param[2] * y;
		oy = param[3] + param[4] * x + param[5] * y;
	}
}

std::ostream& operator<<(std::ostream& ofs, const motionavg::AffineModel& m);

std::ostream& operator<<(std::ostream& ofs, const motionavg::XfmGraph::Node& n);
std::istream& operator>>(std::istream& ifs, motionavg::XfmGraph::Node& n);

std::ostream& operator<<(std::ostream& ofs, const motionavg::PoseGraph& g);
std::ostream& operator<<(std::ostream& ofs, const motionavg::PoseGraph::Edge& e);

std::istream& operator>>(std::istream& ifs, motionavg::PoseGraph& g);
std::istream& operator>>(std::istream& ifs, motionavg::PoseGraph::Edge& e);

std::ostream& operator<<(std::ostream& ofs, const motionavg::BundleGraph& g);
std::ostream& operator<<(std::ostream& ofs, const motionavg::BundleGraph::Edge& e);

std::istream& operator>>(std::istream& ifs, motionavg::BundleGraph& g);
std::istream& operator>>(std::istream& ifs, motionavg::BundleGraph::Edge& e);