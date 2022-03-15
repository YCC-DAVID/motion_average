#include <fstream>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "myutils.h"
#include "mystructs.h"

using namespace std;
using namespace motionavg;

int main(int argc, char** argv)
{
	//int graphtype = atoi(argv[1]);
	fs::path graphfilepath(argv[1]);
	fs::path rootdir = graphfilepath.parent_path();
	string filename = graphfilepath.filename().stem().string();
	fs::path outgraphpath = rootdir / (filename + "_CastPoseGraph.txt");

	int mode = 0;	// convert from pose
	if (argc > 2) mode = 1;	//estimate from points
	/*if (graphtype == 0)
	{
		fs::path outgraphpath = rootdir / (filename + "_castbundle.txt");
	}
	else {
		fs::path outgraphpath = rootdir / (filename + "_castpose.txt");
	}*/

	BundleGraph bgraph;
	ifstream ifs(graphfilepath.string());
	ifs >> bgraph;
	ifs.close();

	PoseGraph pgraph;
	pgraph.basepath = bgraph.basepath;
	pgraph.nodes.resize(bgraph.nodes.size());
	pgraph.edges.resize(bgraph.edges.size());
	// Cast node
	for (size_t i = 0; i < pgraph.nodes.size(); ++i)
	{
		const auto& bnode = bgraph.nodes[i];
		auto& pnode = pgraph.nodes[i];
		pnode.name = bnode.name;
		pnode.path = bnode.path;
		pnode.left = bnode.left;
		pnode.top = bnode.top;
		pnode.right = bnode.right;
		pnode.bottom = bnode.bottom;
		std::copy_n(bnode.poseXfm, 6, pnode.poseXfm);
	}
	// Cast Edge
	if (mode == 0)
	{ // from node pose
		double tmpGeo[6];
		for (size_t i = 0; i < pgraph.edges.size(); ++i)
		{
			const auto& bedge = bgraph.edges[i];
			auto& pedge = pgraph.edges[i];
			pedge.name = bedge.name;
			pedge.source = bedge.source;
			pedge.target = bedge.target;

			const auto& sourcenode = pgraph.nodes[pedge.source];
			const auto& targetnode = pgraph.nodes[pedge.target];
			XfmInv(targetnode.poseXfm, tmpGeo);
			XfmComposite(tmpGeo, sourcenode.poseXfm, pedge.regXfm);
		}
	}
	else { // from tie points
		double tmpGeo[6], tmpGeo2[6], tmpPt[2];
		for (size_t i = 0; i < pgraph.edges.size(); ++i)
		{
			const auto& bedge = bgraph.edges[i];
			auto& pedge = pgraph.edges[i];
			pedge.name = bedge.name;
			pedge.source = bedge.source;
			pedge.target = bedge.target;

			const auto& sourcenode = pgraph.nodes[pedge.source];
			const auto& targetnode = pgraph.nodes[pedge.target];

			size_t num_pt = bedge.tiepoints.size();
			Eigen::MatrixXd A = Eigen::MatrixXd(2 * num_pt, 6),
				b = Eigen::MatrixXd(2*num_pt,1);

			A.setZero();

			XfmInv(targetnode.poseXfm, tmpGeo);
			XfmComposite(tmpGeo, sourcenode.poseXfm, tmpGeo2);

			for (size_t j = 0; j < num_pt; ++j)
			{
				const auto& tp = bedge.tiepoints[j];
				//XfmApply(tmpGeo2, tp.source[0], tp.source[1], tmpPt[0], tmpPt[1]);
				b(2 * j) = tp.target[0];
				b(2 * j + 1) = tp.target[1];
				A(2 * j, 0) = 1;
				A(2 * j, 1) = tp.source[0];
				A(2 * j, 2) = tp.source[1];
				A(2 * j + 1, 3) = 1;
				A(2 * j + 1, 4) = tp.source[0];
				A(2 * j + 1, 5) = tp.source[1];
			}

			
			Eigen::MatrixXd x = A.colPivHouseholderQr().solve(b);
			
			Eigen::MatrixXd cov = (A.transpose() * A).inverse()*num_pt;

			std::copy_n(x.data(), 6, pedge.regXfm);
			
			std::copy_n(cov.data(), 36, pedge.covXfm);
		}
	}

	ofstream ofs(outgraphpath.string());

	ofs << pgraph;
	ofs.close();

	//
	cout << "Write to " << outgraphpath << endl;

	return 0;
}