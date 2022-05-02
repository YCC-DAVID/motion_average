#include "myutils.h"
#include "mystructs.h"

#include <ceres/solver.h>
#include <ceres/problem.h>
#include <ceres/autodiff_cost_function.h>
#include <ceres/local_parameterization.h>
#include <ceres/loss_function.h>

using namespace std;
using namespace motionavg::Affine2D;


struct Reprojection_Error {

	Reprojection_Error(const double* target, const double* source)
	{
		std::copy_n(target, 2, observed_target);
		std::copy_n(source, 2, observed_source);
	}
	Reprojection_Error(const double& target_x, const double& target_y, const double& source_x, const double& source_y)
	{
		observed_target[0] = target_x;
		observed_target[1] = target_y;
		observed_source[0] = source_x;
		observed_source[1] = source_y;
	}

	template <typename T>
	bool operator()(const T* poseTarget,
		const T* poseSource,
		T* residuals) const
	{
		T projTarget[2], projSource[2];

		XfmApply(poseTarget, observed_target[0], observed_target[1], projTarget[0], projTarget[1]);
		XfmApply(poseSource, observed_source[0], observed_source[1], projSource[0], projSource[1]);

		residuals[0] = projTarget[0] - projSource[0];
		residuals[1] = projTarget[1] - projSource[1];
		return true;
	}

	static ceres::CostFunction* Create(const double& target_x, const double& target_y, const double& source_x, const double& source_y)
	{
		return (new ceres::AutoDiffCostFunction<Reprojection_Error, 2, 6, 6>(new Reprojection_Error(target_x, target_y, source_x, source_y)));
	}
	double observed_target[2], observed_source[2];
};

int main(int argc, char** argv)
{
	fs::path input_graph_path = argv[1];
	fs::path inputdir = input_graph_path.parent_path();
	string inputname = input_graph_path.filename().stem().string();
	fs::path output_graph_path = inputdir / (inputname + "_geoaware.txt");
	fs::path output_ceres_report = inputdir / (inputname + "_geoawarereport.txt");

	ifstream ifs(input_graph_path.string());
	motionavg::Affine2D::PoseGraph graph;
	ifs >> graph;
	ifs.close();


	// Declare optimize variables
	std::vector<double> posedata(graph.nodes.size() * 6);
	for (size_t i = 0; i < graph.nodes.size(); ++i)
		std::copy_n(graph.nodes[i].poseXfm, 6, posedata.data() + 6 * i);
	// Declare problem
	ceres::Problem problem;
	for (size_t i = 0; i < graph.edges.size(); ++i)
	{
		const auto& e = graph.edges[i];
		const auto& src = graph.nodes[e.source];
		const auto& tgt = graph.nodes[e.target];

		BBox bb_target, bb_source;
		bb_target.ptMin[0] = tgt.left;
		bb_target.ptMin[1] = tgt.top;
		bb_target.ptMax[0] = tgt.right;
		bb_target.ptMax[1] = tgt.bottom;

		bb_source.ptMin[0] = src.left;
		bb_source.ptMin[1] = src.top;
		bb_source.ptMax[0] = src.right;
		bb_source.ptMax[1] = src.bottom;

		bb_target = bb_target.transform(tgt.poseXfm);
		bb_source = bb_source.transform(src.poseXfm);

		BBox bb_intsec = bb_target.intersect(bb_source);

		double invSrc[6], invTgt[6];
		XfmInv(src.poseXfm, invSrc); XfmInv(tgt.poseXfm, invTgt);
		BBox bb_intsec_prj2src = bb_intsec.transform(invSrc);
		BBox bb_intsec_prj2tgt = bb_intsec_prj2src.transform(e.regXfm);
		BBox bb_intsec_prj2tgt_check = bb_intsec.transform(invTgt);

		{
			double srcPt[2], tgtPt[2];
			srcPt[0] = bb_intsec_prj2src.ptMin[0];
			srcPt[1] = bb_intsec_prj2src.ptMin[1];
			XfmApply(e.regXfm, srcPt[0], srcPt[1], tgtPt[0], tgtPt[1]);
			ceres::CostFunction* cost_function = Reprojection_Error::Create(tgtPt[0], tgtPt[1],
				srcPt[0], srcPt[1]);
			problem.AddResidualBlock(cost_function,
				new ceres::HuberLoss(1),
				posedata.data() + e.target * 6, posedata.data() + e.source * 6);
		}
		{
			double srcPt[2], tgtPt[2];
			srcPt[0] = bb_intsec_prj2src.ptMin[0];
			srcPt[1] = bb_intsec_prj2src.ptMax[1];
			XfmApply(e.regXfm, srcPt[0], srcPt[1], tgtPt[0], tgtPt[1]);
			ceres::CostFunction* cost_function = Reprojection_Error::Create(tgtPt[0], tgtPt[1],
				srcPt[0], srcPt[1]);
			problem.AddResidualBlock(cost_function,
				new ceres::HuberLoss(1),
				posedata.data() + e.target * 6, posedata.data() + e.source * 6);
		}
		{
			double srcPt[2], tgtPt[2];
			srcPt[0] = bb_intsec_prj2src.ptMax[0];
			srcPt[1] = bb_intsec_prj2src.ptMin[1];
			XfmApply(e.regXfm, srcPt[0], srcPt[1], tgtPt[0], tgtPt[1]);
			ceres::CostFunction* cost_function = Reprojection_Error::Create(tgtPt[0], tgtPt[1],
				srcPt[0], srcPt[1]);
			problem.AddResidualBlock(cost_function,
				new ceres::HuberLoss(1),
				posedata.data() + e.target * 6, posedata.data() + e.source * 6);
		}
		{
			double srcPt[2], tgtPt[2];
			srcPt[0] = bb_intsec_prj2src.ptMax[0];
			srcPt[1] = bb_intsec_prj2src.ptMax[1];
			XfmApply(e.regXfm, srcPt[0], srcPt[1], tgtPt[0], tgtPt[1]);
			ceres::CostFunction* cost_function = Reprojection_Error::Create(tgtPt[0], tgtPt[1],
				srcPt[0], srcPt[1]);
			problem.AddResidualBlock(cost_function,
				new ceres::HuberLoss(1),
				posedata.data() + e.target * 6, posedata.data() + e.source * 6);
		}
	}
	// Solve problem
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::ITERATIVE_SCHUR;
	options.minimizer_progress_to_stdout = true;

	bool FREENET = false;
	if (!FREENET)
	{
		// Choose fixed frame
		int fixed_frame_id = 0;

		// find 20211114_013021_ssc18d2_0015_basic_analytic_or
		for (int i = 0; i < graph.nodes.size(); ++i)
			if (graph.nodes[i].name == "20211114_013021_ssc18d2_0015_basic_analytic_or")
			{
				fixed_frame_id = i; break;
			}
		cout << "Found fixed frame Node: " << fixed_frame_id << " Name: " << graph.nodes[fixed_frame_id].name << endl;

		std::vector<int> sspv;
		for (int i = 0; i < 6; ++i)
			sspv.push_back(i);
		auto cssp = new ceres::SubsetParameterization(6, sspv);
		problem.SetParameterization(posedata.data() + 6 * fixed_frame_id, cssp);
	}

	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	std::cout << summary.FullReport() << endl;

	// Write out
	for (size_t i = 0; i < graph.nodes.size(); ++i)
		std::copy_n(posedata.data() + 6 * i, 6, graph.nodes[i].poseXfm);
	std::ofstream ofs(output_graph_path.string());
	ofs << graph;
	ofs.close();

	ofs.open(output_ceres_report.string());
	ofs << summary.FullReport() << endl;
	ofs.close();

	cout << "Write to " << output_graph_path << endl;
}