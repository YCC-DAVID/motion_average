#include <iostream>
#include <fstream>
#include <ceres/problem.h>
#include <ceres/solver.h>
#include <ceres/local_parameterization.h>
#include <ceres/autodiff_cost_function.h>
#include <ceres/loss_function.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>


#include "myutils.h"
#include "mystructs.h"

using namespace std;
using namespace motionavg;

struct AdditiveWithPrior_PoseError {

	AdditiveWithPrior_PoseError(const double* data, const double* cov)
	{
		std::copy_n(data, 6, xfm);

		Eigen::Matrix<double, 6, 6> covMat(cov);
		Eigen::Matrix<double, 6, 6> PMat = covMat.inverse().sqrt();
		std::copy_n(PMat.data(), 36, P);
	}
	/*
	* rel_hat = inv(v0)*v1
	* xfm = observed(inv(v0)*reg*v1)
	* r = v0 - v0_hat
	*/
	template <typename T>
	bool operator()(const T* v0,
		const T* v1,
		T* residuals) const
	{
		T inv0[6];
		XfmInv(v0, inv0);
		T rel01[6];
		XfmComposite(inv0, v1, rel01);

		T res[6];
		for (int i = 0; i < 6; ++i)
			res[i] = rel01[i] - xfm[i];

		for (int i = 0; i < 6; ++i)
		{
			residuals[i] = T(0);
			for (int j = 0; j < 6; ++j)
			{
				residuals[i] += res[j] * P[i * 6 + j];
			}
		}
		return true;
	}

	static ceres::CostFunction* Create(const motionavg::PoseGraph::Edge& e)
	{
		return (new ceres::AutoDiffCostFunction<AdditiveWithPrior_PoseError, 6, 6, 6>(new AdditiveWithPrior_PoseError(e.regXfm, e.covXfm)));
	}
	double xfm[6]; // tx a11 a12 ty a21 a22
	double P[36]; // 6x6 weight matrix
};


//struct Compositional_PoseError {
//
//	Compositional_PoseError(const double* data)
//	{
//		std::copy_n(data, 6, xfm);
//	}
//	/*
//	* rel_hat = inv(v0)*v1
//	* xfm = observed(inv(v0)*reg*v1)
//	* r = inv(rel_hat)*xfm - eye(3)
//	*/
//	template <typename T>
//	bool operator()(const T* v0,
//		const T* v1,
//		T* residuals) const
//	{
//		T inv0[6];
//		XfmInv(v0, inv0);
//		T rel01[6];
//		XfmComposite(inv0, v1, rel01);
//
//		T irel01[6];
//		XfmInv(rel01, irel01);
//		XfmComposite(irel01, xfm, residuals);
//
//		residuals[1] -= T(1.);
//		residuals[5] -= T(1.);
//
//		return true;
//	}
//
//	static ceres::CostFunction* Create(const motionavg::PoseGraph::Edge& e)
//	{
//		return (new ceres::AutoDiffCostFunction<Compositional_PoseError, 6, 6, 6>(new Compositional_PoseError(e.regXfm)));
//	}
//	double xfm[6]; // tx a11 a12 ty a21 a22
//	double P[36]; // 6x6 weight matrix
//};

int main(int argc, char** argv)
{
	fs::path input_graph_path = argv[1];
	fs::path inputdir = input_graph_path.parent_path();
	string inputname = input_graph_path.filename().stem().string();
	fs::path output_graph_path = inputdir / (inputname + "_wpriorceres.txt");
	fs::path output_ceres_report = inputdir / (inputname + "_wpriorceresreport.txt");

	ifstream ifs(input_graph_path.string());
	motionavg::PoseGraph graph;
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
		ceres::CostFunction* cost_function = AdditiveWithPrior_PoseError::Create(e);
		problem.AddResidualBlock(cost_function,
			new ceres::HuberLoss(1),
			posedata.data() + e.target * 6, posedata.data() + e.source * 6);
	}
	// Solve problem
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	options.minimizer_progress_to_stdout = true;
	options.function_tolerance = 1e-14;
	options.gradient_tolerance = 1e-14;
	options.parameter_tolerance = 1e-14;
	options.max_num_iterations = 1000;

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
	return 0;
}