#ifndef _MOTION_AVERAGE
#define _MOTION_AVERAGE
#include <iostream>
#include <vector>

#include "Mathapi.h"

#ifndef MAX
#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b)
#endif

/***************************************************************************
A general class template for motion average algorithm.
This class only deal with numerical computations
--
***************************************************************************/

/*****************************************************************************
basic unit to build the connectivity graph, other ConGUnit can inherit from it
for specific applications
******************************************************************************/
class ConG_Unit {
 public:
  ConG_Unit() {
    LinkWeight = 1;
    RefID = -1;
    OtherID = -1;
    reserved1 = -1;  //  this could customized for different purposes, for example, processed or unprocessed.
    reserved2 = -1;  //  ditto
  }
  ~ConG_Unit() {}

  double LinkWeight;  // mostly kept one.
  int RefID;
  int OtherID;
  int reserved1;
  int reserved2;
  std::vector<double> RelParaVec;    // this vector can store P_ij
  std::vector<double> InvCovMatrix;  // this covariance matrix is problem specific, _LocalSolver defines how this should be used.
};
/***************************************************************************************************
Concept & Intro: The motion average implements a solution that use either direct, or iterative approach to resolve
the BA problem for parameter estimation, to allow featureless estimation (featureless in global estimation). It is especially useful for parametric estimation
of linear problems (parameter linear with respect to observations). The concept of featureless will avoid complex formulation of
multi-ray points to build correlations/constraints across the estimated parameters.

Problem: Given a Graph G(V,E), we estimate for each node V_i, its unknown parameters P_i need to be estimated, in order to achieve the
object function of F(P_i,x_i,X) = 0, where x_i are set of node specific measurements, X is a set of global measurements (unknown, but auxiliary),
as a result of F(P_i,x_i,X) = 0 (when P_i is resolved). So X is essentially used to build the correlation among the nodes.

For Gauss-Newton (GN) Method, X can be eliminated in the normal matrix for each iteration. But, when such correlation can be explicitly built up
as  F_ij(P_i,P_j,x_i,x_j), or more, for example F_ijk(P_i,P_j,P_k,x_i,x_j,x_k) etc. The elimination are not need in the normal matrix for the G-N,
as a result, the design matrix can be built in a more structured way to allow the following possible solutions:

        1) Featureless solution: if F_ij(P_i,P_j,x_i,x_j) = 0 can be further translated to F_ij(P_i,P_j,P_ij,cov(x_i,x_j)) = 0 - possible for linear parameter P_i and P_j
        which are linearly invertible. P_ij encodes the relative mapping between x_i and x_j, often estimated directly between x_i and x_j between two pairs, and cov(x_i,x_j) can even be approximated.
        In this case, even cov(x_i,x_j) can be approximated. So far, it is proven to work with linear parameters (first order polynomial).

        2) Iterative Solution: The F_ij(P_i,P_j,x_i,x_j) constrain allows a good structure to perform a fix-point iterative algorithm, which does not require an explicit matrix solver. For zero order,
        Equivalent to  Gauss-Seidel, the matrix structure of the normal matrix guarantees the convergence (Diagonally Dominant).Extremely memory-friendly,and conceptually has higher numerical
precision , for linear, it does not require For non-linear F_ij(P_i,P_j,x_i,x_j), unproven, should have high chance of convergence,Although the norm of jacobian may not be bounded, A factor n could be
added to the datum to for a contraction mapping to get obtain a Cauchy sequence.

Input:
      Parameterset: std::vector<std::vector<double>> _Node_Para_list, n nodes, for each node the parameter sets are presented as a std::vector<double>
          Graph G(V,E): std::vector<std::vector<ConG_Unit>> _Con_Graph, this allows a sparse structure, you can implement many low-level things by inhert from this ConG_Unit class.

Output: The updated _Node_Para_list


****************************************************************************************************/
class MotionAverage_T {
 public:
  MotionAverage_T();
  ~MotionAverage_T();

  /*************************************************************************************
  Setting necessary parameters - better to do before initialization, add as needed
  **************************************************************************************/
  void Set_Max_Iter_Num(int num_iter) { _num_max_Iter = num_iter; }
  void Set_ConvergenceEpslon(double epslon) { _epslon = epslon; }
  void Set_Reference_Nodes(std::vector<int> Ref_node_list) { _RefNode_List.assign(Ref_node_list.begin(), Ref_node_list.end()); }
  /***************************************************
  Core Computing, functions to drive the computation.
  ****************************************************/

  /***************************************************************************************
  initialize, Node_Para_list should have the initial value, Con_Graph is the pre-built connectivity graph
  this function can be overloaded to introduce more internal variables for LocalSolver.
  ****************************************************************************************/
  bool virtual Initialize(std::vector<std::vector<double>> Node_Para_list, int num_para_per_node, std::vector<std::vector<ConG_Unit>>& Con_Graph);

  /***************************************************************************************
  This virtual function needs to be overload to the specific problem.
  ****************************************************************************************/
  virtual bool DirectSolver_Run() { return true; };

  /****************************************************************************************
  This iteratively runs the inference function to perform the parameter estimation, no overloading
  This largely depends on the _local_solver

  Algorithm:

    While convergence not met & within max_iter
            For each node
                   Compute the relative contribution to this node through the _LocalSolver
            End
                  Check Convergence or max_iter
    End
  ****************************************************************************************/
  bool IterSolver_Run();

  /*********************************************************************************************
  The parallel version of the above, use Jacoby concept.
  theoretically slower but practically faster (due parallelism), but you need to check the synchronization in _LocalSolver.
  **********************************************************************************************/
  bool IterSolver_Run_Parallel();

 protected:
  /***************************************************************************************
  _LocalSolver is a core function that does the following:

          Given a P_i from Node_Para_list, it goes over its neighborhood,and provide an update.

          This implements F_ij(P_i,P_j,x_i,x_j), or F_ijk...n(P_i,P_j,P_k,...,P_n,x_i,x_j,x_k,...,x_n)
      It should be problem specific. All working parameter should be stored through a space within this class.
          This updates the _Node_Para_list_buf, but not the _Node_Para_list.
  ****************************************************************************************/
  virtual bool _LocalSolver(int index_i) { return true; };

  /**************************************************************************************************
  This is an overload function that computes the error metric between iterations.
  index_i is reserved to call out CovMatrix stored in the Inherited variables
  covMatrix is optional to be used to compute the Mahalanobis distance, one may have n covariance matrix.
  ****************************************************************************************************/
  void _Compute_Error_Between_Iterations(std::vector<double> Node_para_1, std::vector<double> Node_para_2, int index_i, double& ThisError);

  /***************************************************************************************************
  This function detect clusters using the connectivity graph, to set up the datum.
  This will fill the _RefNode_List to anchor their value from what is given as the initial.
  ****************************************************************************************************/
  bool _Detecting_Models();

  /*****************************************************************************************************
  This function aggregates the covariance matrix using the pair-wise covariance matrix -- taking the average
  *****************************************************************************************************/
  bool _Aggregate_PairWise_InvCov();

 private:
  void _Normalize_Per_Node_InvCov_Matrix();
  void _update_ref_node_status();

 public:
  std::vector<std::vector<double>> _Node_Para_list;      // The final solution, get access to it.
  std::vector<std::vector<double>> _Node_Para_list_buf;  // This gives a copy of the final solution during iterative processes.

 protected:
  int _num_para_per_node;          // this is to regularize the inconsistences among different initialization
  int _num_max_Iter;               // max number of iteration for BP
  double _Relax_Factor;            // Weight to anchor the neighbors between the reference and the
  std::vector<int> _RefNode_List;  // this is the reference Node that define the datum, its parameter will never change, if it has multiple, the network is not fully connected and has several models.
                                   // if this is empty (not initialized externally),the program will initialize this using the Con_Graph.
  double _epslon;
  std::vector<std::vector<ConG_Unit>>* _Con_Graph;
  std::vector<std::vector<double>> _NodeParaInvCovMatrix;  // this is the Node_level CovarianceMatrix,describing the relative importance of the parameters
  std::vector<unsigned char> _Ref_Node_Status;             // dimension the same as the number of nodes,255 means ref node and 0 means non-ref node, for weighting purpose.
};

#endif