#include "Motion_Average.h"

MotionAverage_T::MotionAverage_T() {
  _num_max_Iter = 50;
  _Node_Para_list.resize(0);
  _RefNode_List.resize(0);
  _epslon = 1e-6;
  _Relax_Factor = 1;
  _num_para_per_node = 0;  // 0 means non intialized.
  _Con_Graph = nullptr;
  _Unary_Graph = nullptr;
}

MotionAverage_T::~MotionAverage_T() {}

bool MotionAverage_T::Initialize(std::vector<std::vector<double>> Node_Para_list, int num_para_per_node, std::vector<std::vector<ConG_Unit>>& Con_Graph) {
  _NodeParaInvCovMatrix.resize(Node_Para_list.size());
  _Node_Para_list.assign(Node_Para_list.begin(), Node_Para_list.end());
  _num_para_per_node = num_para_per_node;
  _Con_Graph = &Con_Graph;

  for (int p = 0; p < _Node_Para_list.size(); p++) {
    if (_Node_Para_list[p].size() < _num_para_per_node) {
      std::cout << "warning: defined per node dimension does not match the initialized node list, initialization failed" << std::endl;
      return false;
    }

  }  // end for p.
  // if the size of each node is bigger than the defined dimension, use the first a few.

  // here initialize per_node CoVarMatrix
  // this can be updated by _LocalSolver
  std::vector<double> temInvCovIdentity(_num_para_per_node * _num_para_per_node, 0);
  for (int p = 0; p < _num_para_per_node; p++) {
    temInvCovIdentity[p * _num_para_per_node + p] = 1;
  }  // end for p.

  for (int p = 0; p < _Node_Para_list.size(); p++) {
    _NodeParaInvCovMatrix[p].assign(temInvCovIdentity.begin(), temInvCovIdentity.end());
  }

  _update_ref_node_status();

  // check and initialize the Con_Graph CovMatrix.

  std::vector<double> temInvCov(_num_para_per_node * _num_para_per_node, 0);

  for (int p = 0; p < _num_para_per_node; p++) temInvCov[p * _num_para_per_node + p] = 1;

  bool to_update_pairwise_Cov = false;
  for (int p1 = 0; p1 < (*_Con_Graph).size(); p1++) {
    for (int p2 = 0; p2 < (*_Con_Graph)[p1].size(); p2++) {
      if ((*_Con_Graph)[p1][p2].InvCovMatrix.size() != _num_para_per_node * _num_para_per_node) {
        // this invalidates the use of the invconvariance matrix.
        to_update_pairwise_Cov = true;

      }  // end if
    }
  }

  // turn this to indentity matrix.
  if (to_update_pairwise_Cov) {
    for (int p1 = 0; p1 < (*_Con_Graph).size(); p1++)
      for (int p2 = 0; p2 < (*_Con_Graph)[p1].size(); p2++) (*_Con_Graph)[p1][p2].InvCovMatrix.assign(temInvCov.begin(), temInvCov.end());
  }

  return true;
}

bool MotionAverage_T::SetUnaryGraph(std::vector<std::vector<UnaryG_Unit>>& Unary_Graph) {
  this->_Unary_Graph = &Unary_Graph;
  return true;
}

#include <fstream>
#include <string>
using namespace std;
void write_nodes(std::string filepath, const std::vector<std::vector<double>>& nodes) {
  std::ofstream ofs(filepath, std::ios::trunc | std::ios::binary);
  size_t num = nodes.size();
  ofs.write(reinterpret_cast<const char*>(&num), sizeof(size_t));
  for (size_t i = 0; i < num; ++i) ofs.write(reinterpret_cast<const char*>(nodes[i].data()), sizeof(double) * nodes[i].size());
  ofs.close();
}

/**************************************************************************************************************************************************
Algorithm:
While convergence not met & within max_iter
For each node
Compute the relative contribution to this node through the _LocalSolver
End
Check Convergence or max_iter
End
Note: this version is not parallelized. A parallel version will be numerically solver by concept but practically faster
****************************************************************************************************************************************************/
bool MotionAverage_T::IterSolver_Run() {
  double relax_factor = 1;                                // this is the relaxation factor for the datum
  int IterCount = 0;                                      // iterator for solver.
  double ThisError = std::numeric_limits<double>::max();  // current error between two lines
  double CumError_last_iter = std::numeric_limits<double>::max();

  if (_RefNode_List.size() <= 0) {
    _Detecting_Models();  // this fills the _RefNode_List.
  }

  std::vector<std::vector<double>> Node_Para_list_last_iter;  // temporary variable to store the last iteration
  Node_Para_list_last_iter.assign(_Node_Para_list.begin(), _Node_Para_list.end());
  _Node_Para_list_buf.assign(_Node_Para_list.begin(), _Node_Para_list.end());

  while (IterCount < _num_max_Iter) {
    for (int p = 0; p < _Node_Para_list.size(); p++) {
      if (std::find(_RefNode_List.begin(), _RefNode_List.end(), p) != _RefNode_List.end()) continue;
      _LocalSolver(p);  // for each node, this overloading function updates _Node_Para_list_buf
      _Node_Para_list[p].assign(_Node_Para_list_buf[p].begin(), _Node_Para_list_buf[p].end());

    }  //

    // write_nodes("N:\\dump_graph\\iter" + to_string(IterCount) + ".bin", _Node_Para_list);

    // check the convergence, update the relaxation factor.
    double CumError = 0;
    double TemThisError = 0;

    for (int p = 0; p < _Node_Para_list.size(); p++) {
      _Compute_Error_Between_Iterations(_Node_Para_list[p], Node_Para_list_last_iter[p], p, TemThisError);
      CumError += pow(TemThisError, 2);
    }

    CumError /= (double)_Node_Para_list.size();

    if (CumError < _epslon) break;

    // here update the relaxation factor.

    if (CumError > CumError_last_iter) {
      // enhance the tie from the reference.
      _Relax_Factor *= 2;
      _Node_Para_list.assign(Node_Para_list_last_iter.begin(), Node_Para_list_last_iter.end());
    } else {
      // relax the tie to let the network fluctuate
      _Relax_Factor = MAX(_Relax_Factor / 2, 1);
      Node_Para_list_last_iter.assign(_Node_Para_list.begin(), _Node_Para_list.end());
      CumError_last_iter = CumError;
    }

    IterCount++;

    std::cout << "Iter " << IterCount << " :" << CumError << std::endl;

  }  // end for while

  return true;
}

/**************************************************************************************************************************************************
Algorithm:
While convergence not met & within max_iter
For each node
Compute the relative contribution to this node through the _LocalSolver
End
Check Convergence or max_iter
End
Note: this is the parallelized version, Jacobi concept, do after all done.
****************************************************************************************************************************************************/
bool MotionAverage_T::IterSolver_Run_Parallel() { return true; }

/***********************************************************************************************************************
Take the average of the covariance matrix for each neighborhood
this is an optional function for pair-wise solvers.
note: this is an approximation.
************************************************************************************************************************/
bool MotionAverage_T::_Aggregate_PairWise_InvCov() {
  // check if the pair-wise cov is valid.

  for (int p1 = 0; p1 < (*_Con_Graph).size(); p1++) {
    for (int p2 = 0; p2 < (*_Con_Graph)[p1].size(); p2++) {
      if ((*_Con_Graph)[p1][p2].InvCovMatrix.size() != _num_para_per_node * _num_para_per_node) {
        std::cout << "Cannot aggregate covMatrix, because the pair_wise Inv_cov is not initialized, Node level InvCov remain identity matrix" << std::endl;
        return false;
      }
    }  // end for p2
  }    // end for p1

  for (int p1 = 0; p1 < _Node_Para_list.size(); p1++) {
    std::vector<double> temMatrixSum(_num_para_per_node * _num_para_per_node, 0);

    //
    for (int p2 = 0; p2 < (*_Con_Graph)[p1].size(); p2++) {
      for (int pp = 0; pp < _num_para_per_node * _num_para_per_node; pp++) {
        temMatrixSum[pp] += (*_Con_Graph)[p1][p2].InvCovMatrix[pp];
      }  // end for pp
    }    //

    double temcumsum = 0;
    for (int pp = 0; pp < _num_para_per_node * _num_para_per_node; pp++) {
      temMatrixSum[pp] /= (double)(*_Con_Graph)[p1].size();
      temcumsum += temMatrixSum[pp];
    }  // end for pp

  }  // end for p1

  _Normalize_Per_Node_InvCov_Matrix();  // make the sum of the inv_cov valuesone;

  return true;
}
/*********************************************************************************************
This function normalizes the InvCov Matrix stored in _NodeParaInvCovMatrix
**********************************************************************************************/
void MotionAverage_T::_Normalize_Per_Node_InvCov_Matrix() {
  for (int p = 0; p < _NodeParaInvCovMatrix.size(); p++) {
    if (_NodeParaInvCovMatrix[p].size() != _num_para_per_node * _num_para_per_node) {
      std::cout << "warning: the invCov is not appropriately filled" << std::endl;
      return;
    }

  }  // end for p

  for (int p = 0; p < _NodeParaInvCovMatrix.size(); p++) {
    double temcumsum = 0;

    for (int pp = 0; pp < _num_para_per_node * _num_para_per_node; pp++) {
      temcumsum += _NodeParaInvCovMatrix[p][pp];
    }  // end for pp

    for (int pp = 0; pp < _num_para_per_node * _num_para_per_node; pp++) {
      _NodeParaInvCovMatrix[p][pp] /= temcumsum;
    }  // end for pp
  }    // end for p

}  //
/***************************************************************************************************
given the _ref_node_list, create an index array indicating whether a node is a reference node or not
this fills _Ref_Node_Status
****************************************************************************************************/
void MotionAverage_T::_update_ref_node_status() {
  if (_RefNode_List.size() <= 0) {
    _Detecting_Models();
  }

  int NumNodes = _Node_Para_list.size();

  _Ref_Node_Status.resize(NumNodes);
  for (int pp = 0; pp < NumNodes; pp++) _Ref_Node_Status[pp] = 0;

  for (int p = 0; p < _RefNode_List.size(); p++) {
    _Ref_Node_Status[_RefNode_List[p]] = 255;
  }
}
/********************************************************************************************
this overloading function is used to compute the errors between iterations.
*********************************************************************************************/
void MotionAverage_T::_Compute_Error_Between_Iterations(std::vector<double> Node_para_1, std::vector<double> Node_para_2, int index_i, double& ThisError) {
  ThisError = 0;

  std::vector<double> tem_row_vector(_num_para_per_node, 0);

  for (int p = 0; p < _num_para_per_node; p++) {
    for (int pp = 0; pp < _num_para_per_node; pp++) {
      tem_row_vector[p] += _NodeParaInvCovMatrix[index_i][pp] * (Node_para_1[pp] - Node_para_2[pp]);
    }

  }  // end for p

  for (int pp = 0; pp < _num_para_per_node; pp++) {
    ThisError += tem_row_vector[pp] * (Node_para_1[pp] - Node_para_2[pp]);
  }  // end for pp
}  //

/*************************************************************************************************************************
Connected component analysis for get clusters of models based on the co-graph
This function will fill the _RefNode_List, it assumes the link of the nodes are successful.
This function needs to be tested.
**************************************************************************************************************************/
bool MotionAverage_T::_Detecting_Models() {
  if (_RefNode_List.size() != 0) {
    std::cout << "Reference Node Already defined" << std::endl;
    return true;
  }

  int num_nodes = _Node_Para_list.size();

  if (num_nodes <= 0) {
    std::cout << "Node list uninitialized" << std::endl;
    return false;
  }

  std::vector<double> tem_Connectivity_Count(num_nodes, 0);             // temporary vector for ranking connectivities
  std::vector<unsigned int> tem_Connectivity_Count_Rank(num_nodes, 0);  // index for ranking
  std::vector<int> tem_visited_ind(num_nodes, -1);                      // indicating if this node belongs to a cluster, value following the cluster ID defined by _Ref_Node.
  std::vector<unsigned char> tem_neighbors_visited(num_nodes, 0);       // indicating if the neighborhood of a node has been all visited.

  // fill tem_Connectivity_Count and tem_Connectivity_Count_Rank
  for (int p = 0; p < num_nodes; p++) {
    tem_Connectivity_Count[p] = (double)(*_Con_Graph)[p].size();
    tem_Connectivity_Count_Rank[p] = p;
  }  // end for p

  sort2(num_nodes, tem_Connectivity_Count.data(), tem_Connectivity_Count_Rank.data());  // from smallest to largest.

  int tem_cluster_ID = 0;
  int tem_candidate_ref_node = num_nodes - 1;
  _RefNode_List.push_back(tem_Connectivity_Count_Rank[tem_candidate_ref_node]);
  tem_visited_ind[tem_Connectivity_Count_Rank[tem_candidate_ref_node]] = tem_cluster_ID;

  std::vector<int> tem_added_nodes;  // record for each cycle the added nodes,a criterion to stop a cluster is that no elements are added in it.
  tem_added_nodes.push_back(tem_candidate_ref_node);

  while (1) {
    int tem_new_clustered_node_count = 0;  // indicating if new nodes are found.
    for (int p1 = 0; p1 < tem_added_nodes.size(); p1++) {
      if (tem_neighbors_visited[tem_added_nodes[p1]]) continue;

      bool is_neighbors_visited = true;

      for (int p2 = 0; p2 < (*_Con_Graph)[tem_added_nodes[p1]].size(); p2++) {
        if (tem_visited_ind[(*_Con_Graph)[tem_added_nodes[p1]][p2].OtherID] == -1) {
          tem_added_nodes.push_back((*_Con_Graph)[tem_added_nodes[p1]][p2].OtherID);
          tem_visited_ind[(*_Con_Graph)[tem_added_nodes[p1]][p2].OtherID] = tem_cluster_ID;
          is_neighbors_visited = false;
          tem_new_clustered_node_count++;

        }  // end if

        if (is_neighbors_visited) tem_neighbors_visited[tem_added_nodes[p1]] = 255;

      }  // end for p2

    }  // end for p1

    if (tem_new_clustered_node_count == 0)  // the current cluster is done
    {
      tem_added_nodes.resize(0);

      tem_candidate_ref_node--;

      while (tem_candidate_ref_node >= 0)  // checking for another reference node which was not visited yet.
      {
        if (tem_visited_ind[tem_Connectivity_Count_Rank[tem_candidate_ref_node]] == -1)  // not visited yet
        {
          tem_added_nodes.push_back(tem_Connectivity_Count_Rank[tem_candidate_ref_node]);
          _RefNode_List.push_back(tem_Connectivity_Count_Rank[tem_candidate_ref_node]);
          tem_cluster_ID++;
          tem_visited_ind[tem_Connectivity_Count_Rank[tem_candidate_ref_node]] = tem_cluster_ID;
          break;
        }
        tem_candidate_ref_node--;
      }

      if (tem_added_nodes.size() <= 0)  // fail to find another reference node.
      {
        break;  // while loop finished and reference nodes are found.
      }         // end if tem_added_nodes.size

    }  // end if no points are added.

  }  // end while (1)

  _update_ref_node_status();

  return true;
}