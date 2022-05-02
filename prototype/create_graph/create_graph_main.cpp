#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>

#include "gdal.h"

#include "myutils.h"
#include "mystructs.h"

#include <mkl.h>
#include <mkl_solvers_ee.h>

using namespace std;
namespace fs = std::filesystem;
using namespace motionavg::Affine2D;

int main(int argc, char** argv)
{
	GDALAllRegister();
	string rootdirpath = argv[1];
	string posegraphfilepath = (fs::path(rootdirpath) / "posegraph.txt").string();
	string bundlegraphfilepath = (fs::path(rootdirpath) / "bundlegraph.txt").string();

	// --------------------------------
	// Collect necessary data from Dr. Qin's output
	cout << "Searching in directory ... [" << rootdirpath << "]"<<endl;
	auto rootdir = fs::path(rootdirpath);

	std::vector<fs::path> imglist;
	std::vector<fs::path> paramlist;
	
	for (const auto & f : fs::recursive_directory_iterator(rootdir))
	{
		if (f.is_directory()) continue;
		if(isImg(f)) imglist.push_back(f);
		else if(isParamFile(f)) paramlist.push_back(f);
	}

	cout << "Found " << imglist.size() << " Images and " << paramlist.size() << " parameters" << endl;
	cout << "Build Graph ..." << endl;


	std::unordered_map<string, fs::path> imgindex;
	for (const auto& p : imglist)
		imgindex[p.filename().stem().string()] = p;
	cout << "Img name index: " << imgindex.size() << endl;
	if (imgindex.size() != imglist.size())
		cout << "Warning: images with same name detected!\n";

	std::unordered_map<string, PairwiseData> paramindex;
	for(const auto& p: paramlist)
	{
		PairwiseData model;
		model.load(p.string());
		paramindex[model.affine.filename] = model;
	}
	cout << "Param name index: " << paramindex.size() << endl;
	if (paramindex.size() != paramindex.size())
		cout << "Warning: transformations with same name detected!\n";
	// --------------------------------------------
	// Build Pose Graph
	PoseGraph pgraph;
	pgraph.basepath = rootdir.string();

	std::unordered_map<string, size_t> imgname_to_nid;
	
	for (const auto& kv : paramindex)
	{
		string name = kv.second.affine.filename;
		string target_img = kv.second.affine.target_file;
		string source_img = kv.second.affine.source_file;

		int target_nid=-1, source_nid=-1;
		if (imgname_to_nid.find(target_img) != imgname_to_nid.end())
		{
			target_nid = imgname_to_nid[target_img];
		}
		else {
			PoseGraph::Node n;
			n.name = target_img;
			n.path = fs::relative(imgindex[target_img], rootdir).string();
			GDALDatasetH pds= GDALOpen(imgindex[target_img].string().c_str(), GA_ReadOnly);
			n.left = 0; n.right = GDALGetRasterXSize(pds);
			n.top = 0; n.bottom = GDALGetRasterYSize(pds);
			double pGeoTrans[6];
			GDALGetGeoTransform(pds, pGeoTrans);
			std::copy_n(pGeoTrans, 6, n.poseXfm);
			GDALClose(pds);
			target_nid = pgraph.insertNode(n);
			imgname_to_nid[target_img] = target_nid;

		}

		if (imgname_to_nid.find(source_img) != imgname_to_nid.end())
		{
			source_nid = imgname_to_nid[source_img];
		}
		else {
			PoseGraph::Node n;
			n.name = source_img;
			n.path = fs::relative(imgindex[source_img], rootdir).string();
			GDALDatasetH pds = GDALOpen(imgindex[source_img].string().c_str(), GA_ReadOnly);
			n.left = 0; n.right = GDALGetRasterXSize(pds);
			n.top = 0; n.bottom = GDALGetRasterYSize(pds);
			double pGeoTrans[6];
			GDALGetGeoTransform(pds, pGeoTrans);
			std::copy_n(pGeoTrans, 6, n.poseXfm);
			GDALClose(pds);

			source_nid = pgraph.insertNode(n);
			imgname_to_nid[source_img] = source_nid;
		}

		PoseGraph::Edge e;
		e.name = name;
		e.target = target_nid;
		e.source = source_nid;
		// GeoParam is a adjustment for Right Pose. The relative position should be:
		const PairwiseData& pairdata = kv.second;
		const AffineModel& affinemodel = pairdata.affine;
		double originSourceGeo[6], resampleSourceGeo[6], composeSourceGeo[6];
		std::copy_n(pgraph.nodes[source_nid].poseXfm, 6, originSourceGeo);
		std::copy_n(originSourceGeo, 6, resampleSourceGeo);
		XfmApply(affinemodel.geoParam, originSourceGeo[0] - affinemodel.pixelParam[0], originSourceGeo[3] + affinemodel.pixelParam[3], resampleSourceGeo[0], resampleSourceGeo[3]);
		XfmComposite(resampleSourceGeo, affinemodel.pixelParam, composeSourceGeo);	// This is the final affine transformation from pixel location to registered UTM

		double invTargetGeo[6];
		double invTargetSourceGeo[6];
		XfmInv(pgraph.nodes[target_nid].poseXfm, invTargetGeo);
		XfmComposite(invTargetGeo, composeSourceGeo, invTargetSourceGeo);
		std::copy_n(invTargetSourceGeo, 6, e.regXfm);

		//double _target_tp[2], _source_tp[2], _reg_source_tp[2], _rep_tp[2], _tp[2];
		//for (auto& tp : pairdata.tiepoints)
		//{
		//	XfmApply(pgraph.nodes[target_nid].poseXfm, tp.target[0], tp.target[1], _target_tp[0], _target_tp[1]);
		//	XfmApply(pgraph.nodes[source_nid].poseXfm, tp.source[0], tp.source[1], _source_tp[0], _source_tp[1]);
		//	XfmApply(composeSourceGeo, tp.source[0], tp.source[1], _reg_source_tp[0], _reg_source_tp[1]);
		//	XfmApply(invTargetSourceGeo, tp.source[0], tp.source[1], _rep_tp[0], _rep_tp[1]);
		//	std::copy_n(tp.target, 2, _tp);
		//}
		
		pgraph.insertEdge(e);

		cout << "Pose Edge: " << e.name << endl;
		cout << "  From " << e.target << " [" << pgraph.nodes[target_nid].name << "] to " << e.source << " [" << pgraph.nodes[source_nid].name << "]" << endl;
	}

	ofstream ofs(posegraphfilepath, ios::binary | std::ios::out);
	ofs << pgraph;
	ofs.close();

	//PoseGraph graph2;
	//ifstream ifs(graphfilepath, ios::binary | std::ios::in);
	//ifs >> graph2;

	//ofs.open("N:\\g.txt");
	//ofs << graph2;
	//ofs.close();

	cout << "Write to " << posegraphfilepath << endl;
	// --------------------------------------------------
	// Build BundleGraph
	BundleGraph bgraph;
	bgraph.basepath = rootdir.string();

	imgname_to_nid.clear();

	for (const auto& kv : paramindex)
	{
		string name = kv.second.affine.filename;
		string target_img = kv.second.affine.target_file;
		string source_img = kv.second.affine.source_file;

		int from_nid = -1, to_nid = -1;
		if (imgname_to_nid.find(target_img) != imgname_to_nid.end())
		{
			from_nid = imgname_to_nid[target_img];
		}
		else {
			BundleGraph::Node n;
			n.name = target_img;
			n.path = fs::relative(imgindex[target_img], rootdir).string();
			GDALDatasetH pds = GDALOpen(imgindex[target_img].string().c_str(), GA_ReadOnly);
			n.left = 0; n.right = GDALGetRasterXSize(pds);
			n.top = 0; n.bottom = GDALGetRasterYSize(pds);
			double pGeoTrans[6];
			GDALGetGeoTransform(pds, pGeoTrans);
			std::copy_n(pGeoTrans, 6, n.poseXfm);
			GDALClose(pds);
			from_nid = bgraph.insertNode(n);
			imgname_to_nid[target_img] = from_nid;
		}

		if (imgname_to_nid.find(source_img) != imgname_to_nid.end())
		{
			to_nid = imgname_to_nid[source_img];
		}
		else {
			BundleGraph::Node n;
			n.name = source_img;
			n.path = fs::relative(imgindex[source_img], rootdir).string();
			GDALDatasetH pds = GDALOpen(imgindex[source_img].string().c_str(), GA_ReadOnly);
			n.left = 0; n.right = GDALGetRasterXSize(pds);
			n.top = 0; n.bottom = GDALGetRasterYSize(pds);
			double pGeoTrans[6];
			GDALGetGeoTransform(pds, pGeoTrans);
			std::copy_n(pGeoTrans, 6, n.poseXfm);
			GDALClose(pds);

			to_nid = bgraph.insertNode(n);
			imgname_to_nid[source_img] = to_nid;
		}

		BundleGraph::Edge e;
		e.name = name;
		e.target = from_nid;
		e.source = to_nid;
		const BundleGraph::Node& fn = bgraph.nodes[from_nid];
		const BundleGraph::Node& tn = bgraph.nodes[to_nid];
		const PairwiseData& pairdata = kv.second;
		e.tiepoints.resize(pairdata.tiepoints.size());
		/*for (size_t i = 0; i < pairdata.tiepoints.size(); ++i)
		{
			const TiePoint& tp_img = pairdata.tiepoints[i];
			TiePoint& tp_geo = e.tiepoints[i];
			XfmApply(fn.poseXfm, tp_img.target[0], tp_img.target[1], tp_geo.target[0], tp_geo.target[1]);
			XfmApply(tn.poseXfm, tp_img.source[0], tp_img.source[1], tp_geo.source[0], tp_geo.source[1]);
		}*/
		std::copy(pairdata.tiepoints.begin(), pairdata.tiepoints.end(), e.tiepoints.begin());

		bgraph.insertEdge(e);

		cout << "Bundle Edge: " << e.name << endl;
		cout << "  From " << e.target <<" to " << e.source << endl;
	}

	ofs.open(bundlegraphfilepath, ios::binary | std::ios::out);
	ofs << bgraph;
	ofs.close();

	//BundleGraph graph2;
	//ifstream ifs(bundlegraphfilepath, ios::binary | std::ios::in);
	//ifs >> graph2;

	//ofs.open("N:\\g.txt");
	//ofs << graph2;
	//ofs.close();
	return 0;
}


/*
* 
* snippet auto paring

		for (const auto& kv : imgindex) {
			string key = kv.first;
			size_t complen = key.size();
			if (key.size() > name.size()) continue;
			if (key.compare(0, complen, name, 0, complen)==0)
			{
				from_img = key;
				break;
			}
		}
		if (from_img.empty())
		{
			cout << "Warning: Edge fromnode " << name << " failed"<<endl;
			continue;
		}
		for (const auto& kv : imgindex) {
			string key = kv.first;
			size_t keylen = key.size();
			size_t complen = name.size() - from_img.size() - 1;
			//if (key.size() < name.size() - from_img.size() - 1) continue;
			auto ret = name.compare(from_img.size() + 1, complen, key.data(), 0, complen);
			if ( ret == 0)
			{
				to_img = key;
				break;
			}
		}
		if (to_img.empty())
		{
			cout << "Warning: Edge tonode " << name << " failed" << endl;
			continue;
		}
*/