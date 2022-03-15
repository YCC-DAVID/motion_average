#include <gdal.h>
#include <cpl_string.h>
#include <gdalwarper.h>
#include <filesystem>
#include <iostream>
#include <fstream>

#include "myutils.h"
#include "mystructs.h"


using namespace std;
namespace fs = std::filesystem;
using namespace motionavg;

bool ensureDir(fs::path p)
{
	if (!fs::exists(p))
		return fs::create_directories(p);
	else
		return true;
}


void GDAL_Affine_TransformImage(GDALDatasetH& pdata1, const GDALDatasetH pdata2, double *srcXfm, double *dstXfm, bool isbilinear_or_nearest =1)
{
	GDALDereferenceDataset(pdata1);
	GDALDereferenceDataset(pdata2);

	//double backupXfm[6];
	//GDALGetGeoTransform(pdata1, backupXfm);

	//GDALSetGeoTransform(pdata1, srcXfm);


	GDALWarpOptions* psWarpOptions = GDALCreateWarpOptions();

	psWarpOptions->hSrcDS = pdata1;
	psWarpOptions->hDstDS = pdata2;
	psWarpOptions->nBandCount = GDALGetRasterCount(pdata1);
	psWarpOptions->panSrcBands =
		(int*)CPLMalloc(sizeof(int) * psWarpOptions->nBandCount);
	for (int p = 0; p < GDALGetRasterCount(pdata1); p++)
		psWarpOptions->panSrcBands[p] = p + 1;

	psWarpOptions->panDstBands =
		(int*)CPLMalloc(sizeof(int) * psWarpOptions->nBandCount);
	for (int p = 0; p < GDALGetRasterCount(pdata1); p++)
		psWarpOptions->panDstBands[p] = p + 1;

	psWarpOptions->pfnProgress = GDALTermProgress;
	if (isbilinear_or_nearest)
		psWarpOptions->eResampleAlg = GRA_Bilinear;
	else
		psWarpOptions->eResampleAlg = GRA_NearestNeighbour;

	GDALSetGeoTransform(pdata2, dstXfm);
	psWarpOptions->pTransformerArg = GDALCreateGenImgProjTransformer3("", srcXfm, "", dstXfm);
	CPLAssert(psWarpOptions->pTransformerArg != NULL);

	psWarpOptions->pfnTransformer = GDALGenImgProjTransform;

	// Initialize and execute the warp operation. 

	GDALWarpOperation oOperation;

	oOperation.Initialize(psWarpOptions);
	oOperation.ChunkAndWarpImage(0, 0,
		GDALGetRasterXSize(pdata2),
		GDALGetRasterYSize(pdata2));

	GDALDestroyGenImgProjTransformer(psWarpOptions->pTransformerArg);
	GDALDestroyWarpOptions(psWarpOptions);

	GDALFlushCache(pdata2);
}

int main(int argc, char** argv)
{
	GDALAllRegister();
	string graphfilepath = argv[1];
	string outdirpath = argv[2];


	XfmGraph* graph;

	ifstream ifs(graphfilepath);
	string identifier;
	std::getline(ifs, identifier);
	ifs.seekg(0);
	if (identifier == "PoseGraph")
	{
		auto* g = new PoseGraph();
		ifs >> *g;
		graph = g;
	}
	else if (identifier == "BundleGraph") {
		auto* g = new BundleGraph();
		ifs >> *g;
		graph = g;
	}
	else {
		cerr << "Unknown format: " << identifier << endl;
		return 1;
	}

	fs::path rootdir(graph->basepath);
	fs::path outdir(outdirpath);
	if (!fs::exists(outdir))
		fs::create_directories(outdir);
	size_t num_nodes = graph->nodes.size();

	GDALDriverH driver = GDALGetDriverByName("GTiff");
	char** papszCreateOptions = NULL;
	papszCreateOptions = CSLSetNameValue(papszCreateOptions, "TFW", "YES");
	for (size_t i = 0; i < num_nodes; ++i)
	{
		XfmGraph::Node& n = graph->nodes[i];
		fs::path inputfile = (rootdir / fs::path(n.path));
		fs::path outputfile = (outdir / fs::path(n.path));

		
		auto srcDS = GDALOpen(inputfile.string().c_str(), GA_ReadOnly);
		double sourceXfm[6];
		GDALGetGeoTransform(srcDS, sourceXfm);
		int srcWidth = GDALGetRasterXSize(srcDS);
		int srcHeight = GDALGetRasterYSize(srcDS);
		

		double globalXfm[6];
		XfmComposite(graph->baseGeo, n.poseXfm, globalXfm);
		BBox bbox;
		bbox.add(0, 0);
		bbox.add(srcWidth, srcHeight);
		BBox globbox = bbox.transform(globalXfm);

		double dstXfm[6], invdstXfm[6];
		std::copy_n(sourceXfm, 6, dstXfm);
		dstXfm[0] = globbox.ptMin[0];
		dstXfm[3] = globbox.ptMax[1];
		XfmInv(dstXfm, invdstXfm);

		BBox rectBbox = globbox.transform(invdstXfm);
		int dstWidth = std::ceil(rectBbox.ptMax[0] - rectBbox.ptMin[0]);
		int dstHeight = std::ceil(rectBbox.ptMax[1] - rectBbox.ptMin[1]);

		if (!ensureDir(outputfile.parent_path()))
		{
			cerr << "Cannot create folder for " << outputfile << endl;
		}
		int numbands = GDALGetRasterCount(srcDS);

		GDALDataType dtype = GDALGetRasterDataType(GDALGetRasterBand(srcDS, 1));

		auto dstDS = GDALCreate(driver, outputfile.string().c_str(), dstWidth, dstHeight, numbands, dtype, 0);
		
		GDAL_Affine_TransformImage(srcDS, dstDS, globalXfm, dstXfm, 1);

		GDALClose(srcDS);
		GDALClose(dstDS);
	}

	delete graph;

	return 0;
}