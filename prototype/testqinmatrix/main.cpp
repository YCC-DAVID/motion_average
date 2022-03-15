#include "../motionaverage/tem_codes.h"
#include "../motionaverage/myutils.h"
#include <fstream>
#include <iomanip>

using namespace std;

int main(int argc, char** argv)
{
	GDALAllRegister();
	fs::path inputpath = fs::path(argv[1]);
	fs::path inputparampath = fs::path(argv[2]);
	
	fs::path outdir = inputparampath.parent_path();
	fs::path barename = inputparampath.filename().stem();
	fs::path outputpath = outdir / (barename.string() +".tif");
	fs::path outputeqpath = outdir / (barename.string() + "_eq.tif");

	rspopt::PairwiseData pairdata;
	pairdata.load(inputparampath.string());
	rspopt::AffineModel& affinemodel = pairdata.affine;
	
	cout << affinemodel << endl;

	double originGeo[6];
	int originWidth, originHeight;
	GDALDatasetH inds = GDALOpen(inputpath.string().c_str(), GA_ReadOnly);
	GDALGetGeoTransform(inds, originGeo);
	originWidth = GDALGetRasterXSize(inds);
	originHeight = GDALGetRasterYSize(inds);
	GDALClose(inds);

	//////////// Compute New TFW
	double updatedGeo[6];//, updatedTfw[6];
	std::copy_n(originGeo, 6, updatedGeo);

	rspopt::XfmApply(affinemodel.geoParam, originGeo[0] - affinemodel.pixelParam[0], originGeo[3] + affinemodel.pixelParam[3], updatedGeo[0], updatedGeo[3]);

	cout << "------------ eq \n";
	double composeGeo[6];// , composeTfw[6];
	//double invPixel[6];
	rspopt::XfmComposite(updatedGeo, affinemodel.pixelParam, composeGeo);

	///////// Create new image
	int updatedWidth, updatedHeight;
	double Ul_cornerx, Ul_cornery;
	Prepare_Affine_Parameters(affinemodel.pixelParam, originWidth, originHeight, updatedWidth, updatedHeight, Ul_cornerx, Ul_cornery);
	PerformAffineImage_2(inputpath.string(), outputpath.string(), updatedWidth, updatedHeight, affinemodel.pixelParam);

	/// Output TFW

	inds = GDALOpen(inputpath.string().c_str(), GA_ReadOnly);

	GDALDriverH driver = GDALGetDriverByName("GTiff");
	GDALDatasetH outeq = GDALCreateCopy(driver, outputeqpath.string().c_str(), inds, FALSE, NULL, NULL, NULL);
	GDALSetGeoTransform(outeq, composeGeo);
	GDALClose(outeq);

	GDALDatasetH out = GDALOpen(outputpath.string().c_str(), GA_Update);
	GDALSetGeoTransform(out, updatedGeo);
	GDALClose(out);

	GDALClose(inds);
	return 0;
}