#include <gdal.h>
#include <gdal_priv.h>
#include <string>
#include "../motionaverage//myutils.h"
#include "../motionaverage/tem_codes.h"

using namespace std;
int main(int argc, char** argv)
{
	GDALAllRegister();
	string inputpath = "J:\\Temp_for_SS\\skysat\\processed\\greedy_single\\20211114_013021_ssc18d1_0010\\basic_analytic_udm2\\111\\raw.tif";
	string outputpath = "J:\\Temp_for_SS\\skysat\\processed\\greedy_single\\20211114_013021_ssc18d1_0010\\basic_analytic_udm2\\111\\warp.tif";
	double geo2tfw[6] = { 0.5, 1., 0., 0.5, 0., 1. }, tfw2geo[6] = { -0.5, 1., 0., -0.5, 0., 1. };
	double originTfw[6] = { 739236.042309402490000, 0.781135141849518,0, 2555789.764514132400000,0,-0.781154811382294 };
	double originGeo[6], updatedGeo[6];
	rspopt::XfmComposite(originTfw, tfw2geo, originGeo);
	std::copy_n(originGeo, 6, updatedGeo);
	double geoParam[6] = { -3184.326911186100900, 0.998476120398305, 0.001686528212597, 9000.293614128604500, 0.000670916989548, 0.996284464652968 };
	double pixelParam[6] = { 4.385083769467428, 0.998476120398306, - 0.001686570680564, 0.742015505967174, - 0.000670900095811, 0.996284464652968 };
	rspopt::XfmApply(geoParam, originGeo[0] - pixelParam[0], originGeo[3] + pixelParam[3], updatedGeo[0], updatedGeo[3]);

	double idenGeo[6] = { 1e-7,1,0,1e-7,0,1 };
	int originWidth, originHeight;
	GDALDatasetH inds = GDALOpen(inputpath.c_str(), GA_Update);
	GDALGetGeoTransform(inds, originGeo);
	GDALSetGeoTransform(inds, geoParam);
	originWidth = GDALGetRasterXSize(inds);
	originHeight = GDALGetRasterYSize(inds);
	GDALClose(inds);

	int updatedWidth, updatedHeight;
	double Ul_cornerx, Ul_cornery;
	Prepare_Affine_Parameters(geoParam, originWidth, originHeight, updatedWidth, updatedHeight, Ul_cornerx, Ul_cornery);
	PerformAffineImage_2(inputpath, outputpath, updatedWidth, updatedHeight, pixelParam);


	inds = GDALOpen(inputpath.c_str(), GA_Update);
	double affineUTM[6];
	rspopt::XfmComposite(updatedGeo, pixelParam, affineUTM);
	GDALSetGeoTransform(inds, affineUTM);
	GDALClose(inds);


	GDALDatasetH ouds = GDALOpen(outputpath.c_str(), GA_Update);
	double rectUTM[6];
	rspopt::XfmComposite(updatedGeo, idenGeo, rectUTM);
	GDALSetGeoTransform(ouds, rectUTM);
	GDALClose(ouds);

	return 0;
}