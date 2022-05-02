#pragma once

#include <gdal.h>
#include <gdal_priv.h>
#include <gdalwarper.h>

#include <iostream>
#include <string>

/*Note: the following is image level affine transformation
for the transformed tfw information, the upper left corner be computed by Affine_GEO, the cell size info follows the reference image
*/

void Affine_TransformImage(GDALDataset* pdata1, GDALDataset* pdata2, double* affpara, bool isbilinear_or_nearest = 1);

void get_max_min(double* X, int numpt, double& maxv, double& minv) {
  if (numpt <= 0) return;
  maxv = minv = X[0];

  for (int p = 1; p < numpt; p++) {
    if (X[p] > maxv)
      maxv = X[p];
    else if (X[p] < minv)
      minv = X[p];
  }

  // end
}

void PerformAffineG(double* affpara, double x, double y, double& x1, double& y1) {
  x1 = affpara[0] + x * affpara[1] + y * affpara[2];
  y1 = affpara[3] + x * affpara[4] + y * affpara[5];
}

/*****************************************************************************************
this prepare the affine transformation parameters, this will not update the affine model
******************************************************************************************/
bool Prepare_Affine_Parameters(double* AffineModelin, int imgwidth, int imgheight, int& imgwidthnew, int& imgheightnew, double& Ul_cornerX, double& Ul_cornery) {
  double corner2[8] = {0};
  corner2[0] = 0;
  corner2[4] = 0;
  corner2[1] = imgwidth - 1;
  corner2[5] = 0;
  corner2[2] = imgwidth - 1;
  corner2[6] = imgheight - 1;
  corner2[3] = 0;
  corner2[7] = imgheight - 1;

  double res_corner2[8] = {0};
  // now get the transformed image
  for (int p = 0; p < 4; p++) {
    PerformAffineG(AffineModelin, corner2[p], corner2[p + 4], res_corner2[p], res_corner2[p + 4]);
  }  // this gives the position after the affine
  // revised codes from the epi estimation
  double maxv, minv;
  double ulcorner_X_2, ulcorner_Y_2;
  get_max_min(res_corner2, 4, maxv, minv);
  Ul_cornerX = minv;
  imgwidthnew = ceil(maxv - minv);
  get_max_min(res_corner2 + 4, 4, maxv, minv);
  Ul_cornery = minv;
  imgheightnew = ceil(maxv - minv);
  // now getting the parameters for transformation
  return true;
}

//////////////////////////////////////////////////////////////////////////
// enforce nearest interpolation.
bool PerformAffineImage_2(std::string imgnamein, std::string imgnameout, int imgwidthnew, int imgheightnew, double* AffineModel) {
  GDALAllRegister();
  GDALDataset* datasetin = (GDALDataset*)GDALOpen(imgnamein.c_str(), GA_Update);
  if (!datasetin) {
    std::cout << "error: image cannot be open: " << imgnamein << std::endl;
    return false;
  }

  int imgwidth = datasetin->GetRasterXSize();
  int imgheight = datasetin->GetRasterYSize();
  int numbands = datasetin->GetRasterCount();
  GDALDataType dtype = datasetin->GetRasterBand(1)->GetRasterDataType();

  GDALDriver* pdrive = (GDALDriver*)GDALGetDriverByName("GTIFF");
  GDALDataset* datasetout = pdrive->Create(imgnameout.c_str(), imgwidthnew, imgheightnew, numbands, dtype, 0);

  // double aff[6];
  // GDALGetGeoTransform(datasetin, aff);

  Affine_TransformImage(datasetin, datasetout, AffineModel, 1);

  // GDALSetGeoTransform(datasetin, aff);

  // done

  GDALClose(datasetin);
  GDALClose(datasetout);

  return true;
}

void Affine_TransformImage(GDALDataset* pdata1, GDALDataset* pdata2, double* affpara, bool isbilinear_or_nearest /*=1*/) {  // affpara has 6 parameters.
  // second dataset containts the necessary information.
  //	Xp = padfTransform[0] + P*padfTransform[1] + L*padfTransform[2];            // new coordinates.
  //	Yp = padfTransform[3] + P*padfTransform[4] + L*padfTransform[5];
  pdata1->Dereference();
  pdata2->Dereference();

  pdata1->SetGeoTransform(affpara);
  GDALWarpOptions* psWarpOptions = GDALCreateWarpOptions();

  psWarpOptions->hSrcDS = pdata1;
  psWarpOptions->hDstDS = pdata2;

  psWarpOptions->nBandCount = pdata1->GetRasterCount();
  psWarpOptions->panSrcBands = (int*)CPLMalloc(sizeof(int) * psWarpOptions->nBandCount);
  for (int p = 0; p < pdata1->GetRasterCount(); p++) psWarpOptions->panSrcBands[p] = p + 1;

  psWarpOptions->panDstBands = (int*)CPLMalloc(sizeof(int) * psWarpOptions->nBandCount);
  for (int p = 0; p < pdata1->GetRasterCount(); p++) psWarpOptions->panDstBands[p] = p + 1;

  psWarpOptions->pfnProgress = GDALTermProgress;
  if (isbilinear_or_nearest)
    psWarpOptions->eResampleAlg = GRA_Bilinear;
  else
    psWarpOptions->eResampleAlg = GRA_NearestNeighbour;

  // psWarpOptions->eResampleAlg = GRA_NearestNeighbour;
  // Establish reprojection transformer.

  auto wtk1 = GDALGetProjectionRef(pdata1);
  auto wtk2 = GDALGetProjectionRef(pdata2);
  CPLAssert(wtk1 != NULL && strlen(wtk1) > 0);
  CPLAssert(wtk2 != NULL && strlen(wtk2) > 0);
  // psWarpOptions->pTransformerArg =
  //	GDALCreateGenImgProjTransformer(pdata1,
  //		NULL,
  //		pdata2,
  //		NULL,
  //		FALSE, 0.0, 1);
  double identaff[6] = {0., 1., 0., 0., 0., 1.};
  GDALSetGeoTransform(pdata2, identaff);
  psWarpOptions->pTransformerArg = GDALCreateGenImgProjTransformer3("", affpara, "", identaff);
  CPLAssert(psWarpOptions->pTransformerArg != NULL);

  psWarpOptions->pfnTransformer = GDALGenImgProjTransform;
  // GDALSetGenImgProjTransformerDstGeoTransform(psWarpOptions->pTransformerArg,affpara);

  // Initialize and execute the warp operation.

  GDALWarpOperation oOperation;

  oOperation.Initialize(psWarpOptions);
  oOperation.ChunkAndWarpImage(0, 0, GDALGetRasterXSize(pdata2), GDALGetRasterYSize(pdata2));

  GDALDestroyGenImgProjTransformer(psWarpOptions->pTransformerArg);
  GDALDestroyWarpOptions(psWarpOptions);

  pdata2->FlushCache();
}
