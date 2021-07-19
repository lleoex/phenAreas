from typing import List, Dict, Any

import numpy as np
from osgeo import gdal, gdalnumeric, ogr, osr
import os

from settings import Settings


def reproject_layer(in_shp_file, target_srs):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    in_data_set = driver.Open(in_shp_file)
    in_layer = in_data_set.GetLayer()
    in_spatial_ref = in_layer.GetSpatialRef()
    # output SpatialReference
    outSpatialRef = target_srs
    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(in_spatial_ref, outSpatialRef)
    # create the output layer
    output_shapefile = '{}_rep.shp'.format(os.path.splitext(in_shp_file)[0])

    if os.path.exists(output_shapefile):
        driver.DeleteDataSource(output_shapefile)
    out_data_set = driver.CreateDataSource(output_shapefile)
    out_layer = out_data_set.CreateLayer("shapes", outSpatialRef, geom_type=ogr.wkbMultiPolygon)
    out_layer.srs = outSpatialRef
    # add fields
    in_layer_defn = in_layer.GetLayerDefn()
    for i in range(0, in_layer_defn.GetFieldCount()):
        field_defn = in_layer_defn.GetFieldDefn(i)
        out_layer.CreateField(field_defn)
    # get the output layer's feature definition
    out_layer_defn = out_layer.GetLayerDefn()
    # loop through the input features
    in_feature = in_layer.GetNextFeature()
    while in_feature:
        # get the input geometry
        geom = in_feature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        out_feature = ogr.Feature(out_layer_defn)
        # set the geometry and attribute
        out_feature.SetGeometry(geom)
        for i in range(0, out_layer_defn.GetFieldCount()):
            out_feature.SetField(out_layer_defn.GetFieldDefn(i).GetNameRef(), in_feature.GetField(i))
        # add the feature to the shapefile
        out_layer.CreateFeature(out_feature)
        # dereference the features and get the next input feature
        out_feature = None
        in_feature = in_layer.GetNextFeature()

    # Save and close the shapefiles
    in_data_set = None
    out_data_set = None
    return output_shapefile


def map_raster(in_arr: np.ndarray, in_gt: List[int], in_srs_wkt: str, region_shapes: Dict[int,str]) -> Dict[int, Any]:
    out_spatial_ref = osr.SpatialReference()
    out_spatial_ref.ImportFromWkt(in_srs_wkt)
    nlons = in_arr.shape[1]
    nlats = in_arr.shape[0]
    minx = in_gt[0]
    maxx = minx + nlons * in_gt[1] + nlats * in_gt[2]
    miny = in_gt[3]
    maxy = miny + nlats * in_gt[5] + nlons * in_gt[4]
    result: Dict[int, Any] = {}

    driver = gdal.GetDriverByName('GTiff')
    ds_out = driver.Create('tmp.tiff', nlons, nlats, 1, gdal.GDT_Float32)
    ds_out.SetGeoTransform(in_gt)
    ds_out.SetProjection(in_srs_wkt)
    ds_out.GetRasterBand(1).WriteArray(in_arr)
    ds_out = None

    for reg_id in region_shapes:
        projected_shape = reproject_layer(region_shapes[reg_id], out_spatial_ref)
        res_filename = 'tmp_map_{0}.tiff'.format(reg_id)
        options = gdal.WarpOptions(dstAlpha=False, multithread=True, cutlineDSName=projected_shape,
                                   xRes=in_gt[1], yRes=in_gt[5],
                                   outputBounds=(minx, miny, maxx, maxy),
                                   # width=src_ds.RasterXSize, height=src_ds.RasterYSize,
                                   # cropToCutline=True,
                                   dstNodata=-9999)  # , dstSRS=outSrsWkt, srcSRS=outSrsWkt)

        gdal.Warp(res_filename, 'tmp.tiff', options=options)
        reg_ds = gdal.Open(res_filename)
        res_arr = gdalnumeric.LoadFile(res_filename)
        src_band = reg_ds.GetRasterBand(1)
        filled_idx = np.where(res_arr[:] != src_band.GetNoDataValue())
        result[reg_id] = filled_idx
    return result
