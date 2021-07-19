import os
import sys
from datetime import datetime

import fiona
from affine import Affine
from shapely.geometry import shape, mapping
import numpy as np
from osgeo import osr, gdal, ogr, gdalnumeric
from osgeo.osr import OAMS_TRADITIONAL_GIS_ORDER
import rasterio
import rasterio.features

from geo import reproject_layer, map_raster
from nc_functions import get_transform, read_netcdf
from render import render_field, render_rgb_field
from settings import Settings


def generate_masked_phenomenas_raster(input_raster_path: str, out_file: str, phen_id: str, var_id: str,
                                      pixel_to_regions) -> str:
    sets = Settings()

    src_array = np.flipud(gdalnumeric.LoadFile(input_raster_path))
    nlons: int = src_array.shape[1]
    nlats: int = src_array.shape[0]
    res_array = np.zeros_like(src_array)

    src_image = gdal.Open(input_raster_path)
    proj = src_image.GetProjection()
    # src_band = src_image.GetRasterBand(1)
    geo_trans = src_image.GetGeoTransform()

    # var->reg->crit

    if phen_id in sets.criterias:
        if isinstance(sets.criterias[phen_id][var_id], dict):
            for reg_id in sets.criterias[phen_id][var_id]:
                crit = sets.criterias[phen_id][var_id][reg_id]
                masked = src_array[pixel_to_regions[reg_id]]
                idxs = np.empty_like(masked)
                if sets.mark_phen_operator[var_id] == 'lt':
                    idxs = np.where(masked < crit)
                elif sets.mark_phen_operator[var_id] == 'gt':
                    idxs = np.where(masked > crit)
                # for id in idxs:
                i = pixel_to_regions[reg_id][0][idxs]
                j = pixel_to_regions[reg_id][1][idxs]
                res_array[i, j] = 1
        else:
            idxs = np.empty_like(src_array)
            crit = sets.criterias[phen_id][var_id]
            if sets.mark_phen_operator[var_id] == 'lt':
                idxs = np.where(src_array < crit)
            elif sets.mark_phen_operator[var_id] == 'gt':
                idxs = np.where(src_array > crit)
            res_array[idxs] = 1
    has_ones = np.any(np.where(res_array == 1))
    # print('has_ones={0}'.format(has_ones))
    # phen_code= sets.phenomenas[phen_id]
    # res_file = '{0}_{1}_masked.tiff'.format(input_raster_path, phen_code)
    res_file = out_file
    if has_ones:
        driver = gdal.GetDriverByName('GTiff')
        ds_out = driver.Create(res_file, nlons, nlats, 1, gdal.GDT_Byte)
        ds_out.SetGeoTransform(geo_trans)
        ds_out.SetProjection(proj)
        rot = np.flipud(res_array)
        ds_out.GetRasterBand(1).WriteArray(rot)
        ds_out = None
    else:
        res_file = ''
    return res_file


def generate_masked_phenomenas_raster_rgb(input_raster_path: str, out_file: str, phen_id: str, var_id: str,
                                          pixel_to_regions) -> str:
    sets = Settings()

    src_array = np.flipud(gdalnumeric.LoadFile(input_raster_path))
    nlons: int = src_array.shape[1]
    nlats: int = src_array.shape[0]
    res_array = np.zeros_like(src_array)

    src_image = gdal.Open(input_raster_path)
    proj = src_image.GetProjection()
    # src_band = src_image.GetRasterBand(1)
    geo_trans = src_image.GetGeoTransform()

    # var->reg->crit

    if phen_id in sets.criterias:
        if isinstance(sets.criterias[phen_id][var_id], dict):
            for reg_id in sets.criterias[phen_id][var_id]:
                crit = sets.criterias[phen_id][var_id][reg_id]
                masked = src_array[pixel_to_regions[reg_id]]
                idxs = np.empty_like(masked)
                if sets.mark_phen_operator[var_id] == 'lt':
                    idxs = np.where(masked < crit)
                elif sets.mark_phen_operator[var_id] == 'gt':
                    idxs = np.where(masked > crit)
                # for id in idxs:
                i = pixel_to_regions[reg_id][0][idxs]
                j = pixel_to_regions[reg_id][1][idxs]
                res_array[i, j] = 1
        else:
            idxs = np.empty_like(src_array)
            crit = sets.criterias[phen_id][var_id]
            if sets.mark_phen_operator[var_id] == 'lt':
                idxs = np.where(src_array < crit)
            elif sets.mark_phen_operator[var_id] == 'gt':
                idxs = np.where(src_array > crit)
            res_array[idxs] = 1
    has_ones = np.any(np.where(res_array == 1))
    # print('has_ones={0}'.format(has_ones))
    # phen_code= sets.phenomenas[phen_id]
    # res_file = '{0}_{1}_masked.tiff'.format(input_raster_path, phen_code)
    res_file = out_file
    if has_ones:
        render_rgb_field(np.flipud(res_array), res_file, phen_id, geo_trans, proj)
    else:
        res_file = ''
    return res_file


def generate_phenomenas_contour_wkt(masked_raster: str) -> str:
    if masked_raster == '':
        return 'MULTIPOLYGON EMPTY'
    src_image = gdal.Open(masked_raster)
    proj = src_image.GetProjection()
    # proj.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER)
    # src_band = src_image.GetRasterBand(1)

    with rasterio.open(masked_raster, 'r') as src:
        image = src.read(1)  # read
        gt1 = src.get_transform()
        features = rasterio.features.shapes(image, transform=Affine.from_gdal(*gt1))
        results = list(
            {'properties': {'val': v}, 'geometry': s}
            for i, (s, v)
            in enumerate(
                features))
        res_file: str = '{0}_poly.shp'.format(masked_raster)
        with fiona.open(
                res_file, 'w',
                driver='ESRI Shapefile',
                crs_wkt=proj,
                schema={'properties': [('val', 'int')],
                        'geometry': 'Polygon'}) as dst:
            for res in results:
                geom = shape(res['geometry'])
                if not geom.is_valid:
                    clean = geom.buffer(0.0)
                    assert clean.is_valid
                    assert clean.geom_type == 'Polygon'
                    geom = clean
                res['geometry'] = mapping(geom)
                dst.write(res)
            dst.close()

    # out_spatial_ref = osr.SpatialReference()
    # out_spatial_ref.ImportFromWkt(proj)
    #
    wgs_srs = osr.SpatialReference()
    wgs_srs.ImportFromEPSG(4326)
    wgs_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER)
    #
    # res_file: str = '{0}_poly.shp'.format(masked_raster)
    driver = ogr.GetDriverByName("ESRI Shapefile")
    # ogr_ds = driver.CreateDataSource(res_file)
    # contour_shp = ogr_ds.CreateLayer('contour', out_spatial_ref)
    #
    # field_defn = ogr.FieldDefn("ID", ogr.OFTInteger)
    # contour_shp.CreateField(field_defn)
    # field_defn = ogr.FieldDefn("val", ogr.OFTInteger)
    # contour_shp.CreateField(field_defn)
    #
    # contour_shp.srs = out_spatial_ref
    # # Generate Contourlines
    # # gdal.ContourGenerate(src_band, 0, 0, [1], 0, 0, contour_shp, 0, 1)
    # gdal.Polygonize(src_band, None, contour_shp, 1, [], callback=None)
    # ogr_ds.Destroy()

    f = reproject_layer(res_file, wgs_srs)
    mul_wkt = build_multipolygon_wkt(f, 0, 1, True)
    driver.DeleteDataSource(res_file)
    driver.DeleteDataSource(f)

    return mul_wkt


def build_multipolygon_wkt(in_shp_file: str, fld_num: int, fld_val: float, save_shape: bool) -> str:
    driver = ogr.GetDriverByName('ESRI Shapefile')

    in_data_set = driver.Open(in_shp_file)
    in_layer = in_data_set.GetLayer()

    multipolygon = ogr.Geometry(ogr.wkbMultiPolygon)

    # loop through the input features
    in_feature = in_layer.GetNextFeature()
    while in_feature:
        # get the input geometry
        geom = in_feature.GetGeometryRef()
        val = in_feature.GetField(fld_num)
        if val == fld_val:
            geom.CloseRings()
            if geom.IsValid():
                multipolygon.AddGeometry(geom)
        in_feature = in_layer.GetNextFeature()

    # create the output layer
    if save_shape:
        in_spatial_ref = in_layer.GetSpatialRef()
        output_shapefile = '{}_mul.shp'.format(os.path.splitext(in_shp_file)[0])
        if os.path.exists(output_shapefile):
            driver.DeleteDataSource(output_shapefile)
        out_data_set = driver.CreateDataSource(output_shapefile)
        out_layer = out_data_set.CreateLayer("shapes", in_spatial_ref, geom_type=ogr.wkbMultiPolygon)

        # add fields
        in_layer_defn = in_layer.GetLayerDefn()
        for i in range(0, in_layer_defn.GetFieldCount()):
            field_defn = in_layer_defn.GetFieldDefn(i)
            out_layer.CreateField(field_defn)

        # get the output layer's feature definition
        out_layer_defn = out_layer.GetLayerDefn()
        out_feature = ogr.Feature(out_layer_defn)
        for i in range(0, out_layer_defn.GetFieldCount()):
            out_feature.SetField(out_layer_defn.GetFieldDefn(i).GetNameRef(), 1)
        out_feature.SetGeometry(multipolygon)
        out_layer.CreateFeature(out_feature)
        out_data_set = None
    # Save and close the shapefiles
    in_data_set = None

    return multipolygon.ExportToWkt()


def process(data_nc_file, data_ini_str, out_dir):
    sets = Settings()
    date_ini = datetime.strptime(data_ini_str, '%Y-%m-%d-%H')
    # xtrm_ds = Dataset(xtrm_nc_file)
    nc_data = read_netcdf(data_nc_file)

    # data_ds = Dataset(data_nc_file)
    gt, projWkt = get_transform(data_nc_file)
    first_varid = next(iter(nc_data))
    first_dt = next(iter(nc_data[first_varid]))

    pix_to_reg = map_raster(nc_data[first_varid][first_dt], gt, projWkt, sets.region_shapes)

    files, files_rgb, res_dir = generate_tiffs(date_ini, gt, nc_data, out_dir, projWkt)

    if not os.path.exists('{0}/phens'.format(res_dir)):
        os.makedirs('{0}/phens'.format(res_dir))

    for phen_id in sets.criterias:
        phen_dir = "{0}/phens/{1}".format(res_dir, phen_id)
        if not os.path.exists(phen_dir):
            os.makedirs(phen_dir)
        phen_dir_rgb = "{0}/{1}".format(res_dir, phen_id)
        if not os.path.exists(phen_dir_rgb):
            os.makedirs(phen_dir_rgb)
        dt_dir_rgb = "{0}/{1}".format(phen_dir_rgb, data_ini_str)
        if not os.path.exists(dt_dir_rgb):
            os.makedirs(dt_dir_rgb)

        for var_id in sets.criterias[phen_id]:
            if var_id in files:
                for dt in files[var_id]:
                    #dt_str = datetime.fromtimestamp(dt.item() / 10**9).strftime('%Y-%m-%d-%H')
                    dt_str = np.datetime_as_string(dt).replace("T", "-")
                    file = files[var_id][dt]
                    if file != '':
                        out_file = '{0}/{1}.tiff'.format(phen_dir, dt_str)
                        marked_raster_file = generate_masked_phenomenas_raster(file, out_file, phen_id, var_id, pix_to_reg)
                        if marked_raster_file != '':
                            out_file_rgb = '{0}/{1}.tiff'.format(dt_dir_rgb, dt_str)
                            generate_masked_phenomenas_raster_rgb(file, out_file_rgb, phen_id, var_id, pix_to_reg)
                        poly_wkt = generate_phenomenas_contour_wkt(marked_raster_file)
                        print('{0};{1};{2}'.format(dt, phen_id, poly_wkt))


def generate_tiffs(date_ini, gt, nc_data, out_dir, projWkt):
    sets = Settings()
    files = {}
    files_rgb = {}
    res_dir = "{0}/{1}".format(sets.tiff_root, out_dir)
    res_dir_1b = "{0}/{1}".format(sets.tiff_root_1b, out_dir)
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)
        # os.makedirs('{0}/vars'.format(res_dir))
    if not os.path.exists(res_dir_1b):
        os.makedirs(res_dir_1b)

    for var_id in nc_data:
        files[var_id] = {}
        files_rgb[var_id] = {}
        var_dir = "{0}/{1}".format(res_dir_1b, var_id)
        var_rgb_dir = "{0}/{1}".format(res_dir, var_id)
        if not os.path.exists(var_dir):
            os.makedirs(var_dir)
        if not os.path.exists(var_rgb_dir):
            os.makedirs(var_rgb_dir)
        dts_dir = "{0}/{1}".format(var_rgb_dir, datetime.strftime(date_ini, '%Y-%m-%d-%H'))
        dts_dir_1b = "{0}/{1}".format(var_dir, datetime.strftime(date_ini, '%Y-%m-%d-%H'))
        if not os.path.exists(dts_dir):
            os.makedirs(dts_dir)
        if not os.path.exists(dts_dir_1b):
            os.makedirs(dts_dir_1b)

        for date, arr in nc_data[var_id].items():
            # start_dir = "{0}/{1}".format(var_dir, data_ini.strftime('%Y-%m-%d-%H00'))
            # if not os.path.exists(start_dir):
            #     os.makedirs(start_dir)
            dt = date.astype('M8[ms]').astype('O')
            dt_str = datetime.strftime(dt, '%Y-%m-%d-%H')
            filename = '{0}/{1}.tiff'.format(dts_dir_1b, dt_str)
            filename_rgb = '{0}/{1}.tiff'.format(dts_dir, dt_str)

            has_any = np.any(np.where(arr != 0))
            if has_any:
                render_field(arr, filename, gt, projWkt)
                render_rgb_field(arr, filename_rgb, var_id, gt, projWkt)
            else:
                filename = ''
                filename_rgb = ''

            files[var_id][date] = filename
            files_rgb[var_id][date] = filename_rgb

        # cut_rasters(files[var_id],  projWkt)
    return files, files_rgb, res_dir


if __name__ == '__main__':

    dom_id = sys.argv[1]
    #dom_id = '2'
    dt_ini_str = sys.argv[2]
    #dt_ini_str = '2021-05-10-12'
    # root_dir = '\\\\10.11.25.55\\NWP_RawData\\WRF_RUMS\\9-3-1'
    #root_dir = 'f:\\wrf_Res'
    root_dir = '/mnt/home/export'
    # in_data_nc_file = 'f:\\wrfout_d03_cut.nc'
    in_data_nc_file = '{0}/{1}/wrfout_d0{2}_cut.nc'.format(root_dir, dt_ini_str, dom_id)

    # xtrm_nc_file = sys.argv[1]
    # data_nc_file = sys.argv[2]
    # dt_ini_str = sys.argv[3]
    out_dir = 'wrf'
    if dom_id == '3':
        out_dir = 'wrf1'
    elif dom_id == '2':
        out_dir = 'wrf3'
    else:
        exit(0)
    process(in_data_nc_file, dt_ini_str, out_dir)
