import os

import matplotlib
import matplotlib.cm
import matplotlib.colors
import matplotlib.pyplot as plt

import numpy as np
from osgeo import gdal

from settings import Settings


def scale(A, min_val, max_val):
    # return (A-np.min(A))/(20 - np.min(A))
    # return (A-np.min(A))/(np.max(A) - np.min(A))
    result = np.zeros_like(A)
    if max_val != min_val:
        result = (A - min_val) / (max_val - min_val)
    return result


def render_field(data: np.ndarray, filename, gt, out_proj_wkt):
    # print("min={0},max={1}".format(min_val, max_val))
    nlons = data.shape[1]
    nlats = data.shape[0]

    driver = gdal.GetDriverByName('GTiff')
    ds_out = driver.Create(filename, nlons, nlats, 1, gdal.GDT_Float32)
    ds_out.SetGeoTransform(gt)
    ds_out.SetProjection(out_proj_wkt)
    ds_out.GetRasterBand(1).WriteArray(data[:])
    ds_out = None


def render_rgb_field(data, filename, varname, gt, out_proj_wkt):
    sets = Settings()
    min_val = np.min(data)
    max_val = np.max(data)
    # print("min={0},max={1}".format(min_val, max_val))
    # clmap = matplotlib.cm.get_cmap('coolwarm',12)

    cols = ((0.0, np.asarray([256, 0, 0, 0]) / 256),
            (0.5, np.asarray([0, 256, 0, 0.87 * 256]) / 256),
            (1.0, np.asarray([0, 0, 256, 1 * 256]) / 256)

            )

    cmp = matplotlib.colors.LinearSegmentedColormap.from_list("test", cols)
    cb_vals = ()
    if varname == "wspdmax":
        min_val = 0
        max_val = 42
        max_val = 42
        cb_vals = (0,2,4,6,8,10,12,14,18,22,26,30,34,38,42)
        wind_cols = (
            (42 / max_val, np.asarray([43, 0, 1, 1 * 256]) / 256),
            (38 / max_val, np.asarray([68, 10, 25, 1 * 256]) / 256),
            (34 / max_val, np.asarray([131, 23, 67, 0.93 * 256]) / 256),
            (30 / max_val, np.asarray([181, 47, 95, 0.85 * 256]) / 256),
            (26 / max_val, np.asarray([227, 103, 87, 0.85 * 256]) / 256),
            (22 / max_val, np.asarray([234, 168, 61, 0.85 * 256]) / 256),
            (18 / max_val, np.asarray([235, 224, 53, 0.85 * 256]) / 256),
            (14 / max_val, np.asarray([153, 220, 69, 0.85 * 256]) / 256),
            (12 / max_val, np.asarray([92, 208, 74, 0.85 * 256]) / 256),
            (10 / max_val, np.asarray([76, 196, 92, 0.85 * 256]) / 256),
            (8 / max_val, np.asarray([73, 178, 158, 0.85 * 256]) / 256),
            (6 / max_val, np.asarray([59, 143, 193, 0.85 * 256]) / 256),
            (4 / max_val, np.asarray([68, 103, 195, 0.85 * 256]) / 256),
            (2 / max_val, np.asarray([83, 82, 180, 0.68 * 256]) / 256),
            (0 / max_val, np.asarray([82, 71, 141, 0.42 * 256]) / 256)
        )
        cmp = matplotlib.colors.LinearSegmentedColormap.from_list("test", wind_cols[::-1])

    if varname == "tmin" or varname == "tmax":
        min_val = -40
        max_val = 50
        tmp_cols = (
            ((50 - min_val) / (max_val - min_val), np.asarray([43, 0, 1, 1 * 256]) / 256),
            ((40 - min_val) / (max_val - min_val), np.asarray([107, 21, 39, 1 * 256]) / 256),
            ((30 - min_val) / (max_val - min_val), np.asarray([190, 48, 102, 0.92 * 256]) / 256),
            ((25 - min_val) / (max_val - min_val), np.asarray([229, 109, 83, 0.92 * 256]) / 256),
            ((20 - min_val) / (max_val - min_val), np.asarray([234, 164, 62, 0.92 * 256]) / 256),
            ((15 - min_val) / (max_val - min_val), np.asarray([235, 215, 53, 0.92 * 256]) / 256),
            ((10 - min_val) / (max_val - min_val), np.asarray([190, 228, 61, 0.92 * 256]) / 256),
            ((5 - min_val) / (max_val - min_val), np.asarray([89, 208, 73, 0.92 * 256]) / 256),
            ((0 - min_val) / (max_val - min_val), np.asarray([75, 182, 152, 0.92 * 256]) / 256),
            ((-5 - min_val) / (max_val - min_val), np.asarray([62, 121, 198, 0.92 * 256]) / 256),
            ((-10 - min_val) / (max_val - min_val), np.asarray([85, 78, 177, 0.92 * 256]) / 256),
            ((-15 - min_val) / (max_val - min_val), np.asarray([36, 24, 106, 0.92 * 256]) / 256),
            ((-20 - min_val) / (max_val - min_val), np.asarray([145, 9, 145, 0.92 * 256]) / 256),
            ((-30 - min_val) / (max_val - min_val), np.asarray([255, 170, 255, 0.92 * 256]) / 256),
            ((-40 - min_val) / (max_val - min_val), np.asarray([238, 238, 238, 0.92 * 256]) / 256)
        )
        #cb_vals = (-40,-30,-20,-15,-10,-5,0,5,10,15,20,25,30,40,50)
        cmp = matplotlib.colors.LinearSegmentedColormap.from_list("test", tmp_cols[::-1])


    if varname == "rain1h" or varname == "rain12h" or varname == "snow1h" or varname == "snow12h":
        min_val = 0
        max_val = 50
        precip_cols = (
            (50 / max_val, np.asarray([84, 16, 41, 1 * 256]) / 256),
            (40 / max_val, np.asarray([147, 23, 78, 1 * 256]) / 256),
            (30 / max_val, np.asarray([190, 48, 102, 0.87 * 256]) / 256),
            (20 / max_val, np.asarray([225, 94, 93, 0.87 * 256]) / 256),
            (15 / max_val, np.asarray([233, 123, 72, 0.87 * 256]) / 256),
            (10 / max_val, np.asarray([234, 164, 62, 0.87 * 256]) / 256),
            (8 / max_val, np.asarray([235, 192, 56, 0.87 * 256]) / 256),
            (6 / max_val, np.asarray([220, 234, 55, 0.87 * 256]) / 256),
            (4 / max_val, np.asarray([149, 219, 70, 0.87 * 256]) / 256),
            (2 / max_val, np.asarray([78, 194, 98, 0.87 * 256]) / 256),
            (1 / max_val, np.asarray([64, 160, 180, 0.87 * 256]) / 256),
            (0.5 / max_val, np.asarray([67, 105, 196, 0.75 * 256]) / 256),
            (0.2 / max_val, np.asarray([85, 78, 177, 0.5 * 256]) / 256),
            (0.1 / max_val, np.asarray([82, 71, 141, 0.38 * 256]) / 256),
            (0, np.asarray([82, 71, 141, 0 * 256]) / 256)
        )
        #cb_vals = (0,0.1,0.2,0.5,1,2,4,6,8,10,15,20,30,40,50)
        cmp = matplotlib.colors.LinearSegmentedColormap.from_list("test", precip_cols[::-1])

    if varname in sets.phenomena_colors:
        min_val = 0
        max_val = 1
        cols = ((0.0, np.asarray([0, 0, 0, 0]) / 256),
                (1.0, np.asarray(sets.phenomena_colors[varname]) / 256)
                )
        # cb_vals = (-40,-30,-20,-15,-10,-5,0,5,10,15,20,25,30,40,50)
        cmp = matplotlib.colors.LinearSegmentedColormap.from_list("test", cols)


    #
    # if varname == "rain":
    #     min_val = 0.1
    #     max_val = 100
    #     clmap = matplotlib.cm.get_cmap('brg',12)
    #
    # if varname == "tmin" or varname == "tmax" or varname == "t2":
    #     min_val = -40
    #     max_val = 50

    # f = open('{0}/minmax.txt'.format(start_dir), "w")
    # f.write('{0};{1}'.format(min, max))
    # f.close()

#    if not os.path.exists('{0}_cbar.png'.format(varname)):
        #plt.figure()
        #mpb = plt.pcolormesh(data[:], cmap=cmp, vmin=min_val, vmax=max_val)
        #cb = plt.colorbar()



        # fig, ax = plt.subplots()
        # cb = plt.colorbar(mpb, ax=ax)
        # #cb.set_ticks(cb_vals)
        # #cb.set_ticklabels(cb_vals)
        # ax.remove()
        #plt.savefig('{0}_cbar.png'.format(varname), bbox_inches='tight')

    # print("min={0},max={1}".format(min_val, max_val))
    nlons = data.shape[1]
    nlats = data.shape[0]

    arr = scale(data[:], min_val, max_val)

    rgba = cmp(arr)
    #html_leg = cmp._repr_html_()

    rgba[data[:] == 0] = 0
    driver = gdal.GetDriverByName('GTiff')
    ds_out = driver.Create(filename, nlons, nlats, 4, gdal.GDT_Byte)
    ds_out.SetGeoTransform(gt)
    ds_out.SetProjection(out_proj_wkt)
    ds_out.GetRasterBand(1).WriteArray(rgba[:, :, 0] * 255)
    ds_out.GetRasterBand(2).WriteArray(rgba[:, :, 1] * 255)
    ds_out.GetRasterBand(3).WriteArray(rgba[:, :, 2] * 255)
    ds_out.GetRasterBand(4).WriteArray(rgba[:, :, 3] * 255)
    ds_out = None
