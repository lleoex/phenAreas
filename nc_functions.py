from datetime import datetime

from netCDF4 import Dataset
import numpy as np
from osgeo import osr

from settings import Settings


def get_time_array(ds: Dataset) -> []:
    times_char = ds.variables['Times']
    times = []
    for t in times_char:
        times.append(get_time(t))
    time_arr = np.array(times)
    return time_arr


def get_time(arr_char) -> datetime:
    time = []
    for c in arr_char:
        time.append(c.decode('UTF-8'))
    # yyyy-MM-dd_HH:mm:ss -> yyyy-MM-dd-HH
    time_str = "".join(time).replace("_", " ")
    t = datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')
    return t


def get_transform(ds_file_path: str):
    nc_ds = Dataset(ds_file_path)
    tl1 = float(nc_ds.getncattr('TRUELAT1'))
    tl2 = float(nc_ds.getncattr('TRUELAT2'))
    clat = float(nc_ds.getncattr('MOAD_CEN_LAT'))
    clon = float(nc_ds.getncattr('STAND_LON'))

    wgsProjWKT = "GEOGCS[\"WGS84 datum, Latitude-Longitude; Degrees\", DATUM[\"WGS_1984\",  \
        SPHEROID[\"World Geodetic System of 1984, GEM 10C\",6378137.0,298.257223563, AUTHORITY[\"EPSG\",\"7030\"]], \
        AUTHORITY[\"EPSG\",\"6326\"]], PRIMEM[\"Greenwich\",0], UNIT[\"degree\",0.0174532925199433], \
        AUTHORITY[\"EPSG\",\"4326\"]]"

    llcProjWKT = 'PROJCS["llc_132_44_wrf",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",' \
                 'SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0]' \
                 ',UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic"],' \
                 'PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",{2}],' \
                 'PARAMETER["Standard_Parallel_1",{0}],PARAMETER["Standard_Parallel_2",{1}],' \
                 'PARAMETER["Scale_Factor",1.0],PARAMETER["Latitude_Of_Origin",{3}],UNIT["Meter",1.0]]'.format(tl1, tl2,
                                                                                                               clon,
                                                                                                               tl1)

    wgs_srs = osr.SpatialReference()
    wgs_srs.ImportFromWkt(wgsProjWKT)

    llc_srs = osr.SpatialReference()
    llc_srs.ImportFromWkt(llcProjWKT)

    wgs_to_llc_transf = osr.CoordinateTransformation(wgs_srs, llc_srs)

    longs = nc_ds.variables['XLONG'][0][:]
    lats = nc_ds.variables['XLAT'][0][:]

    nlats = longs.shape[0]
    nlons = longs.shape[1]

    ul = [float(lats[0, 0]), float(longs[0, 0])]
    ur = [float(lats[0, nlons - 1]), float(longs[0, nlons - 1])]
    lr = [float(lats[nlats - 1, nlons - 1]), float(longs[nlats - 1, nlons - 1])]
    ll = [float(lats[nlats - 1, 0]), float(longs[nlats - 1, 0])]

    ul_llc = wgs_to_llc_transf.TransformPoint(ul[0], ul[1])
    ur_llc = wgs_to_llc_transf.TransformPoint(ur[0], ur[1])
    lr_llc = wgs_to_llc_transf.TransformPoint(lr[0], lr[1])
    ll_llc = wgs_to_llc_transf.TransformPoint(ll[0], ll[1])

    gt1 = (ur_llc[0] - ul_llc[0]) / (nlons - 1)
    gt0 = ul_llc[0] - gt1 / 2.0
    gt2 = (lr_llc[0] - ul_llc[0] + gt1 - gt1 * nlons) / (nlats - 1)
    gt5 = (ll_llc[1] - ul_llc[1]) / (nlats - 1)
    gt3 = ul_llc[1] - gt5 / 2.0
    gt4 = (lr_llc[1] - ul_llc[1] + gt5 - gt5 * nlats) / (nlons - 1)
    gt = (gt0, gt1, gt2, gt3, gt4, gt5)
    return gt, llcProjWKT


def read_netcdf(ds_file_path: str):
    sets = Settings()
    result = {}

    ds = Dataset(ds_file_path)
    time_arr = get_time_array(ds)[1:]
    time_min = np.datetime64(time_arr.min(initial=None))
    time_max = np.datetime64(time_arr.max(initial=None))

    h = time_min.astype(object).hour

    # adjust to day/night begining. but why?
    if 9 <= h < 21:
        dts = np.datetime64(time_min, 'D') + np.timedelta64(9, 'h')
    elif h >= 21:
        dts = np.datetime64(time_min, 'D') + np.timedelta64(21, 'h')
    else:
        dts = np.datetime64(time_min, 'D') - np.timedelta64(3, 'h')

    t = dts  # + np.timedelta64(8,'h')

    t_min_ds = ds.variables['T2MIN'][1:, :, :]
    t_max_ds = ds.variables['T2MAX'][1:, :, :]
    t_ds = ds.variables['T2'][1:, :, :]
    wspd_max_ds = ds.variables['SPDUV10MAX'][1:, :, :]
    precip_ds = ds.variables['RAINC'][:] + ds.variables['RAINNC'][:]
    prate_ds = precip_ds[1:, :, :] - precip_ds[:-1, :, :]
    snow_ds = ds.variables['SNOWNC'][:]
    snowrate_ds = snow_ds[1:, :, :] - snow_ds[:-1, :, :]

    result['tmax'] = {}
    result['tmin'] = {}

    result['rain1h'] = {}
    result['snow1h'] = {}
    result['sr'] = {}
    result['wspdmax'] = {}

    result['rain12h'] = {}
    result['snow12h'] = {}

    # te = t + np.timedelta64(1,'h')
    # for var_id in sets.variables:
    #    result[var_id] = {}
    td_1h = np.timedelta64(1, 'h')
    td_12h = np.timedelta64(12, 'h')
    while t < time_max:
        h_offset_td = np.timedelta64(t - time_min + td_1h, 'h')
        h_offset = int(h_offset_td.astype(object).seconds / 3600)
        # tt = np.logical_and(time_arr > t, time_arr <= (t + td_1h))
        tt = np.logical_and(time_arr > t, time_arr <= (t + td_1h))
        t = t + td_1h
        if tt.max():
            t_min = t_min_ds[tt].min(0) - 273.16
            t_max = t_max_ds[tt].max(0) - 273.16
            precip1h = prate_ds[tt].sum(0)
            snow1h = snowrate_ds[tt].sum(0)
            sr = snow1h / (precip1h + 0.0001)
            solid_idx = np.where(sr > 0.8)
            liq_idx = np.where(sr <= 0.8)
            wspd_max = wspd_max_ds[tt].max(0)

            result['tmin'][t] = t_min
            result['tmax'][t] = t_max
            result['rain1h'][t] = precip1h-snow1h
            result['rain1h'][t][liq_idx] = precip1h[liq_idx]
            result['snow1h'][t] = np.zeros_like(snow1h)
            result['snow1h'][t][solid_idx] = snow1h[solid_idx]
            result['rain1h'][t][np.where(result['rain1h'][t] < 0)] = 0
            result['snow1h'][t][np.where(result['snow1h'][t] < 0)] = 0
            # result['snow1h'][t] = snow1h
            # result['sr'][t] = sr
            result['wspdmax'][t] = wspd_max
    t = dts
    while t < time_max:
        h_offset_td = np.timedelta64(t - time_min + td_12h, 'h')
        h_offset = int(h_offset_td.astype(object).seconds / 3600)
        tt = np.logical_and(time_arr > t, time_arr <= (t + td_12h))
        # tt = np.logical_and(time_arr > t, time_arr <= (t + td_1h))
        t = t + td_1h
        if tt.max():
            precip12h = prate_ds[tt].sum(0)
            snow12h = snowrate_ds[tt].sum(0)
            sr12 = snow12h / (precip12h + 0.0001)
            solid_idx = np.where(sr12 > 0.8)
            liq_idx = np.where(sr12 <= 0.8)
            result['rain12h'][t] = precip12h-snow12h
            result['rain12h'][t][liq_idx] = precip12h[liq_idx]
            result['snow12h'][t] = np.zeros_like(snow12h)
            result['snow12h'][t][solid_idx] = snow12h[solid_idx]
            result['rain12h'][t][np.where(result['rain12h'][t] < 0)] = 0
            result['snow12h'][t][np.where(result['snow12h'][t] < 0)] = 0

    return result
