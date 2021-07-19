class Settings:
    def __init__(self):
        self.shape_root = '/home/gonchukov-lv/src/phenomenaAreas/shp'
        self.tiff_root = '/home/gonchukov-lv/src/phenomenaAreas/tiffs'
        self.tiff_root_1b = '/home/gonchukov-lv/src/phenomenaAreas/tiffs_1b_lead'
        self.regions = {2746: 'east',
                        2747: 'center',
                        2748: 'west',
                        233: 'south',
                        2579: 'prim'}
        self.region_shapes = {
            2746: '{0}/prim_east.shp'.format(self.shape_root),
            2747: '{0}/prim_center.shp'.format(self.shape_root),
            2748: '{0}/prim_west.shp'.format(self.shape_root),
            233: '{0}/prim_south.shp'.format(self.shape_root)
        }
        # self.delta_hours = 1
        # self.variables =   {1516: 'tmin1h', 1517: 'tmax1h', 1515: 'spdmax', 1250: 'precip1h', 1023: 'precipPhase',
        # 16034 : 'sr' }
        self.variables = {'tmin': 1516, 'tmax': 1517, 'wspdmax': 1515, 'precip1h': 1250, 'precipPhase': 1023,
                          'sr': 16034}
        self.vars_to_wrf = {1516: ['T2MIN'], 1517: ['T2MAX'], 1515: ['SPDUV10MAX'], 1250: ['RAINC', 'RAINNC']}

        self.phenomenas = {'very_strong_wind': 11, #+
                           'hurricane': 12,  #+
                           'very_heavy_rain': 64, #+
                           'strong_rainfall': 82,
                           'very_heavy_snowfall': 86, #+
                           'heat': 21,
                           'frost': 25,
                           'strong_wind': 10, #+
                           'heavy_rain': 63,
                           'rainfall': 81,
                           'heavy_snowfall': 85}
        self.phenomena_colors = {
            'very_strong_wind':  [255, 69, 0, 256],  # + #  	orange red
             'hurricane': [255, 69, 0, 256],  # + #  	orange red
             'very_heavy_rain': [0, 0, 139, 256],  # +# dark blue
             'strong_rainfall': [139, 0, 139, 256], #dark magenta
             'very_heavy_snowfall': [0, 128, 0, 256],  # +# green
             'heat': [255,0,0,255], #red
             'frost': [255,0,0,255], #red
             'strong_wind': [255, 165, 0, 256],  # + #orange
             'heavy_rain': [0, 0, 255, 256],  # +#  blue,
             'rainfall': [255,0,255,255], # Magenta
             'heavy_snowfall': [0,255,0,256] # green

        }

        self.mark_phen_operator = {'tmin': 'lt', 'tmax': 'gt', 'wspdmax': 'gt', 'rain1h': 'gt', 'rain12h': 'gt',
                                   'snow12h': 'gt'}
        self.aggregators = {1516: 'min', 1517: 'max', 1515: 'max'}
        # self.corrections = {1516: -293.16, 1517: -243.16, 1515: 15}
        # var->region->value
        self.criterias = {
            'very_strong_wind': {'wspdmax':
                {
                    233: 35,
                    2748: 25,
                    2747: 25,
                    2746: 35
                }
            },
            'frost': {'tmin':
                {
                    233: -35,
                    2748: -40,
                    2747: -43,
                    2746: -35
                }
            },
            'heat': {'tmax':
                {
                    233: 33,
                    2748: 35,
                    2747: 37,
                    2746: 37
                }
            },
            'hurricane': {'wspdmax': 33},
            'strong_wind': {'wspdmax': 15},
            'strong_rainfall': {'rain1h': 30},
            'rainfall': {'rain1h': 15},
            'very_heavy_rain': {'rain12h': 50},
            'very_heavy_snowfall': {'snow12h': 20},
            'heavy_rain': {'rain12h': 15},
            'heavy_snowfall': {'snow12h': 6}
        }
