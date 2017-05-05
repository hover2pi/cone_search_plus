from distutils.extension import Extension

def get_package_data():
    return {'all_sky_target_tool': ['data/*', 'data/radii/*']}
