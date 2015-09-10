# Target specifications as of Jay Anderson's report
jay_lists = {
    'list0_all_lt_k_8.5': {'k_mag': (-5, 8.5)},
    'list1_early_commissioning': {
        'k_mag': (4.5, 5.5),
        'neighbors': (
            {'delta_k': 5, 'r_arcmin': 16},
        )
    },
    'list2_global_alignment': {
        'k_mag': (4.5, 5.5),
        'neighbors': (
            {'delta_k': 7, 'r_arcmin': 2.5},
        )
    },
    'list3_pil_imaging': {
        'k_mag': (6.5, 7.5),
        'neighbors': (
            {'delta_k': 7, 'r_arcmin': 2.5},
        )
    },
    'list4_coarse_phasing': {
        'k_mag': (4.5, 5.5),
        'neighbors': (
            {'delta_k': 7, 'r_arcmin': 2.5},
        )
    },
    'list5_fine_phasing': {
        'k_mag': (8.5, 9.5),
        'neighbors': (
            {'delta_k': 7, 'r_arcmin': 2.5},
        )
    },
    'list6_routine_maintenance': {
        'k_mag': (8.5, 9.5),
        'neighbors': (
            {'delta_k': 7, 'r_arcmin': 0.5},
            {'delta_k': 2, 'r_arcmin': 2.5}
        )
    },
}

# Target specifications as of WFSCOWG May 2015

target_lists = {
    'initial_image_mosaic': {
        'k_mag': (-5, 5), # there are some really bright stars in K (~ -3), make sure we get them
        'neighbors': (
            {'delta_k': 5, 'r_arcmin': 30},
            {'delta_k': 0.5, 'r_arcmin': 45},
        ),
    },
    'early_commissioning': {
        'k_mag': (4.5, 5.5),
        'neighbors': (
            {'delta_k': 5, 'r_arcmin': 16},
        ),
    },
    'global_alignment': {
        'k_mag': (4.5, 5.5),
        'neighbors': (
            {'delta_k': 7, 'r_arcmin': 2.5},
        )
    },
    'coarse_phasing': {
        'k_mag': (8.5, 9.5),
        'neighbors': (
            {'delta_k': 5, 'r_arcmin': 4.027},
        )
    },
    'fine_phasing_routine_maintenance': {
        'k_mag': (8.5, 9.5),
        'neighbors': (
            {'delta_k': 7, 'r_arcmin': 0.5,},
        ),
        # "neighbors in SI FoV must be dimmer than target"
        # NIRCam FoV 2.16' x 2.16'
        'no_brighter_neighbors_r_arcmin': 3.05,
    },
    'mimf_miri': {
        'k_mag': (10, 12),
        'neighbors': (
            {'delta_k': 2.5, 'r_arcmin': 2},
        )
    },
}
