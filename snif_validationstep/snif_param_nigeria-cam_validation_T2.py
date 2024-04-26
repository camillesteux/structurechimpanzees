import sys
import os
os.chdir('~/snif/snif-main-new')
sys.path.append('~/snif/snif-main-new')

from snif import *

gentime = 25
mut = 1.5e-8

for c in [7]:
    for i in [1,2,3,4,5]:
        h_inference_parameters = InferenceParameters(
            data_source = '~/snif/Chimps/ellioti_12-23/inf/231213-172353_c07YN_dA_w050_SQ_Koto_inferred-model-005_T2.ms',
            source_type = SourceType.MSCommand,
            IICR_type = IICRType.T_sim,
            ms_reference_size = 1161.93,
            ms_simulations = int(1e6),
            psmc_mutation_rate = mut,
            psmc_number_of_sequences = 100,
            psmc_length_of_sequences = int(2e6),
            infer_scale = True,
            data_cutoff_bounds = (1e4/gentime, 1e7/gentime),
            data_time_intervals = 64,
            distance_function = ErrorFunction.ApproximatePDF,
            distance_parameter = 0.5,
            distance_max_allowed = 7000,
            distance_computation_interval = (1e4/gentime, 1e7/gentime),
            rounds_per_test_bounds = (25, 50),
            repetitions_per_test = 10,
            number_of_components = c,
            bounds_islands = (2, 30),
            bounds_migrations_rates = [(0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50)],
            bounds_deme_sizes = [(1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1)],
            bounds_event_times = [(1e4/gentime, 1e5/gentime), (1.5e5/gentime, 3e5/gentime), (3.5e5/gentime, 5e5/gentime), (8e5/gentime, 1.5e6/gentime), (2e6/gentime, 7e6/gentime), (2e6/gentime, 7e6/gentime)],
            bounds_effective_size = (10, 10000)
            )

        h_settings = Settings(
            static_library_location = '~/snif/snif-main-new/libs/libsnif.so',
            custom_filename_tag = "".join(['val_c7eb_Kotorep5_12-23_c', str(c), "w050_ind", str(i)]),
            output_directory = '~/snif/Chimps/ellioti_12-23/inf/validation/',
            default_output_dirname = 'snif_res_default'
        )
        infer(inf = h_inference_parameters, settings = h_settings)
