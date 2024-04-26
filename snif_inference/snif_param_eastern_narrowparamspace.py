import sys
import os
os.chdir('~/snif/snif-main')
sys.path.append('~/snif/snif-main')

from snif import *

gentime = 25
mut = 1.5e-8
ind = ["A910_Bwambale", "A911_Kidongo", "A912_Nakuu"]

for i in ind:
    for w in [0.5]:
            h_inference_parameters = InferenceParameters(
                data_source = "".join(['/ifs/igc/folders/PCG/csteux/data/gape-psmc/FINAL_CLEAN_Pan_troglodytes_schweinfurthii-', i, '.psmc']),
                source_type = SourceType.PSMC,
                IICR_type = IICRType.Exact,
                ms_reference_size = 5000,
                ms_simulations = int(1e5),
                psmc_mutation_rate = mut,
                psmc_number_of_sequences = 100,
                psmc_length_of_sequences = int(2e6),
                infer_scale = True,
                data_cutoff_bounds = (1e4/gentime, 2e7/gentime),
                data_time_intervals = 64,
                distance_function = ErrorFunction.ApproximatePDF,
                distance_parameter = w,
                distance_max_allowed = 7000,
                distance_computation_interval = (1e4/gentime, 2e7/gentime),
                rounds_per_test_bounds = (25, 50),
                repetitions_per_test = 10,
                number_of_components = 7,
                bounds_islands = (2, 50),
                bounds_migrations_rates = [(0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50)],
                bounds_deme_sizes = [(1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1)],
                bounds_event_times = [(5e4/gentime, 2e5/gentime), (4e5/gentime, 1.5e6/gentime), (4e5/gentime, 1.5e6/gentime), (2e6/gentime, 6e6/gentime), (2e6/gentime, 6e6/gentime), (7e6/gentime, 1e7/gentime)],
                bounds_effective_size = (10, 10000)
            )

            h_settings = Settings(
                static_library_location = '~/snif/snif-main/libs/libsnif.so',
                custom_filename_tag = "".join(['c7eb2_', i]),
                output_directory = '~/snif/Chimps/schweinfurthii/schwein_all',
                default_output_dirname = '~/snif/default_directory/'
            )

            infer(inf = h_inference_parameters, settings = h_settings)
