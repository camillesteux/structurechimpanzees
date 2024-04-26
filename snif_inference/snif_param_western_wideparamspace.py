import sys
import os
os.chdir('~/snif/snif-main')
sys.path.append('~/snif/snif-main')

from snif import *

gentime = 25
mut = 1.5e-8
clist = [4, 5, 6, 7, 8]

for c in clist:
    for w in [0.5, 1]:
            h_inference_parameters = InferenceParameters(
                data_source = '~/gape-psmc/FINAL_CLEAN_Pan_troglodytes_verus-X00100_Koby.psmc',
                source_type = SourceType.PSMC,
                IICR_type = IICRType.Exact,
                ms_reference_size = 5000,
                ms_simulations = int(1e5),
                psmc_mutation_rate = mut,
                psmc_number_of_sequences = 100,
                psmc_length_of_sequences = int(2e6),
                infer_scale = True,
                data_cutoff_bounds = (1e4/gentime, 1e7/gentime),
                data_time_intervals = 64,
                distance_function = ErrorFunction.ApproximatePDF,
                distance_parameter = w,
                distance_max_allowed = 7000,
                distance_computation_interval = (1e4/gentime, 1e7/gentime),
                rounds_per_test_bounds = (25, 50),
                repetitions_per_test = 5,
                number_of_components = c,
                bounds_islands = (2, 100),
                bounds_migrations_rates = (0.05, 100),
                bounds_deme_sizes = (1, 1),
                bounds_event_times = (1e4/gentime, 1e7/gentime),
                bounds_effective_size = (10, 20000)
            )

            h_settings = Settings(
                static_library_location ='~/snif/snif-main/libs/libsnif.so',
                custom_filename_tag = '',
                output_directory = '~/snif/Chimps/verus/verus_koby',
                default_output_dirname = '~/snif/default_directory/'
            )

            infer(inf = h_inference_parameters, settings = h_settings)
