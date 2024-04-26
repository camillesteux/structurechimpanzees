import sys
import os
os.chdir('~/snif/snif-main-new')
sys.path.append('~/snif/snif-main-new')

from snif import *
from snif2tex import *

parser = argparse.ArgumentParser(description='Param file to run SNIF on simulated PSMC (03-24)')
parser.add_argument("-i", "--individual", help="id of the individual")
args = parser.parse_args()

gentime = 25
mut = 1.5e-8

for c in [7]:
    for w in [0.5]:
        h_inference_parameters = InferenceParameters(
            data_source = ''.join(['~/snif/Chimps/val_psmc/input_data/230323-172956_c07YN_dA_w050_SQ_c7eb2_A911_Kidongo_10x100Mb-', args.individual, '.psmc']),
            source_type = SourceType.PSMC,
            IICR_type = IICRType.Exact,
            ms_reference_size = 723,
            ms_simulations = int(1e6),
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
            number_of_components = c,
            bounds_islands = (2, 50),
            bounds_migrations_rates = [(0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50)],
            bounds_deme_sizes = [(1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1)],
            bounds_event_times = [(5e4/gentime, 2e5/gentime), (4e5/gentime, 1.5e6/gentime), (4e5/gentime, 1.5e6/gentime), (2e6/gentime, 6e6/gentime), (2e6/gentime, 6e6/gentime), (7e6/gentime, 1e7/gentime)],
            bounds_effective_size = (10, 10000)
            )

        h_settings = Settings(
            static_library_location = '~/snif/snif-main-new/libs/libsnif.so',
            custom_filename_tag = ''.join(['schwein_valpsmc_c', str(c), 'eb_w', str(w), '_', args.individual]),
            output_directory = '~/snif/Chimps/val_psmc/val_psmc_res/',
            default_output_dirname = './snif_res_default'
        )

        basename = infer(inf = h_inference_parameters, settings = h_settings)

        config = Configuration(
            SNIF_basename = basename,
            plot_width_cm = 13,
            plot_height_cm = 6,
            IICR_plots_style  = OutputStyle.Full,
            PDF_plots_style = OutputStyle.Excluded,
            CDF_plots_style = OutputStyle.Excluded,
            islands_plot_style = OutputStyle.Excluded,
            Nref_plot_style = OutputStyle.Excluded,
            test_numbers = "all",
            one_file_per_test = False,
            versus_plot_style = OutputStyle.Excluded,
            CG_style = OutputStyle.Full,
            CG_size_history = False,
            CG_source = '~/snif/Chimps/val_psmc/input_data/230323-172956_c07YN_dA_w050_SQ_c7eb2_A911_Kidongo_inferred-model-007.json',
            CG_source_legend = "Target",
            Nref_histograms_bins = 100,
            islands_histograms_bins = 100,
            time_histograms_bins = 100,
            migration_histograms_bins = 100,
            size_histograms_bins = 100,
            scaling_units = TimeScale.Years,
            generation_time = gentime
        )

        TeXify(config)
