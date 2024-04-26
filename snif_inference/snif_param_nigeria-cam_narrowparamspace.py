import argparse
import sys
import os
import numpy as np

# To run snif without having the parameter file in snif-main directory
os.chdir('~/snif/snif-main-new')
sys.path.append('~/snif/snif-main-new')

from snif import *
from snif2tex import *

parser = argparse.ArgumentParser(description='Run SNIF')
#mandatory arguments
parser.add_argument("--tag", help="custom tag for the output files")
parser.add_argument("--psmc", help="path to the psmc curve file")
args = parser.parse_args()

print(args.tag)
print(args.psmc)

# SNIF PARAM
gentime = 25 #generation time
mut = 1.5e-8 #mutation rate

c = 7 # nb of components
w = 0.5 # w distance
islands = (2, 30) # nb of islands
demesize = (10, 10000) # deme size
timewindows = [(1e4/gentime, 1e5/gentime), (1.5e5/gentime, 3e5/gentime), (3.5e5/gentime, 5e5/gentime), (8e5/gentime, 1.5e6/gentime), (2e6/gentime, 6e7/gentime), (2e6/gentime, 7e6/gentime)]
mig = [(0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50), (0.05, 50)]
chsize = [(1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1)]
t = 10 # nb of repetitions
rd = (25, 50) #nb of rounds per repetitions
tag = args.tag
psmc = args.psmc


#INFERENCES
h_inference_parameters = InferenceParameters(
    data_source = psmc,
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
    rounds_per_test_bounds = rd,
    repetitions_per_test = t,
    number_of_components = c,
    bounds_islands = islands,
    bounds_migrations_rates = mig,
    bounds_deme_sizes = chsize,
    bounds_event_times = timewindows,
    bounds_effective_size = demesize
)

h_settings = Settings(
    static_library_location = '~/snif/snif-main-new/libs/libsnif.so',
    custom_filename_tag = tag,
    output_directory = '~/snif/Chimps/ellioti_12-23/',
    default_output_dirname = '~/default_directory'
)

basename = infer(inf = h_inference_parameters, settings = h_settings)
print("Basename: ", basename)

#PLOTS
config = Configuration(
    SNIF_basename = basename,
    plot_width_cm = 13,
    plot_height_cm = 6,
    IICR_plots_style  = OutputStyle.Full,
    PDF_plots_style = OutputStyle.Excluded,
    CDF_plots_style = OutputStyle.Excluded,
    islands_plot_style = OutputStyle.Full,
    Nref_plot_style = OutputStyle.Full,
    test_numbers = "all",
    one_file_per_test = False,
    versus_plot_style = OutputStyle.Full,
    CG_style = OutputStyle.Full,
    CG_size_history = False,
    CG_source = "",
    CG_source_legend = "",
    Nref_histograms_bins = 100,
    islands_histograms_bins = 100,
    time_histograms_bins = 100,
    migration_histograms_bins = 100,
    size_histograms_bins = 100,
    scaling_units = TimeScale.Years,
    generation_time = gentime,
)

TeXify(config)
