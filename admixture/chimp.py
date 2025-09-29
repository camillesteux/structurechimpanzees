#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: rtournebize
"""

import demes, msprime

def samples(ARGS):
    Samples, PopLabels = [], []
    for x in ARGS["samples"]:
        x = x.split(":")
        assert len(x)==4, "Sample definition should have 4 `:`-delimited values, DemeName:nDiploids:TimeYa:Label"
        assert all(c.isalpha for c in x[0]), "For sample definition, you must provide the name of the deme to sample in."
        Samples.append( msprime.SampleSet(  num_samples = int(x[1]),
                                            population = x[0],
                                            time = float(x[2]) / ARGS["generation_time"],
                                            ploidy = 2) )
        PopLabels += [x[3]] * int(x[1]) * 2
    return Samples, PopLabels

def history(Args, Pars):
    Demography = demes.load("chimp.yaml")
    Demography = msprime.Demography.from_demes(Demography)
    Demography.generation_time = Args["generation_time"]

    Samples, PopLabels = samples(Args)
    return Args, Demography, Samples, PopLabels

#___
