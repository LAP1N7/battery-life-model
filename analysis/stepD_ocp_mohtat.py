# -*- coding: utf-8 -*-
"""
OCP pair for the Mohtat (UMBL) graphite / NMC pouch cells.

Source: PyBaMM parameter set Mohtat2020, functions graphite_ocp_PeymanMPM and NMC_ocp_PeymanMPM
  repo pybamm-team/PyBaMM, file packages/pybamm/src/pybamm/input/parameters/lithium_ion/Mohtat2020.py,
  last commit touching the file 43fc19759669a6b6f4a58bbc3107e496afa3da66.
  (The old path src/pybamm/input/parameters/lithium_ion/Mohtat2020.py returns 404.)
Coefficients re-typed here and compared with the fetched file (identical).

Chemistry statement in that file (get_parameter_values docstring): "Parameters for a graphite/NMC532 pouch cell from the paper
Mohtat2020 and references therein". The OCP docstrings only cite "Peyman MPM manuscript (to be submitted)". The earlier
prototype notes NMC111 for the dataset paper's cell table; label mismatch (532 vs 111) is reported in D0_spec.md section 7.
"""
import numpy as np
from stepD_core import OCPPair, register_ocp


def Un_graphite(x):
    x = np.asarray(x, float)
    return (0.063 + 0.8 * np.exp(-75 * (x + 0.001))
            - 0.0120 * np.tanh((x - 0.127) / 0.016)
            - 0.0118 * np.tanh((x - 0.155) / 0.016)
            - 0.0035 * np.tanh((x - 0.220) / 0.020)
            - 0.0095 * np.tanh((x - 0.190) / 0.013)
            - 0.0145 * np.tanh((x - 0.490) / 0.020)
            - 0.0800 * np.tanh((x - 1.030) / 0.055))


def Up_nmc(y):
    y = np.minimum(np.asarray(y, float), 1.2)      # exp term explodes above ~1.2; the cell never goes there
    return (4.3452 - 1.6518 * y + 1.6225 * y ** 2 - 2.0843 * y ** 3 + 3.5146 * y ** 4
            - 2.2166 * y ** 5 - 0.5623e-4 * np.exp(109.451 * y - 100.006))


MOHTAT = register_ocp(OCPPair(
    name="mohtat_peyman",
    Up=Up_nmc, Un=Un_graphite,
    x_range=(-0.05, 1.10), y_range=(-0.05, 1.20),
    source="PyBaMM Mohtat2020.py @43fc19759669 (graphite_ocp_PeymanMPM, NMC_ocp_PeymanMPM)",
    chemistry="graphite / NMC (PyBaMM docstring: NMC532; dataset paper table: NMC111)",
    notes="Fit to the UMBL pouch-cell family itself (Peyman MPM manuscript, unpublished); absolute y offset is a parametrisation artefact.",
))
