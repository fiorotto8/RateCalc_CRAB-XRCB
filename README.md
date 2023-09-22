# Overview

- `Extracted_NuStar_power.csv`
  - Crab nebula Xray spectra from NuStar experiment
  - sampled with webplotdigitizer from [https://iopscience.iop.org/article/10.1088/0004-637X/801/1/66/pdf]
- `FuncPar_XRCB.txt`
  - Parameters of the double power fit of the Xray Cosmic Background
  - from [https://www.sciencedirect.com/science/article/pii/S0370269320304275?ref=pdf_download&fr=RR-2&rr=80a30983af794c4a]
- `HeCF4_6040_NIST.txt`
  - Mass attenuation factors of He/CF4 60/40
  - from [https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html]
- `param.txt`
  - contains detector dimensions, energy scan range and step and gas densities
  - gas densities are from quick google search

The output is a ROOT file with all the TGraph for fluxes and sensitivities. The total integrated flux (in Hz) is printed at the end of the code.

## Requirements

- pandas
- numpy
- argparse
- ROOT
