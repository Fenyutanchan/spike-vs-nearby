# Data Sources

This directory contains the published, binned cosmic-ray lepton fluxes used by
the numerical analysis. The CSV files preserve the units and quoted precision
of their collaboration sources. They are external inputs and should not be
edited in place during fitting; analysis cuts and unit conversions belong in
the data-loading code.

## AMS-02 electron and positron fluxes

- Local files: `ams02_e-.csv` and `ams02_e+.csv`
- Measurement: time-averaged electron and positron fluxes collected from May
  2011 through May 2018.
- Publication: M. Aguilar et al. (AMS Collaboration), *The Alpha Magnetic
  Spectrometer (AMS) on the International Space Station: Part II — Results from
  the First Seven Years*, Physics Reports 894 (2021) 1–116.
- DOI: <https://doi.org/10.1016/j.physrep.2020.09.003>
- Official data page: <https://ams02.space/publications/202102>
- Official positron CSV (Table 1):
  <https://ams02.space/sites/default/files/publication/202102/table-1.csv>
- Official electron CSV (Table 2):
  <https://ams02.space/sites/default/files/publication/202102/table-2.csv>

The files are direct copies of the official AMS-02 CSV tables. Energies are in
GeV and fluxes are in $({\rm m^2\,sr\,s\,GeV})^{-1}$. Statistical and systematic
uncertainties remain in separate columns.

## DAMPE combined electron and positron flux

- Local derived table: `dampe_e±.csv`
- Online primary source: <https://arxiv.org/pdf/1711.10981v1>
- Online source archive: <https://arxiv.org/src/1711.10981v1>
- Measurement: combined $e^-+e^+$ flux from data collected between
  2015-12-27 and 2017-06-08.
- Publication: G. Ambrosi et al. (DAMPE Collaboration), *Direct detection of a
  break in the teraelectronvolt cosmic-ray spectrum of electrons and
  positrons*, Nature 552 (2017) 63–66.
- DOI: <https://doi.org/10.1038/nature24475>
- arXiv record: <https://arxiv.org/abs/1711.10981v1>

`dampe_e±.csv` is a faithful transcription of the 38 published rows in Table 1
of `DAMPE-2017CRE-final.tex`. It covers energy-bin edges from 24.0 to
4570.9 GeV and retains the characteristic-energy uncertainty, acceptance,
event count, background fraction, flux, and separate statistical and
systematic flux uncertainties. The multiplicative powers of ten printed in the
LaTeX table have only been expanded into standard CSV scientific notation; the
background fraction remains in percent as stated in the column header.

The paper states that the published spectrum and its uncertainties are
available in Table 1; event-level DAMPE data are not public and other supporting
data are available from the collaboration on reasonable request.

### Visual check

From the `code/` directory, reproduce the published table directly from the
CSV with:

```sh
julia 'data/ext/dampe_e±_check.jl'
```

The script validates the schema, row count, numerical values, and contiguous
energy bins before printing a formatted character table to the terminal. An
alternative input CSV may be supplied as the first command-line argument.

## Integrity checks

SHA-256 digests of the collaboration-provided inputs:

```text
63d868da6dccd76e1ecef8b36bbcf70a4a748c8369b3b1983510800e94cc4628  ams02_e+.csv
a01ff72cbdd696c6ea28a0d5046c8a0d3774f4a7c1e3291b267fd2489dd70a36  ams02_e-.csv
2a7ecc972e9b8b4b784e00be5195375f639bbada1acb76199d141a46d8f9ebe7  <https://arxiv.org/src/1711.10981v1>
33adb229febbee116332b37b0588c5a23b4c0d7c1f0d037673cfe33491bb6322  `DAMPE-2017CRE-final.tex` in <https://arxiv.org/src/1711.10981v1>
```
