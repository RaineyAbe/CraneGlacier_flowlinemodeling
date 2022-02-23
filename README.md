# CraneGlacier_flowlinemodeling

Code package used to model Crane Glacier response to the 2002 Larsen B ice shelf collapse and future changes in climate

Rainey Aberle and Ellyn Enderlin


![](figures/studyArea.png)

---

### Order of operations (see `workflows/`)
- `0_synthesize_datasets`
- `1_steady-state`
- `2_tune_parameters`
- `3_sensitivity_tests`

---

### Contact

Rainey Aberle (raineyaberle@u.boisestate.edu)

---

### Datasets not included in this repository which are used in the workflow

Velocity maps:  

- NASA ITS\_LIVE (1999-2018) (Gardner et al., 2021)
- ERS-derived (1994-5) (Wuite et al., 2015)

Surface elevation:

- GTOPO30 (~1995) (Gesch et al., 1999)
- NASA OIB L2 (2009-11, 2016-18)

Terminus positions:

- Landsat-derived (Dryak and Enderlin, 2020)

Bed elevation: 

- NASA OIB L1B, manually digitized picks using the CReSIS toolbox (Paden)

---

### References
Dryak, M. C., and Enderlin, E. M. (2020). Analysis of Antarctic Peninsula glacier frontal ablation rates with respect to iceberg melt-inferred variability in ocean conditions. _Journal of Glaciology_. _66_(257): 457-470. https://doi.org/10.1017/JOG.2020.21

Gardner, A. S., M. A. Fahnestock, and Scambos, T. A. (2021). MEaSUREs ITS_LIVE Landsat Image-Pair Glacier and Ice Sheet Surface Velocities: Version 1. https://doi.org/10.5067/IMR9D3PEI28U

Gesch, D. B., Verdin, K. L., and Greenlee, S. K. (1999). New land surface digital elevation model covers the Earth. _Eos Trans. AGU_. _80_(6), 69–70, https://doi.org/10.1029/99EO00050

Wuite, J., Rott, H., Hetzenecker, M., Floricioiu, D., De Rydt, J., Gudmundsson, G. H., Nagler, T., and Kern, M. (2015). Evolution of surface velocities and ice discharge of Larsen B outlet glaciers from 1995 to 2013. _9_(3): 957-969. _The Cryosphere._ https://doi.org/10.5194/tc-9-957-2015

