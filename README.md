# Heavy Neutrino Limits

This package and repo tracks the constraints on the coupling and masses of heavy neutral leptons (HNL). To access it, simply clone this repository.

To use the underlying functions to create your own plots and latex files, you can install HNLimits as you would a standard PyPI package, ``python -m pip install HNLimits``. You can then import it in your Python script ``import HNLimits``.

All limits are kept track in this [![Google Spreadsheets](https://img.shields.io/badge/Google_Sheets-Database-brightgreen.svg)](https://docs.google.com/spreadsheets/d/1TIpmkgOa63-8Sy75qh0YutI5XdRtiClU3aquUdmjqpc/edit?usp=sharing).

* We consider single flavor dominance scenarios, where the HNL mixes predominantly with either the electron, muon, or tau flavor.

* Following the [accompanying paper](https://arxiv.org/abs/2304.06772), we provide limits on the Wilson coefficients of dimesnion-six $\nu$ SMEFT operators as a function of the HNL mass $M\_{N}$.

* The code keeps track of limits in the region 1 MeV $< M_{N} < 100$ GeV.

Additions, comments, or suggestions should be directed to:
* Josu Hernández-García (garcia.josu.hernandez@ttk.elte.hu)
* Matheus Hostert (mhostert@g.harvard.edu)

--- 
**Citation info:**

``arxiv.org/abs/2304.06772``

```
@article{Fernandez-Martinez:2023phj,
    author = "Fern\'andez-Mart\'\i{}nez, Enrique and Gonz\'alez-L\'opez, Manuel and Hern\'andez-Garc\'\i{}a, Josu and Hostert, Matheus and L\'opez-Pav\'on, Jacobo",
    title = "{Effective portals to heavy neutral leptons}",
    eprint = "2304.06772",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "FTUV-23-0303.1224, IFIC/23-09",
    doi = "10.1007/JHEP09(2023)001",
    journal = "JHEP",
    volume = "09",
    pages = "001",
    year = "2023"
}
```


---
## Limits on the mixing

Current compilation of bounds on the mixing of HNLs with a single flavor:


[Electron mixing](https://raw.githubusercontent.com/mhostert/N-SMEFT-Limits/main/plots/mixing/UeN_majorana.png)
![e flavor](https://raw.githubusercontent.com/mhostert/N-SMEFT-Limits/main/plots/mixing/UeN_majorana.png)

[Muon mixing](https://raw.githubusercontent.com/mhostert/N-SMEFT-Limits/main/plots/mixing/UmuN_majorana.png)
![mu flavor](https://raw.githubusercontent.com/mhostert/N-SMEFT-Limits/main/plots/mixing/UmuN_majorana.png)

[Tau mixing](https://raw.githubusercontent.com/mhostert/N-SMEFT-Limits/main/plots/mixing/UtauN_majorana.png)
![tau flavor](https://raw.githubusercontent.com/mhostert/N-SMEFT-Limits/main/plots/mixing/UtauN_majorana.png)


---
## Limits on the dimension-six $\nu$SMEFT operators

| Type                 | Name                                      | Operator                                                                                              | Notebook                               | Figure                                                                                                                                                                                                                                                                                                                                                 |
|----------------------|-------------------------------------------|-------------------------------------------------------------------------------------------------------|----------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Mixing               | $O\_{\rm mixing}$                                        | $\overline{L}\widetilde{H}N$                                                                              | [``0_limits_mixing.ipynb``](https://github.com/mhostert/Heavy-Neutrino-Limits/blob/5299bb1705eb20c4d38e4fb784b0d79805d698b7/0_limits_mixing.ipynb)                     | [electron](https://github.com/mhostert/N-SMEFT-Limits/blob/main/plots/mixing/UeN_majorana.pdf)  <br />  [muon](https://github.com/mhostert/N-SMEFT-Limits/blob/main/plots/mixing/UmuN_majorana.pdf)  <br />  [tau](https://github.com/mhostert/N-SMEFT-Limits/blob/main/plots/mixing/UtauN_majorana.pdf) |
| Higgs-dressed Mixing | $\mathcal{O}_{\rm LNH}^\alpha$            | $\overline{L}\widetilde{H}N (H^\dagger H)$                                                                | [``1_NSMEFT_LHN.ipynb``](https://github.com/mhostert/Heavy-Neutrino-Limits/blob/9319051131e9882ae451405bfc668275c44d5ec8/1_NSMEFT_LHN.ipynb)                 | [electron]()  <br />  [muon]()  <br />  [tau]()                                                                                                                                                                                                                                                          |
| Bosonic Currents     | $\mathcal{O}_{\rm HN}$                    | $\overline{N}\gamma^\mu N (H^\dagger i \overleftrightarrow{D}_\mu H)$                                 | [``2_NSMEFT_bosonic_NC.ipynb``](https://github.com/mhostert/Heavy-Neutrino-Limits/blob/f8abe77c5588d3891425ce6519fb3ac918982ec4/2_NSMEFT_bosonic_NC.ipynb)          | [Bosonic NC]()                                                                                                                                                                                                                                                                                                                                         |
|                      | $\mathcal{O}_{\rm HN\ell}^{\alpha}$       | $\overline{N}\gamma^\mu \ell_\alpha (\widetilde{H}^\dagger i \overleftrightarrow{D}_\mu H)$               | [``3_NSMEFT_bosonic_CC.ipynb``](https://github.com/mhostert/Heavy-Neutrino-Limits/blob/eb39159704ff02b221ff2793d39f9ff7faae5146/3_NSMEFT_bosonic_CC.ipynb)          | [Bosonic CC]()                                                                                                                                                                                                                                                                                                                                         |
| Moments              | $\mathcal{O}_{\rm NB}^{\alpha}$             | $\left(\overline{L}\_\alpha \sigma\_{\mu\nu} N \right)\widetilde{H} B^{\mu \nu }$                        | [``4_NSMEFT_moment_NB.ipynb``](https://github.com/mhostert/Heavy-Neutrino-Limits/blob/eb39159704ff02b221ff2793d39f9ff7faae5146/4_NSMEFT_moment_NB.ipynb)           | [Moment hypercharge]()                                                                                                                                                                                                                                                                                                                                 |
|                      | $\mathcal{O}_{\rm NW}^{\alpha}$             | $\left(\overline{L}\_\alpha \sigma\_{\mu \nu } N\right) \tau^a \widetilde{H} W^{\mu \nu }_{a}$               | [``5_NSMEFT_moment_NW.ipynb``](https://github.com/mhostert/Heavy-Neutrino-Limits/blob/f8abe77c5588d3891425ce6519fb3ac918982ec4/5_NSMEFT_moment_NW.ipynb)           | [Moment W]()                                                                                                                                                                                                                                                                                                                                           |
| 4-fermion NC     | $\mathcal{O}_{\rm ff}$   | $(\overline{f} \gamma^\mu f) (\overline{N} \gamma_\mu N)$ | [``6_NSMEFT_4fermion_NC_ff_LN.ipynb``](https://github.com/mhostert/Heavy-Neutrino-Limits/blob/4f1280cf6d940cf75e5384d54f849bc9c377cac0/6_NSMEFT_4fermion_NC_ff.ipynb)      | [4-fermion ee]()  <br /> [4-fermion uu]()  <br /> [4-fermion dd]()                                                                                                                                                                                                                                                                                                                            |
|                      | $\mathcal{O}_{\rm LN}^\alpha$             | $(\overline{L}\_\alpha \gamma^\mu L\_\alpha) (\overline{N} \gamma_\mu N)$                               | [``6_NSMEFT_4fermion_NC_ff_LN.ipynb``](https://github.com/mhostert/Heavy-Neutrino-Limits/blob/4f1280cf6d940cf75e5384d54f849bc9c377cac0/6_NSMEFT_4fermion_NC_ff.ipynb)      | [4-fermion LN]()                                                                                                                                                                                                                                                                                                                                    |
|                      | $\mathcal{O}_{\rm QN}$                    | $(\overline{Q}\_i \gamma^\mu Q\_i) (\overline{N} \gamma_\mu N)$                | [``7_NSMEFT_4fermion_NC_QN.ipynb``](https://github.com/mhostert/Heavy-Neutrino-Limits/blob/f8abe77c5588d3891425ce6519fb3ac918982ec4/7_NSMEFT_4fermion_NC_QN.ipynb)      | [4-fermion QN]()                                                                                                                                                                                                                                                                                                                                    |
| 4-fermion NC & CC     | $\mathcal{O}_{\rm LNL\ell}^{\alpha\beta}$ | $(\overline{L}\_\alpha N)\epsilon (\overline{L}\_\alpha \ell\_\beta)$                                    | [``8_NSMEFT_4fermion_LNLell.ipynb``](https://github.com/mhostert/Heavy-Neutrino-Limits/blob/f8abe77c5588d3891425ce6519fb3ac918982ec4/8_NSMEFT_4fermion_LNLell.ipynb)  | [4-fermion LNLell NC]()  <br />  [4-fermion LNLell CC]()                                                                                                                                                                                                                                                                                                                                |
| 4-fermion CC               | $\mathcal{O}_{\rm duN\ell}^{\alpha}$      | $\mathcal{Z}\_{ij}^{\rm duN\ell}(\overline{d}\_i \gamma^\mu u\_j) (\overline{N} \gamma_\mu \ell_\alpha)$ | [``9_NSMEFT_4fermion_CC_duNell.ipynb``](https://github.com/mhostert/Heavy-Neutrino-Limits/blob/f8abe77c5588d3891425ce6519fb3ac918982ec4/9_NSMEFT_4fermion_CC_duNell.ipynb) | [4-fermion duNell: electron]() <br /> [4-fermion duNell: muon]() <br /> [4-fermion duNell: tau]() <br />                                                                                                                                                                                                                                                                                                                                |
|                      | $\mathcal{O}_{\rm LNQd}^\alpha$           | $\mathcal{Z}^{\rm LNQd}\_{ij} (\overline{L}\_\alpha N)\epsilon (\overline{Q_i} d\_j)$                    | [``10_NSMEFT_4fermion_CC_LNQd.ipynb``](https://github.com/mhostert/Heavy-Neutrino-Limits/blob/f8abe77c5588d3891425ce6519fb3ac918982ec4/10_NSMEFT_4fermion_CC_LNQd.ipynb)   | [4-fermion LNQd: electron]()  <br /> [4-fermion LNQd: muon]()  <br /> [4-fermion LNQd: tau]()  <br />                                                                                                                                                                                                                                                                                                                                |
|                      | $\mathcal{O}_{\rm QuNL}^\alpha$           | $\mathcal{Z}^{\rm QuNL}\_{ij}(\overline{Q}\_i u\_j)(\overline{N} L\_\alpha)$                              | [``11_NSMEFT_4fermion_CC_QuNL.ipynb``](https://github.com/mhostert/Heavy-Neutrino-Limits/blob/f8abe77c5588d3891425ce6519fb3ac918982ec4/11_NSMEFT_4fermion_CC_QuNL.ipynb)   | [4-fermion QuNL: electron]()   <br /> [4-fermion QuNL: muon]()   <br /> [4-fermion QuNL: tau]()   <br />                                                                                                                                                                                                                                                                                                                                 |
