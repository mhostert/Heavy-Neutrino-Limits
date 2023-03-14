# Heavy Neutrino Limits

This package and repo tracks the constraints on the coupling and masses of heavy neutral leptons (HNL). To access it, simply clone this repository. 

To use the underlying functions to create your own plots and latex files, you can install HNLimits as you would a standard PyPI package, ``python -m pip install HNLimits``. You can tehn import it in your Python script ``import HNLimits``.

* All limits are kept track in this [![Google Spreadsheets](https://img.shields.io/badge/Google_Sheets-Database-brightgreen.svg)](https://docs.google.com/spreadsheets/d/1p_fslIlThKMOThGl4leporUsogq9TmgXwILntUZOscg/edit?usp=sharing)

* We consider single flavor dominance scenarios, where the HNL mixes predominantly with either the electron, muon, or tau flavor. 

* Following the [accompanying paper](www.arxiv.org/abs/XXXXXXX), we provide limits on the Wilson coefficients of dimesnion-six $\nu$ SMEFT operators as a function of the HNL mass $m\_{N}$.

* So far, the code keeps track of limits in the region 1 MeV $< m_{N} < 100$ GeV.

Additions, comments, or suggestions should be directed to:
   * Josu Hernández-García (garcia.josu.hernandez@ttk.elte.hu)
* Matheus Hostert (mhostert@pitp.ca)

--- 
**Citation info:**

```
@article{Fernandez-Martinez2023:TBA, 

}
```

---
## Limits on the dimension-six $\nu$SMEFT operators

| Type                 | Name                                      | Operator                                                                                              | Notebook                               | Figure                                                                                                                                                                                                                                                                                                                                                 |
|----------------------|-------------------------------------------|-------------------------------------------------------------------------------------------------------|----------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Mixing               | $O\_{\rm mixing}$                                        | $\overline{L}\tilde{H}N$                                                                              | ``0_mixing.ipynb``                     | [electron](https://github.com/mhostert/N-SMEFT-Limits/blob/main/plots/mixing/UeN_majorana.pdf)  <br />  [muon](https://github.com/mhostert/N-SMEFT-Limits/blob/main/plots/mixing/UmuN_majorana.pdf)  <br />  [tau](https://github.com/mhostert/N-SMEFT-Limits/blob/main/plots/mixing/UtauN_majorana.pdf) |
| Higgs-dressed Mixing | $\mathcal{O}_{\rm LNH}^\alpha$            | $\overline{L}\tilde{H}N (H^\dagger H)$                                                                | ``1_NSMEFT_LHN.ipynb``                 | [electron]()  <br />  [muon]()  <br />  [tau]()                                                                                                                                                                                                                                                          |
| Bosonic Currents     | $\mathcal{O}_{\rm HN}$                    | $\overline{N}\gamma^\mu N (H^\dagger i \overleftrightarrow{D}_\mu H)$                                 | ``2_NSMEFT_bosonic_NC.ipynb``          | [NC bosonic]()                                                                                                                                                                                                                                                                                                                                         |
|                      | $\mathcal{O}_{\rm HN\ell}^{\alpha}$       | $\overline{N}\gamma^\mu \ell_\alpha (\tilde{H}^\dagger i \overleftrightarrow{D}_\mu H)$               | ``3_NSMEFT_bosonic_CC.ipynb``          | [CC bosonic]()                                                                                                                                                                                                                                                                                                                                         |
| Moments              | $\mathcal{O}_{\rm NB}^{\alpha}$             | $\left(\overline{L}\_\alpha \sigma\_{\mu\nu} N \right)\widetilde{H} B^{\mu \nu }$                        | ``4_NSMEFT_moment_NB.ipynb``           | [Moment hypercharge]()                                                                                                                                                                                                                                                                                                                                 |
|                      | $\mathcal{O}_{\rm NW}^{\alpha}$             | $\left(\overline{L}\_\alpha \sigma\_{\mu \nu } N\right) \tau^a \widetilde{H} W^{\mu \nu }_{a}$               | ``5_NSMEFT_moment_NW.ipynb``           | [Moment W]()                                                                                                                                                                                                                                                                                                                                           |
| Neutral Currents     | $\mathcal{O}_{\rm ff}$                    | $(\overline{f} \gamma^\mu f) (\overline{N} \gamma^\mu N)$                                             | ``6_NSMEFT_4fermion_NC_ff.ipynb``      | [Four fermion ff]()                                                                                                                                                                                                                                                                                                                                    |
|                      | $\mathcal{O}_{\rm LN}^\alpha$             | $(\overline{L}\_\alpha \gamma^\mu L\_\alpha) (\overline{N} \gamma^\mu N)$                               | ``7_NSMEFT_4fermion_NC_LN.ipynb``      | [Four fermion LN]()                                                                                                                                                                                                                                                                                                                                    |
|                      | $\mathcal{O}_{\rm QN}$                    | $\mathcal{Z}^{\rm QN}\_{ij}(\overline{Q}\_i \gamma^\mu Q\_j) (\overline{N} \gamma^\mu N)$                | ``8_NSMEFT_4fermion_NC_QN.ipynb``      | [Four fermion QN]()                                                                                                                                                                                                                                                                                                                                    |
| Charged Currents     | $\mathcal{O}_{\rm LNL\ell}^{\alpha\beta}$ | $(\overline{L}\_\alpha N)\epsilon (\overline{L}\_\alpha \ell\_\beta)$                                    | ``9_NSMEFT_4fermion_CC_LNLell.ipynb``  | [Four fermion LNLell]()                                                                                                                                                                                                                                                                                                                                |
|                      | $\mathcal{O}_{\rm duN\ell}^{\alpha}$      | $\mathcal{Z}\_{ij}^{\rm duN\ell}(\overline{d}\_i \gamma^\mu u\_j) (\overline{N} \gamma^\mu \ell_\alpha)$ | ``10_NSMEFT_4fermion_NC_duNell.ipynb`` | [Four fermion duNell]()                                                                                                                                                                                                                                                                                                                                |
|                      | $\mathcal{O}_{\rm LNQd}^\alpha$           | $\mathcal{Z}^{\rm LNQd}\_{ij} (\overline{L}\_\alpha N)\epsilon (\overline{Q_i} d\_j)$                    | ``11_NSMEFT_4fermion_NC_LNQd.ipynb``   | [Four fermion LNQd]()                                                                                                                                                                                                                                                                                                                                  |
|                      | $\mathcal{O}_{\rm QuNL}^\alpha$           | $\mathcal{Z}^{\rm QuNL}\_{ij}(\overline{Q}\_i u\_j)(\overline{N} L\_\alpha)$                              | ``12_NSMEFT_4fermion_NC_QuNL.ipynb``   | [Four fermion QuNL]()                                                                                                                                                                                                                                                                                                                                  |




---
## Limits on the mixing

Current compilation of bounds on the mixing of HNLs with a single flavor:


[Electron mixing](https://raw.githubusercontent.com/mhostert/N-SMEFT-Limits/main/plots/mixing/UeN_majorana.png)
![e flavor](https://raw.githubusercontent.com/mhostert/N-SMEFT-Limits/main/plots/mixing/UeN_majorana.png)

[Muon mixing](https://raw.githubusercontent.com/mhostert/N-SMEFT-Limits/main/plots/mixing/UmuN_majorana.png)
![mu flavor](https://raw.githubusercontent.com/mhostert/N-SMEFT-Limits/main/plots/mixing/UmuN_majorana.png)

[Tau mixing](https://raw.githubusercontent.com/mhostert/N-SMEFT-Limits/main/plots/mixing/UtauN_majorana.png)
![tau flavor](https://raw.githubusercontent.com/mhostert/N-SMEFT-Limits/main/plots/mixing/UtauN_majorana.png)
