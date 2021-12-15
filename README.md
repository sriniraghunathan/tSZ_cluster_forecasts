# tSZ cluster forecasts for future CMB surveys
 * Paper 1 (2107.10250):
   * Reference: Raghunathan, Whitehorn, Alvarez et al. (2021); [arXiv:2107.10250](https://arxiv.org/abs/2107.10250)
   * Thermal SZ selected cluster catalogue forecasts for [CMB-S4](https://arxiv.org/abs/1907.04473) and [CMB-HD](https://arxiv.org/abs/1906.10134). 
   * Astrophysical and cosmological constraints with tSZ selected clusters.
* Paper 2 (2112.07656):
   * Reference: Raghunathan (2021); [arXiv:2112.07656](https://arxiv.org/abs/2112.07656)
   * Assesing the importance of noise from Thermal SZ signals for cluster surveys and cluster cosmology.
   * Experiments considered: [SPT](https://pole.uchicago.edu/public/Home.html) (SPT-SZ, SPTpol, SPT-3G); [Simons Observatory](https://arxiv.org/abs/1808.07445) (SO-Baseline and SO-Goal); [CMB-S4](https://arxiv.org/abs/1907.04473) (S4-Wide and S4-Ultra deep surveys); and [CMB-HD](https://arxiv.org/abs/1906.10134). 
   * Cosmological constraints with tSZ selected clusters, primary CMB, and tSZ power spectrum.

<!-- # Coming soon. Please check back on July 26, 2021 -->

## Binned cluster counts

## Fisher matrices
* Paper 1 [(2107.10250)](https://arxiv.org/abs/2107.10250):
  * Parameters constrained: Cosmology + <img src="https://render.githubusercontent.com/render/math?math=Y_{\rm SZ}-M"> scaling relation + cluster virialisation model - 16 parameters.
  * Cosmology - 6+2 parameters: <img src="https://render.githubusercontent.com/render/math?math=\Lambda CDM, \sum m_{\nu}, w_{\rm DE} ">
  * <img src="https://render.githubusercontent.com/render/math?math=Y_{\rm SZ}-M"> scaling relation - 6 parameters: <img src="https://render.githubusercontent.com/render/math?math=\alpha_{\rm Y}, \beta_{\rm Y}, \gamma_{\rm Y}, \sigma_{\rm logY}, \alpha_{\sigma}, \gamma_{\sigma}">
  * Cluster virialisation model - 2 parameters:
    * Model 1: <img src="https://render.githubusercontent.com/render/math?math={\rm v}(z) = \eta_{\rm v}(z) (1 - b_{\rm HSE})^{\alpha_{Y}}">
    * Model 2: <img src="https://render.githubusercontent.com/render/math?math={\rm v}(z) = A_{\rm v} {\rm ln}(1">+<img src="https://render.githubusercontent.com/render/math?math=z)"> + <img src="https://render.githubusercontent.com/render/math?math=B_{\rm v}">
     
  * Data products: [2107.10250/data_products/fisher](https://github.com/sriniraghunathan/tSZ_cluster_forecasts/tree/main/2107.10250/data_products/fisher)
  * Script to read Fisher matrices: [2107.10250/read_fisher_mat.ipynb](https://github.com/sriniraghunathan/tSZ_cluster_forecasts/blob/main/2107.10250/read_fisher_mat.ipynb)

* Paper 2 [(2112.07656)](https://arxiv.org/abs/2112.07656):
  * Parameters constrained: Cosmology + <img src="https://render.githubusercontent.com/render/math?math=Y_{\rm SZ}-M"> scaling relation - 15 parameters.
  * Cosmology - 6+2 parameters: <img src="https://render.githubusercontent.com/render/math?math=\Lambda CDM, \sum m_{\nu}, w_{\rm DE} ">
  * <img src="https://render.githubusercontent.com/render/math?math=Y_{\rm SZ}-M"> scaling relation - 6 parameters: <img src="https://render.githubusercontent.com/render/math?math=\alpha_{\rm Y}, \beta_{\rm Y}, \gamma_{\rm Y}, \sigma_{\rm logY}, \alpha_{\sigma}, \gamma_{\sigma}">
  * Data products:
    * ILC weights and residuals: [2112.07656/ilc_weights_residuals](https://github.com/sriniraghunathan/tSZ_cluster_forecasts/tree/main/2112.07656/ilc_weights_residuals)
    * Fisher matrices: [2112.07656/Fisher matrices](https://github.com/sriniraghunathan/tSZ_cluster_forecasts/tree/main/2112.07656/Fisher)
    * Simulation products: [2112.07656/Simulation products](https://github.com/sriniraghunathan/tSZ_cluster_forecasts/tree/main/2112.07656/sim_products)
