# tSZ cluster forecasts for future CMB surveys
 * Paper 1:
  * Reference: Raghunathan, Whitehorn, Alvarez et al. (2021); [arXiv:2107.10250](https://arxiv.org/abs/2107.10250)
  * Thermal SZ selected cluster catalogue forecasts for [CMB-S4](https://arxiv.org/abs/1907.04473) and [CMB-HD](https://arxiv.org/abs/1906.10134). 
  * Astrophysical and cosmological constraints with tSZ selected clusters.
* Paper 2:

<!-- # Coming soon. Please check back on July 26, 2021 -->

## Binned cluster counts

## Fisher matrices
 * Parameters constrained: Cosmology + <img src="https://render.githubusercontent.com/render/math?math=Y_{\rm SZ}-M"> scaling relation + cluster virialisation model - 16 parameters.
   * Cosmology - 6+2 parameters: <img src="https://render.githubusercontent.com/render/math?math=\Lambda CDM, \sum m_{\nu}, w_{\rm DE} ">
   * <img src="https://render.githubusercontent.com/render/math?math=Y_{\rm SZ}-M"> scaling relation - 6 parameters: <img src="https://render.githubusercontent.com/render/math?math=\alpha_{\rm Y}, \beta_{\rm Y}, \gamma_{\rm Y}, \sigma_{\rm logY}, \alpha_{\sigma}, \gamma_{\sigma}">
   * Cluster virialisation model - 2 parameters:
     * Model 1: <img src="https://render.githubusercontent.com/render/math?math={\rm v}(z) = \eta_{\rm v}(z) (1 - b_{\rm HSE})^{\alpha_{Y}}">
     * Model 2: <img src="https://render.githubusercontent.com/render/math?math={\rm v}(z) = A_{\rm v} {\rm ln}(1">+<img src="https://render.githubusercontent.com/render/math?math=z)"> + <img src="https://render.githubusercontent.com/render/math?math=B_{\rm v}">
     
 * Data products: [2107.10250/data_products/fisher](https://github.com/sriniraghunathan/tSZ_cluster_forecasts/tree/main/2107.10250/data_products/fisher)
 * Script to read Fisher matrices: [2107.10250/read_fisher_mat.ipynb](https://github.com/sriniraghunathan/tSZ_cluster_forecasts/blob/main/2107.10250/read_fisher_mat.ipynb)
