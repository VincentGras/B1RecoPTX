# B1RecoPTX

Reconstruction of complex B1 maps for pTx RF coils from pre-saturated TurboFLASH (satTFL) data and optimization of interferometric modes.

## Requirements

- Optimization Toolbox

## Scripts

- **CreateScheme.m** Generates an interferometric encoding scheme from a reference B1+ distribution. The normalized interferometric modes are optimized, and the default nominal flip angles (α = 100°, β = 10°) are written to the scheme .txt file.

- **OptimFA.m** Optimizes the nominal flip angles for a given encoding scheme using the reference B1+ distribution. If enabled, the optimized flip angles replace the default values in the scheme file.

- **TestRecoSphere.m** Evaluates the reconstruction algorithm using simulation data from a spherical phantom. Synthetic satTFL signals with additive Gaussian noise are generated, and the reconstruction algorithm is applied to the simulated data. The resulting B1+ maps are then compared with the ground truth. Available reconstruction algorithms are xfl_proposed_reco (FitM), xfl_simple_reco (standard hybrid method), and xfl_pairwiseCP_reco (standard pairwise method), located in the Functions subfolder. Example encoding schemes are provided in the Schemes subfolder.

- **TestRecoInVivo.m** Demonstrates B1 map reconstruction from in vivo multi-Tx/multi-Rx data. Both transmit (B1+) and receive (B1-) sensitivities are reconstructed. The example scheme is Schemes/InVivo_vcc03.xflparam.txt.

- **SchemeCharacteristics.m** Allows checking the actual flip-angle distributions for a given encoding scheme and RF pulse integrals using a reference B1+ distribution.
