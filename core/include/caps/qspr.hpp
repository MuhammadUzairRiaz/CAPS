// CAPS structure–property models: molecular descriptors, regression models with validation, clustering and
// discriminant analysis, for screening repeat units, curatives and coupling agents in a study table.
//
// Descriptors (hydrogen-suppressed graph unless said; δ_i heavy neighbours, h_i hydrogens):
//   Wiener W = ½ Σ_ij d_ij (Wiener, J. Am. Chem. Soc. 69, 17 (1947)); Zagreb M1 = Σ δ_i², M2 = Σ_bonds δ_i δ_j (Gutman &
//   Trinajstić, Chem. Phys. Lett. 17, 535 (1972)); connectivity mχ = Σ_paths Π (δ)^−½, orders 0–3 (Randić, J. Am. Chem.
//   Soc. 97, 6609 (1975)) and valence mχv with δv = (Z_v − h)/(Z − Z_v − 1) (Kier & Hall, J. Pharm. Sci. 65, 1806 (1976));
//   shape κ1 = A(A−1)²/P1², κ2 = (A−1)(A−2)²/P2², κ3 = (A−1)(A−3)²/P3² (A odd) or (A−3)(A−2)²/P3² (A even) (Kier, Quant.
//   Struct.-Act. Relat. 4, 109 (1985); 5, 1 (1986)); Balaban J = q/(μ+1) Σ_bonds (s_i s_j)^−½ (Chem. Phys. Lett. 89, 399
//   (1982)); path counts; rings, aromatic atoms; H-bond donors (H on N, O) and acceptors (N, O) as Lipinski et al. (Adv.
//   Drug Deliv. Rev. 23, 3 (1997)); rotatable bonds.
//   3D: van der Waals volume and area (Bondi radii), principal moments, radius of gyration; shadow areas on the three
//   principal planes and their ratios to the enclosing rectangles (Rohrbaugh & Jurs, Anal. Chim. Acta 199, 99 (1987));
//   charged partial surface areas PPSA/PNSA/DPSA/FPSA/FNSA/WPSA/WNSA 1–3 from the solvent-accessible areas and the charges
//   (Stanton & Jurs, Anal. Chem. 62, 2323 (1990)).
// Models (y on columns x):
//   MLR, optionally stepwise by AIC (Akaike, IEEE Trans. Autom. Control 19, 716 (1974)) over the columns, their squares
//   and products; PLS1 by NIPALS on autoscaled columns (Wold, Sjöström & Eriksson, Chemom. Intell. Lab. Syst. 58, 109
//   (2001)); GFA, a genetic search over sums of linear, spline ⟨x − t⟩₊ / ⟨t − x⟩₊ and quadratic terms scored by
//   Friedman's lack of fit LOF = LSE / (1 − (c + d p)/M)² (Rogers & Hopfinger, J. Chem. Inf. Comput. Sci. 34, 854 (1994));
//   a neural network of one tanh hidden layer trained by back-propagated gradients (Rumelhart, Hinton & Williams, Nature
//   323, 533 (1986)) with L-BFGS and weight decay.
//   Validation: k-fold or leave-one-out q² = 1 − PRESS/SS (the model's form refitted in each fold), and y-randomisation:
//   the whole fit, selection included, repeated on shuffled y (Rücker, Rücker & Meringer, J. Chem. Inf. Model. 47, 2345
//   (2007)). Outliers: leverage, standardised residuals and Cook's distance (Technometrics 19, 15 (1977)).
// Clustering: hierarchical by the Lance–Williams update (single, complete, average, Ward; Lance & Williams, Comput. J. 9,
//   373 (1967); Ward, J. Am. Stat. Assoc. 58, 236 (1963)) and k-means (Lloyd, IEEE Trans. Inf. Theory 28, 129 (1982));
//   Fisher's linear discriminant analysis (Ann. Eugenics 7, 179 (1936)).
// AlogP and molar refractivity are not given: their atom-type tables are not in CAPS.
#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

struct Descriptor {
  std::string name, unit, what;
  double value = std::numeric_limits<double>::quiet_NaN();
};

// The descriptors of a structure taken as one molecule (every atom; hydrogens explicit). Charges: the structure's own
// when it has them, else Gasteiger–Marsili (charges_used says which); none: the CPSA family is NaN and charges_used
// gives the reason. A repeat unit's [*] attachment points (element 0) are capped with H for the 3D descriptors and the
// charges (volume, areas, shadows, CPSA); mass, Rg and the graph indices use the structure as given.
std::vector<Descriptor> molecular_descriptors(const System& s, std::string* charges_used = nullptr);

// ---- models
struct ModelOptions {
  std::string method = "mlr";      // mlr, pls, gfa, nn
  bool stepwise = false;           // mlr: AIC stepwise over the candidate terms
  bool squares = false, products = false;   // mlr stepwise: add x², x_i x_j as candidates
  int components = 0;              // pls: 0 = chosen by cross-validation (highest q²)
  int gfa_population = 200, gfa_generations = 3000, gfa_max_terms = 6;
  double gfa_smoothing = 1.0;      // d in the LOF
  int nn_hidden = 4, nn_restarts = 5;
  double nn_decay = 1e-3;
  int folds = 0;                   // cross-validation: 0 none, the row count (or more) leave-one-out
  int randomise = 0;               // y-randomisation repeats
  uint64_t seed = 1;
};

struct ModelTerm {
  std::string name;                // "1", "x1", "x1^2", "x1*x2", "<x1-0.5>+"
  double coef = 0, se = std::numeric_limits<double>::quiet_NaN();
};

struct ModelResult {
  std::string method, equation;
  std::vector<ModelTerm> terms;    // the intercept first (mlr, gfa, pls in the original units)
  std::vector<double> fitted, residuals;
  double r2 = 0, r2_adj = 0, rmse = 0, f = std::numeric_limits<double>::quiet_NaN();
  int components = 0;              // pls
  double q2 = std::numeric_limits<double>::quiet_NaN(), rmse_cv = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> predicted_cv;
  std::vector<double> q2_by_components;   // pls
  double lof = std::numeric_limits<double>::quiet_NaN();   // gfa
  std::vector<std::pair<std::string, double>> gfa_models;  // gfa: the best equations and their LOF
  std::vector<double> random_r2, random_q2;                // y-randomisation
  std::vector<double> leverage, std_residual, cooks;       // linear models, per row
  std::vector<std::string> notes;
};

// y and the columns x (each as long as y) with their names.
ModelResult fit_model(const std::vector<double>& y, const std::vector<std::vector<double>>& x, const std::vector<std::string>& names,
                      const ModelOptions& o);

// ---- clustering and discriminant analysis
struct RowClusters {
  std::vector<int> label;                          // per row, 0-based
  std::vector<std::array<double, 3>> merges;       // hierarchical: (a, b, height); rows 0..n−1, merged clusters n, n+1, …
  double within_ss = 0;                            // Σ squared distances to the centroids (scaled columns)
  std::vector<std::vector<double>> centroids;      // in the original units
};
// method: single, complete, average, ward (hierarchical, cut at k clusters) or kmeans; rows are x[r] (each a vector of
// the columns); the columns autoscaled when scale. Ward's heights are the increases of the within-cluster sum of squares.
RowClusters cluster_rows(const std::vector<std::vector<double>>& rows, const std::string& method, int k, bool scale = true, uint64_t seed = 1,
                           int restarts = 10);

struct DiscriminantResult {
  std::vector<std::string> classes;
  std::vector<std::vector<double>> directions;     // discriminant functions (coefficients of the columns), most separating first
  std::vector<double> eigenvalues;
  std::vector<int> predicted, predicted_loo;       // class index per row: training and leave-one-out
  double accuracy = 0, accuracy_loo = 0;
};
DiscriminantResult discriminant(const std::vector<std::vector<double>>& rows, const std::vector<std::string>& cls);

Json model_json(const ModelResult& r);
Json descriptors_json(const std::vector<Descriptor>& d);

}  // namespace caps
