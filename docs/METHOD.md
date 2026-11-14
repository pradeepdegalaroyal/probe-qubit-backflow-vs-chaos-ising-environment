# Method

## Model

**Environment.** Open mixed-field Ising chain of L spins,

    H_E = J sum_i Z_i Z_{i+1} + sum_i hx_i X_i + hz sum_i Z_i,   J = 1,

with transverse fields hx_i = 0.9045 + uniform(-0.1, 0.1). The weak disorder removes the
reflection symmetry and keeps the chain free-fermionic at hz = 0. The longitudinal field hz is
swept over nine values from 0 (integrable) to 0.809 (chaotic).

**Probe.** One qubit coupled by pure dephasing, H = H_E + lambda sigma_z B, where B is a Pauli
operator on the central spin. A probe prepared in |+> keeps the coherence

    c(t) = <psi| exp(+i (H_E - lambda B) t) exp(-i (H_E + lambda B) t) |psi>,

a Loschmidt echo of the environment state psi. It is computed exactly from two dense
diagonalisations (`chaosprobe.probe.loschmidt_coherence`), on a graphics processing unit (GPU)
when one is present.

**Environment state.** Either a Haar-random state or the product state with all spins along +y.
Both have the infinite-temperature energy density for every hz, so a difference between them is
not an energy effect.

## Statistics

- **Information backflow.** Breuer-Laine-Piilo (BLP) measure for the optimal state pair: the sum
  of all increases of |c(t)| on the time grid (0 to 40/J, step 0.1/J).
- **Late-time plateau.** Mean of |c(t)|^2 over the second half of the record.
- **Finite-shot versions.** Ramsey readout of sigma_x and sigma_y with S shots each at 41 time
  points (step 1/J). Naive backflow; thresholded backflow (only the part of each rise above two
  shot-noise standard deviations counts); plateau with the 2/S bias removed.

## Baselines

- **Classical baseline.** Mean adjacent-gap ratio of the full spectrum of H_E from exact
  diagonalisation (0.386 for Poisson statistics, 0.531 for the Gaussian orthogonal ensemble). It
  is the ground truth for "chaotic" and costs a full diagonalisation.
- **Standard quantum baseline.** Out-of-time-order correlator (OTOC)
  F(t) = <psi| W(t) V W(t) V |psi> with W = X on the first spin and V = X on the last spin,
  averaged over the second half of the record.
- **Memoryless reference.** A Markovian dephasing model has monotone |c(t)| and therefore zero
  backflow; any measured backflow is a memory effect or a finite-size or shot-noise artefact.

## Experiments

| Label | Question | Sweep |
|---|---|---|
| E1 | Do the statistics follow the gap ratio? | hz in 9 values, L = 8 and 10, 5 seeds |
| E2 | How do they scale with size? | L = 6 to 12 at hz = 0 and 0.809, 5 seeds (3 at L = 12) |
| E3 | Ablation | coupling X or Z, lambda in 0.2, 0.4, 0.8, state Haar or product |
| E4 | Noise and robustness | shots 30 to 10000, extra probe dephasing 0 to 0.2 J, 40 repetitions per seed |

Metrics: Spearman rank correlation with the gap ratio, contrast (integrable mean divided by
chaotic mean), area under the receiver operating characteristic curve (AUC) between the two end
points, exponential decay rate per site, and the smallest shot budget with AUC of at least 0.95.
All values are reported as mean ± standard deviation over seeds. A seed fixes the disorder
realisation and the random state; the same disorder is used at every hz.

## Validation

`tests/test_core.py` checks the Hamiltonian against explicit Kronecker products, the coherence and
the OTOC against matrix exponentials, the gap ratio against Poisson and random-matrix values, and
the shot model against its large-shot limit.
