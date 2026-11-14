# Research gap

## State of the art

Open-quantum-system theory offers a precise notion of memory. Breuer et al. (2016,
doi:10.1103/RevModPhys.88.021002) define non-Markovian dynamics through the revival of the
trace distance between two system states, the Breuer-Laine-Piilo (BLP) information backflow, and
point out that the memory of a small open system carries information about its environment, so
the system can serve as a sensor of it. Tanimura (2020, doi:10.1063/5.0011599) reviews the hierarchical equations of motion
(HEOM), which treat such memory numerically exactly, but for environments made of harmonic
oscillators. Neither work asks what the memory measured by a probe says when the environment is an
interacting many-body system that is either integrable or quantum chaotic.

On the quantum-information side, Lu and Fang (2026, doi:10.1038/s42005-026-02655-y) estimate
relative entropies on quantum computers and Wang et al. (2026, doi:10.1038/s42005-026-02698-1)
recover quantum information from completely depolarising channels; both treat noise as a channel
with no internal many-body structure. The remaining recent papers in the brief study receivers and
channel models under imperfect channel knowledge, which is the classical counterpart of the same
question: how much can be inferred about a channel from a noisy, finite record.

## Related work table

| Paper | Year | Role | What it does | What it does NOT cover |
|---|---|---|---|---|
| Breuer et al., doi:10.1103/RevModPhys.88.021002 | 2016 | landmark | Defines and quantifies non-Markovian memory; proposes probes of the environment | Interacting many-body environments; integrable versus chaotic baths; shot cost |
| Tran et al., doi:10.1038/s41586-019-0975-z | 2019 | landmark | Observes moiré excitons in a twisted heterobilayer | Any dynamical or information-theoretic diagnostic of the many-body state |
| Tanimura, doi:10.1063/5.0011599 | 2020 | review | Numerically exact non-Markovian dynamics with harmonic baths | Spin-chain baths; link between memory and spectral statistics |
| Lu and Fang, doi:10.1038/s42005-026-02655-y | 2026 | adjacent | Variational estimation of quantum relative entropy and channel capacity | Time-resolved memory; structured environments |
| Wang et al., doi:10.1038/s42005-026-02698-1 | 2026 | adjacent | Transmits quantum information through depolarising channels by superposing orders | Noise that carries information about a many-body bath |
| Ren et al., doi:10.1038/s41377-025-02173-6 | 2026 | adjacent | Entanglement-assisted communication over 20 km fibre | Environment diagnostics; memory effects |
| Devipriya, doi:10.1016/j.aeue.2026.156615 | 2026 | adjacent | Variational-circuit reliability scores for a receiver under imperfect channel knowledge | Exact baselines; quantum environments |
| Kumar et al., doi:10.1186/s13634-026-01327-x | 2026 | adjacent | Learned detector under channel-estimation mismatch | Quantum probes; finite-shot statistics of a quantum readout |
| Pan et al., doi:10.1007/s44163-026-02499-z | 2026 | adjacent | Channel-aware denoising for time-varying channels | Physical model of the channel memory |
| Zhang et al., doi:10.1186/s13634-026-01314-2 | 2026 | adjacent | Resource allocation under stochastic channel models | Quantum channels |
| Gowthaman et al., doi:10.1038/s41598-026-44620-3 | 2026 | adjacent | Fuzzy-rule channel allocation from sensed channel statistics | Quantum channels; exact reference solutions |

## The gap

The landmark definition of memory (information backflow) and the review's exact methods exist, and
the recent work treats noisy channels as objects to be characterised from finite data. None of
them tests the natural many-body claim that follows from the "environment probe" idea: that the
memory seen by a single dephasing qubit distinguishes an integrable from a quantum-chaotic
environment. For a pure-dephasing probe the coherence is exactly a Loschmidt echo of the
environment, so the claim can be tested without approximation at simulation scale.

**Novelty in one sentence.** This study measures, with exact dynamics of 6-12 spin environments,
whether BLP information backflow of one probe qubit tracks the integrable-to-chaotic crossover of
a mixed-field Ising chain, and what it costs in shots and probe dephasing compared with the
late-time coherence plateau and with an out-of-time-order correlator (OTOC).

## Hypotheses (fixed before the experiments were run)

- **H1.** With the environment in a random (Haar) state, backflow anticorrelates with the mean
  gap ratio along the sweep: Spearman rank correlation at most -0.8 at 10 spins.
- **H2.** At the chaotic point, backflow decays exponentially with chain length at a rate of at
  least 0.25 per site; at the integrable point it decays at least three times more slowly.
- **H3.** With at most 1000 shots per quadrature and time point and extra probe dephasing up to
  0.05 J, the measured statistic separates integrable from chaotic environments with an area under
  the receiver operating characteristic curve (AUC) of at least 0.95.
- **H4.** A coupling that breaks the free-fermion structure of the integrable chain (Pauli Z)
  loses the contrast that an integrability-preserving coupling (Pauli X) has.

**Success metric.** Spearman correlation with the gap ratio, AUC between the two end points, and
the shot budget and dephasing rate at which AUC reaches 0.95. A refuted hypothesis is reported as
such.

**Scale.** Exact simulation, 6-12 spins, one probe qubit. No claim is made about larger systems
beyond the measured finite-size trends.
