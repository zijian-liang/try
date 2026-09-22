# Provenance and model conventions

This document consolidates the source notes supplied with the three experiment archives.

## Circuit and decoder sources

- The BB72 supports and seven CNOT layers follow `decoder_setup.py` from [IBM BivariateBicycleCodes](https://github.com/sbravyi/BivariateBicycleCodes), with `l=m=6`, `A=x^3+y+y^2`, and `B=y^3+x+x^2`. The same-shot Pauli failure criterion and two noiseless closing cycles follow its `decoder_run.py`. The associated paper is [High-threshold and low-overhead fault-tolerant quantum memory](https://arxiv.org/abs/2308.07915).
- The weight-7 supports were transcribed from the supplied `weight7_56_12_7_square_edge_lattice.pdf`. On `Z2 x Z14`, `HX=[A|B]`, `HZ=[B^T|A^T]`, `A=1+y+x*y^10`, and `B=1+y^3+x*y^11+x*y^2`. Each sector has rank 22. The supplied static-distance label is 7; circuit distance is a separate question. The reference PDF is not bundled here.
- Circuit sampling uses [Stim](https://github.com/quantumlib/Stim), and joint decoding uses [Tesseract](https://github.com/quantumlib/tesseract-decoder). Each experiment records its pinned dependency requirements.

The supplied IBM license text is retained once in [licenses/IBM_SOURCE_LICENSE.txt](../licenses/IBM_SOURCE_LICENSE.txt).

## Campaign distinctions

The ideal encoded Bell reference qubits let both commuting encoded Bell stabilizers be tracked. They are bookkeeping references, receive no noise, and do not count as physical hardware overhead.

The comparison campaign uses zero idle noise and `R=d`: 7, 6, and 7 noisy rounds for weight-7, BB72, and surface d=7, respectively. Two ideal closing cycles are additional. These conventions differ from the original IBM idle-noise setting; the saved plots should not be described as an exact reproduction of that paper's Figure 3. An upstream BB72 circuit-distance upper bound alone is not a proof of equality.

The schedule scan uses weight-7 only, `R=7`, and physical probabilities `0.003` and `0.0025`. Candidate orders come from a finite earlier schedule search. Historical search bounds do not transfer automatically to the present circuit. The distance gate checks the current schedule and noise convention, and identical labeled detector-error models are filtered; this is not an exhaustive classification of equivalent circuits.

The high-slope experiment selects candidates using historical observations, then records a separate confirmation run. Historical selection counts and independent confirmation counts are not pooled. A two-point slope does not establish an asymptotic exponent or a globally optimal schedule.

## Saved assets

The figures and result snapshots were supplied locally. They are copied without changing scientific values. [The data manifest](../data/manifest.json) records relative source paths, sizes, and SHA256 hashes for these copied assets. Documentation and layout were reorganized separately. No new Monte Carlo experiment or distance proof was run as part of this repository cleanup.
