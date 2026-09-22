# Provenance

The BB72 code and its seven CNOT layers reproduce the uploaded IBM BivariateBicycleCodes `decoder_setup.py` schedules, with l=m=6, A=x^3+y+y^2, B=y^3+x+x^2. The same-shot joint Pauli failure criterion and two noiseless closing cycles follow its `decoder_run.py`. The paper is https://arxiv.org/abs/2308.07915 . The source repository is https://github.com/sbravyi/BivariateBicycleCodes . IBM's published BB72 circuit-distance entry is an upper bound <=6, not a certified equality.

The weight-7 supports were read from the user-supplied `weight7_56_12_7_square_edge_lattice.pdf`: on Z2 x Z14, HX=[A|B], HZ=[B^T|A^T], A=1+y+x*y^10, B=1+y^3+x*y^11+x*y^2. The rank is 22 in each sector. The supplied static-distance label is 7; circuit distance is audited separately.

Joint decoding uses Tesseract: https://github.com/quantumlib/tesseract-decoder . Circuit fault sampling uses Stim: https://github.com/quantumlib/Stim . The simulator's ideal reference qubits permit both commuting encoded Bell stabilizers to be tracked, equivalent to the IBM classical Pauli-frame joint failure test; they are not physical hardware overhead and receive no noise.

Default idle noise is zero by user preference. IBM's original idle noise is p. The parent comparison campaign used R=d: 7, 6, 7 noisy rounds for weight-7, BB72, and surface7. This separate schedule-scan package runs only weight-7, always R=7 plus two ideal closing cycles, at p=0.003 and 0.0025. Do not label its results as an exact replication of Figure 3.

The original selected eight-layer W7 circuit is preserved as w7_baseline. Candidate seed orders come from the earlier finite schedule search, with historical information retained only as search provenance in seed_schedules.json. Every candidate must undergo the new R7 distance gate. Identical labeled detector-error models are filtered; this is not an exhaustive classification of circuit equivalences. Alternative schedules leave check supports, logical labels, number of noisy primitive locations, and the decoder settings unchanged.
