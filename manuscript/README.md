# DeBrota--Stacey, Section III, Conjecture 1 after Theorem 2

This package treats DeBrota--Stacey, Section III, Conjecture 1 following Theorem 2,
Physical Review A **100**, 062327 (2019), arXiv:1907.10999v2.
The main result characterizes the MIC attaining the Haar-averaged output-entropy
bound as a SIC. The note also proves a weighted-two-design version and quantitative
entropy-deficit bounds. For arbitrary finite POVMs, equality characterizes
rank-one weighted projective two-designs, with at least d^2 effects and SICs
at exactly d^2. Further corollaries recover the depolarizing endpoint and the
previously proved entropy-exchange equality case, with an exact deficit formula.
The postscript's Pauli-covariance statement for unbiased rank-one qubit MICs
also has a complete Bloch-vector proof.

An explicit anisotropic unbiased qubit MIC shows that a two-design alone need
not reproduce a Haar mean of the largest output eigenvalue. The postscript's
saturation hypothesis restores the constant spectrum needed for its averaging
step; the example is not a counterexample to its conclusion under that hypothesis.

## Files

- `note.tex`: self-contained mathematical text with its bibliography.
- `metadata.json`: title, creators, license, problem identifier and mathematical scope.
- `verify.py`: finite exact examples and algebraic checks.
- `checker-output.json`: expected finite-check output.
- `references.json`: bibliographic identifiers and mathematical locators.
- `LICENSE.txt`: CC BY 4.0 for the original material.
- `SHA256SUMS`: checksums of the other source-package members.

The PDF is distributed alongside the source archive.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively use a TeX distribution providing the standard packages named in
`note.tex`:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is embedded; BibTeX is unnecessary.

## Finite exact checks

With Python 3.10 or later and SymPy 1.14.0:

```sh
python3 -m pip install sympy==1.14.0
python3 verify.py > checker-output.json
```

The checks construct the one-dimensional SIC, tetrahedral qubit SIC and
nine-vector Hesse qutrit SIC. They test positivity through Hermitian rank-one
projectors, normalization, every Gram entry, linear independence, channel action
on every matrix unit, and the exact projective two-design identity. They also
test 1,092 integer spectra in dimensions 1 through 6 and a full-rank qubit MIC
with a strictly positive entropy deficit. Additional exact checks cover a six-effect
qubit design, its Choi state, and the anisotropic qubit MIC's Gram determinant,
Bloch action and tetrahedral squared output radii. The strict spectral-mean
comparison follows from the written strict-concavity proof, not finite sampling.
The Python-version field in the output
reflects the interpreter used; the mathematical fields are unchanged.

The all-dimensional theorem and all corollaries have complete written proofs in
`note.tex`. In the pinned Lean file linked in the note, `result : claim` expresses
the implication from entropy saturation to existence of a SIC; the file-local
`depolarizing_identity` and `sic_coefficients` express the channel identity and
the structure of the original MIC under their displayed hypotheses.
The finite examples illustrate these statements and are not a substitute for
their proofs.

## Creators and license

Haobo Ma (ChronoAI Pte Ltd, Singapore), followed by Wenlin Zhang (National
University of Singapore; ChronoAI Pte Ltd, Singapore). Original material is
licensed under Creative Commons Attribution 4.0 International, `cc-by-4.0`.
