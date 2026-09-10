# Deprecated

This page combines three deprecated backward-compatibility wrapper functions located in the `deprecated/` directory. New code should not call any of these functions directly.

## PYTHIA2

**Signature:** `out = PYTHIA2(Z, Y, Ybin, Ybest, algolabels, opts)`

**Replacement Guidance:** Use `PYTHIA` with `opts.classifier = 'knn'`.

**Context:** This is a backward-compatibility wrapper; all new code should call `PYTHIA` directly with the desired classifier type set in `opts.classifier`.

## PYTHIAtest

**Signature:** `out = PYTHIAtest(model, Z, Y, Ybin, Ybest, algolabels)`

**Replacement Guidance:** Use `PYTHIA` with a trained model (7-arg eval mode).

**Context:** This exists only for backward compatibility. New code calling `exploreIS` should use `PYTHIA` eval mode directly.

## SIFTED2

**Signature:** `[X, out] = SIFTED2(X, Y, Ybin, featlabels, opts)`

**Replacement Guidance:** Use `SIFTED`.

**Context:** SIFTED2 was renamed to SIFTED, with this thin alias kept in its place. All new code should call `SIFTED` directly.
