# Deprecated Functions

Functions kept only for compatibility with older code

These functions still work but print a warning and forward to their replacement. Do not use them in new code.

| Function | Replacement | Deprecated in |
|---|---|---|
| `PYTHIA2(Z,Y,Ybin,Ybest,algolabels,opts)` | `PYTHIA` with `opts.classifier = 'knn'` | v0.9.0 |
| `PYTHIAtest(model,Z,Y,Ybin,Ybest,algolabels)` | `PYTHIA(Z,Y,Ybin,Ybest,algolabels,opts,model)` (evaluation mode) | v0.9.0 |
| `SIFTED2(X,Y,Ybin,featlabels,opts)` | `SIFTED` (renamed) | v0.9.0 |

## Updating Your Code

```matlab
% Before
out = PYTHIA2(Z, Y, Ybin, Ybest, algolabels, opts);
res = PYTHIAtest(out, Znew, Ynew, Ybinnew, Ybestnew, algolabels);

% After
opts.classifier = 'knn';
out = PYTHIA(Z, Y, Ybin, Ybest, algolabels, opts);
res = PYTHIA(Znew, Ynew, Ybinnew, Ybestnew, algolabels, opts, out);
```

## See Also

`PYTHIA` | `SIFTED` | [What's New](WhatsNew.html)
