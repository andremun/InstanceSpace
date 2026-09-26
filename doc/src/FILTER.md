# FILTER

Remove near-duplicate instances to even out the instance density

<!-- opts: selvars -->

## Syntax

```
[subsetIndex,isDissimilar,isVISA,unif] = FILTER(X,Y,Ybin,opts)
```

## Description

`[subsetIndex,isDissimilar,isVISA,unif] = FILTER(X,Y,Ybin,opts)` marks instances as redundant when they are closer than `opts.mindistance` in feature space to an instance already kept, and they also meet the condition in `opts.type`. Removing them gives a smaller, more evenly spread set of instances, which reduces the bias of densely sampled regions on the projection and the footprints (Alipour et al., 2023).

`InstanceSpace` calls FILTER when `opts.selvars.densityflag` is `true`: once on all features during PRELIM, and again after SIFTED on the selected features only.

## Examples

### Filter the reference data

```matlab
rootdir = 'test/data/example/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile('test/data/metadata.csv', rootdir);
opts.perf = struct('MaxPerf', false, 'AbsPerf', true, 'epsilon', 0.20);
obj = InstanceSpace(rootdir, opts);
obj = obj.build('stages', {'prelim'});
d = obj.model.data;

fopts = struct('mindistance', 0.5, 'type', 'Ftr&Good');
[redundant, ~, ~, unif] = FILTER(d.X, d.Y, d.Ybin, fopts);
fprintf('Kept %d of %d instances, uniformity %.2f\n', sum(~redundant), numel(redundant), unif)
```

### Filter inside the pipeline

```matlab
obj.opts.selvars.densityflag = true;
obj.opts.selvars.mindistance = 0.5;
obj = obj.build();
```

## Input Arguments

### `X` — Feature matrix

*numeric matrix*

Preprocessed features, one row per instance.

### `Y` — Performance matrix

*numeric matrix*

### `Ybin` — Good-performance labels

*logical matrix*

### `opts` — Filter options

*structure*

Normally `obj.opts.selvars`.

#### `opts.mindistance` — Distance threshold

*`0.10` (default) | positive scalar*

Instances closer than this in feature space are candidates for removal.

#### `opts.type` — Removal condition

*`'Ftr&Good'` (default) | `'Ftr'` | `'Ftr&AP'` | `'Ftr&AP&Good'`*

Condition, in addition to feature closeness, for an instance to be removed:

| Value | The instance is removed when |
|---|---|
| `'Ftr'` | it is close in feature space |
| `'Ftr&AP'` | the two performance vectors are also within `sqrt(nalgos/nfeats)*mindistance` |
| `'Ftr&Good'` | every algorithm is also good on both instances |
| `'Ftr&AP&Good'` | both of the above |

## Output Arguments

### `subsetIndex` — Redundant instances

*logical vector*

`true` for each instance marked redundant. Keep `~subsetIndex`.

### `isDissimilar` — Isolated instances

*logical vector*

`false` for each instance that was close to another instance.

### `isVISA` — Close but kept

*logical vector*

`true` for an instance that was close to a kept instance in feature space but did not meet the `opts.type` condition, so it was kept.

### `unif` — Uniformity of the kept set

*scalar*

One minus the coefficient of variation of the nearest-neighbour distances of the kept instances. Values near 1 mean an even spread.

## Version History

### v0.9.0 — Uniformity output

`unif` is returned; before, it was computed and discarded.

## References

- Alipour, H., Muñoz, M.A. & Smith-Miles, K. (2023). Enhanced instance space analysis for the maximum flow problem. *European Journal of Operational Research*, 304(2), 411–428. <https://doi.org/10.1016/j.ejor.2022.04.012>
- Smith-Miles, K. & Muñoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

`ISAsubsetData` | `PRELIM` | `SIFTED` | [Options Reference](OptionsReference.html#opts-selvars)
