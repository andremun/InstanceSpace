# scriptweb

Write colour-scaled data for the MATILDA web platform

## Syntax

```
scriptweb(container,rootdir)
```

## Description

`scriptweb(container,rootdir)` writes the features, performance and number of good algorithms as colour indices, with the colour table, to CSV files in `rootdir`. The [MATILDA](https://matilda.unimelb.edu.au) web platform uses these files to draw its interactive plots. You do not need them for local work.

`InstanceSpace.build` and `InstanceSpace.explore` call scriptweb when both `opts.outputs.csv` and `opts.outputs.web` are `true`.

| File | Contents |
|---|---|
| `color_table.csv` | the colour map |
| `feature_raw_color.csv`, `feature_process_color.csv` | feature colour indices |
| `algorithm_raw_color.csv`, `algorithm_process_color.csv` | performance colour indices, scaled over all algorithms |
| `algorithm_raw_single_color.csv`, `algorithm_process_single_color.csv` | performance colour indices, scaled per algorithm |
| `good_algos_color.csv` | colour index of the number of good algorithms |

## Examples

### Turn on the web output

```matlab
obj.opts.outputs.web = true;
obj = obj.build();
```

## Input Arguments

### `container` — Model or evaluation result

*structure*

### `rootdir` — Output folder

*character vector*

Must exist and end with a file separator.

## See Also

`scriptcsv` | `InstanceSpace`
