# ISAvalidateOpts

Check the type and range of user-supplied options

## Syntax

```
opts = ISAvalidateOpts(opts)
```

## Description

`opts = ISAvalidateOpts(opts)` checks each recognised field that is present in `opts` and raises an error on the first invalid one. Absent fields are not checked; `ISAdefaults` fills them afterwards. `opts` is returned unchanged.

An error here names the field and the expected value, instead of failing later inside a stage. A field set to `[]` counts as present and is rejected. Field names the function does not recognise, including misspellings, are not reported.

`InstanceSpace` calls ISAvalidateOpts when an object is created.

## Examples

### Catch an invalid value early

```matlab
opts.pilot.dims = 4;
try
    ISAvalidateOpts(opts);
catch err
    disp(err.identifier)   % ISA:ISAvalidateOpts:notMember
    disp(err.message)
end
```

## Input Arguments

### `opts` — Options

*structure*

## Output Arguments

### `opts` — Options

*structure*

The input, unchanged.

## Tips

Error identifiers, by the check that failed:

| Identifier | Expected |
|---|---|
| `ISA:ISAvalidateOpts:notStruct` | a structure, for `opts` and each group such as `opts.pilot` |
| `ISA:ISAvalidateOpts:notLogical` | a logical scalar, or 0 or 1 |
| `ISA:ISAvalidateOpts:notFiniteNumericScalar` | a finite real number |
| `ISA:ISAvalidateOpts:notInteger` | a whole number |
| `ISA:ISAvalidateOpts:notPositive` | a positive number (0 is allowed for seeds) |
| `ISA:ISAvalidateOpts:notInUnitRange` | a number in [0, 1] |
| `ISA:ISAvalidateOpts:notMember` | one of a fixed set of values |
| `ISA:ISAvalidateOpts:notText` | a character vector or string |
| `ISA:ISAvalidateOpts:notCellOfText` | a cell array of character vectors |
| `ISA:ISAvalidateOpts:badViewGroups` | a cell array of algorithm index vectors |

## Version History

### v0.9.0 — Introduced

## See Also

`ISAdefaults` | `InstanceSpace` | [Options Reference](OptionsReference.html)
