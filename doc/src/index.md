# Instance Space Analysis Toolbox

Objective algorithm testing through the instance space of a problem

The Instance Space Analysis (ISA) Toolbox relates the structural features of problem instances to the performance of a portfolio of algorithms. It projects every instance into a 2D or 3D *instance space*, shows where each algorithm performs well (its *footprint*), predicts which algorithm to use for a new instance, and estimates where in the space no instances have been observed yet.

The toolbox implements the methodology in Smith-Miles & Muñoz (2023), including the 3D extension (ISA3) of Simpson et al. (2025). It runs in MATLAB R2025a or later and powers the [MATILDA](https://matilda.unimelb.edu.au) web platform.

<div class="figure">
<svg viewBox="0 0 860 150" role="img" aria-labelledby="pipe-title" xmlns="http://www.w3.org/2000/svg">
<title id="pipe-title">The ISA pipeline: INIT and PRELIM prepare the metadata, SIFTED selects features, PILOT projects them, then CLOISTER, PYTHIA and TRACE work in the projected space.</title>
<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>
<g font-family="sans-serif" font-size="13" text-anchor="middle" fill="currentColor" stroke="currentColor">
<g stroke-width="1.2" fill="none">
<rect x="10" y="50" width="110" height="50" rx="6"/><rect x="150" y="50" width="110" height="50" rx="6"/>
<rect x="290" y="50" width="110" height="50" rx="6"/><rect x="430" y="50" width="110" height="50" rx="6"/>
<rect x="590" y="5" width="120" height="40" rx="6"/><rect x="590" y="55" width="120" height="40" rx="6"/>
<rect x="590" y="105" width="120" height="40" rx="6"/><rect x="740" y="55" width="110" height="40" rx="6"/>
<path d="M120,75 H148" marker-end="url(#arr)"/><path d="M260,75 H288" marker-end="url(#arr)"/>
<path d="M400,75 H428" marker-end="url(#arr)"/><path d="M540,75 L588,25" marker-end="url(#arr)"/>
<path d="M540,75 H588" marker-end="url(#arr)"/><path d="M540,75 L588,125" marker-end="url(#arr)"/>
<path d="M710,75 H738" marker-end="url(#arr)"/><path d="M650,95 V103" />
</g>
<g stroke="none">
<text x="65" y="72" font-weight="bold">INIT</text><text x="65" y="89" font-size="11">read metadata</text>
<text x="205" y="72" font-weight="bold">PRELIM</text><text x="205" y="89" font-size="11">label &amp; normalise</text>
<text x="345" y="72" font-weight="bold">SIFTED</text><text x="345" y="89" font-size="11">select features</text>
<text x="485" y="72" font-weight="bold">PILOT</text><text x="485" y="89" font-size="11">project to 2D/3D</text>
<text x="650" y="23" font-weight="bold">CLOISTER</text><text x="650" y="38" font-size="11">empirical bounds</text>
<text x="650" y="73" font-weight="bold">PYTHIA</text><text x="650" y="88" font-size="11">predict algorithms</text>
<text x="650" y="123" font-weight="bold">TRACE</text><text x="650" y="138" font-size="11">footprints</text>
<text x="795" y="73" font-weight="bold">Outputs</text><text x="795" y="88" font-size="11">CSV · PNG · .mat</text>
</g></g></svg>
<p class="caption">The ISA pipeline. <code>InstanceSpace</code> runs the stages in this order; each stage is also a standalone function.</p>
</div>

## Get Started

<div class="tiles">
<div class="tile"><h3><a href="GettingStarted.html">Getting Started</a></h3><p>Build and explore an instance space from the bundled reference data in a few lines.</p></div>
<div class="tile"><h3><a href="InteractiveWalkthrough.html">Stage-by-Stage Walkthrough</a></h3><p>Run each stage separately and inspect what it produces.</p></div>
<div class="tile"><h3><a href="MetadataFormat.html">Metadata File Format</a></h3><p>Prepare <code>metadata.csv</code> for your own problem domain.</p></div>
<div class="tile"><h3><a href="OptionsReference.html">Options Reference</a></h3><p>Every field of <code>opts</code> and <code>options.json</code>, with defaults.</p></div>
</div>

## Functions and Classes

<div class="tiles">
<div class="tile"><h3><a href="InstanceSpace.html">InstanceSpace</a></h3><p>Build, explore, save, load and plot an instance space. The main entry point.</p></div>
<div class="tile"><h3><a href="FunctionList.html#pipeline-stages">Pipeline Stages</a></h3><p>INIT, PRELIM, FILTER, SIFTED, PILOT, PILOTviewpoint, CLOISTER, PYTHIA, TRACE.</p></div>
<div class="tile"><h3><a href="FunctionList.html#utilities">Utilities</a></h3><p>Option defaults and validation, classifier registry, legacy-model migration.</p></div>
<div class="tile"><h3><a href="FunctionList.html#output">Output</a></h3><p>Write CSV files and figures, and restore 3D viewpoints.</p></div>
</div>

## Topics

- [Migrating a Legacy Model](MigratingLegacyModel.html) — use a `model.mat` from a version before v0.9.0.
- [Backward-Compatible Wrappers](buildIS.html) — `buildIS` and `exploreIS`, used by MATILDA.
- [Deprecated Functions](Deprecated.html) — `PYTHIA2`, `PYTHIAtest`, `SIFTED2`.
- [What's New](WhatsNew.html) — release notes.

## Installation

Clone or download the repository, then add it to the MATLAB path from its root folder:

```matlab
cd InstanceSpace
addpath(pwd)     % the class, the wrappers, and info.xml for the Help browser
startup          % adds core/, utils/, output/ and deprecated/
```

Required products: MATLAB R2025a or later with the Statistics and Machine Learning, Optimization, Global Optimization, Parallel Computing, and Financial toolboxes.

With the repository root on the path, this documentation also appears in the MATLAB Help browser under **Supplemental Software**.

## References

- Smith-Miles, K. & Muñoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>
- Simpson, C., Muñoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B. (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis. *Machine Learning*, 114, 240. <https://doi.org/10.1007/s10994-025-06871-5>
