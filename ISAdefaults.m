function opts = ISAdefaults(opts)
% ISAdefaults  Fill in missing opts fields with default values.
%   opts = ISAdefaults(opts) ensures every pipeline function receives a
%   complete options struct, eliminating scattered isfield chains across
%   buildIS, PILOT, TRACE, and CLOISTER. Call once at the buildIS entry
%   point after loading options.json.

% general
if ~isfield(opts, 'general'),           opts.general           = struct; end
if ~isfield(opts.general, 'seed'),      opts.general.seed      = 42;    end
if ~isfield(opts.general, 'verbose'),   opts.general.verbose   = true;  end
% Legacy: top-level opts.parallel.flag/.ncores -> opts.general.parallel/.ncores.
if isfield(opts, 'parallel')
    if isfield(opts.parallel, 'flag') && ~isfield(opts.general, 'parallel')
        opts.general.parallel = opts.parallel.flag;
    end
    if isfield(opts.parallel, 'ncores') && ~isfield(opts.general, 'ncores')
        opts.general.ncores = opts.parallel.ncores;
    end
end
if ~isfield(opts.general, 'parallel'),  opts.general.parallel  = false; end  % conservative default; spec Appendix A default is true
if ~isfield(opts.general, 'ncores'),    opts.general.ncores    = 18;    end

% perf
if ~isfield(opts, 'perf'),                  opts.perf                  = struct; end
if ~isfield(opts.perf, 'MaxPerf'),          opts.perf.MaxPerf          = false;  end
if ~isfield(opts.perf, 'AbsPerf'),          opts.perf.AbsPerf          = false;  end
if ~isfield(opts.perf, 'epsilon'),          opts.perf.epsilon          = 0.05;   end
if ~isfield(opts.perf, 'betaThreshold'),    opts.perf.betaThreshold    = 0.55;   end

% prelim (buildIS-level pre-processing config; iqrMultiplier is also passed
% into PRELIM itself, nanThreshold is read directly by buildIS beforehand)
if ~isfield(opts, 'prelim'),                 opts.prelim                 = struct; end
if ~isfield(opts.prelim, 'iqrMultiplier'),   opts.prelim.iqrMultiplier   = 5;      end
if ~isfield(opts.prelim, 'nanThreshold'),    opts.prelim.nanThreshold    = 0.20;   end

% auto
if ~isfield(opts, 'auto'),              opts.auto              = struct; end
if ~isfield(opts.auto, 'preproc'),      opts.auto.preproc      = true;  end

% bound
if ~isfield(opts, 'bound'),             opts.bound             = struct; end
if ~isfield(opts.bound, 'flag'),        opts.bound.flag        = true;  end

% norm
if ~isfield(opts, 'norm'),              opts.norm              = struct; end
if ~isfield(opts.norm, 'flag'),         opts.norm.flag         = true;  end

% selvars
if ~isfield(opts, 'selvars'),                    opts.selvars                    = struct;     end
if ~isfield(opts.selvars, 'smallscaleflag'),     opts.selvars.smallscaleflag     = false;      end
if ~isfield(opts.selvars, 'smallscale'),         opts.selvars.smallscale         = 0.30;       end
if ~isfield(opts.selvars, 'fileidxflag'),        opts.selvars.fileidxflag        = false;      end
if ~isfield(opts.selvars, 'fileidx'),            opts.selvars.fileidx            = '';         end
if ~isfield(opts.selvars, 'densityflag'),        opts.selvars.densityflag        = false;      end
if ~isfield(opts.selvars, 'mindistance'),        opts.selvars.mindistance        = 0.10;       end
if ~isfield(opts.selvars, 'type'),               opts.selvars.type               = 'Ftr&Good'; end

% sifted
if ~isfield(opts, 'sifted'),                opts.sifted                = struct; end
if ~isfield(opts.sifted, 'flag'),           opts.sifted.flag           = true;   end
if ~isfield(opts.sifted, 'rho'),            opts.sifted.rho            = 0.10;   end
if ~isfield(opts.sifted, 'pval'),           opts.sifted.pval           = 0.05;   end
if ~isfield(opts.sifted, 'K'),              opts.sifted.K              = 10;     end
if ~isfield(opts.sifted, 'MaxIter'),        opts.sifted.MaxIter        = 1000;   end
if ~isfield(opts.sifted, 'Replicates'),     opts.sifted.Replicates     = 100;    end

% pilot
if ~isfield(opts, 'pilot'),             opts.pilot             = struct; end
if ~isfield(opts.pilot, 'analytic'),    opts.pilot.analytic    = false;  end
if ~isfield(opts.pilot, 'ntries'),      opts.pilot.ntries      = 10;     end
% Legacy: opts.pilot.ISA3D (boolean) -> opts.pilot.dims (2|3, spec Appendix A).
if isfield(opts.pilot, 'ISA3D') && ~isfield(opts.pilot, 'dims')
    opts.pilot.dims = 2 + double(logical(opts.pilot.ISA3D));
end
if ~isfield(opts.pilot, 'dims'),        opts.pilot.dims        = 2;      end
if ~isfield(opts.pilot, 'method'),      opts.pilot.method      = 'standard'; end
if ~isfield(opts.pilot, 'alpha'),       opts.pilot.alpha       = 1.0;    end
if ~isfield(opts.pilot, 'viewGroups'),  opts.pilot.viewGroups  = {};     end
if ~isfield(opts.pilot, 'topoWeight'),  opts.pilot.topoWeight  = 0;      end
if ~isfield(opts.pilot, 'verbose'),     opts.pilot.verbose     = opts.general.verbose; end

% cloister
if ~isfield(opts, 'cloister'),              opts.cloister              = struct; end
if ~isfield(opts.cloister, 'pval'),         opts.cloister.pval         = 0.05;   end
% Legacy: opts.cloister.cthres -> opts.cloister.corrThreshold.
if isfield(opts.cloister, 'cthres') && ~isfield(opts.cloister, 'corrThreshold')
    opts.cloister.corrThreshold = opts.cloister.cthres;
end
if ~isfield(opts.cloister, 'corrThreshold'),opts.cloister.corrThreshold= 0.70;   end
if ~isfield(opts.cloister, 'maxFeatures'),  opts.cloister.maxFeatures  = 20;     end

% pythia
if ~isfield(opts, 'pythia'),                   opts.pythia                   = struct;                   end
if ~isfield(opts.pythia, 'flag'),              opts.pythia.flag              = true;                    end
% Legacy: opts.pythia.cvfolds -> opts.pythia.kFold.
if isfield(opts.pythia, 'cvfolds') && ~isfield(opts.pythia, 'kFold')
    opts.pythia.kFold = opts.pythia.cvfolds;
end
if ~isfield(opts.pythia, 'kFold'),             opts.pythia.kFold             = 5;                      end
if ~isfield(opts.pythia, 'tuning'),            opts.pythia.tuning            = 'sobol';                 end
if ~isfield(opts.pythia, 'nTuningIter'),       opts.pythia.nTuningIter       = 20;                     end
if ~isfield(opts.pythia, 'params'),            opts.pythia.params            = [];                     end
if ~isfield(opts.pythia, 'skip'),              opts.pythia.skip              = false;                   end
if ~isfield(opts.pythia, 'ispolykrnl'),        opts.pythia.ispolykrnl        = false;                  end
if ~isfield(opts.pythia, 'useweights'),        opts.pythia.useweights        = false;                  end
if ~isfield(opts.pythia, 'ensembleMethod'),    opts.pythia.ensembleMethod    = 'Bag';                  end
if ~isfield(opts.pythia, 'verbose'),           opts.pythia.verbose           = opts.general.verbose;   end
if ~isfield(opts.pythia, 'seed'),              opts.pythia.seed              = opts.general.seed;      end
% classifier: registry name ('knn' default). Honour legacy useknn flag.
if ~isfield(opts.pythia, 'classifier')
    if isfield(opts.pythia, 'useknn') && ~opts.pythia.useknn
        opts.pythia.classifier = 'svm';
    else
        opts.pythia.classifier = 'knn';
    end
end

% trace
% Note: no nn/prior fields. PYTHIA always runs before TRACE (mandatory coupling,
% spec §4.5), so TRACE never trains its own KNN classifier; a per-algorithm
% opts.trace.nn/prior configuration for that purpose would never be read.
if ~isfield(opts, 'trace'),                     opts.trace                     = struct;    end
if ~isfield(opts.trace, 'method'),              opts.trace.method              = 'trace3';  end
if ~isfield(opts.trace, 'PI'),                  opts.trace.PI                  = 0.6;       end
if ~isfield(opts.trace, 'minInstances'),        opts.trace.minInstances        = 4;         end
if ~isfield(opts.trace, 'minAreaFrac'),         opts.trace.minAreaFrac         = 0.01;      end
if ~isfield(opts.trace, 'contra')
    % Contradiction removal defaults to true only for the legacy method (spec §4.1).
    opts.trace.contra = strcmpi(opts.trace.method, 'legacy');
end

% outputs
if ~isfield(opts, 'outputs'),           opts.outputs           = struct; end
if ~isfield(opts.outputs, 'csv'),       opts.outputs.csv       = true;   end
if ~isfield(opts.outputs, 'png'),       opts.outputs.png       = true;   end
if ~isfield(opts.outputs, 'web'),       opts.outputs.web       = false;  end

end
