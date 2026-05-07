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

% perf
if ~isfield(opts, 'perf'),                  opts.perf                  = struct; end
if ~isfield(opts.perf, 'MaxPerf'),          opts.perf.MaxPerf          = false;  end
if ~isfield(opts.perf, 'AbsPerf'),          opts.perf.AbsPerf          = false;  end
if ~isfield(opts.perf, 'epsilon'),          opts.perf.epsilon          = 0.05;   end
if ~isfield(opts.perf, 'betaThreshold'),    opts.perf.betaThreshold    = 0.55;   end

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
if ~isfield(opts.sifted, 'NTREES'),         opts.sifted.NTREES         = 50;     end
if ~isfield(opts.sifted, 'MaxIter'),        opts.sifted.MaxIter        = 1000;   end
if ~isfield(opts.sifted, 'Replicates'),     opts.sifted.Replicates     = 100;    end

% pilot
if ~isfield(opts, 'pilot'),             opts.pilot             = struct; end
if ~isfield(opts.pilot, 'analytic'),    opts.pilot.analytic    = false;  end
if ~isfield(opts.pilot, 'ntries'),      opts.pilot.ntries      = 10;     end
if ~isfield(opts.pilot, 'ISA3D'),       opts.pilot.ISA3D       = false;  end
if ~isfield(opts.pilot, 'verbose'),     opts.pilot.verbose     = opts.general.verbose; end

% cloister
if ~isfield(opts, 'cloister'),          opts.cloister          = struct; end
if ~isfield(opts.cloister, 'pval'),     opts.cloister.pval     = 0.05;   end
if ~isfield(opts.cloister, 'cthres'),   opts.cloister.cthres   = 0.70;   end

% pythia
if ~isfield(opts, 'pythia'),                opts.pythia                = struct; end
if ~isfield(opts.pythia, 'flag'),           opts.pythia.flag           = true;   end
if ~isfield(opts.pythia, 'useknn'),         opts.pythia.useknn         = false;  end
if ~isfield(opts.pythia, 'cvfolds'),        opts.pythia.cvfolds        = 5;      end
if ~isfield(opts.pythia, 'ispolykrnl'),     opts.pythia.ispolykrnl     = false;  end
if ~isfield(opts.pythia, 'useweights'),     opts.pythia.useweights     = false;  end
if ~isfield(opts.pythia, 'uselibsvm'),      opts.pythia.uselibsvm      = false;  end
if ~isfield(opts.pythia, 'verbose'),        opts.pythia.verbose        = opts.general.verbose; end

% trace
if ~isfield(opts, 'trace'),                     opts.trace                     = struct;    end
if ~isfield(opts.trace, 'method'),              opts.trace.method              = 'trace3';  end
if ~isfield(opts.trace, 'PI'),                  opts.trace.PI                  = 0.6;       end
if ~isfield(opts.trace, 'nn'),                  opts.trace.nn                  = 50;        end  % legacy only
if ~isfield(opts.trace, 'prior'),               opts.trace.prior               = [0.6,0.4]; end  % legacy only
if ~isfield(opts.trace, 'minInstances'),        opts.trace.minInstances        = 4;         end
if ~isfield(opts.trace, 'minAreaFrac'),         opts.trace.minAreaFrac         = 0.01;      end

% outputs
if ~isfield(opts, 'outputs'),           opts.outputs           = struct; end
if ~isfield(opts.outputs, 'csv'),       opts.outputs.csv       = true;   end
if ~isfield(opts.outputs, 'png'),       opts.outputs.png       = true;   end
if ~isfield(opts.outputs, 'web'),       opts.outputs.web       = false;  end

% parallel
if ~isfield(opts, 'parallel'),          opts.parallel          = struct; end
if ~isfield(opts.parallel, 'flag'),     opts.parallel.flag     = false;  end
if ~isfield(opts.parallel, 'ncores'),   opts.parallel.ncores   = 18;     end

end
