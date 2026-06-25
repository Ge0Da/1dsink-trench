function [domains, boundaries, diagnostics] = segmentDomains(s_km, curvature, edgeMask, opts)
%SEGMENTDOMAINS  Spatially stationary domain segmentation from curvature
%
%   [domains, boundaries, diagnostics] = segmentDomains(s_km, curvature, edgeMask)
%   [domains, boundaries, diagnostics] = segmentDomains(s_km, curvature, edgeMask, opts)
%
%   STATIONARITY GUARANTEE
%   ======================
%   Every decision at a point depends ONLY on data within a declared
%   physical neighborhood of that point.  Interior results are identical
%   regardless of total profile extent or trim position.
%
%   Three design principles enforce this:
%     1. DEADBAND computed LOCALLY (running-window MAD, fixed km width)
%     2. CONFIDENCE normalized against LOCAL curvature amplitude and a
%        FIXED physical dwell scale (not global percentiles or profile length)
%     3. BOUNDARY SELECTION has NO EDGE PENALTIES (profile edges are
%        observation limits, not geological features)
%
%   ALGORITHM
%   =========
%   1. Pre-smooth curvature (optional Gaussian, fixed physical width)
%   2. Classify each point as +1 / -1 / 0 using LOCAL deadband
%   3. Find persistent runs (run-length filter, fixed physical length)
%   4. Identify transitions between opposite-sign persistent runs
%   5. Score transitions (local amplitude + fixed dwell scale)
%   6. Select boundaries (greedy, no edge penalties)
%   7. Build domain structure
%   8. Optional stability audit via parameter sweep
%
%   INPUTS
%       s_km        - Along-trench distance (km), monotonic increasing
%       curvature   - Along-trench curvature (m/km^2)
%       edgeMask    - Logical mask, true = valid interior point
%
%   OPTIONS (name-value)
%       deadbandMode              'fixed'|'percentile'|'mad'   default 'mad'
%       deadbandAmplitude         fixed value (m/km^2)         default 0.02
%       deadbandFactor            percentile fraction           default 0.10
%       deadbandSigmas            MAD multiplier                default 1.5
%       deadbandWindow_km         local window half-width (km)  default 100
%       persistenceLength_km      min run length (km)           default 10
%       minDomainLength_km        min boundary spacing (km)     default 20
%       confidenceWeights         [sharp, dwell, amp]           default [0.4 0.3 0.3]
%       dwellScale_km             fixed dwell norm (km)         default 100
%       localAmplitudeRadius_km   local norm radius (km)        default 50
%       maxBoundaries             safety cap                    default 30
%       smoothWindow_km           pre-smooth width (km)         default 5
%       stabilityTest             run sweep audit               default true
%       stabilityPerturbation     +/- fraction                  default 0.30
%       stabilityGrid             levels per param              default 5
%       stabilitySnapRadius_km    match radius (km)             default 15
%       minStability              survival filter               default 0
%       verbose                   print diagnostics             default false
%
%
%   Part of 1DSink Along-Trench Analysis Pipeline
%
%   Author: David Fernández-Blanco
%   Date: 2025
%   Version: 3.4

arguments
    s_km        (:,1) double
    curvature   (:,1) double
    edgeMask    (:,1) logical
    opts.deadbandMode              (1,:) char {mustBeMember(opts.deadbandMode,{'fixed','percentile','mad'})} = 'mad'
    opts.deadbandAmplitude         (1,1) double {mustBeNonnegative} = 0.02
    opts.deadbandFactor            (1,1) double {mustBeInRange(opts.deadbandFactor,0,1)} = 0.10
    opts.deadbandSigmas            (1,1) double {mustBePositive} = 1.5
    opts.deadbandWindow_km         (1,1) double {mustBePositive} = 100
    opts.persistenceLength_km      (1,1) double {mustBePositive} = 10
    opts.minDomainLength_km        (1,1) double {mustBePositive} = 20
    opts.confidenceWeights         (1,3) double {mustBeNonnegative} = [0.4, 0.3, 0.3]
    opts.dwellScale_km             (1,1) double {mustBePositive} = 100
    opts.localAmplitudeRadius_km   (1,1) double {mustBePositive} = 50
    opts.maxBoundaries             (1,1) double {mustBePositive} = 30
    opts.smoothWindow_km           (1,1) double {mustBeNonnegative} = 5
    opts.stabilityTest             (1,1) logical = true
    opts.stabilityPerturbation     (1,1) double = 0.30
    opts.stabilityGrid             (1,1) double = 5
    opts.stabilitySnapRadius_km    (1,1) double = 15
    opts.minStability              (1,1) double = 0
    opts.verbose                   (1,1) logical = false
end

%% ---- initialise outputs ------------------------------------------------
domains     = struct('id',{},'s_start',{},'s_end',{},'length_km',{},...
                     'meanCurv',{},'sign',{},'type',{},'meanRoughness',{});
boundaries  = [];
diagnostics = struct();

%% ---- validate ----------------------------------------------------------
n = numel(s_km);
assert(numel(curvature)==n && numel(edgeMask)==n, ...
       'segmentDomains:size','All inputs must have same length');

if sum(edgeMask) < 10
    warning('segmentDomains:data','Fewer than 10 valid points');
    diagnostics.deadband_median = NaN;
    diagnostics.nRawCrossings   = 0;
    diagnostics.nPersistentRuns = 0;
    return
end

%% ---- Step 0  pre-smooth (local Gaussian, fixed km width) ---------------
ds_km = median(diff(s_km),'omitnan');
c     = curvature(:);

if opts.smoothWindow_km > 0
    sw = max(3, round(opts.smoothWindow_km / ds_km));
    if mod(sw,2)==0, sw = sw+1; end
    c = smoothdata(c,'gaussian',sw);
end
c(~edgeMask) = NaN;

%% ---- Step 1  LOCAL adaptive deadband -----------------------------------
%  Each point's noise threshold depends only on data within
%  +/- deadbandWindow_km/2.  Trimming distant trench has zero effect.

halfW = round(opts.deadbandWindow_km / (2*ds_km));
wSz   = 2*halfW + 1;

switch opts.deadbandMode
    case 'fixed'
        db = opts.deadbandAmplitude * ones(n,1);

    case 'percentile'
        db = NaN(n,1);
        for i = 1:n
            lo_i = max(1,i-halfW);  hi_i = min(n,i+halfW);
            lc   = abs(c(lo_i:hi_i));
            lc   = lc(~isnan(lc));
            if numel(lc)>=5, db(i) = opts.deadbandFactor*prctile(lc,80); end
        end

    case 'mad'
        lMed = movmedian(c, wSz,'omitnan','Endpoints','shrink');
        lDev = abs(c - lMed);
        lMAD = 1.4826 * movmedian(lDev, wSz,'omitnan','Endpoints','shrink');
        db   = opts.deadbandSigmas * lMAD;
end

db = max(db, 1e-6);
db(isnan(db)) = Inf;

diagnostics.deadband_median    = median(db(isfinite(db)));
diagnostics.deadband_local     = db;
diagnostics.deadbandMode       = opts.deadbandMode;
diagnostics.deadbandWindow_km  = opts.deadbandWindow_km;

if opts.verbose
    fdb = db(isfinite(db));
    fprintf('segmentDomains: local deadband (%s, win=%.0f km)  med=%.4f  [%.4f, %.4f]\n',...
            opts.deadbandMode, opts.deadbandWindow_km, median(fdb), min(fdb), max(fdb));
end

%% ---- Step 2  classify with LOCAL deadband ------------------------------
state = zeros(n,1);
state(c >  db) =  1;
state(c < -db) = -1;

rs = sign(c);
diagnostics.nRawCrossings = sum(diff(rs(~isnan(rs)))~=0);

%% ---- Step 3  persistent runs ------------------------------------------
runs = findPersistentRuns(s_km, state, opts.persistenceLength_km);
diagnostics.nPersistentRuns = numel(runs);

if opts.verbose
    fprintf('segmentDomains: %d persistent runs (>= %.1f km)\n',...
            numel(runs), opts.persistenceLength_km);
end

if numel(runs) < 2
    domains = buildSingleDomain(s_km,c,edgeMask);
    diagnostics.allTransitions    = [];
    diagnostics.confidenceScores  = [];
    diagnostics.boundaryStability   = [];
    diagnostics.stabilityNcombos    = 0;
    diagnostics.featureStability    = [];
    diagnostics.featureBoundaryStab = [];
    diagnostics.boundaryRange_km    = [];
    diagnostics.boundaryStd_km      = [];
    return
end

%% ---- Step 4  transitions -----------------------------------------------
trans = identifyTransitions(runs, s_km, c);
diagnostics.allTransitions = trans;

if isempty(trans)
    domains = buildSingleDomain(s_km,c,edgeMask);
    diagnostics.confidenceScores    = [];
    diagnostics.boundaryStability   = [];
    diagnostics.stabilityNcombos    = 0;
    diagnostics.featureStability    = [];
    diagnostics.featureBoundaryStab = [];
    diagnostics.boundaryRange_km    = [];
    diagnostics.boundaryStd_km      = [];
    return
end

%% ---- Step 5  score with LOCAL normalisation ----------------------------
w = opts.confidenceWeights / sum(opts.confidenceWeights);
for i = 1:numel(trans)
    trans(i).confidence = scoreTransition(trans(i), c, s_km, w, ...
        opts.localAmplitudeRadius_km, opts.dwellScale_km);
end

%% ---- Step 6  select boundaries (NO edge penalties) ---------------------
[boundaries, selIdx] = selectBoundaries(trans, ...
    opts.minDomainLength_km, opts.maxBoundaries);
diagnostics.confidenceScores = [trans(selIdx).confidence];

if opts.verbose
    fprintf('segmentDomains: %d boundaries\n', numel(boundaries));
end

%% ---- Step 6b  stability audit ------------------------------------------
if opts.stabilityTest && ~isempty(boundaries)
    nG  = opts.stabilityGrid;
    pR  = opts.stabilityPerturbation;
    snR = opts.stabilitySnapRadius_km;
    mts = linspace(1-pR, 1+pR, nG);
    nC    = nG^3;
    sB    = cell(nC,1);
    sDoms = cell(nC,1);
    idx   = 0;

    for iDB = 1:nG
        db_s = db * mts(iDB);
        st_i = zeros(n,1);
        st_i(c >  db_s) =  1;
        st_i(c < -db_s) = -1;

        for iPL = 1:nG
            r_i = findPersistentRuns(s_km, st_i, opts.persistenceLength_km*mts(iPL));
            if numel(r_i)<2
                for sk = 1:nG
                    idx = idx+1;
                    sB{idx}    = [];
                    sDoms{idx} = buildDomains(s_km, c, edgeMask, []);
                end
                continue
            end
            t_i = identifyTransitions(r_i, s_km, c);
            if ~isempty(t_i)
                for t=1:numel(t_i)
                    t_i(t).confidence = scoreTransition(t_i(t),c,s_km,w,...
                        opts.localAmplitudeRadius_km, opts.dwellScale_km);
                end
            end
            for iMD = 1:nG
                idx = idx+1;
                if ~isempty(t_i)
                    [sB{idx},~] = selectBoundaries(t_i, ...
                        opts.minDomainLength_km*mts(iMD), opts.maxBoundaries);
                else
                    sB{idx} = [];
                end
                sDoms{idx} = buildDomains(s_km, c, edgeMask, sB{idx});
            end
        end
    end

    nB = numel(boundaries);
    surv = zeros(nB,1);
    for b = 1:nB
        ct = 0;
        for cc = 1:nC
            if ~isempty(sB{cc}) && any(abs(sB{cc}-boundaries(b))<=snR)
                ct = ct+1;
            end
        end
        surv(b) = ct/nC;
    end
    diagnostics.boundaryStability     = surv;
    diagnostics.stabilityNcombos      = nC;
    diagnostics.stabilityPerturbation = pR;

    % ---- feature-level stability ----------------------------------------
    % Build reference domains from the main-run boundaries.
    refDoms  = buildDomains(s_km, c, edgeMask, boundaries);
    nDom     = numel(refDoms);

    % domSurv(d,cc) = true when domain d survives perturbation combo cc.
    % A domain survives when >= 2 of 3 probe points land in a perturbed
    % domain of the same type (trap / bypass).
    domSurv = false(nDom, nC);
    for cc = 1:nC
        pd = sDoms{cc};
        if isempty(pd), continue; end
        for d = 1:nDom
            rd       = refDoms(d);
            probe_km = rd.s_start + [0.25, 0.50, 0.75]*(rd.s_end - rd.s_start);
            matchCt  = 0;
            for p = 1:3
                for pd_idx = 1:numel(pd)
                    if probe_km(p) >= pd(pd_idx).s_start && ...
                            probe_km(p) <= pd(pd_idx).s_end
                        if strcmp(pd(pd_idx).type, rd.type)
                            matchCt = matchCt + 1;
                        end
                        break
                    end
                end
            end
            if matchCt >= 2, domSurv(d, cc) = true; end
        end
    end
    featureStability = sum(domSurv, 2) / nC;

    % Feature boundary stability: boundary b (between domain b and b+1)
    % survives when BOTH adjacent domains survive.
    % Also collect matched boundary positions for range / std statistics.
    nB_fbs      = nDom - 1;   % == nB
    featBndStab = zeros(nB_fbs, 1);
    bndRange    = NaN(nB_fbs, 2);
    bndStd      = NaN(nB_fbs, 1);
    for b = 1:nB_fbs
        bothSurv = domSurv(b,:) & domSurv(b+1,:);
        featBndStab(b) = sum(bothSurv) / nC;
        pos = [];
        for cc = 1:nC
            if bothSurv(cc) && ~isempty(sB{cc})
                [minD, minI] = min(abs(sB{cc} - boundaries(b)));
                if minD <= snR
                    pos(end+1) = sB{cc}(minI); %#ok<AGROW>
                end
            end
        end
        if numel(pos) >= 2
            bndRange(b,:) = [min(pos), max(pos)];
            bndStd(b)     = std(pos);
        elseif numel(pos) == 1
            bndRange(b,:) = [pos, pos];
            bndStd(b)     = 0;
        end
    end

    diagnostics.featureStability    = featureStability;
    diagnostics.featureBoundaryStab = featBndStab;
    diagnostics.boundaryRange_km    = bndRange;
    diagnostics.boundaryStd_km      = bndStd;

    if opts.verbose
        for b=1:nB
            lab='robust';
            if surv(b)<0.5, lab='FRAGILE'; end
            if surv(b)<0.2, lab='VERY FRAGILE'; end
            fprintf('  %.1f km  survival=%.0f%%  [%s]\n',boundaries(b),surv(b)*100,lab);
        end
        for d = 1:nDom
            fprintf('  D%d: %.0f--%.0f km  %s  featureStab=%.0f%%\n', ...
                    d, refDoms(d).s_start, refDoms(d).s_end, ...
                    refDoms(d).type, featureStability(d)*100);
        end
    end

    if opts.minStability > 0
        keep = surv >= opts.minStability;
        if any(~keep)
            boundaries = boundaries(keep);
            diagnostics.confidenceScores    = diagnostics.confidenceScores(keep);
            diagnostics.boundaryStability   = surv(keep);
            diagnostics.featureBoundaryStab = featBndStab(keep);
            diagnostics.boundaryRange_km    = bndRange(keep,:);
            diagnostics.boundaryStd_km      = bndStd(keep);
        end
    end
else
    diagnostics.boundaryStability   = [];
    diagnostics.stabilityNcombos    = 0;
    diagnostics.featureStability    = [];
    diagnostics.featureBoundaryStab = [];
    diagnostics.boundaryRange_km    = [];
    diagnostics.boundaryStd_km      = [];
end

%% ---- Step 7  build domains ---------------------------------------------
domains = buildDomains(s_km, c, edgeMask, boundaries);

if opts.verbose
    for d=1:numel(domains)
        fprintf('  D%d: %.1f--%.1f km  %s  k=%.3f\n',...
                d, domains(d).s_start, domains(d).s_end,...
                domains(d).type, domains(d).meanCurv);
    end
end
end

%% ========================================================================
%  HELPER FUNCTIONS  (all spatially local)
%% ========================================================================

function runs = findPersistentRuns(s_km, state, minLen)
    runs = struct('sign',{},'s_start',{},'s_end',{},...
                  'idx_start',{},'idx_end',{},'length_km',{});
    def = state~=0 & ~isnan(state);
    if sum(def)<2, return; end
    dI = find(def);  dS = state(dI);  dK = s_km(dI);
    ch = [1; find(diff(dS)~=0)+1; numel(dS)+1];
    for i = 1:numel(ch)-1
        a = ch(i); b = ch(i+1)-1;
        L = dK(b)-dK(a);
        if L >= minLen
            runs(end+1).sign      = dS(a); %#ok<AGROW>
            runs(end).s_start     = dK(a);
            runs(end).s_end       = dK(b);
            runs(end).idx_start   = dI(a);
            runs(end).idx_end     = dI(b);
            runs(end).length_km   = L;
        end
    end
    if ~isempty(runs), [~,o]=sort([runs.s_start]); runs=runs(o); end
end

function trans = identifyTransitions(runs, s_km, curv)
    trans = struct('s_km',{},'from_sign',{},'to_sign',{},...
                   'from_run',{},'to_run',{},...
                   'sharpness',{},'dwell_before',{},'dwell_after',{},...
                   'amplitude_contrast',{},'confidence',{});
    for i = 1:numel(runs)-1
        r1=runs(i); r2=runs(i+1);
        if r1.sign*r2.sign >= 0, continue; end
        ts = (r1.s_end + r2.s_start)/2;
        ti = find(s_km>=ts,1,'first');
        if isempty(ti)||ti<2||ti>=numel(s_km), sh=0;
        else, ds=s_km(ti)-s_km(ti-1); sh=abs((curv(ti)-curv(ti-1))/max(ds,1e-6));
        end
        ac = abs(mean(curv(r1.idx_start:r1.idx_end),'omitnan') - ...
                 mean(curv(r2.idx_start:r2.idx_end),'omitnan'));
        trans(end+1).s_km               = ts; %#ok<AGROW>
        trans(end).from_sign            = r1.sign;
        trans(end).to_sign              = r2.sign;
        trans(end).from_run             = i;
        trans(end).to_run               = i+1;
        trans(end).sharpness            = sh;
        trans(end).dwell_before         = r1.length_km;
        trans(end).dwell_after          = r2.length_km;
        trans(end).amplitude_contrast   = ac;
        trans(end).confidence           = 0;
    end
end

function conf = scoreTransition(tr, c, s_km, w, ampRad_km, dwellScale_km)
%   LOCAL normalisation only.
    ds = median(diff(s_km),'omitnan');
    ti = find(s_km>=tr.s_km,1,'first');
    rp = round(ampRad_km/ds);
    lo = max(1,ti-rp); hi = min(numel(c),ti+rp);
    lc = c(lo:hi); lc=lc(~isnan(lc));
    if numel(lc)>=5, lr=prctile(lc,90)-prctile(lc,10); else, lr=1; end
    lr = max(lr,1e-6);

    nS = min(1, tr.sharpness / (lr/10));
    nD = min(1, sqrt(tr.dwell_before*tr.dwell_after) / dwellScale_km);
    nA = min(1, tr.amplitude_contrast / lr);

    conf = max(0, min(1, w(1)*nS + w(2)*nD + w(3)*nA));
end

function [bnd, sel] = selectBoundaries(trans, minSp, maxB)
%   NO edge penalties.  Only inter-boundary spacing enforced.
    if isempty(trans), bnd=[]; sel=[]; return; end
    [~,ord] = sort([trans.confidence],'descend');
    bnd=[]; sel=[];
    for i = ord
        cand = trans(i).s_km;
        if isempty(bnd), md=Inf; else, md=min(abs(bnd-cand)); end
        if md >= minSp
            bnd(end+1,1)=cand; sel(end+1,1)=i; %#ok<AGROW>
        end
        if numel(bnd)>=maxB, break; end
    end
    [bnd,ro]=sort(bnd); sel=sel(ro);
end

function doms = buildDomains(s_km, curv, mask, bnd)
    ab = [s_km(1); bnd(:); s_km(end)];
    nD = numel(ab)-1;
    doms = struct('id',{},'s_start',{},'s_end',{},'length_km',{},...
                  'meanCurv',{},'sign',{},'type',{},'meanRoughness',{});
    for d = 1:nD
        in = s_km>=ab(d) & s_km<=ab(d+1) & mask;
        cd = curv(in); cd=cd(~isnan(cd));
        if ~isempty(cd), mc=mean(cd); mr=sqrt(mean(cd.^2));
        else,            mc=0; mr=0; end
        if mc<0, sg=-1; tp='trap'; else, sg=1; tp='bypass'; end
        doms(d).id=d; doms(d).s_start=ab(d); doms(d).s_end=ab(d+1);
        doms(d).length_km=ab(d+1)-ab(d); doms(d).meanCurv=mc;
        doms(d).sign=sg; doms(d).type=tp; doms(d).meanRoughness=mr;
    end
end

function doms = buildSingleDomain(s_km, curv, mask)
    in = mask & ~isnan(curv);
    if any(in), mc=mean(curv(in),'omitnan'); mr=sqrt(mean(curv(in).^2,'omitnan'));
    else,       mc=0; mr=0; end
    if mc<0, sg=-1; tp='trap'; else, sg=1; tp='bypass'; end
    doms = struct('id',1,'s_start',s_km(1),'s_end',s_km(end),...
                  'length_km',s_km(end)-s_km(1),'meanCurv',mc,...
                  'sign',sg,'type',tp,'meanRoughness',mr);
end
