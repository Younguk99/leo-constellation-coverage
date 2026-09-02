%% RUN_SIMULATION  — computes every number reported in Chapter 3
%
%  Produces results.mat, which make_figures.m then reads.
%  Runtime is dominated by the greedy adversary search: expect several minutes.
%
%  METRIC CONVENTION
%    The constellation must cover nine BMOA target areas. A mission fails if
%    the adversary exploits a gap over ANY of them, so the reported gap is the
%    WORST across the nine targets, not an average. The phasing parameter F is
%    likewise chosen to minimise that worst case — i.e. each configuration is
%    evaluated at its own best F, which is the fairest comparison.
%
%  ADVERSARY MODEL
%    Greedy optimal: at every step remove the satellite that maximises the
%    resulting worst-case gap. Do NOT substitute random removal — an opponent
%    holding orbital elements does not choose at random, and random removal
%    understates the degradation.

clear; clc; close all;
cfg = kmw_config();
T = cfg.T;  PL = cfg.P_list;  NMAX = cfg.NMAX;
ntg = size(cfg.targets,1);

fprintf('=== baseline: %d sats, %g km, i=%g deg, incidence %g-%g deg, %d targets ===\n\n', ...
        T, cfg.alt_km, cfg.inc_deg, cfg.inc_min, cfg.inc_max, ntg);

%% ---------- 1. peacetime: best F per P ----------
fprintf('--- peacetime (F chosen to minimise worst-case MRT) ---\n');
fprintf('%4s%5s%5s | %9s%10s | per-target MRT (min)\n','P','s','F*','ART_wc','MRT_wc');
peace = repmat(struct('P',0,'s',0,'F',0,'arts',[],'mrts',[],'A',[]), 1, numel(PL));
for ii = 1:numel(PL)
    P = PL(ii);  bestWC = inf;
    for F = 0:P-1
        A = kmw_accessmask(P, F, cfg);
        arts = zeros(1,ntg);  mrts = zeros(1,ntg);
        for m = 1:ntg
            [~, arts(m), mrts(m)] = kmw_gaps(A(:,:,m), cfg.dt);
        end
        if max(mrts) < bestWC
            bestWC = max(mrts);  bF = F;  bA = A;  bArt = arts;  bMrt = mrts;
        end
    end
    peace(ii).P = P;      peace(ii).s = T/P;   peace(ii).F = bF;
    peace(ii).arts = bArt; peace(ii).mrts = bMrt; peace(ii).A = bA;
    fprintf('%4d%5d%5d | %9.2f%10.1f | %s\n', P, T/P, bF, max(bArt), bestWC, ...
            sprintf('%6.1f', bMrt));
end

fprintf('\n--- Table 4 ---\n%6s%6s%10s%10s%12s\n','P','s','ART','MRT','MRT/ART');
for ii = 1:numel(PL)
    fprintf('%6d%6d%10.2f%10.1f%12.1f\n', peace(ii).P, peace(ii).s, ...
            max(peace(ii).arts), max(peace(ii).mrts), ...
            max(peace(ii).mrts)/max(peace(ii).arts));
end

%% ---------- 2. degradation under greedy attack ----------
fprintf('\n--- degradation: G_max(n), worst BMOA target ---\n');
G = zeros(numel(PL), NMAX+1, 2);          % (config, n, scenario)
scen = {'P','D'};
for ii = 1:numel(PL)
    P = peace(ii).P;  S = peace(ii).s;  A = peace(ii).A;
    plane_of = repelem(1:P, S);
    for sc = 1:2
        alive = true(1,T);
        C = zeros(ntg, size(A,2), 'int16');
        for m = 1:ntg, C(m,:) = sum(A(:,:,m), 1); end
        for step = 0:NMAX
            if step > 0
                if sc == 1                       % plane-concentrated
                    cand = [];
                    for pl = 1:P
                        cand = find(alive & (plane_of == pl));
                        if ~isempty(cand), break; end
                    end
                else                             % distributed
                    cand = find(alive);
                end
                bm = -inf;  bk = cand(1);
                for k = cand
                    w = 0;
                    for m = 1:ntg
                        w = max(w, kmw_maxgap((C(m,:) - int16(A(k,:,m))) > 0, cfg.dt));
                    end
                    if w > bm, bm = w; bk = k; end
                end
                alive(bk) = false;
                for m = 1:ntg, C(m,:) = C(m,:) - int16(A(bk,:,m)); end
            end
            w = 0;
            for m = 1:ntg, w = max(w, kmw_maxgap(C(m,:) > 0, cfg.dt)); end
            G(ii, step+1, sc) = w;
        end
        fprintf('P=%2d s=%2d scen %s: %s\n', P, S, scen{sc}, sprintf('%8.1f', G(ii,:,sc)));
    end
end

%% ---------- 3. N_crit as a function of W ----------
fprintf('\n--- N_crit(W), plane-concentrated attack ---\n');
WL = cfg.W_list;
Ncrit = zeros(numel(PL), numel(WL));
fprintf('%5s%5s |', 'P','s'); fprintf('%8s', compose("W=%d", WL)); fprintf('\n');
for ii = 1:numel(PL)
    for jj = 1:numel(WL)
        idx = find(G(ii,:,1) > WL(jj), 1);
        if isempty(idx)
            Ncrit(ii,jj) = NMAX + 1;      % deadline never exceeded within n <= NMAX
        else
            Ncrit(ii,jj) = idx - 1;       % index 1 corresponds to n = 0
        end
    end
    fprintf('%5d%5d |', peace(ii).P, peace(ii).s);
    fprintf('%8d', Ncrit(ii,:)); fprintf('\n');
end
fprintf('  (a value of %d means the deadline was never exceeded within n <= %d)\n', NMAX+1, NMAX);

%% ---------- 4. daily repeatability ----------
fprintf('\n--- daily maximum gap (worst target of each configuration) ---\n');
daily = nan(numel(PL), cfg.days);
for ii = 1:numel(PL)
    A = peace(ii).A;
    [~, mi] = max(peace(ii).mrts);
    u = any(A(:,:,mi), 1);
    d = diff(int8(u));  st = find(d==-1);  en = find(d==1);
    en = en(en > st(1));  st = st(1:numel(en));
    gm = (en-st)*cfg.dt/60;  th = st*cfg.dt/3600;
    for dd = 1:cfg.days
        sel = th >= (dd-1)*24 & th < dd*24;
        if any(sel), daily(ii,dd) = max(gm(sel)); end
    end
    sp = 100*(max(daily(ii,:))-min(daily(ii,:)))/min(daily(ii,:));
    fprintf('P=%2d (target #%d): %s  spread %.1f%%\n', peace(ii).P, mi, ...
            sprintf('%7.1f', daily(ii,:)), sp);
end

%% ---------- 5. single-satellite clustering ----------
c1 = cfg;  c1.T = 1;  c1.days = 6;
A1 = kmw_accessmask(1, 0, c1);
u1 = squeeze(A1(1,:,1));
t1 = (0:c1.dt:(c1.days*86400-c1.dt))/3600;
d = diff(int8(u1));  acc1 = t1(find(d==1)+1);
iv1 = diff(acc1);
fprintf('\n--- single satellite ---\n');
fprintf('accesses: %d in %d days (%.2f/day)\n', numel(acc1), c1.days, numel(acc1)/c1.days);
fprintf('short intervals ~%.2f h, long intervals ~%.2f h\n', ...
        median(iv1(iv1<3)), median(iv1(iv1>3)));

%% ---------- 6. geometry constants ----------
[~, info] = kmw_accessmask(PL(1), 0, cfg);
fprintf('\n--- geometry ---\n');
fprintf('Keplerian period %.2f min, nodal period %.2f min\n', ...
        info.period_min, info.nodal_period_min);
fprintf('nodal day %.3f h -> %.4f revs/day  (15.0 confirms the gamma = 15/1 repeat)\n', ...
        info.nodal_day_h, info.revs_per_day);
fprintf('ground-track spacing %.3f deg  (360/15 = 24.000)\n', info.track_walk_deg);
fprintf('annulus %.2f-%.2f deg = %.0f-%.0f km\n', info.lam_in_deg, info.lam_out_deg, ...
        info.lam_in_deg*111.2, info.lam_out_deg*111.2);
fprintf('access window +-%.2f deg longitude; spacing/window = %.2f\n', ...
        info.halfwin_deg, info.track_walk_deg/(2*info.halfwin_deg));

save('results.mat','cfg','peace','G','Ncrit','daily','info','u1','t1','iv1','acc1');
fprintf('\nsaved results.mat — now run make_figures\n');
