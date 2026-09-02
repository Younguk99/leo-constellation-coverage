function [gaps, ART, MRT, Pfail] = kmw_gaps(A, dt, W_s)
%KMW_GAPS  Observation-gap statistics and mission-window failure probability.
%
%   [gaps, ART, MRT, Pfail] = KMW_GAPS(A, dt, W_s)
%
%   INPUTS
%     A    : logical (nsat x nt) access mask, or (1 x nt) already-collapsed.
%     dt   : time step (s)
%     W_s  : mission window length W (s). Optional; if omitted Pfail = NaN.
%
%   OUTPUTS
%     gaps  : vector of observation gap durations (s)
%     ART   : mean gap  (min)   -- "average revisit time"
%     MRT   : max  gap  (min)   -- "maximum revisit time" = G_max
%     Pfail : P(no valid access anywhere inside a random window of length W)
%
%   Pfail is computed exactly, not by Monte Carlo. For a gap of length g, the
%   set of window start times that see no access has measure max(0, g - W).
%   Hence  Pfail = sum_i max(0, g_i - W) / T_total.
%
%   IMPORTANT: gaps are measured strictly BETWEEN two access intervals, so the
%   leading and trailing partial gaps are discarded. Run long enough that this
%   is immaterial (>= ~15 days for a 30-satellite LEO constellation).

    if nargin < 3, W_s = NaN; end

    anyacc = any(A, 1);
    d = diff(int8(anyacc));
    starts = find(d == -1) + 1;     % access -> no access
    ends   = find(d ==  1) + 1;     % no access -> access

    if isempty(starts) || isempty(ends)
        gaps = []; ART = NaN; MRT = NaN; Pfail = NaN;
        if all(anyacc), ART = 0; MRT = 0; Pfail = 0; end
        return
    end

    ends   = ends(ends > starts(1));
    starts = starts(1:numel(ends));
    gaps   = (ends - starts) * dt;

    ART = mean(gaps) / 60;
    MRT = max(gaps)  / 60;

    if isnan(W_s)
        Pfail = NaN;
    else
        total_s = numel(anyacc) * dt;
        Pfail = sum(max(0, gaps - W_s)) / total_s;
    end
end
