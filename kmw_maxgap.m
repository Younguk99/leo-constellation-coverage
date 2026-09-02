function g = kmw_maxgap(u, dt)
%KMW_MAXGAP  Largest observation gap in a boolean access series, in minutes.
%
%   Faster than kmw_gaps because it returns only the maximum. Used inside the
%   greedy adversary loop, which calls it tens of thousands of times.
%
%   Gaps are measured strictly BETWEEN two access intervals, so leading and
%   trailing partial gaps are discarded.

    d  = diff(int8(u));
    st = find(d == -1);
    en = find(d ==  1);
    if isempty(st) || isempty(en)
        g = 0;                       % never loses access (or never has it)
        if ~any(u), g = NaN; end
        return
    end
    en = en(en > st(1));
    if isempty(en), g = 0; return; end
    st = st(1:numel(en));
    g  = max(en - st) * dt / 60;
end
