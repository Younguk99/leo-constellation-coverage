function cfg = kmw_config()
%KMW_CONFIG  Baseline analysis conditions.
%
%   Conditions follow Lee et al. (2023) and Kim (2025) so that results are
%   directly comparable with the two studies closest to the actual programme.
%   Change ONLY here; every other script calls this.

cfg.T       = 40;        % total satellites
cfg.alt_km  = 490;       % gamma = 15/1 repeat ground track
cfg.inc_deg = 43;        % target latitude + 3-5 deg
cfg.inc_min = 15;        % SAR incidence angle, min  -> elevation 75
cfg.inc_max = 35;        % SAR incidence angle, max  -> elevation 55
cfg.look    = 'both';    % 'both' | 'right' | 'left'   <-- STATE THIS IN THE PAPER
cfg.days    = 12;
cfg.dt      = 20;        % s

% BMOA targets, Lee et al. (2023) Table 2   [lat lon]
cfg.targets = [38.4 126.4; 38.7 126.9; 38.9 127.5; 39.9 126.5; 40.3 127.3;
               40.8 128.5; 39.9 125.3; 40.6 126.5; 41.4 126.9];

cfg.P_list  = [2 4 5 8 10 20];   % plane counts (must divide T)
cfg.NMAX    = 8;                 % intercepts to simulate
cfg.W_list  = [30 45 60 90 120 180 240];   % target exposure deadlines (min)
cfg.W_main  = 120;               % deadline used for the Pareto figure

cfg.fs = 9;  cfg.lw = 1.3;       % figure style
end
