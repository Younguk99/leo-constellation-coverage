function [A, info] = kmw_accessmask(P, F, cfg)
%KMW_ACCESSMASK  Access mask of a Walker-Delta constellation to ground targets.
%
%   [A, info] = KMW_ACCESSMASK(P, F, cfg)
%
%   Walker-Delta i: T/P/F, circular orbit, J2 secular propagation of RAAN and
%   argument of latitude.
%
%   OUTPUT
%     A    : logical (T x nt x ntgt).  A(k,j,m) = sat k sees target m at t(j)
%     info : .t, .period_min, .plane_of, .a_km, .track_walk_deg,
%            .lam_out_deg, .lam_in_deg, .halfwin_deg
%
%   GEOMETRY
%     For a spherical Earth the SAR incidence angle at the target equals the
%     satellite zenith angle there, so incidence = 90 - elevation.
%     Incidence 15-35 deg <=> elevation 55-75 deg.

    MU = 398600.4418;  RE = 6378.137;  J2 = 1.08262668e-3;  WE = 7.2921150e-5;

    T = cfg.T;
    if mod(T,P) ~= 0, error('T must be divisible by P.'); end
    S   = T/P;
    a   = RE + cfg.alt_km;
    inc = deg2rad(cfg.inc_deg);
    n   = sqrt(MU/a^3);

    kJ2      = J2 * (RE/a)^2;
    raan_dot = -1.5 * kJ2 * n * cos(inc);
    u_dot    = n * (1 + 0.75*kJ2*(3*cos(inc)^2 - 1)) ...
                 + 0.75 * kJ2 * n * (5*cos(inc)^2 - 1);

    t   = 0 : cfg.dt : (cfg.days*86400 - cfg.dt);
    nt  = numel(t);
    ntg = size(cfg.targets,1);

    el_min = 90 - cfg.inc_max;      % 55
    el_max = 90 - cfg.inc_min;      % 75
    s_lo = sin(deg2rad(el_min));
    s_hi = sin(deg2rad(el_max));

    % target positions in ECI (3 x nt x ntgt)
    RT = zeros(3, nt, ntg);
    for m = 1:ntg
        lat = deg2rad(cfg.targets(m,1));
        lon = deg2rad(cfg.targets(m,2)) + WE*t;
        RT(1,:,m) = RE*cos(lat)*cos(lon);
        RT(2,:,m) = RE*cos(lat)*sin(lon);
        RT(3,:,m) = RE*sin(lat);
    end

    A = false(T, nt, ntg);
    plane_of = zeros(1,T);
    k = 0;
    for p = 0:P-1
        raan = deg2rad(p*360/P) + raan_dot*t;
        cr = cos(raan);  sr = sin(raan);
        for s = 0:S-1
            k = k + 1;  plane_of(k) = p + 1;
            u  = deg2rad(s*360/S + F*p*360/T) + u_dot*t;
            cu = cos(u); su = sin(u);
            rx = a*( cu.*cr - su*cos(inc).*sr );
            ry = a*( cu.*sr + su*cos(inc).*cr );
            rz = a*( su*sin(inc) );

            if ~strcmpi(cfg.look,'both')
                dx = a*(-su.*cr - cu*cos(inc).*sr);
                dy = a*(-su.*sr + cu*cos(inc).*cr);
                dz = a*( cu*sin(inc) );
                hx = ry.*dz - rz.*dy;  hy = rz.*dx - rx.*dz;  hz = rx.*dy - ry.*dx;
            end

            for m = 1:ntg
                d0 = rx - RT(1,:,m);  d1 = ry - RT(2,:,m);  d2 = rz - RT(3,:,m);
                rr = sqrt(d0.^2 + d1.^2 + d2.^2);
                se = (d0.*RT(1,:,m) + d1.*RT(2,:,m) + d2.*RT(3,:,m)) ./ (RE*rr);
                ok = (se >= s_lo) & (se <= s_hi);
                if ~strcmpi(cfg.look,'both')
                    side = hx.*RT(1,:,m) + hy.*RT(2,:,m) + hz.*RT(3,:,m);
                    if strcmpi(cfg.look,'right'), ok = ok & (side>0);
                    else,                         ok = ok & (side<0); end
                end
                A(k,:,m) = ok;
            end
        end
    end

    lam_out = rad2deg(acos(RE/a*cos(deg2rad(el_min)))) - el_min;
    lam_in  = rad2deg(acos(RE/a*cos(deg2rad(el_max)))) - el_max;
    lat0    = mean(cfg.targets(:,1));

    % Ground-track spacing must be referenced to the NODAL period and the
    % NODAL day. Using the Keplerian period against the sidereal day omits
    % J2 nodal regression and understates the spacing (23.67 vs 24.00 deg
    % for the baseline orbit, which satisfies a gamma = 15/1 repeat).
    T_nodal = 2*pi / u_dot;                 % s, ascending-node to ascending-node
    D_nodal = 2*pi / (WE - raan_dot);       % s, Earth rotation rel. to the plane

    info.t                = t;
    info.period_min       = 2*pi/n/60;          % Keplerian
    info.nodal_period_min = T_nodal/60;
    info.nodal_day_h      = D_nodal/3600;
    info.revs_per_day     = D_nodal/T_nodal;    % ~15 for the baseline
    info.plane_of         = plane_of;
    info.a_km             = a;
    info.track_walk_deg   = 360*T_nodal/D_nodal;
    info.lam_out_deg      = lam_out;
    info.lam_in_deg       = lam_in;
    info.halfwin_deg      = lam_out*111.2/(111.2*cosd(lat0));
end
