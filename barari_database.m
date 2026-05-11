function db = barari_database()
% BARARI_DATABASE  Enriched bearing capacity & group efficiency database
%
% Comprehensive database from:
%   [1] Barari et al. (2023) - Table 4, Table 5, Eqs. 20-22
%       Int J Numer Anal Methods Geomech. 2023;47:2872-2898
%   [2] Kim et al. (2014) - Centrifuge test data (T1-T4)
%       J Geotech Geoenviron Eng. 140(5):04014008
%   [3] Byrne & Houlsby (2000) - 1g tests, very dense sand
%   [4] Ibsen et al. (2014,2015) - 1g tests, Dr~80%
%   [5] Gottardi & Houlsby (1999) - Surface footing tests
%
% Usage:
%   db = barari_database();
%   db.SD_range          % S/D values: 0.75 to 5.0
%   db.LD_values         % L/D values: [0.25, 0.5, 0.75, 1.0]
%   db.NcH_tri(i,:)      % N_cH for tripod at L/D(i), all S/D
%   db.EcH_tetra(i,:)    % E_cH for tetrapod at L/D(i), all S/D
%   db.NcM_tri(i,:)      % N_cM for tripod at L/D(i), all S/D
%   db.EcM_tetra(i,:)    % E_cM for tetrapod at L/D(i), all S/D
%   db.table5            % Vertical capacities (V_peak, V_t)
%   db.kim2014           % Kim et al. centrifuge data
%   db.cross_ref         % Cross-reference data from other papers

%% ===================================================================
%  S/D and L/D parameter space
%  ===================================================================
db.SD_range  = 0.75 : 0.25 : 5.0;   % 18 points
db.LD_values = [0.25, 0.50, 0.75, 1.00];
n_SD = length(db.SD_range);
n_LD = length(db.LD_values);

%% ===================================================================
%  TABLE 4 COEFFICIENTS - Barari et al. (2023)
%  Gaussian: N_cH(S/D) = a1*exp(-((S/D-b1)/c1)^2) + a2*exp(-((S/D-b2)/c2)^2)
%  Quadratic: N_cM(S/D) = p1*(S/D)^2 + p2*(S/D) + p3
%  ===================================================================

% --- N_cH coefficients [a1, b1, c1, a2, b2, c2] ---
%           L/D=0.25           L/D=0.50           L/D=0.75           L/D=1.00
NcH_tri_coeff = [
    1.37, 4.16, 5.00, 0.26, 1.40, 0.69;   % L/D=0.25
    2.56, 4.73, 4.42, 0.56, 2.01, 0.99;   % L/D=0.50
    3.62, 4.72, 1.80, 1.95, 1.98, 1.76;   % L/D=0.75
    5.40, 5.04, 3.67, 0.48, 3.98, 0.88;   % L/D=1.00
];

NcH_tetra_coeff = [
    1.18, 5.63, 5.18, 0.52, 0.92, 3.13;   % L/D=0.25
    2.53, 4.66, 2.98, 1.32, 1.64, 1.43;   % L/D=0.50
    3.74, 5.58, 2.30, 2.26, 2.32, 1.90;   % L/D=0.75
    9.70, 10.25, 6.25, 2.71, 3.15, 1.86;  % L/D=1.00
];

% --- E_cH coefficients [a1, b1, c1, a2, b2, c2] ---
EcH_tri_coeff = [
    1.47, 4.16, 5.00, 0.28, 1.40, 0.69;   % L/D=0.25
    2.04, 4.73, 4.42, 0.44, 2.01, 0.99;   % L/D=0.50
    2.73, 4.72, 1.80, 1.53, 1.98, 1.76;   % L/D=0.75
    3.92, 5.04, 3.67, 0.35, 3.98, 0.88;   % L/D=1.00
];

EcH_tetra_coeff = [
    1.46, 5.63, 5.18, 0.64, 0.92, 3.13;   % L/D=0.25
    2.31, 4.66, 2.98, 1.21, 1.64, 1.43;   % L/D=0.50
    3.08, 5.58, 3.00, 1.86, 2.32, 1.90;   % L/D=0.75
    13.21, 13.33, 7.38, 2.33, 3.18, 1.91; % L/D=1.00
];

% --- N_cM coefficients [p1, p2, p3] (quadratic: p1*x^2 + p2*x + p3) ---
NcM_tri_coeff = [
    -0.02, 0.42, 0.09;   % L/D=0.25
    -0.03, 0.51, 0.16;   % L/D=0.50
    -0.04, 0.65, 0.33;   % L/D=0.75
    -0.048, 0.75, 0.62;  % L/D=1.00
];

NcM_tetra_coeff = [
    -0.03, 0.63, 0.05;   % L/D=0.25
    -0.026, 0.62, 0.19;  % L/D=0.50
    -0.05, 0.89, 0.27;   % L/D=0.75
    0.05, 1.00, 0.50;    % L/D=1.00
];

% --- E_cM coefficients [p1, p2, p3] ---
EcM_tri_coeff = [
    -0.08, 1.42, 0.33;   % L/D=0.25
    -0.065, 1.06, 0.34;  % L/D=0.50
    -0.04, 0.74, 0.37;   % L/D=0.75
    -0.032, 0.50, 0.41;  % L/D=1.00
];

EcM_tetra_coeff = [
    -0.12, 2.11, 0.19;   % L/D=0.25
    -0.06, 1.53, 0.47;   % L/D=0.50
    -0.06, 1.12, 0.34;   % L/D=0.75
    -0.03, 0.76, 0.38;   % L/D=1.00
];

%% ===================================================================
%  GENERATE FULL DATABASE - evaluate equations at all S/D points
%  ===================================================================

% Pre-allocate: rows = L/D index, cols = S/D index
db.NcH_tri    = zeros(n_LD, n_SD);
db.NcH_tetra  = zeros(n_LD, n_SD);
db.EcH_tri    = zeros(n_LD, n_SD);
db.EcH_tetra  = zeros(n_LD, n_SD);
db.NcM_tri    = zeros(n_LD, n_SD);
db.NcM_tetra  = zeros(n_LD, n_SD);
db.EcM_tri    = zeros(n_LD, n_SD);
db.EcM_tetra  = zeros(n_LD, n_SD);

for i = 1:n_LD
    for j = 1:n_SD
        SD = db.SD_range(j);

        % Gaussian: f(x) = a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2)
        db.NcH_tri(i,j)   = gauss2(SD, NcH_tri_coeff(i,:));
        db.NcH_tetra(i,j) = gauss2(SD, NcH_tetra_coeff(i,:));
        db.EcH_tri(i,j)   = gauss2(SD, EcH_tri_coeff(i,:));
        db.EcH_tetra(i,j) = gauss2(SD, EcH_tetra_coeff(i,:));

        % Quadratic: f(x) = p1*x^2 + p2*x + p3
        db.NcM_tri(i,j)   = quad_eval(SD, NcM_tri_coeff(i,:));
        db.NcM_tetra(i,j) = quad_eval(SD, NcM_tetra_coeff(i,:));
        db.EcM_tri(i,j)   = quad_eval(SD, EcM_tri_coeff(i,:));
        db.EcM_tetra(i,j) = quad_eval(SD, EcM_tetra_coeff(i,:));
    end
end

% Store raw coefficients for reference
db.coefficients.NcH_tri   = NcH_tri_coeff;
db.coefficients.NcH_tetra = NcH_tetra_coeff;
db.coefficients.EcH_tri   = EcH_tri_coeff;
db.coefficients.EcH_tetra = EcH_tetra_coeff;
db.coefficients.NcM_tri   = NcM_tri_coeff;
db.coefficients.NcM_tetra = NcM_tetra_coeff;
db.coefficients.EcM_tri   = EcM_tri_coeff;
db.coefficients.EcM_tetra = EcM_tetra_coeff;

%% ===================================================================
%  TABLE 5 - Vertical compression & tensile capacities
%  Tripod, D=6.5m, S/D=3.0, Saemangeum sand
%  ===================================================================
db.table5.LD       = [0.25, 0.50, 0.75, 1.00];
db.table5.V_peak   = [49.1, 55.3, 62.4, 69.4];   % [MN] compression capacity
db.table5.V_t      = [-0.49, -1.98, -4.46, -7.93]; % [MN] tensile capacity
db.table5.D        = 6.5;    % [m] bucket diameter
db.table5.SD       = 3.0;    % spacing ratio
db.table5.soil     = 'Saemangeum silty sand, Dr=75%';
db.table5.source   = 'Barari et al. (2023) Table 5';

%% ===================================================================
%  FAILURE ENVELOPE PARAMETERS - Barari et al. (2023) Figures 26-27
%  Tripod, S/D=3.0, V/V_peak ~ 0.5
%  ===================================================================
db.failure_envelope.LD   = [0.25, 0.50, 0.75, 1.00];
% Uniaxial capacity parameters (linear fits from Figure 27):
%   h0 = alpha_H * (L/D) + beta_H   (normalized horizontal capacity)
%   m0 = alpha_M * (L/D) + beta_M   (normalized moment capacity)
db.failure_envelope.h0_fit = [1.94, 0.20];  % [slope, intercept] for h0 vs L/D
db.failure_envelope.m0_fit = [1.72, 0.41];  % [slope, intercept] for m0 vs L/D
db.failure_envelope.h0 = 1.94 * db.failure_envelope.LD + 0.20;
db.failure_envelope.m0 = 1.72 * db.failure_envelope.LD + 0.41;
db.failure_envelope.beta1 = [0.90, 0.85, 0.80, 0.75]; % shape parameter (decreasing with L/D)
db.failure_envelope.beta2 = 1.0;  % fixed
db.failure_envelope.eccentricity_a = -0.2;  % yield surface rotation
db.failure_envelope.source = 'Barari et al. (2023) Figures 26-27, Eqs. 1-3';

%% ===================================================================
%  HARDENING LAW PARAMETERS - Barari et al. (2023)
%  ===================================================================
db.hardening.LD    = [0.25, 0.50, 0.75, 1.00];
db.hardening.K0    = [2.18, 1.99, 1.88, 1.85];       % initial stiffness ratio
db.hardening.K_stiff = [14000, 12700, 12000, 11400];  % [kPa] stiffness
db.hardening.source = 'Barari et al. (2023) hardening law calibration';

%% ===================================================================
%  SOIL PARAMETERS - Barari et al. (2023) Table 3
%  HSsmall model for Saemangeum sand
%  ===================================================================
db.soil_barari.gamma_prime   = 9.21;     % [kN/m3] buoyant unit weight
db.soil_barari.phi_peak      = 39.18;    % [deg] peak friction angle
db.soil_barari.psi           = 9.18;     % [deg] dilatancy angle
db.soil_barari.phi_crit      = 30.0;     % [deg] critical state friction
db.soil_barari.c_prime       = 1.0;      % [kPa] effective cohesion
db.soil_barari.Dr            = 0.75;     % relative density
db.soil_barari.delta         = 22.6;     % [deg] interface friction
db.soil_barari.E50_ref       = 16235;    % [kPa] secant stiffness at sigma_ref
db.soil_barari.Eur_ref       = 48705;    % [kPa] unloading-reloading stiffness
db.soil_barari.G0_ref        = 43387;    % [kPa] small-strain shear modulus
db.soil_barari.gamma_07      = 3.44e-4;  % small-strain threshold
db.soil_barari.nu_prime      = 0.27;     % Poisson's ratio
db.soil_barari.m_power       = 0.80;     % stress-dependency power
db.soil_barari.sigma_ref     = 100;      % [kPa] reference pressure
db.soil_barari.source        = 'Barari et al. (2023) Table 3 - Saemangeum SM sand';

%% ===================================================================
%  KIM ET AL. (2014) - Centrifuge test data at 70g
%  ===================================================================
db.kim2014.source = 'Kim et al. (2014) J Geotech Geoenviron Eng. 140(5):04014008';

% Foundation geometry (prototype scale)
db.kim2014.tripod.D        = 6.5;    % [m]
db.kim2014.tripod.L        = 8.0;    % [m]
db.kim2014.tripod.LD       = 1.23;
db.kim2014.tripod.S        = 26.85;  % [m] center-to-center
db.kim2014.tripod.SD       = 4.13;
db.kim2014.tripod.eccentricity = 33; % [m] load height above seabed

db.kim2014.monopod.D       = 15.5;   % [m]
db.kim2014.monopod.L       = 10.5;   % [m]
db.kim2014.monopod.LD      = 0.68;

% Soil properties (Yellow Sea site)
db.kim2014.soil.SM_Dr      = [0.65, 0.70]; % range
db.kim2014.soil.SM_Dr_mean = 0.70;
db.kim2014.soil.SM_SPT     = [17, 39];     % range
db.kim2014.soil.ML_Dr      = [0.35, 0.65]; % range
db.kim2014.soil.ML_Dr_mean = 0.60;
db.kim2014.soil.ML_SPT     = [10, 31];     % range
db.kim2014.soil.phi_SM     = 33;           % [deg] friction angle

% Test results
db.kim2014.tests.names = {'T1_monopod', 'T2_tripod', 'T3_tripod_cyclic', 'T4_reinforced'};
db.kim2014.tests.Dr_SM = [72.4, 82.7, 75.8, 78.2]; % [%]
db.kim2014.tests.Dr_ML = [58.4, 58.4, 66.3, 66.3]; % [%]
db.kim2014.tests.type  = {'monotonic', 'monotonic', 'cyclic', 'monotonic'};

% Moment-rotation key points
db.kim2014.T1_monopod.M_yield    = 198;  % [MN.m]
db.kim2014.T1_monopod.theta_yield = 3.1; % [deg]

db.kim2014.T2_tripod.M_yield     = 93;   % [MN.m]
db.kim2014.T2_tripod.theta_yield = 0.6;  % [deg]

db.kim2014.T4_reinforced.M_yield     = 232;  % [MN.m]
db.kim2014.T4_reinforced.theta_yield = 0.65; % [deg]

% Pullout capacity
db.kim2014.pullout.T2_measured   = 2.53;  % [MN]
db.kim2014.pullout.T2_K08        = 3.69;  % [MN] estimated K=0.8
db.kim2014.pullout.T2_K05        = 2.31;  % [MN] estimated K=0.5
db.kim2014.pullout.T4_measured   = 6.61;  % [MN]
db.kim2014.pullout.delta_interface = 22;   % [deg] = 2/3*phi'

% Cyclic loading results (T3)
db.kim2014.cyclic.M_over_My = [0.39, 0.33, 0.48, 0.62, 0.70, 1.34, 0.62, 2.35, 2.51];
db.kim2014.cyclic.D_max_mm  = [0.70, 0.74, 1.23, 1.23, 1.85, 4.71, 1.23, 12.8, 50.0];
db.kim2014.cyclic.labels    = {'S1','S2','S3','S4','S5','S6','S7','S8','S9_monotonic'};

%% ===================================================================
%  CROSS-REFERENCE: Byrne & Houlsby (2000)
%  1-g tests on skirted foundations in very dense sand
%  ===================================================================
db.cross_ref.byrne_houlsby.source = 'Byrne & Houlsby (2000) - 1g tests, very dense sand';
db.cross_ref.byrne_houlsby.LD_tested = [0, 0.16, 0.33, 0.66];
db.cross_ref.byrne_houlsby.soil      = 'Very dense dry sand';
db.cross_ref.byrne_houlsby.note      = 'Beta factors needed for parabolic yield surface along V axis';

%% ===================================================================
%  CROSS-REFERENCE: Ibsen et al. (2014, 2015)
%  1-g tests, Dr~80%, five embedment ratios
%  ===================================================================
db.cross_ref.ibsen.source    = 'Ibsen et al. (2014,2015) - 1g tests, Dr~80%';
db.cross_ref.ibsen.LD_tested = [0, 0.25, 0.50, 0.75, 1.00];
db.cross_ref.ibsen.Dr        = 0.80;
db.cross_ref.ibsen.K_e       = 15000;  % [kPa] elastic unloading-reloading stiffness
db.cross_ref.ibsen.K_p       = 6000;   % [kPa] plastic stiffness
db.cross_ref.ibsen.note      = 'Foundation-soil stiffness for macro-element models';

%% ===================================================================
%  CROSS-REFERENCE: Gottardi & Houlsby (1999)
%  Surface footing tests on dense Leighton-Buzzard sand
%  ===================================================================
db.cross_ref.gottardi_houlsby.source = 'Gottardi & Houlsby (1999) - 29 loading tests';
db.cross_ref.gottardi_houlsby.soil   = 'Dry dense Leighton-Buzzard sand (14/25)';
db.cross_ref.gottardi_houlsby.Dr     = 0.75;
db.cross_ref.gottardi_houlsby.phi_triax  = 42.3;  % [deg]
db.cross_ref.gottardi_houlsby.phi_ps     = 47.8;  % [deg] plane strain
db.cross_ref.gottardi_houlsby.h0     = 0.12;  % normalized horizontal capacity (surface)
db.cross_ref.gottardi_houlsby.m0     = 0.09;  % normalized moment capacity (surface)
db.cross_ref.gottardi_houlsby.a      = -0.22; % eccentricity parameter
db.cross_ref.gottardi_houlsby.note   = 'Benchmark for surface footing (L/D=0) failure envelope';

%% ===================================================================
%  CROSS-REFERENCE: Bordon et al. - Stiffness correction factors
%  ===================================================================
db.cross_ref.bordon.source = 'Bordon et al. - Stiffness correction for polygonal foundation groups';
db.cross_ref.bordon.note   = 'K_group = N*K_f with interaction correction; >150 push-over analyses';

%% ===================================================================
%  VERTICAL BEARING CAPACITY FORMULA
%  Normalized: V_peak / (gamma_prime * D^3) = 1 + alpha * (L/D)
%  ===================================================================
db.Vpeak_formula.barari_alpha     = 2.9;   % for Saemangeum sand
db.Vpeak_formula.ahidashti_alpha  = 5.99;  % for post-liquefied soils
db.Vpeak_formula.source = 'Barari: V/(gD3)=1+2.9(L/D); Ahidashti: V/(gD3)=1+5.99(L/D)';

%% ===================================================================
%  HELPER: Lookup functions
%  ===================================================================
% Store function handles for easy interpolation
db.lookup.NcH  = @(config, LD, SD) interp_db(db, 'NcH',  config, LD, SD);
db.lookup.EcH  = @(config, LD, SD) interp_db(db, 'EcH',  config, LD, SD);
db.lookup.NcM  = @(config, LD, SD) interp_db(db, 'NcM',  config, LD, SD);
db.lookup.EcM  = @(config, LD, SD) interp_db(db, 'EcM',  config, LD, SD);
db.lookup.V_peak = @(LD) interp1(db.table5.LD, db.table5.V_peak, LD, 'pchip');
db.lookup.V_t    = @(LD) interp1(db.table5.LD, db.table5.V_t,    LD, 'pchip');

end

%% ===================================================================
%  LOCAL HELPER FUNCTIONS
%  ===================================================================

function y = gauss2(x, c)
% Two-term Gaussian: a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2)
a1 = c(1); b1 = c(2); c1 = c(3);
a2 = c(4); b2 = c(5); c2 = c(6);
y = a1 * exp(-((x - b1) / c1)^2) + a2 * exp(-((x - b2) / c2)^2);
end

function y = quad_eval(x, c)
% Quadratic: p1*x^2 + p2*x + p3
y = c(1)*x^2 + c(2)*x + c(3);
end

function val = interp_db(db, quantity, config, LD, SD)
% INTERP_DB  2D interpolation from the enriched database
%   config: 'tri' or 'tetra'
%   LD: embedment ratio (interpolated between 0.25-1.0)
%   SD: spacing ratio (interpolated between 0.75-5.0)

field = [quantity '_' config];
data = db.(field);  % n_LD x n_SD matrix

% Clamp to valid range
LD = max(0.25, min(1.00, LD));
SD = max(0.75, min(5.00, SD));

% 2D interpolation using interp2
[SD_grid, LD_grid] = meshgrid(db.SD_range, db.LD_values);
val = interp2(SD_grid, LD_grid, data, SD, LD, 'spline');
end
