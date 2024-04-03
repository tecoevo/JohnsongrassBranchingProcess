%% Generation of data from deterministic model (seed bank strength)
% Saving deterministic dynamics of a Johnsongrass population treated 
% one year with the herbicide followed by then years of no treatment
% depending on the seed bank strength. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters: 
% Field size:
A = 10^4;
% Number of years:
n_years = 11;

% Initial seedbank density: 
dens_seeds = 10;
% Initial rhizome density: 
dens_rhizomes = 1;

% 1 x n_years vector of herbicide application. Each entry corresponds 
% to one season and is a logical value stating whether the herbicide is
% applied. 
herb = [ones(1, 1), zeros(1, 10)];

% Proportion of self-pollination: 
p_self = 0.95;
% Fitness cost on seed production associated with resiance:
c = 0.3;
% Factor reducing the fitness cost of RW type relative to RR type:
k_c = 0.5;
% Factor reducing the herbicide efficiency of RW type relative to WW type:
k_h = 0.5;
% Proportion of seed germination:
g = 0.01:0.001:0.47;
% Natural yearly seed mortality in the seedbank:
d_B = 0.48 .* g ./ (1 - g) * (1 - 0.3)/0.3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3 X length(p_germ) array with genotype type frequencies (WW, RW, RR)
% in the seed bank under different seed bank strength:
Bank = zeros(3, length(g));
% 3 X length(p_germ) array with genotype type frequencies (WW, RW, RR)
% in the plants under different seed bank strength:
Plants = zeros(3, length(g));


% Initial population composition:
% Read table with genotype frewuencies at eqilibrium
T = readtable('../Data/Table_standing_variants.txt');
% Initial fraction of the RR type in seeds and plants:
RR = T.RR(round(T.Cost,4) == c & round(T.pSelf,4) == p_self);
% Initial fraction of the RW type in seeds and plants:
RW = T.RW(round(T.Cost,4) == c & round(T.pSelf,4) == p_self);

% Initial seedbank:
% Absolute genotype frequencies (WW, RW, RR) in the initial seed bank:
S0= dens_seeds * A * [1-RR-RW; RW; RR];
% Initial rhizomes:
% Absolute genotype frequencies (WW, RW, RR) in the initial rhizomes:
R0 = dens_rhizomes * A * [1-RR-RW; RW; RR];
% Plant density in presecing season:
dens0 = dens_rhizomes / 0.65;


% Loop over all parameter sets
for i1 = 1:length(g)

% gives the dynamics:
%   P: matrix of absolute genotype frequencies in plants
%   R: matrix of absolute genotype frequencies in rhizomes
%   SB: matrix of absolute genotype frequencies in seed bank
%   P_dens: vector of plant densities
[P, R, SB, P_dens] = Dynamics_deterministic(A, p_self, S0, R0, ...
    dens0, herb, n_years, c, k_c, k_h, g(i1), d_B(i1));

% Genotype frequencies expected in the long run without weed control:
Bank(:, i1) = (1-g(i1)) * SB(:, n_years)/A;
Plants(:, i1) = P(:, n_years)/A;
end

% Create a table
T = table;
% Assign columns to table
T.g = g';
T.dB = d_B';
T.WWseedDensity = Bank(1, :)';
T.RWseedDensity = Bank(2, :)';
T.RRseedDensity = Bank(3, :)';
T.WWplantDensity = Plants(1, :)';
T.RWplantDensity = Plants(2, :)';
T.RRplantDensity = Plants(3, :)';
% Write table to text file 
writetable(T, 'Table_standing_variants_seedbank_strength_1herb_10no');