%% Generation of data for standing genetic variation (asexual reproduction)
% Saving final genotype frequencys in Johnsongrass after 1000 years without 
% control measures from deterministic dynamics as expected standing genetic 
% variation for different propportions of asexual reproduction.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters: 
% Field size:
A = 10^4;
% Number of years:
n_years = 1000;

% Initial seedbank density: 
dens_seeds = 10;
% Initial rhizome density: 
dens_rhizomes = 1;

% Proportion of selfpollination: 
p_self = 0.95;
% Fitness cost on seed production associated with resiance:
c = 0.3;
% Factor reducing the fitness cost of RS type relative to RR type:
k_c = 0.5;

% Maximum number of rhizome buds produced per plant:
b = 140;
% Maximum number of seeds produced per plant:
f = 13000;
% Proportion of seedgermination: 
g = 0.3;
% Proportion of bud sprouting:
g_Z = 0.2;
% Rhizome winter mortality: 
d_Z = 0.35;
% Loss and natural mortality of fresh seeds over winter:
d_S = 0.94;
% Natural yearly seed mortality in the seedbank:
d_B = 0.48;

% Rate plant offspring generated via asexual vs sexual reptoduction (no
% control):
rate = ((1-d_Z) * b * g_Z ) / (f * (1-d_S) * g/(1 - (1-d_B) * (1-g)));

% Number of rhizome buds produced per plant:
b = (0:0.001:2) * 0.9 * 140;
% Number of seeds produced per plant:
f = (1 + (1 - (0:0.001:2)) * rate) * 0.9 * 13000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% length(b) vector with WW type frequencies
% arrived in the long run without weed control:
WW = zeros(length(b), 1);
% length(b) vector with RW type frequencies
% arrived in the long run without weed control:
RW = zeros(length(b), 1);
% length(b) vector with RR type frequencies
% arrived in the long run without weed control:
RR = zeros(length(b), 1);

% Loop over all parameter sets
for i = 1:length(b) 

% Initial seedbank:
% Absolute genotype frequencies (WW, RW, RR) in the initial seed bank:
S0 = dens_seeds * A * [1; 0; 0];
% Initial rhizomes:
% Absolute genotype frequencies (WW, RW, RR) in the initial rhizomes:
R0 = dens_rhizomes * A * [1; 0; 0];
% Plant density in presecing season:
dens0 = dens_rhizomes / 0.65;

% Gives the dynamics:
%   P: matrix of absolute genotype frequencies in plants
P = DeterministicDynamics(A, p_self, f(i), b(i), S0, R0, dens0, ...
    n_years, c, k_c);

% Genotype frequencies expected in the long run without weed control:
WW(i) = P(1, end)/sum(P(:, end));
RW(i) = P(2, end)/sum(P(:, end));
RR(i) = P(3, end)/sum(P(:, end));

end

% Create a table
T = table;
% Assign columns to table
T.AsexualSexual = (b./f)';
T.WW = WW;
T.RW = RW;
T.RR = RR;
% Write table to text file 
writetable(T, 'Table_standing_variants_asexual');