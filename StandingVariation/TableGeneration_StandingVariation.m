%% Generation of data for standing genetic variation
% Saving final genotype frequencys in Johnsongrass after 1000 years without 
% control measures from deterministic dynamicsas expected standing genetic 
% variation for different levels of self-pollination and varying fitness 
% costs associated with resistance.

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

% Fecundity, i.e. number of seeds produced per plant:
f = 13000; 
% Number of rhizome buds produced per plant:
b = 140;

% Proportion of selfpollination: 
p_self = 0:0.001:1;
% Fitness cost on seed production associated with resiance:
c = 0.001:0.001:0.4;
% Factor reducing the fitness cost of RS type relative to RR type:
k_c = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% length(p_self) x length(c) array with WW type frequencies
% arrived in the long run without weed control:
WW = zeros(length(p_self), length(c));
% length(p_self) x length(c) array with RW type frequencies
% arrived in the long run without weed control:
RW = zeros(length(p_self), length(c));
% length(p_self) x length(c) array with RR type frequencies
% arrived in the long run without weed control:
RR = zeros(length(p_self), length(c));

% Loop over all parameter sets
for i = 1:length(c) 

% Initial seedbank:
% Absolute genotype frequencies (WW, RW, RR) in the initial seed bank:
S0 = dens_seeds * A * [1; 0; 0];
% Initial rhizomes:
% Absolute genotype frequencies (WW, RW, RR) in the initial rhizomes:
R0 = dens_rhizomes * A * [1; 0; 0];
% Plant density in presecing season:
dens0 = dens_rhizomes / 0.65;


for j = 1:length(p_self)

% Gives the dynamics:
% Matrix of absolute genotype frequencies in plants
P = DeterministicDynamics(A, p_self(j), f, b, S0, R0, dens0, n_years, ...
    c(i), k_c);

% Genotype frequencies expected in the long run without weed control:
WW(j, i) = P(1, end)/sum(P(:, end));
RW(j, i) = P(2, end)/sum(P(:, end));
RR(j, i) = P(3, end)/sum(P(:, end));

end
end

% Create a table
T = table;
% Assign columns to table
T.pSelf = repmat(p_self', length(c), 1);
T.Cost = reshape(repmat(c, length(p_self), 1), length(p_self)*length(c), 1);
T.WW = reshape(WW(:, :), length(p_self)*length(c), 1);
T.RW = reshape(RW(:, :), length(p_self)*length(c), 1);
T.RR = reshape(RR(:, :), length(p_self)*length(c), 1);
% Write table to text file 
writetable(T, 'Table_standing_variants');