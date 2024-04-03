function [P] = ...
    DeterministicDynamics(A, p_self, f, b, S0, R0, dens0, n_years, c, k_c)
% Dynamics gives the deterministic dynamics of a Johnsongrass
% population over n_year years without control measures.
%
%   Input: 
%   A: field size
%   p_self: proportion of selfpollination
%   f: maximum number of seeds produced per plant
%   b: maximum number of rhizome buds produced per plant
%   S0: 3 x 1 vector of absolute genotype frequencies in initial seeds
%   R0: 3 x 1 vector of absolute genotype frequencies in initial rhizomes
%   dens0: plant density in preceding season
%   n_years: number of years
%   c: fitness cost on seed production associated with resiance 
%   k_c: factor reducing the fitness cost of RW type relative to RR type
%
%   Output:
%   P: matrix of absolute genotype frequencies in plants
%   R: matrix of absolute genotype frequencies in rhizomes
%   SB: matrix of absolute genotype frequencies in seed bank
%   P_dens: vector of plant densities

%% Simulation:
% Setting parameters:

% Ecological:
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

% Inverse of highest possible plant density:
m = 1/220;
% Area needed for a plant to produce f seeds and b buds:
a = 0.1;

% Evolutionary:
% Mutation rate:
mu = 10^(-8);


% Array of inheritance matrices. With the row j column k entry of the
% i-th matrix giving the fraction of type i seeds produced by a type j
% plant pollinated by type k pollen. (1 corresponding to genotype WW, 2 
% corresponding to genotype RW, 3 corresponding to genotype RR) 
MI = [1 0.5 0; 0.5 0.25 0; 0 0 0];
MI(:, :, 2) = [0 0.5 1; 0.5 0.5 0.5; 1 0.5 0];
MI(:, :, 3) = [0 0 0; 0 0.25 0.5; 0 0.5 1];
% Add Mutation:
M = (1 - mu)^2 * MI(:, :, 1)+ mu * (1 - mu) * MI(:, :, 2) + ...
    mu^2 * MI(:, :, 3);
M(:, :, 2) = 2 * mu * (1 - mu) * MI(:, :, 1) + ...
    ((1 - mu)^2 + mu^2) * MI(:, :, 2) + 2 * mu * (1 - mu) * MI(:, :, 3);
M(:, :, 3) = mu^2 * MI(:, :, 1) + mu * (1 - mu) * MI(:, :, 2) + ...
    (1 - mu)^2 * MI(:, :, 3);

% 3 x (n_years+1) array of genotype frequencies in the seed bank. Each  
% column corresponds to one season. Row 1 contains the numbers of WW seeds  
% at season start. Row 1 contains the numbers of RW seeds. Row 1 contains 
% the numbers of RR seeds. 
SB = zeros(3, n_years+1);

% 3 x n_years array of genotype frequencies in the plants. Each column 
% corresponds to one season. Row 1 contains the numbers of WW plants at the 
% end of the season. Row 1 contains the numbers of RW plants. Row 1 
% contains the numbers of RR plants. 
P = zeros(3, n_years);

% 3 x n_years array of genotype frequencies in the rhizomes. Each column 
% corresponds to one season. Row 1 contains the numbers of WW rhizomess at
% the end of the season. Row 1 contains the numbers of RW rhizomes. Row 1 
% contains the numbers of RR rhizomes. 
R = zeros(3, n_years+1);

% Initial seedbank:
% Absolute genotype frequencies (WW, RW, RR) in the initial seed bank:
SB(:, 1) = S0;

% Initial rhizomes:
% Absolute genotype frequencies (WW, RW, RR) in the initial rhizomes:
R(:, 1) = R0;

% Number of buds on initial rhizomes
b_t = b * (1 + a * dens0) ^(-1);

% Vector of absolute genotype frequencies (WW, RW, RR) in the produced
% seeds: 
S = zeros(3, 1);

% Loop over seasons:
for t = 1:n_years
    
    % Plants emerging from rhizomes and seeds during the season t:
    P(:, t) = g * SB(:, t) + g_Z * b_t * R(:, t);
    
    % Self-thinning:
    P(:, t) = P(:, t) * (1 + m * sum(P(:, t)) / A) ^(-1);
    
    % Density dependant reproduction:
    f_t = f * (1 + a * sum(P(:, t)) / A) ^(-1);
    b_t = b * (1 + a * sum(P(:, t)) / A) ^(-1);


    % Seeds produced during the season t:
    for i = 1:3
      S(i) = f_t * P(:, t)' * diag([1 1-k_c*c 1-c]) * ...
         ((1 - p_self) / sum(P(:, t)) * M(:, :, i) * P(:, t) + ...
         p_self * diag(M(:, :, i)));
    end

    % Seed bank present at the beginning of next season t+1:
    SB(:, t+1) = (1 - d_B) * (1 - g) * SB(:, t) + (1 - d_S) * S;
    
    % Rhizomes present at the beginning of next season t+1:
    R(:, t+1) = (1 - d_Z) * P(:, t);
end
end