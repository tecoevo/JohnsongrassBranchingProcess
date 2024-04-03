%% Generation of data from MGWP (impact of initial density)
% Escape and extinction dynamics of Johnsongrass populations modeled as  
% multiype Galton-Watson process with density dependent reproduction
% depending on the initial density.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters: 
% Number of replicates 
n = 10^2;
% Number of population replicates 
n_rep = 10^3;

% Field size:
A = 10^4;
% Number of years:
n_years = 500;

% Initial seedbank density: 
dens_seeds = (1:0.5:5) * 80;
% Initial plant density: 
dens_plants = (1:0.5:5) * 1;

% 1 x n_years vector of herbicide application. Each entry corresponds 
% to one season and is a logical value stating whether the herbicide is
% applied. 
herb = ones(1, n_years);

% Proportion of selfpollination: 
p_self = 0.95;
% Fitness cost on seed production associated with resiance:
c = 0.3;
% Factor reducing the fitness cost of RW type relative to RR type:
k_c = 0.5;
% Factor reducing the herbicide efficiency of RW type relative to WW type:
k_h = 0.5; 
% Rhizome winter mortality: 
d_Z = 0.35;
% Natural yearly seed mortality in the seedbank:
d_B = 0.48;
% Loss and natural mortality of fresh seeds over winter:
d_S = 0.94;
% Proportion of seed germination: 
g = 0.3;
% Proportion of bud sprouting (no tillage):
g_Z = 0.2;

% Number of rhizome buds produced per plant:
b = 0.93 * 140;
% Number of seeds produced per plant:
f = 0.93 * 13000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% Initial population composition:
% Read table with genotype frequencies at eqilibrium
T = readtable('../Data/Table_standing_variants.txt');
% Initial fraction of the RR type in seeds and plants:
RR = T.RR(round(T.Cost,4) == c & round(T.pSelf,4) == p_self);
% Initial fraction of the RW type in seeds and plants:
RW = T.RW(round(T.Cost,4) == c & round(T.pSelf,4) == p_self);

% n*n_rep x length(b) array with logical values stating whether the
% population escaped from control:
Escaped = zeros(n*n_rep, length(dens_seeds));
% n*n_rep x length(b) array with logical values stating whether the
% population escaped from control: the
% population went extinct:
Extinct = zeros(n*n_rep, length(dens_seeds));
% n*n_rep x length(b) array with logical values stating whether 
% RW plants established:
RWplant = zeros(n*n_rep, length(dens_seeds));
% n*n_rep x length(b) array with logical values stating whether 
% RR plants established:
RRplant = zeros(n*n_rep, length(dens_seeds));
% n*n_rep x length(b) array with times till a resistant plant establishes 
% on the field and rescues the population:
timeEscaped = NaN(n*n_rep, length(dens_seeds));
% n*n_rep x length(b) array with times till extinction in populations 
% going extinct:
timeExtinct = NaN(n*n_rep, length(dens_seeds));
% n*n_rep x length(b) array with times till RW plant establishes in 
% escaping population:
timeRWplant = NaN(n*n_rep, length(dens_seeds));
% n*n_rep x length(b) array with times time till RR plant establishes 
% in escaping population:
timeRRplant = NaN(n*n_rep, length(dens_seeds));

% n x length(b) array with proportions of simulated populations escaping 
% control:
pEscape = zeros(n, length(dens_seeds));
% n x length(b) array with proportions of simulated populations going
% extinct:
pExtinct = zeros(n, length(dens_seeds));
% n x length(b) array with proportions of simulated populations with 
% RW plants:
pRWplant = zeros(n, length(dens_seeds));
% n x length(b) array with proportions of simulated populations with 
% RR plants:
pRRplant = zeros(n, length(dens_seeds));
% n x length(b) array with average times till a resistant plant establishes 
% on the field and rescues the population:
avgTimeEscape = zeros(n, length(dens_seeds));
% n x length(b) array with average times till extinction in populations 
% going extinct:
avgTimeExtinct = zeros(n, length(dens_seeds));
% n x length(b) array with average times till RW plant establishes in 
% escaping population:
avgTimeRWplant = zeros(n, length(dens_seeds));
% n x length(b) array with average times time till RR plant establishes 
% in escaping population:
avgTimeRRplant = zeros(n, length(dens_seeds));
    

% Loop over asexual vs sexual reproduction rates
for l = 1:length(dens_seeds)

% Replicates:
for j = 1:n
for i = 1:n_rep

    % Initial seedbank:
    % Absolute genotype frequencies (WW, RW, RR) in the initial seed bank:
    S0 = poissrnd(dens_seeds(l) * A * [1-RR-RW; RW; RR]);
    % Initial plants:
    % Absolute genotype frequencies (WW, RW, RR) in the initial plants:
    P0 = poissrnd(dens_plants(l) * A * [1-RR-RW; RW; RR]);
 

    % gives the dynamics:
    %   P: matrix of absolute genotype frequencies in plants
    %   escape: logical value stating whether the population escaped from
    %   control and started to regrow
    %   extinct: logical value stating whether the population went extinct
    %   t_extinct: year in which the population went extinct
    [P, ~, ~, escape, extinct, t_extinct] = ...
        Dynamics_DensityDependance(A, p_self, S0, ...
        P0, herb, n_years, c, k_c, k_h, b, f, d_Z, d_B);

    % Save simulation results of the current run:
    Escaped((j-1)*n_rep+i, l) = escape;
    Extinct((j-1)*n_rep+i, l) = extinct;
    timeExtinct((j-1)*n_rep+i, l) = t_extinct;
    temp = find(P(2,:)>0, 1);
    if ~isempty(temp)
        timeRWplant((j-1)*n_rep+i, l) = temp - 1;
    end
    temp = find(P(3,:)>0, 1);
    if ~isempty(temp)
        timeRRplant((j-1)*n_rep+i, l) = temp - 1;
    end
    if escape
        timeEscaped((j-1)*n_rep+i, l) = ...
         min(timeRWplant((j-1)*n_rep+i, l), timeRRplant((j-1)*n_rep+i, l));
    end
end

% Logical values stating whether RW plants established:
RWplant((j-1)*n_rep+1:j*n_rep, l) = ...
    ~isnan(timeRWplant((j-1)*n_rep+1:j*n_rep, l));
% Logical values stating whether RR plants established:
RRplant((j-1)*n_rep+1:j*n_rep, l) = ...
    ~isnan(timeRRplant((j-1)*n_rep+1:j*n_rep, l));

% Proportion of simulated populations escaped from control
pEscape(j, l) = sum(Escaped((j-1)*n_rep+1:j*n_rep, l)) / n_rep;
% Proportion of simulated populations went extinct
pExtinct(j, l) = sum(Extinct((j-1)*n_rep+1:j*n_rep, l)) / n_rep;
% Proportion of simulated populations with RW plants
pRWplant(j, l) = sum(~isnan(RWplant((j-1)*n_rep+1:j*n_rep, l))) / n_rep;
% Proportion of simulated populations with RR plants
pRRplant(j, l) = sum(~isnan(RRplant((j-1)*n_rep+1:j*n_rep, l))) / n_rep;
% Average time till resistant plants establish and rescue the population
avgTimeEscape(j, l) = ...
    sum(timeEscaped((j-1)*n_rep+1:j*n_rep, l),"omitnan") / ...
    sum(Escaped(:, 1));
% Average time till extinction in populations going extinct
avgTimeExtinct(j, l) = ...
    sum(timeExtinct((j-1)*n_rep+1:j*n_rep, l),"omitnan") / ...
    sum(Extinct(:, 3));
% Average time till RW plant establishes in escaping population
avgTimeRWplant(j, l) = ...
    sum(timeRWplant((j-1)*n_rep+1:j*n_rep, l),"omitnan") / ...
    sum(~isnan(timeRWplant((j-1)*n_rep+1:j*n_rep, l)));
% Average time till RR plant establishes in escaping population
avgTimeRRplant(j, l) = ...
    sum(timeRRplant((j-1)*n_rep+1:j*n_rep, l),"omitnan") / ...
    sum(~isnan(timeRRplant((j-1)*n_rep+1:j*n_rep, l)));
end

end

% Create a table for the averages over n_rep simulated populations
T1 = table;
% Assign columns to table
T1.SeedDensity = reshape(repmat(dens_seeds, n, 1), ...
    n*length(dens_seeds), 1);
T1.PlantDensity = reshape(repmat(dens_plants, n, 1), ...
    n*length(dens_seeds), 1);
T1.Run = repmat((1:n)', length(dens_seeds), 1);
T1.pEscape = reshape(pEscape, n*length(dens_seeds), 1);
T1.timeEscape = reshape(avgTimeEscape, n*length(dens_seeds), 1);
T1.pExtinct = reshape(pExtinct, n*length(dens_seeds), 1);
T1.timeExtinct = reshape(avgTimeExtinct, n*length(dens_seeds), 1);
T1.pRWplant = reshape(pRWplant, n*length(dens_seeds), 1);
T1.timeRWplant = reshape(avgTimeRWplant, n*length(dens_seeds), 1);
T1.pRRplant = reshape(pRRplant, n*length(dens_seeds), 1);
T1.timeRRplant = reshape(avgTimeRRplant, n*length(dens_seeds), 1);
% Write table to text file 
writetable(T1, 'Table_MGWP_Initial_Density_Avg_DensityDependance93');

% Create a table for the individual runs
T2 = table;
% Assign columns to table
T2.SeedDensity = reshape(repmat(dens_seeds, n*n_rep, 1), ...
    n*n_rep*length(dens_seeds), 1);
T2.PlantDensity = reshape(repmat(dens_plants, n*n_rep, 1), ...
    n*n_rep*length(dens_seeds), 1);
T2.Run = repmat((1:n*n_rep)', length(dens_seeds), 1);
T2.Escaped = reshape(Escaped, n*n_rep*length(dens_seeds), 1);
T2.timeEscaped = reshape(timeEscaped, n*n_rep*length(dens_seeds), 1);
T2.Extinct = reshape(Extinct, n*n_rep*length(dens_seeds), 1);
T2.timeExtinct = reshape(timeExtinct, n*n_rep*length(dens_seeds), 1);
T2.RWplant = reshape(RWplant, n*n_rep*length(dens_seeds), 1);
T2.timeRWplant = reshape(timeRWplant, n*n_rep*length(dens_seeds), 1);
T2.RRplant = reshape(RRplant, n*n_rep*length(dens_seeds), 1);
T2.timeRRplant = reshape(timeRRplant, n*n_rep*length(dens_seeds), 1);
% Write table to text file 
writetable(T2, 'Table_MGWP_Initial_Density_DensityDependance93');