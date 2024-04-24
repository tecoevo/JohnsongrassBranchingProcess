%% Generation of data from MGWP (impact of seed bank strength)
% Escape and extinction dynamics of Johnsongrass populations pre-treated 
% with the herbicide modeled as multiype Galton-Watson process depending 
% on the seed bank strength. 

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
dens_seeds = 80;
% Initial plant density: 
dens_plants = 1;

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
% Number of rhizome buds produced per plant:
b = 0.93*140;
% Number of seeds produced per plant:
f = 0.93*13000; 
% Rhizome winter mortality (no tillage): 
d_Z = 0.35;
% Proportion of seed germination:
g = 0.05:0.05:0.45;
% Natural yearly seed mortality in the seedbank:
d_B = 0.48 .* g ./ (1 - g) * (1 - 0.3)/0.3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% Initial population composition:
% Read table with genotype frequencies after pre-treatment
T = readtable('../Data/Table_standing_variants_seedbank_strength_1herb_10no.txt');

% n x length(g) array with proportions of simulated populations escaping 
% control:
pEscape = zeros(n, length(g));
% n x length(g) array with proportions of simulated populations going
% extinct:
pExtinct = zeros(n, length(g));
% n x length(g) array with proportions of simulated populations with 
% RW plants:
pRWplant = zeros(n, length(g));
% n x length(g) array with proportions of simulated populations with 
% RR plants:
pRRplant = zeros(n, length(g));
% n x length(g) array with average times till a resistant plant 
% establishes on the field and rescues the population:
avgTimeEscape = NaN(n, length(g));
% n x length(g) array with average times till extinction in populations 
% going extinct:
avgTimeExtinct = NaN(n, length(g));
% n x length(g) array with average times till RW plant establishes in 
% escaping population:
avgTimeRWplant = NaN(n, length(g));
% n x length(g) array with average times time till RR plant establishes 
% in escaping population:
avgTimeRRplant = NaN(n, length(g));

% (n_years+1) x length(g) array with probabilities of the first 
% resistant plant appearing in a given year:
pResistantPlant = zeros(n_years+1, length(g));


% Loop over seed bank mortalities 
for l = 1:length(g)

% Initial population composition:
% Initial frequency of WW seeds:
WWseedFrequency = T.WWseedFrequency(round(T.g,4) == round(g(l),4));
% Initial frequency of RW seeds:
RWseedFrequency = T.RWseedFrequency(round(T.g,4) == round(g(l),4));
% Initial frequency of RR seeds:
RRseedFrequency = T.RRseedFrequency(round(T.g,4) == round(g(l),4));
% Initial frequency of WW plants:
WWplantFrequency = T.WWplantFrequency(round(T.g,4) == round(g(l),4));
% Initial frequency of RW plants:
RWplantFrequency = T.RWplantFrequency(round(T.g,4) == round(g(l),4));
% Initial frequency of RR plants:
RRplantFrequency = T.RRplantFrequency(round(T.g,4) == round(g(l),4));

% Replicates:
for j = 1:n
    
% n_rep x 6 array with simulation results.
% Column 1 contains logical values stating whether the population escaped 
% control.
% Column 2 contains the year in which a resistant plant establishes and 
% rescues the opulation.
% Column 3 contains logical values stating whether population went extinct.
% Column 4 contains the year in which the population went extinct.
% Column 5 contains the year of first RW plant surviving till reproduction.
% Column 6 contains the year of first RR plant surviving till reproduction.
Sim = NaN(n_rep, 6);

for i = 1:n_rep

    % Initial seedbank:
    % Absolute genotype frequencies (WW, RW, RR) in the initial seed bank:
    S0 = poissrnd(dens_seeds * A * ...
        [WWseedFrequency; RWseedFrequency; RRseedFrequency]);
    % Initial plants:
    % Absolute genotype frequencies (WW, RW, RR) in the initial plants:
    P0 = poissrnd(dens_plants * A * ...
        [WWplantFrequency; RWplantFrequency; RRplantFrequency]);
 

    % gives the dynamics:
    %   P: matrix of absolute genotype frequencies in plants
    %   escape: logical value stating whether the population escaped from
    %   control and started to regrow
    %   extinct: logical value stating whether the population went extinct
    %   t_extinct: year in which the population went extinct
    [P, ~, ~, escape, extinct, t_extinct] = Dynamics(A, p_self, S0, P0, ...
        herb, n_years, c, k_c, k_h, b, f, d_Z, d_B(l), g(l));

    % Save simulation results of the current run:
    Sim(i, 1) = escape;
    Sim(i, 3) = extinct;
    Sim(i, 4) = t_extinct;
    temp = find(P(2,:)>0, 1);
    if ~isempty(temp)
        Sim(i, 5) = temp - 1; % Time of first RW plant 
    end
    temp = find(P(3,:)>0, 1);
    if ~isempty(temp)
        Sim(i, 6) = temp - 1; % Time of first RR plant 
    end
    if escape
        Sim(i, 2) = min(Sim(i, 5), Sim(i, 6)); % Time first resistant plant 
    end
    if (Sim(i, 5) <= Sim(i, 6)) || isnan(Sim(i, 6))
        if ~isnan(Sim(i, 5))
            pResistantPlant(Sim(i, 5)+1, l) = ...
            pResistantPlant(Sim(i, 5)+1, l) + 1; 
        end
    else
        pResistantPlant(Sim(i, 6)+1, l) = ...
        pResistantPlant(Sim(i, 6)+1, l) + 1;
    end
end

% Proportion of simulated populations escaped from control
pEscape(j, l) = sum(Sim(:, 1)) / n_rep;
% Proportion of simulated populations went extinct
pExtinct(j, l) = sum(Sim(:, 3)) / n_rep;
% Proportion of simulated populations with RW plants
pRWplant(j, l) = sum(~isnan(Sim(:, 5))) / n_rep;
% Proportion of simulated populations with RR plants
pRRplant(j, l) = sum(~isnan(Sim(:, 6))) / n_rep;
% Average time till resistant plants establish and rescue the population
avgTimeEscape(j, l) = sum(Sim(:, 2),"omitnan")/sum(Sim(:, 1));
% Average time till extinction in populations going extinct
avgTimeExtinct(j, l) = sum(Sim(:, 4),"omitnan")/sum(Sim(:, 3));
% Average time till RW plant establishes in escaping population
avgTimeRWplant(j, l) = sum(Sim(:, 5),"omitnan")/sum(~isnan(Sim(:, 5)));
% Average time till RR plant establishes in escaping population
avgTimeRRplant(j, l) = sum(Sim(:, 6),"omitnan")/sum(~isnan(Sim(:, 6)));
end

end

% Create a table with probabilities and average times of resistant plants
% appearing, escape and extinction
T1 = table;
% Assign columns to table
T1.g = reshape(repmat(g, n, 1), n*length(g), 1);
T1.dB = reshape(repmat(d_B, n, 1), n*length(g), 1);
T1.Run = repmat((1:n)', length(g), 1);
T1.pEscape = reshape(pEscape, n*length(g), 1);
T1.timeEscape = reshape(avgTimeEscape, n*length(g), 1);
T1.pExtinct = reshape(pExtinct, n*length(g), 1);
T1.timeExtinct = reshape(avgTimeExtinct, n*length(g), 1);
T1.pRWplant = reshape(pRWplant, n*length(g), 1);
T1.timeRWplant = reshape(avgTimeRWplant, n*length(g), 1);
T1.pRRplant = reshape(pRRplant, n*length(g), 1);
T1.timeRRplant = reshape(avgTimeRRplant, n*length(g), 1);
% Write table to text file 
writetable(T1, 'Table_MGWP_SeedBankStrength');

% Create a table with waiting time distribution til first resistant plant
T2 = table;
% Assign columns to table
T2.g = reshape(repmat(g, (n_years+1), 1), (n_years+1)*length(g), 1);
T2.Year = repmat((0:n_years)', length(g), 1);
T2.pResistantPlant = reshape(pResistantPlant/(n_rep*n), (n_years+1)*length(g), 1);
% Write table to text file 
writetable(T2, 'Table_MGWP_WaitingTime_SeedBankStrength');