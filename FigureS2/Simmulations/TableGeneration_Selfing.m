%% Generation of data from MGWP (impact of self-pollination)
% Escape and extinction dynamics of Johnsongrass populations modeled as  
% multiype Galton-Watson process depending on the proportion of 
% self-pollination. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters: 
% Number of population replicates: 
n_rep = 10^5;

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
p_self = 0:0.5:1;
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
% Rhizome winter mortality: 
d_Z = 0.35;
% Proportion of seed germination:
g = 0.3;
% Natural yearly seed mortality in the seedbank:
d_B = 0.48;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% n_rep x length(p_self) array with logical values stating whether the
% population escaped from control:
Escaped = zeros(n_rep, length(p_self));
% n_rep x length(p_self) array with logical values stating whether the
% population escaped from control: the
% population went extinct:
Extinct = zeros(n_rep, length(p_self));
% n_rep x length(p_self) array with logical values stating whether 
% RW plants established:
RWplant = zeros(n_rep, length(p_self));
% n_rep x length(p_self) array with logical values stating whether 
% RR plants established:
RRplant = zeros(n_rep, length(p_self));
% n_rep x length(p_self) array with times till a resistant plant 
% establishes on the field and rescues the population:
timeEscaped = NaN(n_rep, length(p_self));
% n_rep x length(p_self) array with times till extinction in populations 
% going extinct:
timeExtinct = NaN(n_rep, length(p_self));
% n_rep x length(p_self) array with times till RW plant establishes in 
% escaping population:
timeRWplant = NaN(n_rep, length(p_self));
% n_rep x length(p_self) array with times till RR plant establishes 
% in escaping population:
timeRRplant = NaN(n_rep, length(p_self));
% n_rep x length(p_self) array with time of the first RW plant establishing 
% conditioned on no resistant plant established before:
timeRWplantFirst = NaN(n_rep, length(p_self));
% n_rep x length(p_self) array with time of the first RR plant establishing
% conditioned on no resistant plant established before:
timeRRplantFirst = NaN(n_rep, length(p_self));

% (n_years+1) x length(p_self) array with probabilities of the first 
% resistant plant appearing in a given year:
pResistantPlant = zeros(n_years+1, length(p_self));


% Read table with genotype frequencies at eqilibrium
T1 = readtable('../Data/Table_standing_variants.txt');

% Loop over self-pollination rates
for l = 1:length(p_self)

% Initial population composition:
% Initial fraction of the RR type in seeds and plants
RR = T1.RR(round(T1.Cost,4) == c & round(T1.pSelf,4) == round(p_self(l),4));
% Initial fraction of the RW type in seeds and plants
RW = T1.RW(round(T1.Cost,4) == c & round(T1.pSelf,4) == round(p_self(l),4));


% Replicates:
for i = 1:n_rep

    % Initial seedbank:
    % Absolute genotype frequencies (WW, RW, RR) in the initial seed bank:
    S0 = poissrnd(dens_seeds * A * [1-RR-RW; RW; RR]);
    % Initial plants:
    % Absolute genotype frequencies (WW, RW, RR) in the initial plants:
    P0 = poissrnd(dens_plants * A * [1-RR-RW; RW; RR]);
 

    % gives the dynamics:
    %   P: matrix of absolute genotype frequencies in plants
    %   escape: logical value stating whether the population escaped from
    %   control and started to regrow
    %   extinct: logical value stating whether the population went extinct
    %   t_extinct: year in which the population went extinct
    [P, ~, ~, escape, extinct, t_extinct] = Dynamics(A, p_self(l), S0, ...
        P0, herb, n_years, c, k_c, k_h, b, f, d_Z, d_B, g);

    % Save simulation results of the current run:
    Escaped(i, l) = escape;
    Extinct(i, l) = extinct;
    timeExtinct(i, l) = t_extinct;
    temp = find(P(2,:)>0, 1);
    if ~isempty(temp)
        timeRWplant(i, l) = temp - 1;
    end
    temp = find(P(3,:)>0, 1);
    if ~isempty(temp)
        timeRRplant(i, l) = temp - 1;
    end
    if escape
        timeEscaped(i, l) = min(timeRWplant(i, l), timeRRplant(i, l));
    end
    if (timeRWplant(i, l) < timeRRplant(i, l)) || isnan(timeRRplant(i, l))
        timeRWplantFirst(i, l) = timeRWplant(i, l);
        if ~isnan(timeRWplant(i, l))
            pResistantPlant(timeRWplant(i, l)+1, l) = ...
            pResistantPlant(timeRWplant(i, l)+1, l) + 1;
        end
    elseif timeRWplant(i, l) == timeRRplant(i, l)
        timeRWplantFirst(i, l) = timeRWplant(i, l);
        timeRRplantFirst(i, l) = timeRRplant(i, l);
        pResistantPlant(timeRWplant(i, l)+1, l) = ...
        pResistantPlant(timeRWplant(i, l)+1, l) + 1;
    else
        timeRRplantFirst(i, l) = timeRRplant(i, l);
        pResistantPlant(timeRRplant(i, l)+1, l) = ...
        pResistantPlant(timeRRplant(i, l)+1, l) + 1;
    end
end

% Logical vales stating whether RW plants established:
RWplant(:, l) = ~isnan(timeRWplant(:, l));
% Logical vales stating whether RR plants established:
RRplant(:, l) = ~isnan(timeRRplant(:, l));

end

% Create a table with times of resistant plants appearing, escape and
% extinction
T1 = table;
% Assign columns to table
T1.pSelf = reshape(repmat(p_self, n_rep, 1), n_rep*length(p_self), 1);
T1.Run = repmat((1:n_rep)', length(p_self), 1);
T1.Escaped = reshape(Escaped, n_rep*length(p_self), 1);
T1.timeEscaped = reshape(timeEscaped, n_rep*length(p_self), 1);
T1.Extinct = reshape(Extinct, n_rep*length(p_self), 1);
T1.timeExtinct = reshape(timeExtinct, n_rep*length(p_self), 1);
T1.RWplant = reshape(RWplant, n_rep*length(p_self), 1);
T1.timeRWplant = reshape(timeRWplant, n_rep*length(p_self), 1);
T1.RRplant = reshape(RRplant, n_rep*length(p_self), 1);
T1.timeRRplant = reshape(timeRRplant, n_rep*length(p_self), 1);
T1.timeRWplantFirst = reshape(timeRWplantFirst, n_rep*length(p_self), 1);
T1.timeRRplantFirst = reshape(timeRRplantFirst, n_rep*length(p_self), 1);
% Write table to text file 
writetable(T1, 'Table_MGWP_Selfing_n5');

% Create a table with waiting time distribution til first resistant plant
T2 = table;
% Assign columns to table
T2.pSelf = reshape(repmat(p_self, (n_years+1), 1), ...
    (n_years+1)*length(p_self), 1);
T2.Year = repmat((0:n_years)', length(p_self), 1);
T2.pResistantPlant = reshape(pResistantPlant/n_rep, ...
    (n_years+1)*length(p_self), 1);
% Write table to text file 
writetable(T2, 'Table_MGWP_WaitingTime_Selfing_n5');