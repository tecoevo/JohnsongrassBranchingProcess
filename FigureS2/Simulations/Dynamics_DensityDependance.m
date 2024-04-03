function [P, B, P_dens, escape, extinct, t_extinct] = ...
Dynamics_DensityDependance(A, p_self, S0, P0, herbs, n_years, c, k_c, ...
k_h, b, f, d_Z, d_B)
% Dynamics_DensityDependance gives the dynamics, potential extinction time 
% or escape from control of a Johnsongrass population modelled as a 
% multitype Galton-Watson process with intraspecific competition (density 
% dependant reproduction) over a maximum of n_years years depending on 
% herbicde application, start population and different ecological 
% parameters.
%
%   Input: 
%   A: field size
%   p_self: proportion of selfpollination
%   S0: 3 x 1 vector of absolute genotype frequencies in initial seeds
%   P0: 3 x 1 vector of absolute genotype frequencies in initial plants
%   herbs: 1 x n_years vektor of herbicide application. Each entry 
%   corresponds to one season and is a logical value stating whether 
%   the herbicide is applied
%   n_years: number of years
%   c: fitness cost on seed production associated with resistance 
%   k_c: factor reducing the fitness cost of RS type relative to RR type
%   k_h: factor reducing the herbicide efficiency of RW type relative to 
%        WW type
%   b: number of rhizome buds produced per plant
%   f: number of seeds produced per plant
%   d_Z: rhizome winter mortality 
%   d_B: natural yearly seed mortality in the seedbank
%
%   Output:
%   P: matrix of absolute genotype frequencies in plants
%   B: matrix of absolute genotype frequencies in seed bank
%   P_dens: vector of plant densities
%   escape: logical value stating whether the population escaped from
%   control and started to regrow
%   extinct: logical value stating whether the population went extinct
%   t_extinct: year in which the population went extinct

%% Simulation:
% Setting parameters:

% Ecological:
% Proportion of seed germination: 
g = 0.3;
% Proportion of bud sprouting:
g_Z = 0.2;

% Fitness cost on seed production of the different genotypes:
cr = [0, k_c*c, c];

% Loss and natural mortality of fresh seeds over winter:
d_S = 0.94;

% Area needed for a plant to produce f seeds and b buds:
a = 0.1;

% Evolutionary:
% Mutation rate:
mu = 10^(-8);


% Antropogenic:
% Herbicide efficacy: 
% Seedlings: 
E_L = 0.998;
% Tillers: 
E_T = 0.985;


% 3x3x3 array of inheritance matrices. With the row j column k entry of 
% matrix i giving the fraction of type i seeds produced by a type j
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


% Initialisation:

% 3 x (n_years+1) array of genotype frequencies in the seed bank. Each  
% column corresponds to one season. Row 1 contains the numbers of WW seeds  
% at season start. Row 1 contains the numbers of RW seeds. Row 1 contains 
% the numbers of RR seeds. 
B = NaN(3, n_years+1);

% 3 x n_years array of genotype frequencies in the plants. Each column 
% corresponds to one season. Row 1 contains the numbers of WW plants at the 
% end of the season. Row 1 contains the numbers of RW plants. Row 1 
% contains the numbers of RR plants. 
P = NaN(3, n_years+1);

% 1 x n_years vektor of plant densities. Every entry corresponds to one
% season.
P_dens = NaN(1, n_years+1);

% Initial seed bank:
% Absolute genotype frequencies (WW, RW, RR) in the initial seed bank:
B(:, 1) = S0;

% Initial plants:
% Absolute genotype frequencies (WW, RW, RR) in the initial plants:
P(:, 1) = P0;
% Initial plant density:
P_dens(1) = sum(P0) / A;

% Logical value stating whether the population escaped from control and 
% started to regrow:
escape = 0;

% Logical value stating whether the population went extinct:
extinct = 0;

% Year in which the Population went extinct:
t_extinct = NaN;

% Loop over seasons:
for t = 1:n_years
    % Herbicide efficiency on the different genotypes:
    h_L = herbs(t)*[E_L, (1-k_h)*E_L, 0];
    h_T = herbs(t)*[E_T, (1-k_h)*E_T, 0];

    % Paramters of the offspring distribution:
    % 3 x 3 array with the means of the Poisson distributed numbers of 
    % plant offspring produced by one plant via sexual reproduction. The 
    % rows correspond to the parents type and the columns to the offprings
    % type. 
    lambda1 = ...
        f*(1+a*P_dens(t))^(-1)*(p_self* ...
        [diag(M(:,:,1)),diag(M(:,:,2)),diag(M(:,:,3))].*...
        (1 - repmat(cr',1,3)) + ...
        [M(:,1,1),M(:,1,2),M(:,1,3)].*((1 - repmat(cr',1,3)) + ...
        [0,0,0;1,1,1;1,1,1]*(1 - cr(1)))*(1 - p_self))*...
        (1 - d_S)*g.*(1 - repmat(h_L,3,1));

    % 3 x 3 array with the means of the Poisson distributed numbers of 
    % seed offspring produced by one plant via sexual reproduction. The 
    % rows correspond to the parents type and the columns to the offprings
    % type.
    lambda2 = f*(1+a*P_dens(t))^(-1)*(p_self* ...
        [diag(M(:,:,1)),diag(M(:,:,2)),diag(M(:,:,3))].*...
        (1 - repmat(cr',1,3)) + ...
        [M(:,1,1),M(:,1,2),M(:,1,3)].*((1 - repmat(cr',1,3)) + ...
        [0,0,0;1,1,1;1,1,1]*(1 - cr(1)))*(1 - p_self))*(1 - d_S)*(1-g);

    % 1 x 3 vector with the means of the Poisson distributed numbers of 
    % plant offspring produced by one plant via asexual reproduction. The 
    % rows correspond to the type.
    lambda3 = b*(1+a*P_dens(t))^(-1)*(1 - d_Z)*g_Z*(1 - h_T);

    % 3 x 3 array with the probabilities of the multinomial distributed 
    % numbers of plant and seed offspring produced by one seed. The rows 
    % correspond to the type. The first column gives the probability of 
    % producung a plant, the second of staying dormant as a seed and the 
    % third is the probability that the seed dies.
    p = [((1 - d_B)*g*(1 - h_L))', repmat((1 - d_B)*(1 - g),3,1), ...
        (d_B + (1 - d_B)*g*h_L)'];


    % Plants emerging from rhizomes and new seeds during the season t:
    P(:, t+1) = sum(poissrnd(P(:, t).*(lambda1+diag(lambda3))))';
    
    % New seeds staying dormant during season t:
    B(:, t+1) = sum(poissrnd(P(:, t).*lambda2))';

    % Fate of seeds in the seed bank:
    temp = multinomrand(B(:, t), p);

    % Plants emerging from old seeds in the seed bank during the season t:
    P(:, t+1) = P(:, t+1) + temp(:, 1);

    % Old seeds in the seed bank staying dormant during season t:
    B(:, t+1) = B(:, t+1) + temp(:, 2);

    % Plant density at the end of season t:
    P_dens(t+1) = sum(P(:, t+1)) / A;

    % Stop simulation if plant density is >10 plants/m^2 and resistant 
    % plants are present or plant and seed density is 0
    if (P_dens(t+1) > 10) && (sum(P(2:end, t+1))>0)
        escape = 1;
        return
    elseif (P_dens(t+1) + sum(B(:, t+1))) == 0
        extinct = 1;
        t_extinct = t;
        P(:, t+1:end) = 0;
        B(:, t+1:end) = 0;
        P_dens(t+1) = 0;
        return
    end

end
end