% Algorithm creating multionomial random numbers. Implementation of the
% conditional distribution method. The multinomial distribution is
% decomposed into marginal and conditional binomial distributions.
function [r] = multinomrand(n, p, varargin)
% Generates  random numbers from binomial distributions for a given number 
% of trials n and  a multinomial probability vector p of length k, which 
% must sum to one.
% The output r is a 1 x k vector containing the counts for each multinomial
% category.
% Input:
% n: number of trials
% p: multinomial probabilty vector
% (m: number of replicates)
% Output: 
% r: vector containing the counts for each multinomial category

% r = multinomrand(n, p)
% Input: n is a scalar specifying the number of trials and p is a vector 
% of length k containing the probabilitys of the k categories. p must sum
% up to one.
% Output: a 1 x k vector containing the counts for each multinomial
% category.

% r = multinomrand(n, p, m)
% Input: n is a scalar specifying the number of trials and p is a vector 
% of length k containing the probabilitys of the k categories. p must sum
% up to one. m is a scalar specifying the number of random vectors generated. 
% Output: a m x k matrix where each row is a multinomial outcome containing 
% the counts for each multinomial category.

% r = multinomrand(N, P)
% Input: N is a scalar or m x 1 vector specifying the number of trials. P
% is m x k matrix where each row contains a set of multinomial probabilties
% and must sum up to one.
% Output: a m x k matrix where each row is a multinomial outcome, which has 
% been generated using the corresponding rows of N and P, containing 
% the counts for each multinomial category.

if isnumeric(n) && isnumeric(p) && isnumeric([varargin{1:end}])
    size_n = size(n);
    size_p = size(p);
    if size_p(2) == 1
        p = p';
        size_p = size(p);
    end
    if nargin == 2 
        if size_n(2) == 1
            if isequal(size_n(1), size_p(1)) 
                r = zeros(size_p);
                for i = 1:size_p(1)
                    r(i,:) = multinomrand_single(n(i), p(i,:));
                end
            elseif size_n(1) == 1
                r = zeros(size_p);
                for i = 1:size_p(1)
                    r(i,:) = multinomrand_single(n, p(i,:));
                end
            elseif size_p(1) == 1
                r = zeros(size_n(1),size_p(2));
                for i = 1:size_n(1)
                    r(i,:) = multinomrand_single(n(i), p(1,:));
                end
            else 
                error('Error calling binomrand(n, p). n and p must have the same number of rows or one input argument must have only one row.')
            end
        else 
            error('Error calling multinomrand(n, p) n must be skalar or a vector of size m x 1.')
        end
    
    elseif nargin == 3
        if isscalar(n) && size_p(1) == 1 && isscalar(varargin{1}) 
            r = zeros(varargin{1}, size_p(2));
            for i = 1:varargin{1}
                r(i,:) = multinomrand_single(n, p(1,:));
            end
        else
            error('Error calling binomrand(n, p, sz). n and m must be scalars and p a vector.')
        end
    else
        error('Error calling binomrand. There must be two or three input arguments')
    end
else
    error('Error calling binomrand. Input arguments must be numeric.')
end
end


function [r] = multinomrand_single(n, p)
% Generates  random numbers from multinomial distributions for a given number 
% of trials n and  a multinomial probability vector p of length k, which 
% must sum to one.
% The output r is a 1 x k vector containing the counts for each multinomial
% category.
% Input:
% n: number of trials
% p: multinomial probabilty vector
% Output: 
% r: vector containing the counts for each multinomial category
k = size(p, 2);
p = [0, p];
pp = cumsum(p);
r = zeros(1, k+1);
temp = 0;
for i = 1:k 
    temp = temp + r(i);
    r(i+1) = binomrand_single(n - temp, p(i+1)/(1-pp(i)));
end
r(1) = [];
end

% Algorithm creating bionomial random numbers. Implementation of the bnldev
% algorithm in NUMERICAL RECIPES IN FORTRAN 90 p. 155-156.
function [k] = binomrand_single(n, pp)
% Input:
% n: number of trials
% pp: succes probabilty

% Output: 
% k: number of successes in n bernoulli trials with success
%    probability p 

p = min(pp, 1-pp); % Choose smaller probabilty, B(n-k, 1-p) equals B(n, p)
m = n * p; % Mean of the binomial distribution 

% If n < 25, use direct method:
if n < 25
    % Generate n U(0, 1) distributed random numbers and count which are < p
    k = sum(rand(n, 1) < p); 

% If less than one succes is expected out of 25 or more trials, the 
% distribution is approximately Poisson, use direct Poisson method:
elseif m < 1
    g = exp(-m);
    t = 1;
    for j = 0:n
        t = t * rand;
        if t < g 
            break
        end
    end
    k = j;

% If n > 25 and more than one succes is expected, use the rejection method
% with a Lorentzian comparison function:
else
    % Compute some quantities that are used repeatedly outside the loop
    g = gammaln(n+1);
    plog = log(p);
    qlog = log(1-p);
    sq = sqrt(2 * m * (1-p));

    t = -1;
    while rand > t
        em = -1;
        while em < 0 || em >= n+1
            y = tan(pi * rand); % deviate from a Lorentzian comparison function
            em = sq * y + m; % y, shifted and scaled
        end
        em = floor(em); % Trick for integer-valued distributions
        % The ratio of the desired distribution to the comparison function; 
        % we accept or reject by comparing it to another uniform deviate.
        t = 1.2 * sq * (1 + y^2) * exp(g - gammaln(em + 1) - ...
            gammaln(n - em +1) + em * plog + (n - em) * qlog);
    end
    k = em;
end

if p ~= pp
    k = n - k; % Undo symmetry transformation
end
end