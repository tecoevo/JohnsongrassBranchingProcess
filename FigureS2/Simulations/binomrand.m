% Algorithm creating bionomial random numbers. Implementation of the bnldev
% algorithm in NUMERICAL RECIPES IN FORTRAN 90 p. 155-156.
function [k] = binomrand(n, p, varargin)
% Generates  random numbers from the binomial distribution for a given 
% number of trials n and  probability of success for each trial p.
% Input:
% n: number of trials
% p: succes probabilty
% (sz: vector specifying the output size)
% (sz1, ..., szN: scalars specifying the output size)
% Output: 
% k: number of successes in n bernoulli trials with success
%    probability p

% k = binomrand(n, p)
% Input: n and p can either be vectors, matrices, or multidimensional 
% arrays of the same size or at least one argument must be scalar.
% Output: a vector, matrix, or multidimensional array k of binomial 
% distributed random numbers of the same size as the non scalar argument 
% n or p is returned.

% k = binomrand(n, p, sz)
% Input: n and p  must be scalars. Sz must be a vector. 
% Output: a vector, matrix, or multidimensional array k of random numbers 
% from the binomial distribution B(n,p) of size sz is returned.

% k = binomrand(n, p, sz1, ..., szN)
% Input: n, p and sz1, ... szN must be scalars.
% Output: a vector, matrix, or multidimensional array k of random numbers 
% from the binomial distribution B(n,p) of size s1 x ... x szN is returned.

if isnumeric(n) && isnumeric(p) && isnumeric([varargin{1:end}])
    if nargin == 2 
        size_n = size(n);
        size_p = size(p);
        if isequal(size_n, size_p) 
            k = zeros(size_p);
            for i = 1:prod(size_p)
                k(i) = binomrand_single(n(i), p(i));
            end
        elseif size_n == 1
            k = zeros(size_p);
            for i = 1:prod(size_p)
                k(i) = binomrand_single(n, p(i));
            end
        elseif size_p == 1
            k = zeros(size_n);
            for i = 1:prod(size_n)
                k(i) = binomrand_single(n(i), p);
            end
        else 
            error('Error calling binomrand(n, p). n and p must be vectors, matrices, or multidimensional arrays of the same size or one argument must be a scalar.')
        end
    
    elseif nargin == 3
        if isscalar(n) && isscalar(p) && isvector(varargin{1}) 
            k = zeros(varargin{1});
            for i = 1:prod(varargin{1})
                k(i) = binomrand_single(n, p);
            end
        else
            error('Error calling binomrand(n, p, sz). n and p must be scalars and sz a vector.')
        end
    else
        if isscalar(n) && isscalar(p) && isvector([varargin{1:end}]) 
            k = zeros([varargin{1:end}]);
            for i = 1:prod([varargin{1:end}])
                k(i) = binomrand_single(n, p);
            end
        else
            error('Error calling binomrand(n, p, s1, ..., szN). n, p, sz1, ..., szN must be scalars.')
        end
    end
else
    error('Error using binomrand. Input arguments must be numeric.')
end
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