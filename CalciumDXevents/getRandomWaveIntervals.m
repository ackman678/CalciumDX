function randWaveIntervals = getRandomWaveIntervals(Nwaves,mu,sigma,maxt,deltaspacing,constrainNwaves)
%getRandomWaveIntervals
%can be switched between random normal and random gamma distributions. Default is rand gamma distr.
%The relevant params for the gamma distribution (instead of mu and sigma) are the Shape, A and Scale, B params. Calc from gamfit in matlab or fitdistr in R
%Uncomment lines 35-39 if you want to make sure to squeeze the no. of waves into maxt
%
%Nwaves-- no. of wave intervals you want to simulate
%mu-- mean for wave interval distribution
%sigma-- std dev for wave interval distribution
%maxt-- max time for the simulated dataset, same units as mu and sigma. If the simulated Nwaves dataset can't fit within maxt, the while loop will iteratively trim down mu and sigma until Nwaves fits within maxt
%deltaspacing-- refractory period for waves. Default is 10 secs.
%
%randWaveIntervals = getRandomWaveIntervals(10, myParams.waveinterval_mu, myParams.waveinterval_sigma, myParams.maxt);
%
%James B. Ackman 2011-11-16
if nargin < 6 || isempty(constrainNwaves); constrainNwaves = 1; end;
if nargin < 5 || isempty(deltaspacing); deltaspacing = 10; end;
if nargin < 4 || isempty(maxt); maxt = Nwaves*mu; end;
if nargin < 3 || isempty(sigma); sigma = 46.3379; end;
if nargin < 2 || isempty(mu); mu = 65.15658; end;

if mu < deltaspacing
error('mu is less than deltaspacing, check input params!')
end

%       # shape         rate    <---------Calc gamma params in R on waveisiOns.s dataset--------scale = 1/rate, matlab requires scale, B
%  # 2.517740215   0.038645565 
% # (0.117240497) (0.001989708)
%R = gamrnd(A,B) % A is shape parameter, B is scale parameter

A = 2.517740215;
B = 1/0.038645565;
trimlength_shape = (A*0.01);
%randWaveIntervals = mu + sigma.*randn(Nwaves,1);
randWaveIntervals = deltaspacing + gamrnd(A,B,[Nwaves 1]);
answer = sum(randWaveIntervals);
if constrainNwaves > 0
while answer > maxt+deltaspacing
	A = A-trimlength_shape;
	randWaveIntervals = deltaspacing + gamrnd(A-trimlength_shape,B,[Nwaves 1]);
	answer = sum(randWaveIntervals);
end
end
%answer = find(randWaveIntervals < deltaspacing);
%trimlength_mu = (mu*0.05);
%trimlength_sigma = (sigma*0.10);
%while length(answer) > 0
%	newIntervals = mu + sigma.*randn(length(answer),1);
%	randWaveIntervals(answer) = newIntervals;
%%	if sum(abs(randWaveIntervals)) > maxt+deltaspacing
%%		mu = mu - trimlength_mu; %trim time until the distribution of intervals fits within max time
%%		sigma = sigma - trimlength_sigma; %trim time until the distribution of intervals fits within max time
%%		answer = 1:length(Nwaves);
%%	else
%%		answer = find(randWaveIntervals < deltaspacing);
%%	end
%	answer = find(randWaveIntervals < deltaspacing);
%end
end