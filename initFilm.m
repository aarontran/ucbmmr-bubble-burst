% Returns a column vector representing discretized initial film thickness
% Uses the procedure of Appendix B of Savva's 2007 thesis
% Scaled to give a dimensional result, with non-unit thickness
% Even spacing, with N+1 points
function h = initFilm(alpha,ho,L,dx)
	
	% h runs from -1/2 to 1/2
	% x runs from 0 to 1/2 (and beyond)
	
	% For us...
	% h runs from -ho to + ho
	% x runs from 0 to ho
	
	% Discretization:
	% i runs from 0 to N (extends roughly to ho/dx
	% (to convert i to dimensional distance: x = i*dx)
	
	N = round(L/dx);
	
	x = transpose(((1:N+1) - 1)*dx);
	xS = (1/2) * x / ho;
	
	y = 1/2 - alpha - xS .+ (1/2) * sqrt( (1+2*alpha)^2 + 4.*xS.*(xS + 2*alpha - 1) );
	hS = 1/2 * sqrt(1 - y.^2);
	h = 2*ho*hS;
end