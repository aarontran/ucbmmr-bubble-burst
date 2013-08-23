% Attempt to solve evolution equations
% for straight retraction of a free planar film
% using an explicit finite difference scheme
% (see Erneux/Davis, Brenner/Gueyffier, Savva/Bush)

% Not even going to touch the non-Newtonian case
% until this is working first

% Aaron Tran
% August 22, 2013
% ====================================

clear
clf

% Material constants
% ------------------
mu = 100; % Dynamic viscosity (Pa s)
rho = 1000; % Fluid density (kg/m^3)
gamma = 0.073; % Surface tension (N/m)


% Initial conditions, discretization
% ----------------------------------
% index 1 located at x=0; N+1 located at x=L

Ho = 0.01; % Initial film thickness
L = 0.5; % Initial film length (fixed computational domain, does not evolve/remesh)

dt = 0.0001; % Time step (s)
dx = 0.0001; % Spatial step along film length (m)
N = round(L/dx) % Number of intervals in grid (please set dx, L so that this is an integer)
iter = 2 % Number of iterations

% Using initial profile of Appendix B, Savva 2007 (thesis)
alpha = 1/20;
h = initFilm(alpha,Ho,L,dx); % Initial film thicknesses

u = zeros(N+1,1); % Initial velocity field


hnew = zeros(N+1,1);
unew = zeros(N+1,1);

kappa = zeros(N+1,1);
hx = zeros(N+1,1);
hxx = zeros(N+1,1);


% Time evolution
% --------------

% Super inefficient right now, but simple / easier to follow
% First compute kappa everywhere
% Then compute new u,h fields

% Computation @ edges uses more complex stencils
% Stencils from http://reference.wolfram.com/mathematica/tutorial/NDSolvePDE.html

edgeIndex = 2; % Marks first index where h != 0

for i = 1:iter
	
	% Careful at edges, things blow up... hx goes to infinity, I dunno what happens to kappa...
	
	% Calculate derivatives of h, kappa
	for j = edgeIndex:N+1
		if j == edgeIndex
			hx(j) = [-25 48 -36 16 -3] * h(j:j+4) /(12*dx); % Error ~ (1/5)dx^4
			hxx(j) = [45 -154 214 -156 61 -10] * h(j:j+5) / (12*dx^2); % Error ~ (137/180)dx^4
		elseif j == N+1
			hx(j) = [3 -16 36 -48 25] * h(j-4:j) /(12*dx); % Error ~ (1/5)dx^4... should be forced to 0, but is about 0 anyways...
			hxx(j) = [-10 61 -156 214 -154 45] * h(j-5:j) / (12*dx^2); % Error ~ (137/180)dx^4
		else
			hx(j) = (h(j+1)-h(j-1))/(2*dx); % Error ~ (1/6)dx^2
			hxx(j) = (h(j+1) - 2*h(j) + h(j-1))/(dx^2); % Error ~ (1/12)dx^2
		end
		
		kappa(j) = (hxx(j)/4) / (1 + (hx(j))^2 /4)^(3/2);
	end
	
	% Calculate velocity, film thickness at edges
	uhxEdge = [-25 48 -36 16 -3] * (u(edgeIndex:edgeIndex+4) .* h(edgeIndex:edgeIndex+4)) / (12*dx);
	hnew(edgeIndex) = h(edgeIndex) - dt * uhxEdge;
	
	uxEdge = [-25 48 -36 16 -3] * u(edgeIndex:edgeIndex+4) / (12*dx);
	uxxEdge = [45 -154 214 -156 61 -10] * u(edgeIndex:edgeIndex+5) / (12*dx^2);
	kappaxEdge = [-3 4 -1] * kappa(edgeIndex:edgeIndex+2) / (2*dx);
	unew(edgeIndex) = 4*mu/rho*(uxxEdge + hx(edgeIndex)*uxEdge/h(edgeIndex)) + 2*gamma/rho * kappaxEdge;
	% this blows up: hx(edgeIndex)*uxEdge/h(edgeIndex)
	% hxEdge goes to infinity, but h(edgeIndex) goes to zero...
	
	% Calculate velocity, film thickness
	for j = edgeIndex+1:N
		ux = (u(j+1)-u(j-1))/(2*dx);
		uxx = (u(j+1) - 2*u(j) + u(j-1))/(dx^2);
		
		kappax = (kappa(j+1) - kappa(j-1))/(2*dx);
		% Alternate calculation: write out all derivatives as
		% 1/8 hx (hxx)^2 / (1 + 1/4 (hx)^2)^(5/2) + 1/4 hxxx / (1 + 1/4 (hx)^2)^(3/2)
		
		hnew(j) = h(j) - dt*(u(j+1)*h(j+1) - u(j-1)*h(j-1))/(2*dx);
		% hnew(j) = h(j) - dt/(2*dx)*(hx*u(j) + ux*h(j)); % Alternate formula, not sure which is better
		unew(j) = 4*mu/rho*(uxx + hx(j)*ux/h(j)) + 2*gamma/rho * kappax;
		
		if hnew(j) < Ho/100
			hnew(j) = 0;
			unew(j) = 0;
			edgeIndex = j+1; % SUPER AD HOC, assumes no holes...
		end
	end
	
	hnew(N+1) = hnew(N); % Neumann BC's
	unew(N+1) = unew(N); % Neumann BC's
	
	h = hnew;
    u = unew;
end

subplot(3,1,1)
plot(1:N+1,h);
xlabel("x");
ylabel("height");
subplot(3,1,2)
plot(1:N+1,u);
xlabel("x");
ylabel("veloc");
subplot(3,1,3)
plot(1:N+1,kappa);
xlabel("x");
ylabel("kappa");
