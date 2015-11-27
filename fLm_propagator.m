%% fLm propagator from A. Bovet et al. NF 2014
% A. Bovet 22.11.13
% Apply the fLm propagator to any initial conditions

clear x n0 n t t0
%% test case initial distribution
% spatial range
x=[-30:.25:30];

% initial distribution n_0(x)
x0=0;
n0=normpdf(x,x0,2);

% time range
t=[10 20 30];

figure
plot(x,n0,'k')
%% center and symmetrize the domain for the convolution to work
x0=sum(n0.*x)/sum(n0);

x=x-x0;
if abs(x(1))>abs(x(end))
    [~,ind]=min(abs(x+x(end)));
    x=x(ind:end);
    n0=n0(ind:end);
else
    [~,ind]=min(abs(x+x(1)));
    x=x(1:ind);
    n0=n0(1:ind);
end
%%
% tranport exponents
alpha=0.94667;
beta=0.56075;
gam=2*beta/alpha
H=gam/2

K=0.86365; % diffusivity (fractional) in the units used in the optimizer [L^alpha/T^beta]

clear gamma
% scale factor of the Levy dist
sigma=K^(1/alpha)*gamma((beta-1)/alpha+1);
%sigma=1;

skew=-0.0034855;

% factor
c=beta^(1/alpha)*gamma((beta-1)/alpha+1);

V=0.0; % be careful, not ok when the distribution reaches the spatial boundaries

%% convolution with initial condition

% n(x,t)
n_c=zeros(length(x),length(t));

dx=mean(diff(x));

parfor i=1:length(t)
    n_c(:,i)=conv(n0,c/((t(i))^(beta/alpha))*stblpdf(c*(x-V*t(i))/((t(i))^(beta/alpha)),alpha,skew,sigma,0),'same')*dx;
end

hold on
for i=1:length(t)
    plot(x,n_c(:,i),'g');
end

%% propagate with dirac initial condition

x=[-30:.01:30];
x0=0;

% time range
t=[1:9];

n_c=zeros(length(x),length(t));
parfor i=1:length(t)
    n_c(:,i)=c/((t(i))^(beta/alpha))*stblpdf(c*(x-x0)/((t(i))^(beta/alpha)),alpha,0,sigma,0);
end

clear plt h ha
h=figure('Units', 'pixels', ...
    'Position', [100 100 450 400],'color','w');

ha=axes;
hold on

plt=plot(x,n_c,'k')

blue=[52 55 150]/255;
red=[170 20 20]/255;
set(plt,'color',blue)
set(plt,'Linewidth',1)
set(ha,'Box','on')
set(ha,'Fontsize',14,'Fontname','Helvetica','Linewidth',1.2)
set(txt,'Fontsize',14,'Fontname','Helvetica')

%% variance
clear vars
for i=1:length(t)
    vars(i)=sum(n_c(:,i)'.*(x.^2))/sum(n_c(:,i));
end
%%
figure;
plot(t,vars)
figure
loglog(t,vars)

%% s-th moment
clear mom mom0
s=2;

threshold=0.001;
for i=1:length(t)
    mom(i) = fractional_moment(x,n_c(:,i),s,threshold);
end
mom0=fractional_moment(x,n0,s,threshold);

%%
figure;
plot(t-t0,mom-mom0)
figure
loglog(t-t0,mom-mom0)

