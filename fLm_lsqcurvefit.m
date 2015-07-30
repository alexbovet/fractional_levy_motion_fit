function [a_f,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    fLm_lsqcurvefit(a0,x,t,n0,ndata,fitskew)
% CALL : [a_f,resnorm,residual,exitflag,output,lambda,jacobian] = fLm_lsqcurvefit(a0,x,t,n0,ndata,fitskew)
% this function convolves the flm propagator (see Bovet et al. NF 2014) with an
% initial distribution.
%   A. Bovet 28.11.2013
%   a(1) : alpha
%   a(2) : beta
%   a(3) : K
%   a(4) : skewness
%
%   x       : spatial 1D grid   (size N)
%   t       : temporal points (t>0) (size M)
%   n0      : initial distribution at t=0 (size N)
%   ndata   : profiles of the distribution to fit (size N x M)
%   fitskew : also fit the skewness? default=1

if nargin < 6
    fitskew=1;
end

% prepare data for least-square minimisation
[X,T]=meshgrid(x,t);
xdata(1,:,:)=X;
xdata(2,:,:)=T;

% bounds of parameters

if fitskew
    skew_min=-1;
    skew_max=1;
else
    skew_min=-0.001;
    skew_max=0.001;
end
amin=[0.001 0.001 0.001 skew_min];
amax=[2 2 10 skew_max];

% options of the optimizer
options = optimset('Display','iter-detailed','PlotFcns', @plotdist);


[a_f,resnorm,residual,exitflag,output,lambda,jacobian]= ...
    lsqcurvefit(@fLm_convolution,a0,xdata,ndata,amin,amax,options);


% nested objective function
    function n = fLm_convolution(a0,xdata)
        % convolution with initial condition
        
        xx=squeeze(xdata(1,1,:));
        tt=squeeze(xdata(2,:,1));
        
        % tranport exponents
        alpha=a0(1);
        beta=a0(2);
        K=a0(3); % diffusivity (fractional) in the units used in the optimizer [L^alpha/T^beta]
        sigma=abs(K^(1/alpha)*gamma((beta-1)/alpha+1));
        
        % factor
        c=beta^(1/alpha)*gamma((beta-1)/alpha+1);
        
        %V=a0(4);
        skew=a0(4);
        
        n=zeros(length(xx),length(tt));
        
        dx=mean(diff(xx));
        
        parfor i=1:length(tt)
            n(:,i)=conv(n0,c/((tt(i))^(beta/alpha))*stblpdf(c*xx/((tt(i))^(beta/alpha)),alpha,skew,sigma,0),'same')*dx;
        end
        
    end

% plot function
    function stop = plotdist(a_c, optimValues, state)
        
        stop = false;
        
        switch state
%             case 'init'
%                 hold on
%                 plot(x,ndata,'b')
            case 'iter'
                plot(x,ndata,'b')
                hold on
                plot(x,fLm_convolution(a_c,xdata),'r')
                hold off
                title(['iteration : ' num2str(optimValues.iteration)])
                disp(['---- alpha = ' num2str(a_c(1)) ', beta = ' num2str(a_c(2)) ', K = ' num2str(a_c(3)) ', Skew = ' num2str(a_c(4))])
            case 'done'
%                 hold off
            otherwise
        end
    end

end


