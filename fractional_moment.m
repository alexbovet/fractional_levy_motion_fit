function [mom, mean] = fractional_moment(x,y,s,threshold)
%FRACTIONAL_MOMENT compute moment of order s, where s can be real
%
%   call : [mom, mean] = fractional_moment(x,y,s,threshold)
%
%   mom = <|x-mean|^s>

if nargin == 3
    threshold = 0;
end

if length(x) ~= length(y)
    disp('Error : vectors must have the same length');
else
    
    ind=find(y>threshold*max(y));
    y=y(ind);
    x=x(ind);
    
    y=reshape(y,1,length(y));
    x=reshape(x,1,length(x));
    
    % normalize distribution
    y=y/nansum(y);
    
    
    mean=nansum(y.*x);
    
    mom=nansum(y.*(abs(x-mean).^s));
    
end