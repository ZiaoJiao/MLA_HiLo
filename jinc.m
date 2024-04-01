function[out]=jinc(x);
%
% jinc function
%
%J1(2*pi*x)/x
%
mask=(x~=0); %locate non-zero elements of x, take on a value of 1 for
             %any element where x is non-zero
out=pi*ones(size(x)); % initialize output with pi (value for x=0)
out(mask)=besselj(1,2*pi*x(mask))./(x(mask));
end
