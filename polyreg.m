function [ sol ] = polyreg( x,y,n)
%   polyreg
%   This function fits a polynomial of order n to a given data points (x,y)
%   using regression techniques
%   ----------------------------------------------
%   inputs
%   x                       x points
%   y                       f(x)
%   n                       polynomial order; 1 for linear
%   ----------------------------------------------
%   output
%   sol.constants       	constants that need to be found
%   sol.fn                  curve fittin-g function
%   sol.syx                 the standard erro of the estimate
%   sol.sr                  the best fit criterion
%   sol.r2                  the coefficient of determination; 1 is the best value
%   sol.std                 standard deviation
%   ----------------------------------------------
%   Example 
%   x=[1,4,2,4];
%   y=[.9 3 10 15];
%   n=1;
%   sol = polyreg(x,y,n);
%   ----------------------------------------------
%
%   All copyrights goes to Mohammad Al-Fetyani
%   University of Jordan

if size(x,2) > size(x,1)
    x = x';
end

if size(y,2) > size(y,1)
    y = y';
end
nData = size(x,1); % number of data
nVars = size(x,2); % number of unknown

% Number of unknowns initialize
A=zeros(n^nVars,nVars);

%% all possible combination of x
for i=1:nVars
    A(:,i)=floor(mod((1:n^nVars)/n^(i-1),n));
end
A=[A;eye(nVars)*n];
A=A(sum(A,2)<=n,:);

%% construct array to solve
legend=cell(size(A,1),1);
X=zeros(size(x,1),size(A,1));
for i=1:size(A,1)
    currentPos=find(A(i,:));
    currentLegend='';
    expression = '';
    for j=1:length(currentPos)
        if j==1
            currentLegend=[currentLegend,'x',num2str(currentPos(j))];
            expression=[expression,'x(:,',num2str(currentPos(j)),')'];
            if A(i,currentPos(j)) > 1
                currentLegend=[currentLegend,'.^',num2str(A(i,currentPos(j)))];
                expression=[expression,'.^',num2str(A(i,currentPos(j)))];
            end
        else
            currentLegend=[currentLegend,'.*x',num2str(currentPos(j))];
            expression=[expression,'.*x(:,',num2str(currentPos(j)),')'];
            if A(i,currentPos(j)) > 1
                currentLegend=[currentLegend,'.^',num2str(A(i,currentPos(j)))];
                expression=[expression,'.^',num2str(A(i,currentPos(j)))];
            end
        end
    end
    if isempty(currentPos)
        legend{i,1} = '1';
        X(:,i)=ones(size(x,1),1);
    else
        legend{i,1}=currentLegend;
        X(:,i)=eval(expression);
    end
end
%% solution
a=(X'*X);
b=(X'*y);
sol.ab=[a,b];
c = a\b;
sol.legends = legend;
sol.constants = c;

%% make function

var='@(';
for i=1:size(A,2)
    if i==1
        var=[var,'x',num2str(i)];
    else
        var=[var,',x',num2str(i)];
    end
end
var =[var,')'];
fn='';
for i=1:length(legend)
    if i ==1 
        fn = [fn,num2str(c(i)),'.*',legend{i}];
    else
        fn = [fn,'+',num2str(c(i)),'.*',legend{i}];
    end
end

%% outputs
sol.fn =str2func([var,fn]);
xs = num2cell(x);
ymodel = zeros(size(xs,1),1);
for i=1:size(xs,1)
    ymodel(i) = sol.fn(xs{i,:});
end
sol.sr_sum = sum((y-ymodel).^2);
sol.sr_vector = (y-ymodel).^2;
sol.syx= sqrt(sol.sr_sum/(nData-nVars-1));
avgy = mean(y);
sol.st_vector=(y-avgy).^2;
sol.st_sum=sum((y-avgy).^2);
sol.r2=(sol.st_sum-sol.sr_sum)/sol.st_sum;
sol.std = sqrt(sol.st_sum/(nData-1));
end

