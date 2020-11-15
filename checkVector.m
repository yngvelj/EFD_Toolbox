%% checkVector- check the length of a vector by finding the index where the
%%              vector repeats itself/jumps down in value
%% ========================================================================
%   vec     -   vector to be checked
%   index   -   The length of the vector part that is repeated
%% ========================================================================

function index = checkVector(vec)
index = 0;      %Initialize index
for i = 2:length(vec)
   if((vec(i))<(vec(i-1)))
       index = i-1;
       return;
   end
end
end