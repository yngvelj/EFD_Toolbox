%% calculateVorticity   -   takes a PIV data struct and calculates the 
%%                          vorticity of the flow
%% ========================================================================
% data  -   data struct from PIV
% omega -   output vortex matrix
% -------------------------------------------------------------------------
function omega = calculateVorticity(data)

omega = cell(length(data),1);

% Calculate grid size
dx = data(1).x(2)-data(1).x(1);
dy = data(1).y(2)-data(1).y(1);

% Go through data to calculate finite differences and vorticity

for t=1:length(data)    %Go through timesteps
    omega{t} = zeros(data(1).yl,data(1).xl);
    for i=2:data(1).xl-1  
        for j=2:data(1).yl-1
            
            % Use central differences
            dudy = (data(t).velXArray(j+1,i)-...
                data(t).velXArray(j-1,i))/(2*dy);
            dvdx = (data(t).velYArray(j,i+1)-...
                data(t).velYArray(j,i-1))/(2*dx);
            omega{t}(j,i) = dvdx-dudy;
          
        end
    end
end

%Fill edges with forward and backwards differences

% for t=1:length(data)    %Go through timesteps
%     %Bottom row
%             dudy(1,:) = (data(t).velXArray(2,:)-... %forward
%                 data(t).velXArray(1,:))/(dy);
%             dvdx(1,2:end-1) = (data(t).velYArray(1,3:end)-...%central
%                 data(t).velYArray(1,1:end-2))/(2*dx);
%             dvdx(1,1) = (data(t).velYArray(1,2)-...
%                 data(t).velYArray(1,1))/(dx); %forward
%             dvdx(1,end) = (data(t).velYArray(1,end)-...
%                 data(t).velYArray(1,end-1))/(dx); %Backwards
%             
%             %Top row
%             
%             dudy(end,:) = (data(t).velXArray(end,:)-... %backwards
%                 data(t).velXArray(end-1,:))/(dy);
%             dvdx(end,2:end-1) = (data(t).velYArray(end,3:end)-...%central
%                 data(t).velYArray(end-1,1:end-2))/(2*dx);
%             dvdx(end,1) = (data(t).velYArray(end,2)-...
%                 data(t).velYArray(end,1))/(dx); %forward
%             dvdx(end,end) = (data(t).velYArray(end,end)-...
%                 data(t).velYArray(end,end-1))/(dx); %Backwards
%             
%             %Left side
%             dudy(
%             
%             
%           
%             
%             omega(j,i) = dvdx-dudy;
%         
% end


end