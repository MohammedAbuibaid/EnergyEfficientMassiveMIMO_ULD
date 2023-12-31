function [Points_inside_cell] = TestPoints_MonteCarlo_HexCell(Rmax,Rmin,TestPoints,plott)

%Generate UE locations randomly with uniform distribution inside the cells

rng('shuffle');
UELocations = zeros(TestPoints,1); %Define matrix for storing UE locations
nbrToGenerate = TestPoints;        %Number of UE locations left to generate
notFinished = true(TestPoints,1);  %Indices of the UE locations that are left to generate

%Iterate the generation of UE locations until all of them are inside the
%hexagonal cell
while nbrToGenerate>0
    %Generate new UE locations uniformly at random between circles of radius dmax and dmin
    UELocations(notFinished,1) = sqrt(rand(nbrToGenerate,1)*(Rmax^2-Rmin^2)+ Rmin^2 ) .* exp(1i*2*pi*rand(nbrToGenerate,1));
    %Check which UEs that are inside a hexagonal and declare as finished    
    Angles = angle(UELocations);
    Distances = abs(UELocations);
    %Symmetry allows us to rotate all angles to lie in the area 0, pi/3
    Angles_modulus = mod(Angles,pi/3);
    %Extract the Cartesian coordinates for the rotated points
    x = Distances .* cos(Angles_modulus);
    y = Distances .* sin(Angles_modulus);
    %Check if the points are in the hexagon, in an area limited by three lines
    finished = (x < Rmax - y/sqrt(3));
    
    %Update which UEs that are left to generate
    notFinished = (finished==false);
    %Update how many UEs that are left to generate
    nbrToGenerate = sum(notFinished);
end

Points_inside_cell = UELocations;

if plott == true
    figure % Plotting Serving Cell
    plot(UELocations,'.','linewidth',1);
    title('UE distribution withn a cell')
    xlabel('Distance m') % x-axis label
    ylabel('Distance m') % y-axis label
end

end
