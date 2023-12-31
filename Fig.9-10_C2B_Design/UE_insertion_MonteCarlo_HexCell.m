function [UELocations] = UE_insertion_MonteCarlo_HexCell(TestPoints,UE_pdf,Ri,plott)

Sky_Blue = [0.5294 0.8078 0.9804];
Light_Salmon = [1.0000 0.6250 0.4766];
LightGreen = [0.5625 0.9297 0.5625];
ColorMAT = [Sky_Blue;Light_Salmon;LightGreen];

NumUE = round(TestPoints*UE_pdf);

UELocStr = cell(length(NumUE),1);

for a=1:length(NumUE)-1
    rng('shuffle');
    UELocStr{a} = sqrt(rand(NumUE(a),1)*(Ri(a+1)^2-Ri(a)^2)+ Ri(a)^2 ) .* exp(1i*2*pi*rand(NumUE(a),1));
end


rng('shuffle');
UELoc = [];
nbrToGenerate = NumUE(end);
%disp(['nbrToGenerate = ' num2str(nbrToGenerate)])
notFinished = true(NumUE(end),1);  %Indices of the UE locations that are left to generate
while nbrToGenerate>0
    UELoc(notFinished,1) = sqrt(rand(nbrToGenerate,1)*(Ri(end)^2-Ri(length(NumUE))^2)+ Ri(length(NumUE))^2 ) .* exp(1i*2*pi*rand(nbrToGenerate,1));
    Angles = angle(UELoc);
    Distances = abs(UELoc);
    Angles_modulus = mod(Angles,pi/3);
    x = Distances .* cos(Angles_modulus);
    y = Distances .* sin(Angles_modulus);
    finished = (x < (Ri(end) - y/sqrt(3)));
    notFinished = (finished==false);
    nbrToGenerate = sum(notFinished);
    %disp(['nbrToGenerate = ' num2str(nbrToGenerate)])
end

UELocStr{length(NumUE)} = UELoc;
UELocations = cell2mat(UELocStr);

if plott == true
    figure
    hold on
    for a=1:length(NumUE)
    plot(real(UELocStr{a}),imag(UELocStr{a}),'.','linewidth',1,'color',ColorMAT(mod(a-1,3)+1,:));
    end
    title('UE distribution withn a cell')
    xlabel('Distance m')
    ylabel('Distance m')
    axis(500*[-1 1 -1 1])
    
end

end

