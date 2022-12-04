clc;
close all;
clear all;

z1 = -1:0.01:1;
z2 = -1:0.01:1;
[F1,F2] = meshgrid(z1,z2); 
Function = (F2-F1).^4 + (12.*F1.*F2) - F1 + F2 - 3;
contour(F1, F2, Function, 20); 
pause(0.01);
hold on;

%particle bounds and counts
particleCount = 10;
dimension = 2;
particleLowerBound = -1;
particleUpperBound = 1;
velocityLowerBound = -1;
velocityUpperBound = 1;
iterations = 500;
%forming a iteration x 1 matrix
upperValues = zeros(iterations,1);
averageValues = zeros(iterations,1);
lowerValues = zeros(iterations,1);

%bounding the position of the particle
particlePosition = particleLowerBound + (rand(particleCount, dimension).*(particleUpperBound - particleLowerBound));
%bounding the velocity of the particle
particleVelocity = velocityLowerBound + (rand(particleCount, dimension).*(velocityUpperBound - velocityLowerBound));

%initial Position of the particle
particleBest = particlePosition;

%matrix to store particle's best
particleBestValue = zeros(particleCount,1);

for i=1:particleCount
    fx(i) = (particlePosition(i,2)-particlePosition(i,1)).^4 + 12*particlePosition(i,1)*particlePosition(i,2) - particlePosition(i,1) + particlePosition(i,2)-3;
    particleBestValue(i) = fx(i);
end
particleBestValue;
randomParticlePlot1=plot3(particlePosition(:,1),particlePosition(:,2),particleBestValue,'.k','markersize',10);

%initial Global Best
[globalBestValue , globalBestIndex] = min(particleBestValue);
globalBest = particlePosition(globalBestIndex,:);
%gbestParticlePlot=plot3(globalBest(1),globalBest(2),particleBestValue','*k','markersize',10);

%initializing velocity of the particle
%cognitive coefficient close to 2 --- c1
%social coefficient close to 2 --- c2
% Initial constant as 0.73 --- w
w = 0.729;
c_1 = 2;
c_2 = 2;
for v = 1:iterations
    for p = 1:particleCount
        particleVelocity(p,:) = w*particleVelocity(p,:) + (c_1*rand()*(particleBest(p,:) - particlePosition(p,:))) +(c_2*rand()*(globalBest - particlePosition(p,:)));
    end
    
    particleVelocity;
    
    %new positions w.r.t. velocity of the particle
    
    particlePosition = particlePosition + particleVelocity;
    
    for i=1:particleCount
        fx(i) = (particlePosition(i,2)-particlePosition(i,1)).^4 + 12*particlePosition(i,1)*particlePosition(i,2) - particlePosition(i,1) + particlePosition(i,2)-3;
        particleNewValue(i) = fx(i);
    
        if particleNewValue(i)< particleBestValue(i)
            particleBestValue(i) = particleNewValue(i);
            particleBest(i,:) = particlePosition(i,:);
        end
    end
    particleBestValue;
    
    %Updating Global Best
    [globalBestValue , globalBestIndex] = min(particleBestValue);
    globalBest = particlePosition(globalBestIndex,:);

    set(randomParticlePlot1,'xdata',particlePosition(:,1),'ydata',particlePosition(:,2),'zdata',particleNewValue);
    drawnow
    pause(0.01);
    upperValues(v) = max(particleBestValue);
    averageValues(v) = mean(particleBestValue);
    lowerValues(v) = min(particleBestValue);
end
globalBest
globalBestValue
figure;
t = 1:iterations;
plot(t, upperValues, 'o', t, averageValues, 'x', t, lowerValues, '*');
hold on;
plot(t, [upperValues averageValues lowerValues]);
hold off;
legend('upper values', 'mean values', 'lower values');
xlabel('Iterations'); ylabel('Objective function value');
