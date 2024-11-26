%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GERCHBERG-SAXTON (GS) ITERATIVE PHASE RETRIEVAL ALGORITHM WITH MAXIMUM ENTROPY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is based on the Gerchberg-Saxton iterative phase retrieval 
% algorithm as described by R. W. Gerchberg and W. O. Saxton in 
% "A practical algorithm for the determination of phase from image and 
% diffraction plane pictures" Optik 35(2), pages 237 - 246 (1972) 
% Citation: Tatiana Latychevskaia "Iterative phase retrieval in coherent diffractive imaging:
% practical issues", Applied Optics 57(25), 7187 - 7197 (2018)

close all
clear all

Iterations = 2000; % number of iterations
N = 256;           % number of pixels
p = 0.1;           % time to pause between updates

% Regularization parameter for MEM
lambda = 0.01; % Controls the influence of the entropy term

% Reading diffraction pattern created by "a_simulate_DP.m"
fid = fopen('a_dp.bin', 'r');
dp = fread(fid, [N, N], 'real*4');
fclose(fid);
imagesc(rot90(log(dp + 10))), colorbar;
colormap(gray);
dp_amplitude = sqrt(dp);

% Reading amplitude distribution in the sample plane
fid = fopen('a_sample_amplitude.bin', 'r');
sample_amplitude = fread(fid, [N, N], 'real*4');
fclose(fid);

% Initial guess for the complex-valued wave distribution in the detector plane
phase = (2*rand(N,N) - 1)*pi;
field_detector = dp_amplitude.*exp(i*phase);% This is what improves with iterations

% GS iterative loop
for ii = 1:Iterations
    fprintf('Iteration: %d\n', ii)

    % Update sample distribution in the object plane
    field_object = IFT2Dc(field_detector);
    sample_phase_updated = angle(field_object);
    sample_amplitude_updated = abs(field_object);
    sample_phase_updated = sample_phase_updated.*sample_amplitude_updated;
      
 % Visualization
    % Updated phase distribution of the sample
    subplot(1, 2, 1);
    imshow(flipud(rot90(sample_phase_updated)), []);
    title('Updated Phase Distribution of the Sample (rad)');
    xlabel('x / px');
    ylabel('y / px');
    axis on;
    set(gca, 'YDir', 'normal');
    colormap('gray');
    colorbar;

    % Updated amplitude distribution of the sample
    subplot(1, 2, 2);
    imshow(flipud(rot90(sample_amplitude_updated)), []);
    title('Sample Amplitude with MEM Regularization');
    xlabel('x / px');
    ylabel('y / px');
    axis on;
    set(gca, 'YDir', 'normal');
    colormap('gray');
    colorbar;

    pause(p);
    sample_updated = sample_amplitude.*exp(i*sample_phase_updated);


    % Apply Maximum Entropy Regularization in the object domain
    % Normalize the amplitude to avoid division by zero
    amplitude_norm = sample_amplitude_updated / max(sample_amplitude_updated(:));
    entropy_term = -amplitude_norm .* log(amplitude_norm + eps); % Entropy
    sample_amplitude_mem = sample_amplitude_updated + lambda * entropy_term;

    % Combine updated amplitude and phase
    sample_updated = sample_amplitude_mem .* exp(1i * sample_phase_updated);

    % Update the complex-valued wave distribution in the detector plane
    field_detector_updated = FT2Dc(sample_updated);

    % Replace the updated amplitude with the measured amplitude
    field_detector = dp_amplitude .* exp(1i * angle(field_detector_updated));

end

% Show the final reconstructed phase and amplitude of the sample
figure();

subplot(1,3,1);
imshow(flipud(rot90(sample_phase_updated)), []);
title('Final Phase Distribution of the Sample (rad)');
xlabel('x / px');
ylabel('y / px');
axis on;
set(gca, 'YDir', 'normal');
colormap('gray');
colorbar;

subplot(1,3,2);
imshow(flipud(rot90(sample_amplitude_mem)), []);
title('Final Amplitude Distribution of the Sample (a.u.)');
xlabel('x / px');
ylabel('y / px');
axis on;
set(gca, 'YDir', 'normal');
colormap('gray');
colorbar;


fid_phase = fopen('a_sample_phase.bin', 'r');
phase = fread(fid_phase, [N, N], 'real*4');   
fclose(fid_phase);
% showing phase distribution of the sample
subplot(1,3,3);
imshow(flipud(rot90(phase)), []);
title('simulated Phase distribution of the sample / rad')
xlabel({'x / px'})
ylabel({'y / px'})
axis on
set(gca,'YDir','normal')
colormap('gray')
colorbar;
