%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GERCHBERG-SAXTON (GS) ITERATIVE PHASE RETRIEVAL ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Citation for parts of this algorithm: Tatiana Latychevskaia "Iterative phase retrieval in coherent diffractive imaging:
% practical issues", Applied Optics 57(25), 7187 - 7197 (2018)

close all
clear all

Iterations = 100;  % number of iterations
N = 256;           % number of pixels
p = 0.01;          % time to pause, otherwise images are not shown

% reading diffraction pattern created by "a_simulate_DP.m"
    fid = fopen('a_dp.bin', 'r');
    dp = fread(fid, [N, N], 'real*4');
    fclose(fid);   
    imagesc(rot90(log(dp + 10))), colorbar;
    colormap(gray);
    dp_amplitude = sqrt(dp);

% reading in amplitude distribution in the sample plane: This is known
    fid = fopen('a_sample_amplitude.bin', 'r');
    sample_amplitude = fread(fid, [N, N], 'real*4');
    fclose(fid);  
 
% reading in the true (simulated) phase of the sample
    fid_phase = fopen('a_sample_phase.bin', 'r');
    true_phase = fread(fid_phase, [N, N], 'real*4');   
    fclose(fid_phase);

%Initializing entropy metric
    entropy_object = zeros(Iterations, 1); % Entropy of the object

% Guess for the complex-valued wave distribution in the detector plane (random)
    phase = (2*rand(N,N) - 1)*pi;
    field_detector = dp_amplitude.*exp(i*phase);% This is what improves with iterations

% GS iterative loop
for ii = 1:Iterations 
      fprintf('Iteration: %d\n', ii) 

    % updating complex-valued sample distribution in the object plane
      sample_phase_updated = angle(IFT2Dc(field_detector));
      sample_amplitude_updated = abs(IFT2Dc(field_detector));
      sample_phase_updated = sample_phase_updated.*sample_amplitude_updated;

      % showing updated phase distribution of the sample
      subplot(1,2,1);
      imshow(flipud(rot90(sample_phase_updated)), [0, max(sample_phase_updated(:))]);
      title('Updated Phase Distribution')
      xlabel({'x / px'})
      ylabel({'y / px'})
      axis on
      set(gca,'YDir','normal')
      colormap('gray')
      colorbar;

      % showing updated amplitude distribution of the sample
      subplot(1,2,2);
      imshow(flipud(rot90(sample_amplitude_updated)), []);
      title('Sample Amplitude')
      xlabel({'x'})
      ylabel({'y'})
      axis on
      set(gca,'YDir','normal')
      colormap('gray')
      colorbar;
      
      pause(p);
    % updating the sample distribution in the object plane
      sample_updated = sample_amplitude.*exp(i*sample_phase_updated);
      
    % updating complex-valued wavefront distribution in the detector plane
      field_detector_updated = FT2Dc(sample_updated);

    % replacing updated amplitude for the measured amplitude
      field_detector = dp_amplitude.*exp(i*angle(field_detector_updated));

    % Entropy of the object
      normalized_amplitude = sample_amplitude_updated / sum(sample_amplitude_updated(:));
      entropy_object(ii) = -sum(normalized_amplitude(:) .* log(normalized_amplitude(:) + eps));

end

% showing reconstructed phase of the sample
      figure();
      subplot(1,3,1);
      imshow(flipud(rot90(sample_phase_updated)), []);
      title('Sample Phase Distribution at Iteration: ',num2str(Iterations));
      xlabel({'x / px'})
      ylabel({'y / px'})
      axis on
      set(gca,'YDir','normal')
      colormap('gray')
      colorbar; 
% showing reconstructed amplitude of the sample
      subplot(1,3,2);
      imshow(flipud(rot90(sample_amplitude_updated)), []);
      title('Sample Amplitude Distribution')
      xlabel({'x'})
      ylabel({'y'})
      axis on
      set(gca,'YDir','normal')
      colormap('gray')
      colorbar; 
      
  % showing the real phase distribution of the sample
      subplot(1,3,3);
      imshow(flipud(rot90(true_phase)), []);
      title('True Phase Distribution')
      xlabel({'x / px'})
      ylabel({'y / px'})
      axis on
      set(gca,'YDir','normal')
      colormap('gray')
      colorbar;

      figure()
      plot(1:Iterations, entropy_object, 'LineWidth', 2);
      title('Entropy of the Object vs Iteration');
      xlabel('Iteration');
      ylabel('Entropy');
      grid on;

