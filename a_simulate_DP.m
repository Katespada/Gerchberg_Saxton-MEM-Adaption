%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION OF DIFFRCATION PATTERN
% of a sample with complex-valued transmission function
% amplitude distribution: 
% a cat cartoon, t=1 inside the cat contour and t=0 outside the cat contour
% phase distribution: 
% in form of a word "cat" with maximal phase shift of 0.1 rad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Citation for this code/algorithm or any of its parts: Tatiana Latychevskaia
% "Iterative phase retrieval in coherent diffractive imaging: practical issues",
% Applied Optics 57(25), 7187 - 7197 (2018)

close all
clear all

N = 256;                % number of pixels

% reading sample distribution
    fid_amplitude = fopen('a_sample_amplitude.bin', 'r');
    amplitude = fread(fid_amplitude, [N, N], 'real*4');
    fclose(fid_amplitude);
      % showing amplitude distribution of the sample
      figure, imshow(flipud(rot90(amplitude)), []);
      title('Amplitude distribution of the sample / a.u.')
      xlabel({'x / px'})
      ylabel({'y / px'})
      axis on
      set(gca,'YDir','normal')
      colormap('gray')
      colorbar; 
    
    fid_phase = fopen('a_sample_phase.bin', 'r');
    phase = fread(fid_phase, [N, N], 'real*4');   
    fclose(fid_phase);
      % showing phase distribution of the sample
      figure, imshow(flipud(rot90(phase)), []);
      title('Phase distribution of the sample / rad')
      xlabel({'x / px'})
      ylabel({'y / px'})
      axis on
      set(gca,'YDir','normal')
      colormap('gray')
      colorbar; 
   
    sample = amplitude.*exp(i*phase); 
   
dp = abs(FT2Dc(sample)).^2;

% saving simulated diffraction pattern as bin-file
        fid = fopen(['a_dp.bin'], 'w');
        fwrite(fid, dp, 'real*4');
        fclose(fid);         
% saving simulated diffraction pattern as jpg-file
        p = log(dp+10);
        p = (p - min(min(p)))/(max(max(p)) - min(min(p)));
        filename1 = char(strcat('a_dp.jpg'));
        imwrite(p, filename1);