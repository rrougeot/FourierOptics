%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET PARAMETERS

    % Define the coronagraph system
        Coronagraph = 3;
    % 0: Raw imaging system without occulters A+B+C+D
	% 1: External coronagraph O+A+B
	% 2: Classic Lyot coronagraph A+B+C+D
	% 3: Externally occulted Lyot coronagraph O+A+O'+C+D
    
    % Diffraction from the external occulter to be loaded
	    Psi_0_filename = './*.mat'; % [string]
	% Matlab variable name shall be Psi_0
        Psi_0_radius_filename = './*.mat'; % [string]
	% Matlab variable name shall be Psi_0_radius
	
    % CDF curve to be loaded, for adaptive solar sampling
	    CDF_filename = './*.mat'; % [string]
    % Matlab variable name shall be CDF

    % Define planes to save
        To_save = [1, 1, 1];
    % To_save = [plane B or O', plane C, plane D]
        
    % Optical parameters    
        lambda = 550; %[nm]
        z = 144.348; % [m]
        Rsun = 0.00465421; % [rad]

        R = ; % [mm]
        f = ; % [mm]
        n = ; % refractive index

        Rp = 25; % [mm]
        Rm_1 = 0; % [mm] 
        Rm_2 = 1.629; % [mm]
		Rm_3 = 0; % [mm]
        Rl = 0.97; % [Rp unit]
    
    % Numerical parameters    
        N = 2^13; % Size of the arrays

    % Sampling in plane A
		Sampling_plane_A_type = 1;
        % 0: fixed spatial sampling s_a in plane A
        % 1: smart sampling in plane A
		s_a = 15; % [micron]

    % Normalisation
        Norm_plane = 'D'; % Norm the Sun in plane 'D' or plane 'B'
        Isun = 1; % 2.08e20 * 0.1 * (2.81/3600/180*pi)^2 * 19.6;
   
    % Misalignment & Off-pointing of the Formation    
        Delta_pos = [0,0]; % [mm]
        Delta_ang = [0,0]; % [arcsec]
        
    % Optical aberrations   
        Inc_aber = [0, 0];
	% to include optical aberrations in planes [A, C]
        % Sigma_ab = [tilt_x, tilt_y, defocus, ast_x, ast_y, coma_x, coma_y, spherical]
        Sigma_ab_plane_A = [0, 0, 0, 0, 0, 0, 0, 0]; % [lambda unit] for plane A
		Sigma_ab_plane_C = [0, 0, 0, 0, 0, 0, 0, 0]; % [lambda unit] for plane C
    
    % Roughness scattering    
        Inc_scatt = [0, 0, 0]; % to include scattering in planes [A, B or O', C] 
        % PSD_param = [A [nm2.mm] amplitude, B [mm] correlation width, C slope]
        PSD_param_plane_A = [3.75e-1, 1.79e-1, 0.82]; % ABC parameters for plane A
        PSD_param_plane_BOp = [0, 0, 0]; % ABC parameters for plane B or O'
		PSD_param_plane_C = [0, 0, 0]; % ABC parameters for plane C
        seed_plane_A = 1939; % for the random process in plane A
        seed_plane_B = 1988; % for the random process in plane B or O'
        seed_plane_C = 2017; % for the random process in plane C
    
    % Source definition   
        Source_type = 1; 
        % 0: one single point source (a,b)
        % 1: integration over the full solar disc
		% 2: vignetted region 
	
	% Point source definition for Source_type == 0
        Point_source = [0,0]; % [Rsun unit]
		
    % Solar sampling type for Source_type == 1
        Solar_sampling_type = 0; 
        % 0: regular even radial sampling + circular integration
        % 1: adaptive radial sampling + circular integration
        % 2: adaptive two dimension sampling
        Nr_sun = 1000;
        Ntheta_sun = 10000;
        g = 1; 
		
	% Points source line definition for Source_type == 3
        Point_source_lim = 1; % [Rsun unit]
        Nv = 256; % Size of the vignetted zone
	   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  