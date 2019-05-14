%% PROPAGATION IN THE EXTERNALLY OCCULTED LYOT SOLAR CORONAGRAPH 
% ASPIICS study case
%
% Ref: Rougeot R., et al. 2017, A&A, 590, A2
%      Rougeot R. & Aime, C. 2018, A&A
%      Rougeot R., et al. 2019, A&A, accepted

function [Results] = Coronagraph_model(Parameters,Filename_results)

%% Description
% Compute the full propagation through the coronagraphic system
% based on parameters defined in Parameters Matlab structure
%
%
% -Inputs: Parameters, Matlab structure which contains:
%   - Coronagraph: 0, 1 or 2
%     0: propagation in the raw imaging system without occulters
%     1: propagation through the classic Lyot coronagraph
%     2: propagation through externally occulted Lyot coronagraph
%   - Psi_0_filename: location of 1D diffracted wave front in plane A [string]
%   - Psi_0_radius_filename: location of 1D radius in plane A [string]
%   - CDF_filename: location of CDF curve for adaptive sampling [string]
%   - To_save: 1x3 boolean vectors to save plane B/O', C and D
%   - lambda: wavelength [nm]
%   - f: focal of the telescope [mm]
%   - z: distance external occulter - pupil [m]
%   - Rsun: angular radius of the Sun [rad]
%   - R: radius of the external occulter [mm]
%   - n: refractive index
%   - Rp: radius of the entrance pupil in plane A [mm]
%   - Rm_1: radius of the hole of the internal occulter in plane B/O' [mm]
%   - Rm_2: radius of the hole inside the internal occulter in plane B/O' [mm]
%   - Rm_3: radius of the aperture in plane B/O' [mm]
%   - Rl: radius of the Lyot stop [Rp unit]
%   - N: size of the arrays
%   - Sampling_plane_A_type: 0 or 1
%     0: fixed spatial sampling s_a in plane A
%     1: smart sampling in plane A
%   - s_a: spatial sampling in plane A if type 0 [micron]
%   - Norm_plane: string 'B' or 'D'
%     'B': Results normed to the solar brightness in plane B
%     'D': Results normed to the solar brightness in plane D (Lyot stop
%     reduction)
%   - Isun: Brightness at the centre of the solar disk image
%   - Delta_pos: 1x2 vector with lateral misalignement occulter/pupil [mm]
%   - Delta_ang: 1x2 vector with off-pointing of the formation [arcsec]
%   - Inc_aber: 1x2 bolean to include aberrations in plane A and C
%     0: no optical aberration
%     1: include optical aberrations
%   - Sigma_ab_plane_A: 1x8 vector with RMS level of the aberrations in plane A [lambda unit]
%   - Sigma_ab_plane_C: 1x8 vector with RMS level of the aberrations in plane C [lambda unit]
%     [tilt_x, tilt_y, defocus, ast_x, ast_y, coma_x, coma_y, spherical]
%   - Inc_scatt: 1x3 bolean to include scattering in plane A, B/O' and C
%     0: no roughness scattering
%     1: include roughness scattering
%   - PSD_param_plane_A: 1x3 vectors with parameters of PSD in plane A
%   - PSD_param_plane_BOp: 1x3 vectors with parameters of PSD in plane B/O'
%   - PSD_param_plane_C: 1x3 vectors with parameters of PSD in plane C
%     [A [nm2.mm] amplitude, B [mm] correlation width, C slope]
%   - seed_plane_A: int, to set the seed for the random process in plane A
%   - seed_plane_BOp: int, to set the seed for the random process in plane B/O'
%   - seed_plane_C: int, to set the seed for the random process in plane C
%   - Source_type: 0, 1 or 2 
%     0: study case with one single point (a,b)
%     1: integration over the full solar disc
%	  2: line source for vignetting study
%   - Point_source: 1x2 vector angular of the point source [Rsun unit]
%	- Solar_sampling_type: 0, 1 or 2 
%     0: regular even radial sampling + circular integration
%     1: adaptive radial sampling + circular integration
%     2: adaptive two dimension sampling
%   - Nr_sun: number of points along the segment ]0,1], excluding 0
%   - Ntheta_sun: number of rotation for the circular integration, for type 0 or 1
%     	OR number of points in the most external ring, for type 2
%   - g: >=1, defines the slope to decrease the number of points per ring
%   - Point_source_lim: starting point for the linear source [Rsun unit]
%   - Nv: Size of the vignetted zone
%   - Filename_results: the name (with .mat) to save the results
%	
%
% Outputs:
% - Results: Matlab structure containing all the relevant outputs

% Load the parameters from the structure
Parameters

save('Temp.mat','-struct','Parameters')
load('Temp.mat')
delete('Temp.mat');

if nargin<2
    Filename_results='';
end

Results = struct();

%% Units in meters/rad

    lambda = lambda/10^9;
    R = R/1000;
    f = f/1000;
    Rp = Rp/1000;
    s_a = s_a/10^6;
    Rm_1 = Rm_1/1000;
    Rm_2 = Rm_2/1000;
	Rm_3 = Rm_3/1000;

    Delta_pos = Delta_pos/1000;
    Delta_ang = Delta_ang/3600/180*pi;

    PSD_param_plane_A(1) = PSD_param_plane_A(1) * 10^-21;
    PSD_param_plane_A(2) = PSD_param_plane_A(2) * 10^-3;
    PSD_param_plane_BOp(1) = PSD_param_plane_BOp(1) * 10^-21;
    PSD_param_plane_BOp(2) = PSD_param_plane_BOp(2) * 10^-3;
    PSD_param_plane_C(1) = PSD_param_plane_C(1) * 10^-21;
    PSD_param_plane_C(2) = PSD_param_plane_C(2) * 10^-3;

%% Variables and arrays

    disp('Initialization')

	% Sampling in plane A
    if(Sampling_plane_A_type)
		if(Coronagraph == 0 || Coronagraph == 1)
			s_a = sqrt( (Rp * lambda)/(N * Rsun)); % [m]
		elseif(Coronagraph == 2)
			s_a = sqrt( (Rp * lambda * z)/(N * R)); % [m]
		end
    end
	
	% Sampling in plane B/O'
    if(Coronagraph == 0 || Coronagraph == 1)
    	s_bop = lambda * f / s_a / N; % [m]
    elseif(Coronagraph == 2)
		z1 = f*z/(z-f); % [m]
    	s_bop = lambda * z1 / s_a / N; % [m]
    end
	
	% Sampling in plane D
    s_sun = lambda / Rsun / s_a / N; % [/Rsun] 

	% Vector and arrays
    v = ((-1/2):(1/N):(1/2 - 1/N))*N; % axis vector
    [x,y] = meshgrid(v,v); % pixel x,y meshgrid
    r = sqrt(x.^2 + y.^2); % pixel radius meshgrid  
    vp = (1:1:(N/2-2)); % vector for circular integration

	% Coronagraph elements
    P = (r*s_a) <= Rp; % Pupil
    M = ((r*s_bop) < Rm_1) + ((r*s_bop) > Rm_2); % Lyot mask / internal occulter
	if(Rm_3 > Rm_2)
		M = M.*((r*s_bop) < Rm_3);
	end
    L = (r*s_a) <= (Rp*Rl); % Lyot Stop
    
	% Numerical normalization factors
    Ap = sum(P(:)); % pixel area of the pupil
    if(strcmp(Norm_plane,'D'))
        G1 = 1/Rl^2; % norm to Lyot stop throughput
    elseif(strcmp(Norm_plane,'B'))
        G1 = 1;
    end
    
	% Fresnel propagator
    if(Coronagraph == 2)
        F = Fresnel(lambda, -z, r, s_a);
    end
 
    % Load diffraction complex amplitude
    if(Coronagraph == 2)
        load(Psi_0_filename);
        load(Psi_0_radius_filename);
    end
    
    % Load CDF (ASPIICS, sharp-edged disk, current parameters)
        load(CDF_filename)

    % Arrays for intensities
	if (Source_type == 0 || Source_type == 1)
		I_BOp_2D = zeros(N,N);
		I_C_2D = zeros(N,N);
		I_D_2D = zeros(N,N);
	elseif (Source_type == 2)
		H = zeros(Nv, Nv, Nv); % First index is the point source
	end
	
	
 %% Source definition
    
    if (Source_type == 0) % One single point source
	
        a = Point_source(1);
        b = Point_source(2);
        Ntot = 1;
		
		LB = 1;
		dS = 1;
        G2 = 1/Ap;
        		
	elseif (Source_type == 1) % Solar disc
	
        [a, b, dS, Ntot] = Adaptive_sampling(Solar_sampling_type, Nr_sun, Ntheta_sun, g, CDF);
        
		LB = Limb_darkening(a, b);
        G2 = 1/Ap * (1/s_sun)^2; % Normalization factor
        
    elseif (Source_type == 2) % Vignetted region
        
		a = (0:1:(Nv-1))*s_sun + Point_source_lim;
        b = zeros(size(a));
		Ntot = Nv;
		
		G2 = 1/Ap;
                
		Zone_x = ((N/2+1)-Nv/2) : ((N/2+1)+(Nv/2-1));
		Zone_y = ((N/2+1) + floor(1/s_sun)) + (0:(Nv-1));
        
    end

    % Additional vector to compute CDF (comment/uncomment)
    % I_PS = zeros(size(a));
    
	
%% Perturbations
	
	% Aberrations in plane A
    if(Inc_aber(1))
         Wab_a = Aberrations(Sigma_ab_plane_A, s_a, Rp, x, y);
    else Wab_a = 1;
    end

	% Aberrations in plane C
    if(Inc_aber(2))
         Wab_c = Aberrations(Sigma_ab_plane_C, s_a, Rp, x, y);
    else Wab_c = 1;
    end
     
	% Roughness in plane A	 
    if(Inc_scatt(1))
         [Wr_a, h_a] = Scattering(PSD_param_plane_A, seed_plane_A, lambda, n, N, s_a); 
         sigma_r_a = sum( h_a(:).^2 ) / N^2;
         sigma_r_a = sqrt(sigma_r_a);                
    else
         Wr_a = 1;
         h_a = 0;
         sigma_r_a = 0;
    end

	% Roughness in plane B/O'
    if(Inc_scatt(2))
         [Wr_bop, h_bop] = Scattering(PSD_param_plane_BOp, seed_plane_BOp, lambda, n, N, s_bop); 
         sigma_r_bop = sum( h_bop(:).^2 ) / N^2;
         sigma_r_bop = sqrt(sigma_r_bop);    
    else
         Wr_bop = 1;
         h_bop = 0;
         sigma_r_bop = 0;
    end

	% Roughness in plane C
    if(Inc_scatt(3))
         [Wr_c, h_c] = Scattering(PSD_param_plane_C, seed_plane_C, lambda, n, N, s_a); 
         sigma_r_c = sum( h_c(:).^2 ) / N^2;
         sigma_r_c = sqrt(sigma_r_c);                
    else
         Wr_c = 1;
         h_c = 0;
         sigma_r_c = 0;
    end
    
	% Off-poiting
    dalpha = Delta_ang(1);
    dbeta = Delta_ang(2);
        
	% Misalignments
    dx = Delta_pos(1);
    dy = Delta_pos(2);
         

%% Integration over the sources

    disp('Start propagation')
	
    for i = 1:Ntot
        
		% Angular coordinates in radian
        alpha = a(i)*Rsun;
        beta = b(i)*Rsun;
        
		% Tilt of the wavefront
        T = Tilt(lambda, alpha, beta, x, y, s_a);
        
		% Diffracted wavefront
        if (Coronagraph == 2)
             Wd = Diffraction_interp(alpha+dalpha, beta+dbeta, x, y, dx, dy, z, s_a, Psi_0_radius, Psi_0); 
        else Wd = 1;   
        end
        
		% Wave propagation
        if (Coronagraph == 0)
            Psi_A = P .* Wab_a .* Wr_a .* T;
            Psi_BOp = Fourier( Psi_A, N);
            Psi_C = Fourier( Psi_BOp .* Wr_bop, N);
            Psi_D = Fourier( L .* Psi_C .* Wab_c .* Wr_c, N);
		elseif(Coronagraph == 1)
            Psi_A = P .* Wab_a .* Wr_a .* T;
            Psi_BOp = Fourier( Psi_A, N);
            Psi_C = Fourier( M .* Psi_BOp .* Wr_bop, N);
            Psi_D = Fourier( L .* Psi_C .* Wab_c .* Wr_c, N);
        elseif (Coronagraph == 2)
            Psi_A = P .* Wd .* Wab_a .* Wr_a .* T;
            Psi_BOp = Fourier( Psi_A .* F, N);
            Psi_C = Fourier( M .* Psi_BOp .* Wr_bop, N) .* conj(F);
            Psi_D = Fourier( L .* Psi_C .* Wab_c .* Wr_c, N);
        end
        
		% Elementary intensity computation
		if(Source_type == 0 || Source_type == 1)
		
        if (Coronagraph == 0)
            I_BOp_2D = I_BOp_2D + dS(i) * G2 * LB(i) * abs(Psi_BOp).^2;
            I_C_2D = I_C_2D + dS(i) * G2 * LB(i) * abs(Psi_C).^2;
            I_D_2D = I_D_2D + dS(i) * G2 * G1 * LB(i) * abs(Psi_D).^2;
        elseif (Coronagraph == 1)
            I_BOp_2D = I_BOp_2D + dS(i) * G2 * LB(i) * abs(Psi_BOp).^2;
            I_C_2D = I_C_2D + dS(i) * G2 * LB(i) * abs(Psi_C).^2;
            I_D_2D = I_D_2D + dS(i) * G2 * G1 * LB(i) * abs(Psi_D).^2;
        elseif (Coronagraph == 2)
            I_BOp_2D = I_BOp_2D + dS(i) * G2 * LB(i) * abs(Psi_BOp).^2;
            I_C_2D = I_C_2D + dS(i) * G2 * LB(i) * abs(Psi_C).^2;
            I_D_2D = I_D_2D + dS(i) * G2 * G1 * LB(i) * abs(Psi_D).^2;
        end
		
		elseif(Source_type == 2)
			H(i, :,:) = G2 * G1 * abs(Psi_D(Zone_x, Zone_y)).^2;
		end
		
        % Additional line to compute CDF (comment/uncomment)
        % I_PS(i) = sum( abs(Psi_D(:)).^2 );
        
        disp(strcat(num2str(i/Ntot *100),'%'));

    end

    disp('End propagation')
    
%% Integration and observed intensities

	% Single point source: 2D images
    if (Source_type == 0) 
		if(To_save(1))
		if (Coronagraph == 0 || Coronagraph == 1)
		Results.I_B = I_BOp_2D;	
		elseif (Coronagraph == 2)
        Results.I_Op = I_BOp_2D;
        end
		end
			
        if(To_save(2))
        Results.I_C = I_C_2D;
        end
			
        if(To_save(3))
        Results.I_D = I_D_2D;
        end		
	end
	
	% Solar disc	 
    if (Source_type == 1)
	
		% Radial intensity
		if (Solar_sampling_type == 0 || Solar_sampling_type == 1)
            I_BOp_1D = Circular_Int(v,v,I_BOp_2D,Ntheta_sun) * Isun; 
            I_C_1D = Circular_Int(v,v,I_C_2D,Ntheta_sun) * Isun; 
            I_D_1D = Circular_Int(v,v,I_D_2D,Ntheta_sun) * Isun;
            
			if(To_save(1))
			if (Coronagraph == 0 || Coronagraph == 1)
			Results.I_B = I_BOp_1D;	
			elseif (Coronagraph == 2)
            Results.I_Op = I_BOp_1D;
            end
			end
			
            if(To_save(2))
            Results.I_C = I_C_1D;
            end
			
            if(To_save(3))
            Results.I_D = I_D_1D;
            end
		
		% 2D intensity
		elseif (Solar_sampling_type == 2)  
        
		    I_BOp_2D = I_BOp_2D + flipud(I_BOp_2D);
            I_BOp_2D = I_BOp_2D * Isun;

            I_C_2D = I_C_2D + flipud(I_C_2D);
            I_C_2D = I_C_2D * Isun;

            I_D_2D = I_D_2D + flipud(I_D_2D);
            I_D_2D = I_D_2D * Isun;
        		
            if(To_save(1))
			if (Coronagraph == 0 || Coronagraph == 1)
			Results.I_B = I_BOp_2D;	
			elseif (Coronagraph == 2)
            Results.I_Op = I_BOp_2D;
            end
			end
			
            if(To_save(2))
            Results.I_C = I_C_2D;
            end
			
            if(To_save(3))
            Results.I_D = I_D_2D;
            end
         end
    end
   
    % Image in the vignetting zone
    if (Source_type == 2)
		Results.H = H;
		Results.Zone_x = Zone_x;
		Results.Zone_y = Zone_y;
    end
 
    
%% Useful Outputs
	
    Results.v = v;
    Results.vp = vp;
    Results.s_a = s_a;
    Results.s_sun = s_sun;
    Results.s_bop = s_bop;

    if(Coronagraph == 2)
		Results.z1 = z1;
    end

    if(Solar_sampling_type == 2)
		Results.Ntot = Ntot;
    end

    if(Inc_scatt(1))
		Results.h_a = h_a;
		Results.sigma_r_a = sigma_r_a;
    end

    if(Inc_scatt(2))
		Results.h_bop = h_bop;
		Results.sigma_r_bop = sigma_r_bop;
    end

    if(Inc_scatt(3))
		Results.h_c = h_c;
		Results.sigma_r_c = sigma_r_c;
    end
    
if Filename_results
    save(Filename_results,'-v7.3','-struct','Results');
end
  
end
    
