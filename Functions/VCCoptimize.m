%%%%%%% Encoding Matrices %%%%%%%%
%%%%%%%  saturat modes: X = [0 XSAT] := [O SAT_] %%%%%%%%%
%%%%%%%  readout modes: Y = [Y0 YSAT] := [READ_ SAT_] %%%%%%%%%


classdef VCCoptimize < handle
    
    properties
        
        READ_ = eye(8)/sqrt(8);
        SAT_ = eye(8)/sqrt(8); 
        Vmax = 165; % V
        b1 = []; % rad/s/V, Nc x Nvoxel matrix (unit uT/Volt) rad/V/s

        maxcond_read = 2; 
        maxcond_sat = 2; 

        FAsat_nom = 110;
        FAread_nom = 10;

    end
    
    methods
        
        function p = VCCoptimize(varargin)
        
            if (nargin  == 0)
                return;
            end            
            
            
            for j = 1:2:numel(varargin)      
                par = varargin{j};
                value = varargin{j+1};           
                p.(par) =value;
            end   
        
        end


        function q = new(p)
            q = VCCoptimize_ND(p);          
        end
            
      
        function set_b1(p, b1)
            assert (size(b1, 1) == p.getNc, 'bad shape')
            p.b1 = b1;
        end


        function setREAD(p, val)
            p.READ_ = val;
        end
        
        function setSAT(p, val)
            p.SAT_ = val;
        end

        function set_maxcond_read(p, val)
            p.maxcond_read = val;
        end

        function set_maxcond_sat(p, val)
            p.maxcond_sat = val;
        end

        function set_FAsat_nom(p, val)
            p.FAsat_nom = val;
        end

        function set_FAread_nom(p, val)
            p.FAread_nom = val;
        end


        function val=getI(p)
            val = 1/pi * p.Vmax *1e-3;
        end

        function val = getNc(p) % number of transmit channels
            val = size(p.SAT_,1);
         end
         
         function val = getNsat(p)
            val = size(p.SAT_,2);
        end

        function val = getREAD(p)
            val = p.READ_;
        end
        
        function val = getSAT(p)
            val = p.SAT_;
        end


        %%%%%%%%%%%%%%%%% optimization part %%%%%%%%%%%%%%%%%%
       function val = profiles_read(p, i)
            
            if nargin < 2 
                i = [];
            end
            y = p.getREAD;
            if ~isempty(i)
                y = y(:, i);
            end
            val = y' * p.b1;
            
        end
        

        function val = profiles_sat(p, i)
            
            if nargin < 2
                i = [];
            end
            x = p.getSAT;
            if ~isempty(i)
                x = x(:, i);
            end
            val = x' * p.b1;
            
        end


        function y = optimvec_read(p)
            yc = reshape(p.getREAD(), p.getNc * p.getNc, 1);
            y = cat(1, real(yc), imag(yc));
            
        end
        
        function x = optimvec_sat(p)
            
            xc = reshape(p.getSAT(), p.getNc * p.getNsat, 1);
            x = cat(1, real(xc), imag(xc));
            
        end
        
        
        function p = update_with_optimvec_read(p, y)
            n = p.getNc * p.getNc;
            yc = y(1:n) + 1j * y(1+n:2*n);
            yc = reshape(yc, p.getNc, p.getNc);
            yc = yc ./ vecnorm(yc, 2, 1);
            p.READ_(:, :) = yc;
            
        end
        
        function p = update_with_optimvec_sat(p, x)
            n = p.getNc * p.getNsat;
            xc = x(1:n) + 1j * x(1+n:2*n);
            xc = reshape(xc, p.getNc, p.getNsat);
            xc = xc ./ vecnorm(xc, 2, 1);
            p.SAT_(:, :) = xc;
        end


        
        function C = objfun_read(p) 
            C = local_metric(abs(p.profiles_read()));
        end

        
        function C = objfun_sat(p) 
            C = local_metric(abs(p.profiles_sat()));
        end


        function [C, Ceq] = constrfun_cond_read(p)
            C = cond(p.getREAD()) - p.maxcond_read;
            Ceq = [];
        end

        function [C, Ceq] = constrfun_cond_sat(p)
            C = cond(p.getSAT()) - p.maxcond_sat;
            Ceq = [];
        end
        
        
        
        function p = optimize_read(p)
            
            opt = optimoptions(@fmincon);
            opt.Display = 'iter';
            opt.Algorithm = 'sqp';
            opt.MaxFunctionEvaluations = 1000000;
            opt.MaxIter = 1000;
                        
            if p.getNc <= 0
                return;
            end
            
            y0 = p.optimvec_read();
            
            fun = @(y) p.update_with_optimvec_read(y).objfun_read();
            nlcon = [];
            if p.maxcond_read < inf
                nlcon = @(y) p.update_with_optimvec_read(y).constrfun_cond_read();
            end
            y = fmincon(fun, y0, [], [], [], [], [], [], nlcon, opt);
            p.update_with_optimvec_read(y);
            
        end

        
        function p = optimize_sat(p)
            
            opt = optimoptions(@fmincon);
            opt.Display = 'iter';
            opt.Algorithm = 'sqp';
            opt.MaxFunctionEvaluations = 1000000;
            opt.MaxIter = 1000;
           
            if ~(p.getNsat > 0 && p.getNc > 0)
                return;
            end
            
            x0 = p.optimvec_sat();

            fun = @(x) p.update_with_optimvec_sat(x).objfun_sat();
            nlcon = [];
            if p.maxcond_sat < inf
                nlcon = @(x) p.update_with_optimvec_sat(x).constrfun_cond_sat();
            end
            x = fmincon(fun, x0, [], [], [], [], [], [], nlcon, opt);
            p.update_with_optimvec_sat(x);
            
        end 


        %%%%%%%%%%%%%%%%% writting part %%%%%%%%%%%%%%%%%%
       

        function dist_FA_norm=get_FAdist_norm(p, MatFA) % unitless transmit profile of encoding matrices
            dist_FA_norm=abs(MatFA'*p.b1)*p.getI();
        end


        function MatAlpha = getMSAT_actuel(p) % normalized sat modes + zeros for non-sat cycles
            sat=p.getSAT();
            scaleA=max(abs(sat(:))); 
            MatAlpha=p.getSAT()/scaleA; 
            MatAlpha = cat(2, zeros(size(p.READ_)), MatAlpha); 
        end

        
        function MatBeta = getMREAD_actuel(p) % rescale to the same average value for two matrices

            read=p.getREAD();
            sat=p.getSAT();

            scaleB1=max(abs(read(:))); MatBeta1 = p.getREAD()/scaleB1;
            scaleB2=max(abs(sat(:))); MatBeta2 = p.getSAT()/scaleB2;                 
    
            dist_beta1_norm=p.get_FAdist_norm(MatBeta1);
            mean_beta1=mean(dist_beta1_norm(:));
    
            dist_beta2_norm=p.get_FAdist_norm(MatBeta2);
            mean_beta2=mean(dist_beta2_norm(:));
    
            %rescale one scheme to another
            if mean_beta1>mean_beta2
                MatBeta1=MatBeta1/(mean_beta1/mean_beta2);
            else
                MatBeta2=MatBeta2/(mean_beta2/mean_beta1);
            end

            MatBeta=cat(2, MatBeta1, MatBeta2);

        end


        function M = Alpha(p) % изменить имя
            M=p.get_FAdist_norm(p.getMSAT_actuel())*p.FAsat_nom; % V*s
        end


        function M = Beta(p)
            M=p.get_FAdist_norm(p.getMREAD_actuel())*p.FAread_nom;
        end
 

        function out = expo_txt(p, fname)
            
            if nargin < 2
                fname = 'XFLparameters.txt';
            end


            Alpha=p.Alpha();
            Beta=p.Beta();
            
            s=newline;
            s = [s, sprintf('# Nsat = %d\n', p.getNsat)];
            s = [s, sprintf('# Cond_read = %.3f\n', cond(p.READ_))];
            s = [s, sprintf('# Cond_sat = %.3f\n', cond(p.SAT_))];


            s = [s, sprintf(strcat("# mean alpha for each mode = ", repmat('%.3f ', 1, size(Alpha,1))), mean(Alpha,2))];
            s = [s, sprintf("# mean alpha averaged over modes =  %.3f", mean(Alpha(any(Alpha,2),:),'all'))];
            s = [s, sprintf("# min max alpha over modes =  %.3f", min(max(Alpha,[],1)))];

            s = [s, sprintf(strcat("# mean beta angles for each mode = ", repmat('%.3f ', 1, size(Beta,1))), mean(Beta,2))];
            s = [s, sprintf("# mean beta angles averaged over modes =  %.3f", mean(Beta(:)))];
            s = [s, sprintf("# min max beta over modes =  %.3f", min(max(Beta,[],1)))];

            

            XFLexpoxflparam(fname, p.getMSAT_actuel(), p.getMREAD_actuel(), p.Vmax, p.FAsat_nom, p.FAread_nom);


            if (~isempty(fname))
   
                if (isempty(regexp(fname, '\.xflparam.txt$', 'once')))
                    filename = [fname, '.xflparam.txt'];
                end

                fid = fopen(filename, 'at');
                
                if fid == -1
                    error('Failed to open file.');
                end
                
                fprintf(fid, '%s\n', s);
                
                fclose(fid);

            end
            
            out = fname;
        end
        

    end  

end


function C = local_metric(profile)

    m = size(profile, 1);
    n = size(profile, 2);
    P = profile ./ mean(profile, 2);
    C = norm(vecnorm(P - 1, 2, 1)) * m^(-1/2) * n^(-1/2);

end


