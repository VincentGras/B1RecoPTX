function [s,fname] = XFLexpoxflparam(fname, sat, fl, Vref, alpha_nominal_deg, beta_nominal_deg)

    
s = '';
if nargin > 3 && Vref > 0
    s = [s, sprintf('#Vref = %.6f V\n', Vref)];
end

if nargin > 4 && alpha_nominal_deg > 0
    s = [s, sprintf('#alpha_nominal = %.1f deg\n', alpha_nominal_deg)];
end

if nargin > 5 && beta_nominal_deg > 0
    s = [s, sprintf('#beta_nominal = %.1f deg\n', beta_nominal_deg)];
end

% go through cycles with no saturation
for i = 1:size(sat, 2)
    
	if (all(sat(:,i)==0))
		s = [s, sprintf('RELATIVE_RFEXC %s\n', num2str(get_magn_phase(fl(:,i)), '%10.3f'))];
	end

end

% then go through the cycles with a saturation pulse
for i = 1:size(sat, 2)
    
	if (any(sat(:,i)~=0,1))
		s = [s,sprintf('SATURATE_RFSAT %s\n', num2str(get_magn_phase(sat(:,i)), '%10.3f')),...
			sprintf('SATURATE_RFEXC %s\n', num2str(get_magn_phase(fl(:,i)), '%10.3f'))];
	end
end

if (~isempty(fname))
   
   if (isempty(regexp(fname, '\.xflparam.txt$', 'once')))
       fname = [fname, '.xflparam.txt'];
   end
   f = fopen(fname, 'w');
   
   fprintf(f, ['#XFLparameters\n', s]);
   fclose(f);
else
   fprintf(['#XFLparameters\n', s]);
end


function mp = get_magn_phase(z)

    mp = zeros(1, numel(z)*2);
    mp(1:2:end) = abs(z);
    mp(2:2:end) = mod(angle(z), 2*pi)*180/pi;
    