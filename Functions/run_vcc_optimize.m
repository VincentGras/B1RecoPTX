function run_vcc_optimize(b1, varargin)

par.maxcond_read = 2; 
par.maxcond_sat = 2; 
par.cp = [];
par.Q = [];
par.FAsat_nom = 100;
par.FAread_nom = 10;
par.Nsat = size(b1, 1);
Nc = size(b1, 1);

arg=readoptg(par, varargin{:});


cp = arg.cp;
if isempty(cp)
    cp = ones(Nc, 1);
end

Q = arg.Q;
if isempty(Q)
    [Q, ~] = qr(cp / norm(cp));
end

p = VCCoptimize('maxcond_read', arg.maxcond_read, 'maxcond_sat', arg.maxcond_sat, 'b1', b1, 'FAsat_nom', arg.FAsat_nom, 'FAread_nom', arg.FAread_nom);
p.setREAD(Q); 
p.setSAT(Q(:, 1:arg.Nsat));  

p.optimize_read(); % optimize Y(0)
p.optimize_sat();  % optimize Y(SAT)=X(SAT)

fname = p.expo_txt(sprintf('Sphere_vcc%02d_Acond%d_Bcond%d', arg.Nsat, arg.maxcond_sat, arg.maxcond_read)); % save encoding schemes
