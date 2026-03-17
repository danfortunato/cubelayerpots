function test_reordering

% Test cooordinate reordering

ncheb = 5;
nref = 1;
nside = 2^nref;

% Get the cube surface nodes used by the layer potential operators
[x_pot, y_pot, z_pot] = cubepts(ncheb, nside);
xyz_pot = [x_pot(:) y_pot(:) z_pot(:)];

% Get the cube surface nodes from Hari's code
rr = refel(3, ncheb+1);
root = block(rr, [-0.5 -0.5 -0.5], 1);
root.split(nref);
[x_hari, y_hari, z_hari] = root.getCoords();
xyz_hari = [x_hari(:) y_hari(:) z_hari(:)];

% Reorder each
[HfC, CfH] = cube2hari(ncheb, nside);
xyz_hari_comp = HfC*xyz_pot;
xyz_pot_comp  = CfH*xyz_hari;

err_hfc = norm(xyz_hari_comp(:) - xyz_hari(:), inf);
err_cfh = norm(xyz_pot_comp(:)  - xyz_pot(:),  inf);

fprintf('Error in coordinates (HfC) = %g\n', err_hfc);
fprintf('Error in coordinates (CfH) = %g\n', err_cfh);

end
