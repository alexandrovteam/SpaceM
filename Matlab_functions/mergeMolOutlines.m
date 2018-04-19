function dummy = mergeMolOutlines(outlines_p, mol_p, save_p)
outlines = imread(outlines_p);
mol = imread(mol_p);
C = imfuse(outlines, mol, 'blend','Scaling','joint');
imwrite(C, save_p)
dummy=0;