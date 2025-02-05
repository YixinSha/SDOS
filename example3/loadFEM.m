% ------------------------------------------------------------------------
% Extract FEM matrices and mesh information from COMSOL 
% ------------------------------------------------------------------------
function [mtrE, mtrK, no2xyz] = loadFEM                                                                                                                             
% Returns:                                                             
%    mtrE   = mass matrix
%    mtrE   = stiffness matrix
%    no2xyz = x and y coordinates of the nodes
% Notes: 
%    Save COMSOL model as an m-file and copy it below

% Begin of COMSOL script
import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.param.set('a', '1[m]');
model.param.set('d0', '0.615*a');
model.param.set('d1', '0.615*a');
model.param.set('d2', '0.11*a');
model.param.set('kx', 'pi/a/2');
model.param.set('ky', '0');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol1').set('table', {'-a/4' '-3*a/sqrt(3)/4'; '-a/4+a' '-3*a/sqrt(3)/4'; '-a/4+a-a/2' '-3*a/sqrt(3)/4+a*sqrt(3)/2'; '-a/4-a/2' '-3*a/sqrt(3)/4+a*sqrt(3)/2'});
model.component('comp1').geom('geom1').create('pol2', 'Polygon');
model.component('comp1').geom('geom1').feature('pol2').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol2').set('table', {'(d1+2*d2)/2' '-(d1+2*d2)/2/sqrt(3)'; '-(d1+2*d2)/2' '-(d1+2*d2)/2/sqrt(3)'; '0' '(d1+2*d2)/sqrt(3)'});
model.component('comp1').geom('geom1').create('pol4', 'Polygon');
model.component('comp1').geom('geom1').feature('pol4').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol4').set('table', {'(d1+2*d2)/2' '-(d1+2*d2)/2/sqrt(3)'; '(d1+2*d2)/2-d2' '-(d1+2*d2)/2/sqrt(3)'; '(d1+2*d2)/2-d2/2' '-(d1+2*d2)/2/sqrt(3)+d2*sqrt(3)/2'});
model.component('comp1').geom('geom1').create('pol5', 'Polygon');
model.component('comp1').geom('geom1').feature('pol5').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol5').set('table', {'-(d1+2*d2)/2' '-(d1+2*d2)/2/sqrt(3)'; '-(d1+2*d2)/2+d2' '-(d1+2*d2)/2/sqrt(3)'; '-(d1+2*d2)/2+d2/2' '-(d1+2*d2)/2/sqrt(3)+d2*sqrt(3)/2'});
model.component('comp1').geom('geom1').create('pol6', 'Polygon');
model.component('comp1').geom('geom1').feature('pol6').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol6').set('table', {'0' '(d1+2*d2)/sqrt(3)'; '0-d2/2' '(d1+2*d2)/sqrt(3)-d2*sqrt(3)/2'; '0+d2/2' '(d1+2*d2)/sqrt(3)-d2*sqrt(3)/2'});
model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'pol2'});
model.component('comp1').geom('geom1').feature('dif1').selection('input2').set({'pol4' 'pol5' 'pol6'});
model.component('comp1').geom('geom1').create('rot1', 'Rotate');
model.component('comp1').geom('geom1').feature('rot1').set('keep', true);
model.component('comp1').geom('geom1').feature('rot1').setIndex('rot', '180', 0);
model.component('comp1').geom('geom1').feature('rot1').selection('input').set({'dif1'});
model.component('comp1').geom('geom1').create('mov1', 'Move');
model.component('comp1').geom('geom1').feature('mov1').set('displx', 'a/2');
model.component('comp1').geom('geom1').feature('mov1').set('disply', '-a/2*sqrt(3)');
model.component('comp1').geom('geom1').feature('mov1').selection('input').set({'rot1'});
model.component('comp1').geom('geom1').create('mov2', 'Move');
model.component('comp1').geom('geom1').feature('mov2').set('keep', true);
model.component('comp1').geom('geom1').feature('mov2').set('displx', 'a/2');
model.component('comp1').geom('geom1').feature('mov2').set('disply', '-a/2*sqrt(3)');
model.component('comp1').geom('geom1').feature('mov2').selection('input').set({'pol1'});
model.component('comp1').geom('geom1').create('pol3', 'Polygon');
model.component('comp1').geom('geom1').feature('pol3').active(false);
model.component('comp1').geom('geom1').feature('pol3').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol3').set('table', {'a/2' 'a/2/sqrt(3)';  ...
'0' 'a/sqrt(3)';  ...
'-a/2' 'a/2/sqrt(3)';  ...
'-a/2' '-a/2/sqrt(3)';  ...
'0' '-a/sqrt(3)';  ...
'a/2' '-a/2/sqrt(3)'});
model.component('comp1').geom('geom1').create('arr1', 'Array');
model.component('comp1').geom('geom1').feature('arr1').set('type', 'linear');
model.component('comp1').geom('geom1').feature('arr1').set('linearsize', 2);
model.component('comp1').geom('geom1').feature('arr1').set('displ', {'-a/2' 'a/2*sqrt(3)'});
model.component('comp1').geom('geom1').feature('arr1').selection('input').set({'dif1' 'pol1'});
model.component('comp1').geom('geom1').create('arr2', 'Array');
model.component('comp1').geom('geom1').feature('arr2').set('type', 'linear');
model.component('comp1').geom('geom1').feature('arr2').set('linearsize', 2);
model.component('comp1').geom('geom1').feature('arr2').set('displ', {'a/2' '-a/2*sqrt(3)'});
model.component('comp1').geom('geom1').feature('arr2').selection('input').set({'mov1' 'mov2'});
model.component('comp1').geom('geom1').feature('fin').set('repairtoltype', 'relative');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').physics.create('emw', 'ElectromagneticWaves', 'geom1');
model.component('comp1').physics('emw').create('wee2', 'WaveEquationElectric', 2);
model.component('comp1').physics('emw').feature('wee2').selection.set([2 4 6 8]);
model.component('comp1').physics('emw').create('pc1', 'PeriodicCondition', 1);
model.component('comp1').physics('emw').feature('pc1').selection.set([1 6 13 15 22 24 32 36]);
model.component('comp1').physics('emw').create('pmc1', 'PerfectMagneticConductor', 1);
model.component('comp1').physics('emw').feature('pmc1').selection.set([2 31]);

model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').create('cpe1', 'CopyEdge');
model.component('comp1').mesh('mesh1').create('edg2', 'Edge');
model.component('comp1').mesh('mesh1').create('cpe2', 'CopyEdge');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').create('copy1', 'Copy');
model.component('comp1').mesh('mesh1').create('edg3', 'Edge');
model.component('comp1').mesh('mesh1').create('cpe3', 'CopyEdge');
model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').create('copy2', 'Copy');
model.component('comp1').mesh('mesh1').feature('edg1').selection.set([13]);
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').selection.set([1 6 13 22]);
model.component('comp1').mesh('mesh1').feature('edg2').selection.set([14]);
model.component('comp1').mesh('mesh1').feature('edg2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([5 6]);
model.component('comp1').mesh('mesh1').feature('edg3').selection.set([6]);
model.component('comp1').mesh('mesh1').feature('edg3').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri2').selection.set([3 4]);

model.component('comp1').view('view1').axis.set('xmin', -2.8775582313537598);
model.component('comp1').view('view1').axis.set('xmax', 3.4806787967681885);
model.component('comp1').view('view1').axis.set('ymin', -4.182127475738525);
model.component('comp1').view('view1').axis.set('ymax', 2.417126178741455);

model.component('comp1').physics('emw').prop('ShapeProperty').set('order_electricfield', 1);
model.component('comp1').physics('emw').prop('components').set('components', 'outofplane');
model.component('comp1').physics('emw').feature('wee1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee1').set('mur_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee1').set('sigma_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('mur_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('mur', [13; 0; 0; 0; 13; 0; 0; 0; 13]);
model.component('comp1').physics('emw').feature('wee2').set('sigma_mat', 'userdef');
model.component('comp1').physics('emw').feature('pc1').set('PeriodicType', 'Floquet');
model.component('comp1').physics('emw').feature('pc1').set('kFloquet', {'kx'; 'ky'; '0'});

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 3);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', 0.075);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmin', 3.46E-4);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hnarrow', 2);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('source').set([13]);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('destination').set([32]);
model.component('comp1').mesh('mesh1').feature('edg2').feature('size1').set('hauto', 2);
model.component('comp1').mesh('mesh1').feature('cpe2').selection('source').set([14]);
model.component('comp1').mesh('mesh1').feature('cpe2').selection('destination').set([2 7 23 31]);
model.component('comp1').mesh('mesh1').feature('copy1').selection('source').set([5 6]);
model.component('comp1').mesh('mesh1').feature('copy1').selection('destination').set([7 8]);
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('hmax', 0.07);
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('hmin', 3.46E-4);
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('hnarrow', 2);
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('hnarrowactive', true);
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('hgrad', 1.1);
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('hgradactive', true);
model.component('comp1').mesh('mesh1').feature('cpe3').selection('source').set([6]);
model.component('comp1').mesh('mesh1').feature('cpe3').selection('destination').set([24]);
model.component('comp1').mesh('mesh1').feature('copy2').selection('source').set([3 4]);
model.component('comp1').mesh('mesh1').feature('copy2').selection('destination').set([1 2]);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('eig', 'Eigenfrequency');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('e1', 'Eigenvalue');
model.sol('sol1').feature('e1').create('d1', 'Direct');
model.sol.create('sol2');
model.sol('sol2').study('std1');
model.sol('sol2').label([native2unicode(hex2dec({'53' 'c2'}), 'unicode')  native2unicode(hex2dec({'65' '70'}), 'unicode')  native2unicode(hex2dec({'53' '16'}), 'unicode')  native2unicode(hex2dec({'89' 'e3'}), 'unicode') ' 1']);

model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup1D');
model.result.create('pg3', 'PlotGroup2D');
model.result('pg1').set('data', 'dset2');
model.result('pg1').create('surf1', 'Surface');
model.result('pg2').create('glob1', 'Global');
model.result('pg3').create('surf1', 'Surface');

model.study('std1').feature('eig').set('eigmethod', 'region');
model.study('std1').feature('eig').set('eigunit', 'Hz');
model.study('std1').feature('eig').set('appnreigs', 10);
model.study('std1').feature('eig').set('maxnreigs', 10);
model.study('std1').feature('eig').set('eiglr', '1*3*10^8');
model.study('std1').feature('eig').set('ngen', 5);

model.sol('sol1').attach('std1');
model.sol('sol1').feature('e1').set('transform', 'eigenfrequency');
model.sol('sol1').feature('e1').set('eigmethod', 'region');
model.sol('sol1').feature('e1').set('appnreigs', 10);
model.sol('sol1').feature('e1').set('maxnreigs', 10);
model.sol('sol1').feature('e1').set('eiglr', '1*3*10^8');
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('e1').feature('d1').label([native2unicode(hex2dec({'5e' 'fa'}), 'unicode')  native2unicode(hex2dec({'8b' 'ae'}), 'unicode')  native2unicode(hex2dec({'76' '84'}), 'unicode')  native2unicode(hex2dec({'76' 'f4'}), 'unicode')  native2unicode(hex2dec({'63' 'a5'}), 'unicode')  native2unicode(hex2dec({'6c' '42'}), 'unicode')  native2unicode(hex2dec({'89' 'e3'}), 'unicode')  native2unicode(hex2dec({'56' '68'}), 'unicode') ' (emw)']);
model.sol('sol1').feature('e1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').runAll;

model.result('pg1').label([native2unicode(hex2dec({'75' '35'}), 'unicode')  native2unicode(hex2dec({'57' '3a'}), 'unicode') ' (emw)']);
model.result('pg1').set('solrepresentation', 'solutioninfo');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature('surf1').label([native2unicode(hex2dec({'88' '68'}), 'unicode')  native2unicode(hex2dec({'97' '62'}), 'unicode') ]);
model.result('pg1').feature('surf1').set('colortable', 'RainbowLight');
model.result('pg1').feature('surf1').set('smooth', 'internal');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg2').set('showlegends', false);
model.result('pg2').feature('glob1').set('expr', {'freq*a/3/10^8'});
model.result('pg2').feature('glob1').set('unit', {'m/s'});
model.result('pg2').feature('glob1').set('descr', {''});
model.result('pg2').feature('glob1').set('xdata', 'expr');
model.result('pg2').feature('glob1').set('xdataexpr', 'kx');
model.result('pg2').feature('glob1').set('xdataunit', '1/m');
model.result('pg2').feature('glob1').set('xdatadescr', '');
model.result('pg2').feature('glob1').set('linestyle', 'none');
model.result('pg2').feature('glob1').set('linecolor', 'black');
model.result('pg2').feature('glob1').set('linemarker', 'point');
model.result('pg2').feature('glob1').set('markerpos', 'datapoints');
model.result('pg3').label([native2unicode(hex2dec({'75' '35'}), 'unicode')  native2unicode(hex2dec({'57' '3a'}), 'unicode') ' (emw) 1']);
model.result('pg3').set('looplevel', [5]);
model.result('pg3').set('frametype', 'spatial');
model.result('pg3').feature('surf1').label([native2unicode(hex2dec({'88' '68'}), 'unicode')  native2unicode(hex2dec({'97' '62'}), 'unicode') ]);
model.result('pg3').feature('surf1').set('colortable', 'RainbowLight');
model.result('pg3').feature('surf1').set('smooth', 'internal');
model.result('pg3').feature('surf1').set('resolution', 'normal');
% End of COMSOL script

% Extract mass and stiffness matrices
mtrInfo = mphmatrix(model,'sol1','out',{'K','E','Kc','Ec'},'initmethod','init');
mtrE = (mtrInfo.E);
mtrK = (mtrInfo.K);

% Extract mesh information
meshInfo = mphxmeshinfo(model);
no2xyz = meshInfo.nodes.coords;

end