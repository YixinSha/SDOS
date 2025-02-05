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
model.param.set('kx', 'pi/a/2');
model.param.set('ky', '0');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol1').set('table', {'-a/4' '-3*a/sqrt(3)/4'; '-a/4+a' '-3*a/sqrt(3)/4'; '-a/4+a-a/2' '-3*a/sqrt(3)/4+a*sqrt(3)/2'; '-a/4-a/2' '-3*a/sqrt(3)/4+a*sqrt(3)/2'});
model.component('comp1').geom('geom1').create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('c1').set('pos', {'-a/4' '-3*a/sqrt(3)/4+a/sqrt(3)'});
model.component('comp1').geom('geom1').feature('c1').set('r', '0.2*a');
model.component('comp1').geom('geom1').create('c2', 'Circle');
model.component('comp1').geom('geom1').feature('c2').set('pos', {'-a/4+a-a/2' '-3*a/sqrt(3)/4+a*sqrt(3)/2-a/sqrt(3)'});
model.component('comp1').geom('geom1').feature('c2').set('r', '0.2*a');
model.component('comp1').geom('geom1').create('arr1', 'Array');
model.component('comp1').geom('geom1').feature('arr1').set('type', 'linear');
model.component('comp1').geom('geom1').feature('arr1').set('linearsize', 2);
model.component('comp1').geom('geom1').feature('arr1').set('displ', {'-a/2' '+a/2*sqrt(3)'});
model.component('comp1').geom('geom1').feature('arr1').selection('input').set({'c1' 'c2' 'pol1'});
model.component('comp1').geom('geom1').create('arr2', 'Array');
model.component('comp1').geom('geom1').feature('arr2').set('type', 'linear');
model.component('comp1').geom('geom1').feature('arr2').set('linearsize', 3);
model.component('comp1').geom('geom1').feature('arr2').set('displ', {'a/2' '-a/2*sqrt(3)'});
model.component('comp1').geom('geom1').feature('arr2').selection('input').set({'arr1(1,3)'});
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').physics.create('emw', 'ElectromagneticWaves', 'geom1');
model.component('comp1').physics('emw').create('wee2', 'WaveEquationElectric', 2);
model.component('comp1').physics('emw').feature('wee2').selection.set([1 2 3 4]);
model.component('comp1').physics('emw').create('pc1', 'PeriodicCondition', 1);
model.component('comp1').physics('emw').feature('pc1').selection.set([1 3 5 7 8 10 12 13]);

model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').create('cpe1', 'CopyEdge');
model.component('comp1').mesh('mesh1').create('edg2', 'Edge');
model.component('comp1').mesh('mesh1').create('cpe2', 'CopyEdge');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').create('copy1', 'Copy');
model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').mesh('mesh1').create('copy2', 'Copy');
model.component('comp1').mesh('mesh1').feature('edg1').selection.set([1 3 5 8]);
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg2').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('edg2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([2 6 8]);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('map1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('map1').selection.set([3]);

model.component('comp1').view('view1').axis.set('xmin', -1.6610170602798462);
model.component('comp1').view('view1').axis.set('xmax', 2.240598201751709);
model.component('comp1').view('view1').axis.set('ymin', -2.77394962310791);
model.component('comp1').view('view1').axis.set('ymax', 1.7224241495132446);

model.component('comp1').physics('emw').prop('ShapeProperty').set('order_electricfield', 1);
model.component('comp1').physics('emw').prop('components').set('components', 'outofplane');
model.component('comp1').physics('emw').feature('wee1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee1').set('epsilonr', [15.26; 0; 0; 0; 15.26; 0; 0; 0; 15.26]);
model.component('comp1').physics('emw').feature('wee1').set('mur_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee1').set('mur', {'0.7997'; '0.7154i'; '0'; '-0.7154i'; '0.7997'; '0'; '0'; '0'; '1'});
model.component('comp1').physics('emw').feature('wee1').set('sigma_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('mur_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('sigma_mat', 'userdef');
model.component('comp1').physics('emw').feature('pc1').set('PeriodicType', 'Floquet');
model.component('comp1').physics('emw').feature('pc1').set('kFloquet', {'kx'; 'ky'; '0'});

model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hauto', 4);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('source').set([1 3 5 8]);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('destination').set([7 10 12 13]);
model.component('comp1').mesh('mesh1').feature('edg2').feature('size1').set('hauto', 4);
model.component('comp1').mesh('mesh1').feature('cpe2').selection('source').set([2]);
model.component('comp1').mesh('mesh1').feature('cpe2').selection('destination').set([4 6 9 11]);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 7);
model.component('comp1').mesh('mesh1').feature('copy1').selection('source').set([2 6 8]);
model.component('comp1').mesh('mesh1').feature('copy1').selection('destination').set([1 5 7]);
model.component('comp1').mesh('mesh1').feature('copy2').selection('source').set([3]);
model.component('comp1').mesh('mesh1').feature('copy2').selection('destination').set([4]);
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

model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').create('surf1', 'Surface');

model.study('std1').feature('eig').set('ngen', 5);

model.sol('sol1').attach('std1');
model.sol('sol1').feature('e1').set('transform', 'eigenfrequency');
model.sol('sol1').feature('e1').set('shift', '1[GHz]');
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('e1').feature('d1').label([native2unicode(hex2dec({'5e' 'fa'}), 'unicode')  native2unicode(hex2dec({'8b' 'ae'}), 'unicode')  native2unicode(hex2dec({'76' '84'}), 'unicode')  native2unicode(hex2dec({'76' 'f4'}), 'unicode')  native2unicode(hex2dec({'63' 'a5'}), 'unicode')  native2unicode(hex2dec({'6c' '42'}), 'unicode')  native2unicode(hex2dec({'89' 'e3'}), 'unicode')  native2unicode(hex2dec({'56' '68'}), 'unicode') ' (emw)']);
model.sol('sol1').feature('e1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').runAll;

model.result('pg1').label([native2unicode(hex2dec({'75' '35'}), 'unicode')  native2unicode(hex2dec({'57' '3a'}), 'unicode') ' (emw)']);
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature('surf1').label([native2unicode(hex2dec({'88' '68'}), 'unicode')  native2unicode(hex2dec({'97' '62'}), 'unicode') ]);
model.result('pg1').feature('surf1').set('colortable', 'RainbowLight');
model.result('pg1').feature('surf1').set('smooth', 'internal');
model.result('pg1').feature('surf1').set('resolution', 'normal');

model.result('pg1').run;
% End of COMSOL script

% Extract mass and stiffness matrices
mtrInfo = mphmatrix(model,'sol1','out',{'K','E','Kc','Ec'},'initmethod','init');
mtrE = (mtrInfo.E);
mtrK = (mtrInfo.K);

% Extract mesh information
meshInfo = mphxmeshinfo(model);
no2xyz = meshInfo.nodes.coords;

end