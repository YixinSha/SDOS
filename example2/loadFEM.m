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

model.param.set('a', '1 [m]');
model.param.set('r', '2/13.3 [m]');
model.param.set('kx', '0');
model.param.set('ky', '0');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').useConstrDim(true);
model.component('comp1').geom('geom1').create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('c1').set('r', 'r');
model.component('comp1').geom('geom1').create('sq1', 'Square');
model.component('comp1').geom('geom1').feature('sq1').set('base', 'center');
model.component('comp1').geom('geom1').feature('sq1').set('size', 'a');
model.component('comp1').geom('geom1').create('arr1', 'Array');
model.component('comp1').geom('geom1').feature('arr1').set('type', 'linear');
model.component('comp1').geom('geom1').feature('arr1').set('linearsize', 3);
model.component('comp1').geom('geom1').feature('arr1').set('displ', {'0' '-a'});
model.component('comp1').geom('geom1').feature('arr1').selection('input').set({'c1' 'sq1'});
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('pos', {'0' '3*a/4'});
model.component('comp1').geom('geom1').feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').feature('r1').set('size', {'a' 'a/2'});
model.component('comp1').geom('geom1').create('arr2', 'Array');
model.component('comp1').geom('geom1').feature('arr2').set('type', 'linear');
model.component('comp1').geom('geom1').feature('arr2').set('linearsize', 3);
model.component('comp1').geom('geom1').feature('arr2').set('displ', {'a' '0'});
model.component('comp1').geom('geom1').feature('arr2').selection('input').set({'arr1' 'r1'});
model.component('comp1').geom('geom1').create('ige1', 'IgnoreEdges');
model.component('comp1').geom('geom1').feature('ige1').selection('input').set('fin(1)', [8 10 12 14 16 17 19 21 23 25 26]);
model.component('comp1').geom('geom1').run;

model.component('comp1').physics.create('emw', 'ElectromagneticWaves', 'geom1');
model.component('comp1').physics('emw').create('wee2', 'WaveEquationElectric', 2);
model.component('comp1').physics('emw').feature('wee2').selection.set([4 5 6 7 8 9 10 11 12]);
model.component('comp1').physics('emw').create('wee3', 'WaveEquationElectric', 2);
model.component('comp1').physics('emw').feature('wee3').selection.set([9]);
model.component('comp1').physics('emw').create('pc1', 'PeriodicCondition', 1);
model.component('comp1').physics('emw').feature('pc1').selection.set([1 3 5 8 9 10]);

model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').create('cpe1', 'CopyEdge');
model.component('comp1').mesh('mesh1').create('edg2', 'Edge');
model.component('comp1').mesh('mesh1').create('cpe2', 'CopyEdge');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').create('copy1', 'Copy');
model.component('comp1').mesh('mesh1').feature('edg1').selection.set([1 3 5]);
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg2').selection.set([7]);
model.component('comp1').mesh('mesh1').feature('edg2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([3 6 9 12]);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri2').selection.set([2 5 8 11]);
model.component('comp1').mesh('mesh1').feature('ftri2').create('size1', 'Size');

model.component('comp1').view('view1').axis.set('xmin', -2.3361051082611084);
model.component('comp1').view('view1').axis.set('xmax', 3.926668882369995);
model.component('comp1').view('view1').axis.set('ymin', -2.604316473007202);
model.component('comp1').view('view1').axis.set('ymax', 1.0706831216812134);

model.component('comp1').physics('emw').prop('ShapeProperty').set('order_electricfield', 1);
model.component('comp1').physics('emw').prop('components').set('components', 'outofplane');
model.component('comp1').physics('emw').prop('BackgroundField').set('ktMax', '(2*(sqrt(2*log(10))))/emw.w0');
model.component('comp1').physics('emw').feature('wee1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee1').set('mur_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee1').set('sigma_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee1').set('minput_frequency_src', 'userdef');
model.component('comp1').physics('emw').feature('wee1').set('minput_frequency', 'root.freq');
model.component('comp1').physics('emw').feature('dcont1').set('pairDisconnect', true);
model.component('comp1').physics('emw').feature('dcont1').label([native2unicode(hex2dec({'8f' 'de'}), 'unicode')  native2unicode(hex2dec({'7e' 'ed'}), 'unicode')  native2unicode(hex2dec({'60' '27'}), 'unicode') ]);
model.component('comp1').physics('emw').feature('wee2').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('epsilonr', [15; 0; 0; 0; 15; 0; 0; 0; 15]);
model.component('comp1').physics('emw').feature('wee2').set('mur_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('mur', {'0.83054 '; '0.42371i'; '0'; '-0.42371i'; '0.83054 '; '0'; '0'; '0'; '1'});
model.component('comp1').physics('emw').feature('wee2').set('sigma_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('minput_frequency_src', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('minput_frequency', 'root.freq');
model.component('comp1').physics('emw').feature('wee3').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee3').set('epsilonr', [25; 0; 0; 0; 25; 0; 0; 0; 25]);
model.component('comp1').physics('emw').feature('wee3').set('mur_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee3').set('mur', {'0.83054 '; '0.42371i'; '0'; '-0.42371i'; '0.83054 '; '0'; '0'; '0'; '1'});
model.component('comp1').physics('emw').feature('wee3').set('sigma_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee3').set('minput_frequency_src', 'userdef');
model.component('comp1').physics('emw').feature('wee3').set('minput_frequency', 'root.freq');
model.component('comp1').physics('emw').feature('pc1').set('PeriodicType', 'Floquet');
model.component('comp1').physics('emw').feature('pc1').set('kFloquet', {'kx'; 'ky'; '0'});

model.component('comp1').mesh('mesh1').feature('edg1').set('smoothmaxiter', 8);
model.component('comp1').mesh('mesh1').feature('edg1').set('smoothmaxdepth', 8);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hauto', 3);
model.component('comp1').mesh('mesh1').feature('cpe1').set('smoothmaxiter', 8);
model.component('comp1').mesh('mesh1').feature('cpe1').set('smoothmaxdepth', 8);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('source').set([1 3 5]);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('destination').set([8 9 10]);
model.component('comp1').mesh('mesh1').feature('edg2').set('smoothmaxiter', 8);
model.component('comp1').mesh('mesh1').feature('edg2').set('smoothmaxdepth', 8);
model.component('comp1').mesh('mesh1').feature('edg2').feature('size1').set('hauto', 3);
model.component('comp1').mesh('mesh1').feature('cpe2').set('smoothmaxiter', 8);
model.component('comp1').mesh('mesh1').feature('cpe2').set('smoothmaxdepth', 8);
model.component('comp1').mesh('mesh1').feature('cpe2').selection('source').set([7]);
model.component('comp1').mesh('mesh1').feature('cpe2').selection('destination').set([2 4 6]);
model.component('comp1').mesh('mesh1').feature('ftri1').set('smoothmaxiter', 8);
model.component('comp1').mesh('mesh1').feature('ftri1').set('smoothmaxdepth', 8);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 6);
model.component('comp1').mesh('mesh1').feature('ftri2').set('smoothmaxiter', 8);
model.component('comp1').mesh('mesh1').feature('ftri2').set('smoothmaxdepth', 8);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hauto', 6);
model.component('comp1').mesh('mesh1').feature('copy1').label([native2unicode(hex2dec({'52' '6f'}), 'unicode')  native2unicode(hex2dec({'67' '2c'}), 'unicode') ' 1']);
model.component('comp1').mesh('mesh1').feature('copy1').set('smoothmaxiter', 8);
model.component('comp1').mesh('mesh1').feature('copy1').set('smoothmaxdepth', 8);
model.component('comp1').mesh('mesh1').feature('copy1').selection('source').set([2 5 8 11]);
model.component('comp1').mesh('mesh1').feature('copy1').selection('destination').set([1 4 7 10]);
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
model.sol('sol1').feature('e1').create('i1', 'Iterative');
model.sol('sol1').feature('e1').feature('i1').create('sv1', 'SORVector');

model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').create('surf1', 'Surface');

model.study('std1').feature('eig').set('neigs', 1);
model.study('std1').feature('eig').set('neigsactive', true);
model.study('std1').feature('eig').set('eigunit', 'Hz');
model.study('std1').feature('eig').set('shift', '1.8066E8');
model.study('std1').feature('eig').set('ngen', 1);
model.study('std1').feature('eig').set('ngenactive', false);

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label([native2unicode(hex2dec({'7f' '16'}), 'unicode')  native2unicode(hex2dec({'8b' 'd1'}), 'unicode')  native2unicode(hex2dec({'65' 'b9'}), 'unicode')  native2unicode(hex2dec({'7a' '0b'}), 'unicode') ': ' native2unicode(hex2dec({'72' '79'}), 'unicode')  native2unicode(hex2dec({'5f' '81'}), 'unicode')  native2unicode(hex2dec({'98' '91'}), 'unicode')  native2unicode(hex2dec({'73' '87'}), 'unicode') ]);
model.sol('sol1').feature('v1').label([native2unicode(hex2dec({'56' 'e0'}), 'unicode')  native2unicode(hex2dec({'53' 'd8'}), 'unicode')  native2unicode(hex2dec({'91' 'cf'}), 'unicode') ' 1.1']);
model.sol('sol1').feature('e1').label([native2unicode(hex2dec({'72' '79'}), 'unicode')  native2unicode(hex2dec({'5f' '81'}), 'unicode')  native2unicode(hex2dec({'50' '3c'}), 'unicode')  native2unicode(hex2dec({'6c' '42'}), 'unicode')  native2unicode(hex2dec({'89' 'e3'}), 'unicode')  native2unicode(hex2dec({'56' '68'}), 'unicode') ' 1.1']);
model.sol('sol1').feature('e1').set('transform', 'eigenfrequency');
model.sol('sol1').feature('e1').set('neigs', 1);
model.sol('sol1').feature('e1').set('shift', '1.8066E8');
model.sol('sol1').feature('e1').feature('dDef').label([native2unicode(hex2dec({'76' 'f4'}), 'unicode')  native2unicode(hex2dec({'63' 'a5'}), 'unicode') ' 2']);
model.sol('sol1').feature('e1').feature('aDef').label([native2unicode(hex2dec({'9a' 'd8'}), 'unicode')  native2unicode(hex2dec({'7e' 'a7'}), 'unicode') ' 1']);
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('e1').feature('d1').active(true);
model.sol('sol1').feature('e1').feature('d1').label([native2unicode(hex2dec({'5e' 'fa'}), 'unicode')  native2unicode(hex2dec({'8b' 'ae'}), 'unicode')  native2unicode(hex2dec({'76' '84'}), 'unicode')  native2unicode(hex2dec({'76' 'f4'}), 'unicode')  native2unicode(hex2dec({'63' 'a5'}), 'unicode')  native2unicode(hex2dec({'6c' '42'}), 'unicode')  native2unicode(hex2dec({'89' 'e3'}), 'unicode')  native2unicode(hex2dec({'56' '68'}), 'unicode') ' (emw)']);
model.sol('sol1').feature('e1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('e1').feature('i1').label([native2unicode(hex2dec({'5e' 'fa'}), 'unicode')  native2unicode(hex2dec({'8b' 'ae'}), 'unicode')  native2unicode(hex2dec({'76' '84'}), 'unicode')  native2unicode(hex2dec({'8f' 'ed'}), 'unicode')  native2unicode(hex2dec({'4e' 'e3'}), 'unicode')  native2unicode(hex2dec({'6c' '42'}), 'unicode')  native2unicode(hex2dec({'89' 'e3'}), 'unicode')  native2unicode(hex2dec({'56' '68'}), 'unicode') ' (emw)']);
model.sol('sol1').feature('e1').feature('i1').set('linsolver', 'tfqmr');
model.sol('sol1').feature('e1').feature('i1').set('prefuntype', 'right');
model.sol('sol1').feature('e1').feature('i1').feature('ilDef').label([native2unicode(hex2dec({'4e' '0d'}), 'unicode')  native2unicode(hex2dec({'5b' '8c'}), 'unicode')  native2unicode(hex2dec({'51' '68'}), 'unicode') ' LU ' native2unicode(hex2dec({'52' '06'}), 'unicode')  native2unicode(hex2dec({'89' 'e3'}), 'unicode') ' 1']);
model.sol('sol1').feature('e1').feature('i1').feature('sv1').label('SOR Vector 1.1');
model.sol('sol1').feature('e1').feature('i1').feature('sv1').set('iter', 1);
model.sol('sol1').feature('e1').feature('i1').feature('sv1').set('seconditer', 3);
model.sol('sol1').runAll;

model.result('pg1').label([native2unicode(hex2dec({'75' '35'}), 'unicode')  native2unicode(hex2dec({'57' '3a'}), 'unicode') ' (emw)']);
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature('surf1').label([native2unicode(hex2dec({'88' '68'}), 'unicode')  native2unicode(hex2dec({'97' '62'}), 'unicode') ]);
model.result('pg1').feature('surf1').set('colortable', 'RainbowLightClassic');
model.result('pg1').feature('surf1').set('smooth', 'internal');
model.result('pg1').feature('surf1').set('resolution', 'normal');
% End of COMSOL script

% Extract mass and stiffness matrices
mtrInfo = mphmatrix(model,'sol1','out',{'K','E','Kc','Ec'},'initmethod','init');
mtrE = (mtrInfo.E);
mtrK = (mtrInfo.K);

% Extract mesh information
meshInfo = mphxmeshinfo(model);
no2xyz = meshInfo.nodes.coords;

end