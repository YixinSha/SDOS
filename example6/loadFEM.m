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

model.param.set('a', '1 [cm]');
model.param.set('r1', '0.45 [cm]');
model.param.set('r2', '0.30 [cm]');
model.param.set('kx', '-pi/a');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.result.evaluationGroup.create('std1EvgFrq', 'EvaluationGroup');
model.result.evaluationGroup('std1EvgFrq').create('gev1', 'EvalGlobal');

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').feature('r1').set('size', {'3*a' 'sqrt(3)*a'});
model.component('comp1').geom('geom1').create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('c1').set('pos', {'-1*a' '0'});
model.component('comp1').geom('geom1').feature('c1').set('r', 'r2');
model.component('comp1').geom('geom1').create('c2', 'Circle');
model.component('comp1').geom('geom1').feature('c2').set('pos', {'1*a' '0'});
model.component('comp1').geom('geom1').feature('c2').set('r', 'r2');
model.component('comp1').geom('geom1').create('c3', 'Circle');
model.component('comp1').geom('geom1').feature('c3').set('pos', {'-0.5*a' 'sqrt(3)*a/2'});
model.component('comp1').geom('geom1').feature('c3').set('r', 'r1');
model.component('comp1').geom('geom1').create('c4', 'Circle');
model.component('comp1').geom('geom1').feature('c4').set('pos', {'0.5*a' 'sqrt(3)*a/2'});
model.component('comp1').geom('geom1').feature('c4').set('r', 'r1');
model.component('comp1').geom('geom1').create('c5', 'Circle');
model.component('comp1').geom('geom1').feature('c5').set('pos', {'0.5*a' '-sqrt(3)*a/2'});
model.component('comp1').geom('geom1').feature('c5').set('r', 'r2');
model.component('comp1').geom('geom1').create('c6', 'Circle');
model.component('comp1').geom('geom1').feature('c6').set('pos', {'-0.5*a' '-sqrt(3)*a/2'});
model.component('comp1').geom('geom1').feature('c6').set('r', 'r2');
model.component('comp1').geom('geom1').create('copy1', 'Copy');
model.component('comp1').geom('geom1').feature('copy1').set('disply', '-sqrt(3)*a');
model.component('comp1').geom('geom1').feature('copy1').selection('input').set({'r1'});
model.component('comp1').geom('geom1').create('copy2', 'Copy');
model.component('comp1').geom('geom1').feature('copy2').set('disply', 'sqrt(3)*a');
model.component('comp1').geom('geom1').feature('copy2').selection('input').set({'r1'});
model.component('comp1').geom('geom1').create('copy3', 'Copy');
model.component('comp1').geom('geom1').feature('copy3').set('disply', '-sqrt(3)*a');
model.component('comp1').geom('geom1').feature('copy3').selection('input').set({'c1' 'c2' 'c5' 'c6'});
model.component('comp1').geom('geom1').create('c7', 'Circle');
model.component('comp1').geom('geom1').feature('c7').set('pos', {'-0.5*a' 'sqrt(3)*a/2+sqrt(3)*a'});
model.component('comp1').geom('geom1').feature('c7').set('r', 'r1');
model.component('comp1').geom('geom1').create('c8', 'Circle');
model.component('comp1').geom('geom1').feature('c8').set('pos', {'0.5*a' 'sqrt(3)*a/2+sqrt(3)*a'});
model.component('comp1').geom('geom1').feature('c8').set('r', 'r1');
model.component('comp1').geom('geom1').create('c9', 'Circle');
model.component('comp1').geom('geom1').feature('c9').set('pos', {'-1*a' 'sqrt(3)*a'});
model.component('comp1').geom('geom1').feature('c9').set('r', 'r1');
model.component('comp1').geom('geom1').create('c10', 'Circle');
model.component('comp1').geom('geom1').feature('c10').set('pos', {'a' 'sqrt(3)*a'});
model.component('comp1').geom('geom1').feature('c10').set('r', 'r1');
model.component('comp1').geom('geom1').create('c11', 'Circle');
model.component('comp1').geom('geom1').feature('c11').set('pos', {'-0.5*a' 'sqrt(3)*a/2'});
model.component('comp1').geom('geom1').feature('c11').set('r', 'r1');
model.component('comp1').geom('geom1').create('c12', 'Circle');
model.component('comp1').geom('geom1').feature('c12').set('pos', {'0.5*a' 'sqrt(3)*a/2'});
model.component('comp1').geom('geom1').feature('c12').set('r', 'r1');
model.component('comp1').geom('geom1').create('c13', 'Circle');
model.component('comp1').geom('geom1').feature('c13').set('pos', {'0.5*a' '-sqrt(3)*a/2'});
model.component('comp1').geom('geom1').feature('c13').set('r', 'r2');
model.component('comp1').geom('geom1').create('c14', 'Circle');
model.component('comp1').geom('geom1').feature('c14').set('pos', {'-0.5*a' '-sqrt(3)*a/2'});
model.component('comp1').geom('geom1').feature('c14').set('r', 'r2');
model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'r1'});
model.component('comp1').geom('geom1').feature('dif1').selection('input2').set({'c1' 'c11' 'c12' 'c13' 'c14' 'c2'});
model.component('comp1').geom('geom1').create('dif2', 'Difference');
model.component('comp1').geom('geom1').feature('dif2').selection('input').set({'copy2'});
model.component('comp1').geom('geom1').feature('dif2').selection('input2').set({'c10' 'c3' 'c4' 'c7' 'c8' 'c9'});
model.component('comp1').geom('geom1').create('dif3', 'Difference');
model.component('comp1').geom('geom1').feature('dif3').selection('input').set({'copy1'});
model.component('comp1').geom('geom1').feature('dif3').selection('input2').set({'c5' 'c6' 'copy3'});
model.component('comp1').geom('geom1').create('arr1', 'Array');
model.component('comp1').geom('geom1').feature('arr1').set('type', 'linear');
model.component('comp1').geom('geom1').feature('arr1').set('linearsize', 2);
model.component('comp1').geom('geom1').feature('arr1').set('displ', {'0' '-sqrt(3)*a'});
model.component('comp1').geom('geom1').feature('arr1').selection('input').set({'dif3'});
model.component('comp1').geom('geom1').create('arr2', 'Array');
model.component('comp1').geom('geom1').feature('arr2').set('type', 'linear');
model.component('comp1').geom('geom1').feature('arr2').set('linearsize', 2);
model.component('comp1').geom('geom1').feature('arr2').set('displ', {'0' 'sqrt(3)*a'});
model.component('comp1').geom('geom1').feature('arr2').selection('input').set({'dif2'});
model.component('comp1').geom('geom1').create('pare1', 'PartitionEdges');
model.component('comp1').geom('geom1').feature('pare1').setIndex('param', '0.25', 0);
model.component('comp1').geom('geom1').feature('pare1').selection('edge').set('arr2(1)', 14);
model.component('comp1').geom('geom1').create('pare2', 'PartitionEdges');
model.component('comp1').geom('geom1').feature('pare2').setIndex('param', '0.25', 0);
model.component('comp1').geom('geom1').feature('pare2').selection('edge').set('arr2(2)', 14);
model.component('comp1').geom('geom1').create('pare3', 'PartitionEdges');
model.component('comp1').geom('geom1').feature('pare3').setIndex('param', '0.25', 0);
model.component('comp1').geom('geom1').feature('pare3').selection('edge').set('arr1(1)', 13);
model.component('comp1').geom('geom1').create('pare4', 'PartitionEdges');
model.component('comp1').geom('geom1').feature('pare4').setIndex('param', '0.25', 0);
model.component('comp1').geom('geom1').feature('pare4').selection('edge').set('arr1(2)', 13);
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').physics.create('acpr', 'PressureAcoustics', 'geom1');
model.component('comp1').physics('acpr').create('pc1', 'PeriodicCondition', 1);
model.component('comp1').physics('acpr').feature('pc1').selection.set([1 3 5 7 9 24 25 26 27 28]);

model.component('comp1').mesh('mesh1').create('edg4', 'Edge');
model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').create('cpe1', 'CopyEdge');
model.component('comp1').mesh('mesh1').create('edg2', 'Edge');
model.component('comp1').mesh('mesh1').create('cpe2', 'CopyEdge');
model.component('comp1').mesh('mesh1').create('edg3', 'Edge');
model.component('comp1').mesh('mesh1').create('cpe3', 'CopyEdge');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri3', 'FreeTri');
model.component('comp1').mesh('mesh1').create('cpd1', 'CopyDomain');
model.component('comp1').mesh('mesh1').create('cpd2', 'CopyDomain');
model.component('comp1').mesh('mesh1').feature('edg4').selection.set([7]);
model.component('comp1').mesh('mesh1').feature('edg4').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').selection.set([3 5]);
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').selection.set([1 3 5 7 9]);
model.component('comp1').mesh('mesh1').feature('edg2').selection.set([8 15 21]);
model.component('comp1').mesh('mesh1').feature('edg2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg2').feature('size1').selection.set([8 10 15 16 21 22]);
model.component('comp1').mesh('mesh1').feature('edg3').selection.set([6 14 20]);
model.component('comp1').mesh('mesh1').feature('edg3').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').selection.set([4 6 13 14 19 20]);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri2').selection.set([4]);
model.component('comp1').mesh('mesh1').feature('ftri2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri3').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri3').selection.set([3]);
model.component('comp1').mesh('mesh1').feature('ftri3').create('size1', 'Size');

model.component('comp1').view('view1').axis.set('xmin', -0.038602784276008606);
model.component('comp1').view('view1').axis.set('xmax', 0.047568291425704956);
model.component('comp1').view('view1').axis.set('ymin', -0.045163318514823914);
model.component('comp1').view('view1').axis.set('ymax', 0.046307966113090515);

model.component('comp1').physics('acpr').prop('ShapeProperty').set('order_pressure', 1);
model.component('comp1').physics('acpr').feature('fpam1').set('rho_mat', 'userdef');
model.component('comp1').physics('acpr').feature('fpam1').set('rho', 1.25);
model.component('comp1').physics('acpr').feature('fpam1').set('c_mat', 'userdef');
model.component('comp1').physics('acpr').feature('fpam1').set('c', 343);
model.component('comp1').physics('acpr').feature('pc1').set('PeriodicType', 'Floquet');
model.component('comp1').physics('acpr').feature('pc1').set('kFloquet', {'kx'; '0'; '0'});

model.component('comp1').mesh('mesh1').feature('edg4').label([native2unicode(hex2dec({'8f' 'b9'}), 'unicode') ' 0']);
model.component('comp1').mesh('mesh1').feature('edg4').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg4').feature('size1').set('hmax', 'a/10');
model.component('comp1').mesh('mesh1').feature('edg4').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hauto', 2);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', 'a/5');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('cpe1').label([native2unicode(hex2dec({'59' '0d'}), 'unicode')  native2unicode(hex2dec({'52' '36'}), 'unicode')  native2unicode(hex2dec({'8f' 'b9'}), 'unicode') ' 01']);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('source').set([3 5 7]);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('destination').set([25 26 27]);
model.component('comp1').mesh('mesh1').feature('edg2').feature('size1').set('hauto', 2);
model.component('comp1').mesh('mesh1').feature('edg2').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg2').feature('size1').set('hmax', 'a/6');
model.component('comp1').mesh('mesh1').feature('edg2').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('cpe2').selection('source').set([8 15 21]);
model.component('comp1').mesh('mesh1').feature('cpe2').selection('destination').set([10 16 22]);
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('hauto', 3);
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('hmax', 'a/5');
model.component('comp1').mesh('mesh1').feature('edg3').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('cpe3').selection('source').set([6 14 20]);
model.component('comp1').mesh('mesh1').feature('cpe3').selection('destination').set([4 13 19]);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 7);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hauto', 7);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hauto', 7);
model.component('comp1').mesh('mesh1').feature('cpd1').selection('source').set([4]);
model.component('comp1').mesh('mesh1').feature('cpd1').selection('destination').set([5]);
model.component('comp1').mesh('mesh1').feature('cpd2').selection('source').set([2]);
model.component('comp1').mesh('mesh1').feature('cpd2').selection('destination').set([1]);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('eig', 'Eigenfrequency');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('e1', 'Eigenvalue');

model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup2D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg2').create('surf1', 'Surface');
model.result('pg2').feature('surf1').set('expr', 'acpr.Lp');

model.study('std1').feature('eig').set('ngen', 5);

model.sol('sol1').attach('std1');
model.sol('sol1').feature('e1').set('transform', 'eigenfrequency');
model.sol('sol1').feature('e1').set('shift', '100[Hz]');
model.sol('sol1').feature('e1').set('eigref', '100');
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').runAll;

model.result.evaluationGroup('std1EvgFrq').label([native2unicode(hex2dec({'72' '79'}), 'unicode')  native2unicode(hex2dec({'5f' '81'}), 'unicode')  native2unicode(hex2dec({'98' '91'}), 'unicode')  native2unicode(hex2dec({'73' '87'}), 'unicode') ' (' native2unicode(hex2dec({'78' '14'}), 'unicode')  native2unicode(hex2dec({'7a' '76'}), 'unicode') ' 1)']);
model.result.evaluationGroup('std1EvgFrq').set('data', 'dset1');
model.result.evaluationGroup('std1EvgFrq').set('looplevelinput', {'all'});
model.result.evaluationGroup('std1EvgFrq').feature('gev1').set('expr', {'freq*2*pi' 'imag(freq)/abs(freq)' 'abs(freq)/imag(freq)/2'});
model.result.evaluationGroup('std1EvgFrq').feature('gev1').set('unit', {'rad/s' '1' '1'});
model.result.evaluationGroup('std1EvgFrq').feature('gev1').set('descr', {[native2unicode(hex2dec({'89' 'd2'}), 'unicode')  native2unicode(hex2dec({'98' '91'}), 'unicode')  native2unicode(hex2dec({'73' '87'}), 'unicode') ] [native2unicode(hex2dec({'96' '3b'}), 'unicode')  native2unicode(hex2dec({'5c' '3c'}), 'unicode')  native2unicode(hex2dec({'6b' 'd4'}), 'unicode') ] [native2unicode(hex2dec({'54' 'c1'}), 'unicode')  native2unicode(hex2dec({'8d' '28'}), 'unicode')  native2unicode(hex2dec({'56' 'e0'}), 'unicode')  native2unicode(hex2dec({'5b' '50'}), 'unicode') ]});
model.result.evaluationGroup('std1EvgFrq').run;
model.result('pg1').label([native2unicode(hex2dec({'58' 'f0'}), 'unicode')  native2unicode(hex2dec({'53' '8b'}), 'unicode') ' (acpr)']);
model.result('pg1').feature('surf1').set('colortable', 'Wave');
model.result('pg1').feature('surf1').set('colortablesym', true);
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg2').label([native2unicode(hex2dec({'58' 'f0'}), 'unicode')  native2unicode(hex2dec({'53' '8b'}), 'unicode')  native2unicode(hex2dec({'7e' 'a7'}), 'unicode') ' (acpr)']);
model.result('pg2').feature('surf1').set('resolution', 'normal');
% End of COMSOL script

% Extract mass and stiffness matrices
mtrInfo = mphmatrix(model,'sol1','out',{'K','E','Kc','Ec'},'initmethod','init');
mtrE = (mtrInfo.E);
mtrK = (mtrInfo.K);

% Extract mesh information
meshInfo = mphxmeshinfo(model);
no2xyz = meshInfo.nodes.coords;

end