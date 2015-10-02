function tests = UT02_LoadingModels()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testEquilibriumLoading(a)
m = LoadModelMassAction('Equilibrium.txt');
a.verifyEqual(m.Name, 'Equilibrium');
end

function testSimpleLoading(a)
m = LoadModelMassAction('Simple.txt');
a.verifyEqual(m.Name, 'Simple')
verifyDerivatives(a, m);
end

%% More comprehensive advanced model loading
function testBasicSBMLLoading(a)
opts = [];
opts.UseNames = true;
m = LoadModelSbmlAnalytic('test.xml', opts);
m = FinalizeModel(m);

a.verifyEqual(m.nv, 1);
a.verifyEqual(m.nk, 3);
a.verifyEqual(m.ns, 0);
a.verifyEqual(m.nu, 3);
a.verifyEqual(m.nx, 1);
a.verifyEqual(m.nr, 3);
a.verifyEqual(m.nz, 0);
verifyDerivatives(a, m);
end

function testEnzymeSBMLLoading(a)
m = LoadModelSbmlAnalytic('enzyme-catalysis-basic.xml');

m = AddOutput(m, 'complex', '"E:S"');
m = AddOutput(m, 'product', '"S#P"');
m = AddOutput(m, 'modified_product', '1.5*"S#P"');
m = AddOutput(m, 'random_y', '"E:S" + sqrt("S")');

a.verifyEqual(m.ny, 4);
a.verifyEqual(m.Outputs(1).Name, 'complex');
a.verifyEqual(m.Outputs(1).Expression, '"E:S"');

m = FinalizeModel(m);

% Random values
t = 10*rand;
x = rand(4,1);
u = [];
% Expected output values
y1 = x(2);
y2 = x(4);
y3 = 1.5*x(4);
y4 = x(2) + sqrt(x(3));
yExpected = [y1 y2 y3 y4]';

a.verifyEqual(m.y(t,x,u), yExpected);

verifyDerivatives(a, m);
end

function testSimpleAnalytic(a)
m = simple_analytic_model();
verifyDerivatives(a, m);
end

function testSimpleMassActionSBMLLoading(a)
warning('off', 'symbolic2massaction:repeatedSpeciesNames'); % suppress warning for repeated species C
m = LoadModelSbmlMassAction('simple_massaction.xml');
warning('on', 'symbolic2massaction:repeatedSpeciesNames'); % reenable warning for continued session

m = FinalizeModel(m);

verifyDerivatives(a, m);
end

function testSimpleMassActionAsAnalyticSBMLLoading(a)
opts = [];
opts.EvaluateExternalFunctions = true; % simple_massaction has x^2 terms, and power needs to be evaluated
m = LoadModelSbmlAnalytic('simple_massaction.xml');
m = FinalizeModel(m, opts);

verifyDerivatives(a, m);
end

function testSimulateSimpleAnalytic(a)
m = simple_analytic_model();

verifyDerivatives(a, m);
end

function testSimpleMassActionSimBiologyLoading(a)
load('simple_massaction_simbio_model.mat')

warning('off', 'symbolic2massaction:repeatedSpeciesNames'); % suppress warning for repeated species C
m = LoadModelSimBioMassAction(simbiomodel);
warning('on', 'symbolic2massaction:repeatedSpeciesNames'); % reenable warning for continued session

m = FinalizeModel(m);

verifyDerivatives(a, m);
end

function verifyDerivatives(a, m)
% Make all the values about the same size
x0 = rand(m.nx,1) + 1;
u0 = rand(m.nu,1) + 1;
k0 = rand(m.nk,1) + 1;
s0 = rand(m.ns,1) + 1;
m = m.Update(k0);

verifyClose(a, x0, @(x)m.f(0,x,u0), @(x)m.dfdx(0,x,u0))
verifyClose(a, u0, @(u)m.f(0,x0,u), @(u)m.dfdu(0,x0,u))
verifyClose(a, k0, @(k)k_wrapper(k, 'f'), @(k)k_wrapper(k, 'dfdk'))

verifyClose(a, x0, @(x)m.dfdx(0,x,u0), @(x)m.d2fdx2(0,x,u0))
verifyClose(a, u0, @(u)m.dfdu(0,x0,u), @(u)m.d2fdu2(0,x0,u))
verifyClose(a, k0, @(k)k_wrapper(k, 'dfdk'), @(k)k_wrapper(k, 'd2fdk2'))

verifyClose(a, u0, @(u)m.dfdx(0,x0,u), @(u)m.d2fdudx(0,x0,u))
verifyClose(a, x0, @(x)m.dfdu(0,x,u0), @(x)m.d2fdxdu(0,x,u0))
verifyClose(a, k0, @(k)k_wrapper(k, 'dfdx'), @(k)k_wrapper(k, 'd2fdkdx'))

verifyClose(a, x0, @(x)m.dfdk(0,x,u0), @(x)m.d2fdxdk(0,x,u0))
verifyClose(a, k0, @(k)k_wrapper(k, 'dfdu'), @(k)k_wrapper(k, 'd2fdkdu'))
verifyClose(a, u0, @(u)m.dfdk(0,x0,u), @(u)m.d2fdudk(0,x0,u))

verifyClose(a, s0, @(s)m.x0(s), @(s)m.dx0ds(s))
verifyClose(a, k0, @(k)k_wrapper_x0(k, 'x0'), @(k)k_wrapper_x0(k, 'dx0dk'))

verifyClose(a, s0, @(s)m.dx0dk(s), @(s)m.d2x0dsdk(s))
verifyClose(a, k0, @(k)k_wrapper_x0(k, 'dx0dk'), @(k)k_wrapper_x0(k, 'd2x0dk2'))

verifyClose(a, s0, @(s)m.dx0ds(s), @(s)m.d2x0ds2(s))
verifyClose(a, k0, @(k)k_wrapper_x0(k, 'dx0ds'), @(k)k_wrapper_x0(k, 'd2x0dkds'))

    function val = k_wrapper(k, member)
        m_temp = m.Update(k);
        func = m_temp.(member);
        val = func(0, x0, u0);
    end

    function val = k_wrapper_x0(k, member)
        m_temp = m.Update(k);
        func = m_temp.(member);
        val = func(s0);
    end
end

function verifyClose(a, x0, f, dfdx)
[~, dfdx_finite, dfdx_analytic] = fdiff(x0, f, dfdx);
a.verifyEqual(sparse(dfdx_finite), dfdx_analytic, 'RelTol', 0.001)
end
