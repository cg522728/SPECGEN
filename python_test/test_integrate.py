import unittest
from random import Random
import integrate as an
import tube_atten as ta
import xraylib as xrl
from math import exp

def square( arg1):
    total = arg1**2
    return total

def integrate(xmin, xmax, nstep):
    xstep = (xmax-xmin)/nstep
    ymin = square(xmin)/2
    ymax = square(xmax)/2
    sum = float(0)
    for x in range(1,nstep):
        xval = xmin + x*xstep
        yval = square(xval)
        sum = sum +yval
    out = xstep*(ymin+sum+ymax)
    return out

def tube_atten(E, Z, D):
    expon = -xrl.CS_Total(Z, E)*D*1e-4*xrl.ElementDensity(Z)
    ta = exp(expon)
    return ta
        
def integrate_build_testcase(n):
    XDOMAIN = (0.0,100.0)
    NDOMAIN = (1,1000)
    def random_sample(n):
        random = Random(("integrate",n))
        for i in xrange(n):
            x1 = random.uniform(*XDOMAIN)
            x2 = random.uniform(*XDOMAIN)
            nstep = random.randint(*NDOMAIN)
            if x1<x2:
                xmin = x1
                xmax = x2
            if x2<x1:
                xmin = x2
                xmax = x1
            yield xmin, xmax, nstep
    def make_test_method(test_point):
        def test_integrate(self):
            result_i = integrate(*test_point)
            result_f = an.mathchg.integrate(square, *test_point)
            self.assertAlmostEqual(result_f,result_i)
        return test_integrate
    count = 0
    dico = {}
    for pt in random_sample(n):
        testname = "test_func_to_test_%i" % count
        dico[testname] = make_test_method(pt)
        count += 1
    return type("IntegrateTestCase",(unittest.TestCase,),dico)

def tube_atten_build_testcase(n):
    DDOMAIN = (1,3000)
    EDOMAIN = (0,150)
    def random_sample(n):
        random = Random(("tube_atten",n))
        for i in xrange(n):
            Z = random.randint(1,92)
            D = random.uniform(*DDOMAIN)
            E = random.uniform(*EDOMAIN) 
            yield E, Z, D
    def make_test_method(test_point):
        def test_tube_atten(self):
            result_i = tube_atten(*test_point)
            result_f = ta.anode.tube_atten(*test_point)
            self.assertAlmostEqual(result_f,result_i,places=3)
        return test_tube_atten
    count = 0
    dico = {}
    for pt in random_sample(n):
        testname = "test_tube_atten_%i" % count
        dico[testname] = make_test_method(pt)
        count += 1
    return type("TubeAttenTestCase",(unittest.TestCase,),dico)

IntegrateTestCase = integrate_build_testcase(1024)
TubeAttenTestCase = tube_atten_build_testcase(1024)

if __name__ == '__main__':
    unittest.main()
