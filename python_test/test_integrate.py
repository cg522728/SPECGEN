import unittest
import integrate as mintegrate
import random
XDOMAIN = (0.0,100.0)
NDOMAIN = (1,1000)

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

def random_sample(n):
    for i in xrange(n):
        x1 = random.uniform(*XDOMAIN)
        x2 = random.uniform(*XDOMAIN)
        nstep = int(random.uniform(*NDOMAIN))
        if x1<x2:
            xmin = x1
            xmax = x2
        if x2<x1:
            xmin = x2
            xmax = x1
        yield xmin, xmax, nstep

def build_testcase(n):
    def make_test_method(test_point):
        def a_test(self):
            result_i = integrate(*test_point)
            result_f = mintegrate.mathchg.integrate(square, *test_point)
            self.assertAlmostEqual(result_f,result_i)
        return a_test
    count = 0
    dico = {}
    for pt in random_sample(n):
        testname = "test_func_to_test_%i" % count
        dico[testname] = make_test_method(pt)
        count += 1
    return type("TestFuncTestCase",(unittest.TestCase,),dico)

TestFuncTestCase = build_testcase(1024)

if __name__ == '__main__':
    unittest.main()