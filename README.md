# Surface Equivalence 


## Introduction

Surface-Equivalence is a library for computing projective isomorphisms between surfaces, if such isomorphisms exists.

This library depends on [SageMath](https://SageMath.org) libraries and uses [Maple](https://www.maplesoft.com) for Groebner basis computations.

## Installation

* Install Sage from [SageMath](https://SageMath.org).
We assume that `sage` is accessible from your commandline interface.

* Install [Maple](https://www.maplesoft.com).
We assume for some examples that `maple` is accessible from your commandline interface.

* Install the `surface_equivalence` package: 
```    
sage -pip install surface_equivalence
```    
If you do not have root access use the following command instead:
```    
sage -pip install --user surface_equivalence
```    

* We advice to upgrade the `surface_equivalence` package regularly:
```
sage -pip install --upgrade surface_equivalence
```
 
* To execute some [usecases](https://github.com/niels-lubbes/moebius_aut/blob/master/surface_equivalence/src/surface_equivalence/__main__.py) type:
```    
sage -python -m surface_equivalence
```

* For showing which files were installed 
or for uninstalling the `surface_equivalence` package, 
use one of the following commands:
```
sage -pip show --files surface_equivalence
sage -pip uninstall surface_equivalence
```


## Examples

For running the examples below, either copy paste the code into the Sage interface or run them as a Python module:

    sage -python -m my_module_name.py

See [this file](https://github.com/niels-lubbes/surface_equivalence/blob/master/surface_equivalence/src/surface_equivalence/__main__.py) 
for more example usecases. 
See the [source code](https://github.com/niels-lubbes/surface_equivalence/blob/master/surface_equivalence/src/surface_equivalence)
the io-specification of each function.
The [test functions](https://github.com/niels-lubbes/surface_equivalence/blob/master/surface_equivalence/src/tests)
might be informative for how to call each function.


### Example 1: Projective automorphisms of Roman surfaces

We start by importing the required libraries and initialize parameters.

```python
from surface_equivalence.class_se_ring import ring
from surface_equivalence.class_se_ring import SERing
from surface_equivalence.sage_interface import sage_gcd
from surface_equivalence.sage_interface import sage_matrix
from surface_equivalence.sage_interface import sage_vector
from surface_equivalence.sage_interface import sage_SR
from surface_equivalence.sage_interface import sage_solve
from surface_equivalence.sage_interface import sage_ideal
from surface_equivalence.sage_interface import sage_QQ
from surface_equivalence.sage_interface import sage_identity_matrix

y = ring( 'y0,y1,y2,y3' )
x = ring( 'x0,x1,x2' )
c = ring( 'c0,c1,c2,c3,c4,c5,c6,c7,c8' )
```

We denote the parametrizations of the two [Roman surfaces](https://en.wikipedia.org/wiki/Roman_surface) by `f` and `g`.

```python
f = ring( '[x0^2+x1^2+x2^2,x0*x1,x0*x2,x1*x2]' )
g = f
```

Our goal is to recover matrices `U` that define the projective automorphisms of the Roman surface.
For this purpose we first compute coefficient matrices.

```python
# compatible reparametrizations are linear and indexed by c
r = {x[0]:c[0] * y[0] + c[1] * y[1] + c[2] * y[2],
     x[1]:c[3] * y[0] + c[4] * y[1] + c[5] * y[2],
     x[2]:c[6] * y[0] + c[7] * y[1] + c[8] * y[2]}

# compute kernel and coefficient matrix of f
Mf = SERing.get_matrix_P2( f )
Kf = Mf.right_kernel_matrix().T
assert ( Mf * Kf ).is_zero()

# compute the coefficient matrix of g composed with r
gr = [ comp.subs( r ) for comp in g ]
assert sage_gcd( gr ) == 1
assert SERing.get_degree( gr, 'y0,y1,y2' ) == SERing.get_degree( f )
Mgr = SERing.get_matrix_P2( gr, 'y0,y1,y2' )

# output
print( 'f    =', f )
print( 'g    =', g )
print( 'r    =', r )
print( 'Mf   =', Mf.dimensions(), list( Mf ), SERing.get_mon_P2( 2 ) )
print( 'Kf.T =', Kf.T.dimensions(), list( Kf.T ) )
print( 'Mgr  =', Mgr.dimensions(), list( Mgr ), '\n' + str( Mgr ) )
```

Output:
```
f    = [x0^2 + x1^2 + x2^2, x0*x1, x0*x2, x1*x2] 
g    = [x0^2 + x1^2 + x2^2, x0*x1, x0*x2, x1*x2] 
r    = {x2: c6*y0 + c7*y1 + c8*y2, x1: c3*y0 + c4*y1 + c5*y2, x0: c0*y0 + c1*y1 + c2*y2} 
Mf   = (4, 6) [(1, 0, 0, 1, 0, 1), (0, 1, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0), (0, 0, 0, 0, 1, 0)] [x0^2, x0*x1, x0*x2, x1^2, x1*x2, x2^2] 
Kf.T = (2, 6) [(1, 0, 0, 0, 0, -1), (0, 0, 0, 1, 0, -1)] 
Mgr  = (4, 6) [(c0^2 + c3^2 + c6^2, 2*c0*c1 + 2*c3*c4 + 2*c6*c7, 2*c0*c2 + 2*c3*c5 + 2*c6*c8, c1^2 + c4^2 + c7^2, 2*c1*c2 + 2*c4*c5 + 2*c7*c8, c2^2 + c5^2 + c8^2), (c0*c3, c1*c3 + c0*c4, c2*c3 + c0*c5, c1*c4, c2*c4 + c1*c5, c2*c5), (c0*c6, c1*c6 + c0*c7, c2*c6 + c0*c8, c1*c7, c2*c7 + c1*c8, c2*c8), (c3*c6, c4*c6 + c3*c7, c5*c6 + c3*c8, c4*c7, c5*c7 + c4*c8, c5*c8)] 
[         c0^2 + c3^2 + c6^2 2*c0*c1 + 2*c3*c4 + 2*c6*c7 2*c0*c2 + 2*c3*c5 + 2*c6*c8          c1^2 + c4^2 + c7^2 2*c1*c2 + 2*c4*c5 + 2*c7*c8          c2^2 + c5^2 + c8^2]
[                      c0*c3               c1*c3 + c0*c4               c2*c3 + c0*c5                       c1*c4               c2*c4 + c1*c5                       c2*c5]
[                      c0*c6               c1*c6 + c0*c7               c2*c6 + c0*c8                       c1*c7               c2*c7 + c1*c8                       c2*c8]
[                      c3*c6               c4*c6 + c3*c7               c5*c6 + c3*c8                       c4*c7               c5*c7 + c4*c8                       c5*c8] 
 ```   

We compute all solutions for `c` such that `Mgr*Kf==0`. This takes a few seconds.

```python
ec_lst = ( Mgr * Kf ).list() + [  sage_matrix( SERing.R, 3, 3, c ).det() * ring( 't' ) - 1 ]
pc_lst = sage_ideal( ec_lst ).elimination_ideal( ring( 't' ) ).primary_decomposition()
sol_lst = []
for pc in pc_lst:
    s_lst = list( reversed( sorted( pc.gens() ) ) )
    s_dct = ring( sage_solve( [sage_SR( comp ) for comp in s_lst], [sage_SR( comp ) for comp in c], solution_dict = True )[0] )
    sol_lst += [s_dct]
    print( s_lst, '-->', s_dct )
```

Output:
```
[c0, c1 - c6, c2, c3, c4, c5 - c6, c7, c8] --> {c3: 0, c8: 0, c4: 0, c5: r1, c0: 0, c1: r1, c6: r1, c7: 0, c2: 0} 
[c0, c1 + c6, c2, c3, c4, c5 - c6, c7, c8] --> {c3: 0, c8: 0, c4: 0, c5: r2, c0: 0, c1: -r2, c6: r2, c7: 0, c2: 0} 
[c0, c1 - c6, c2, c3, c4, c5 + c6, c7, c8] --> {c3: 0, c8: 0, c4: 0, c5: r3, c0: 0, c1: -r3, c6: -r3, c7: 0, c2: 0} 
[c0, c1 + c6, c2, c3, c4, c5 + c6, c7, c8] --> {c3: 0, c8: 0, c4: 0, c5: r4, c0: 0, c1: r4, c6: -r4, c7: 0, c2: 0} 
[c0 - c7, c1, c2, c3, c4, c5 - c7, c6, c8] --> {c3: 0, c8: 0, c4: 0, c5: r5, c0: r5, c1: 0, c6: 0, c7: r5, c2: 0} 
[c0 + c7, c1, c2, c3, c4, c5 - c7, c6, c8] --> {c3: 0, c8: 0, c4: 0, c5: r6, c0: -r6, c1: 0, c6: 0, c7: r6, c2: 0} 
[c0 + c7, c1, c2, c3, c4, c5 + c7, c6, c8] --> {c3: 0, c8: 0, c4: 0, c5: r7, c0: r7, c1: 0, c6: 0, c7: -r7, c2: 0} 
[c0 - c7, c1, c2, c3, c4, c5 + c7, c6, c8] --> {c3: 0, c8: 0, c4: 0, c5: r8, c0: -r8, c1: 0, c6: 0, c7: -r8, c2: 0} 
[c0 - c8, c1, c2, c3, c4 - c8, c5, c6, c7] --> {c3: 0, c8: r9, c4: r9, c5: 0, c0: r9, c1: 0, c6: 0, c7: 0, c2: 0} 
[c0 + c8, c1, c2, c3, c4 - c8, c5, c6, c7] --> {c3: 0, c8: r10, c4: r10, c5: 0, c0: -r10, c1: 0, c6: 0, c7: 0, c2: 0} 
[c0 + c8, c1, c2, c3, c4 + c8, c5, c6, c7] --> {c3: 0, c8: r11, c4: -r11, c5: 0, c0: -r11, c1: 0, c6: 0, c7: 0, c2: 0} 
[c0 - c8, c1, c2, c3, c4 + c8, c5, c6, c7] --> {c3: 0, c8: r12, c4: -r12, c5: 0, c0: r12, c1: 0, c6: 0, c7: 0, c2: 0} 
[c0, c1, c2 - c6, c3, c4 - c6, c5, c7, c8] --> {c3: 0, c8: 0, c4: r13, c5: 0, c0: 0, c1: 0, c6: r13, c7: 0, c2: r13} 
[c0, c1, c2 + c6, c3, c4 - c6, c5, c7, c8] --> {c3: 0, c8: 0, c4: -r14, c5: 0, c0: 0, c1: 0, c6: -r14, c7: 0, c2: r14} 
[c0, c1 - c8, c2, c3 - c8, c4, c5, c6, c7] --> {c3: r15, c8: r15, c4: 0, c5: 0, c0: 0, c1: r15, c6: 0, c7: 0, c2: 0} 
[c0, c1 - c8, c2, c3 + c8, c4, c5, c6, c7] --> {c3: -r16, c8: r16, c4: 0, c5: 0, c0: 0, c1: r16, c6: 0, c7: 0, c2: 0} 
[c0, c1 + c8, c2, c3 - c8, c4, c5, c6, c7] --> {c3: r17, c8: r17, c4: 0, c5: 0, c0: 0, c1: -r17, c6: 0, c7: 0, c2: 0} 
[c0, c1 + c8, c2, c3 + c8, c4, c5, c6, c7] --> {c3: -r18, c8: r18, c4: 0, c5: 0, c0: 0, c1: -r18, c6: 0, c7: 0, c2: 0} 
[c0, c1, c2 - c7, c3 - c7, c4, c5, c6, c8] --> {c3: r19, c8: 0, c4: 0, c5: 0, c0: 0, c1: 0, c6: 0, c7: r19, c2: r19} 
[c0, c1, c2 - c7, c3 + c7, c4, c5, c6, c8] --> {c3: r20, c8: 0, c4: 0, c5: 0, c0: 0, c1: 0, c6: 0, c7: -r20, c2: -r20} 
[c0, c1, c2 + c7, c3 - c7, c4, c5, c6, c8] --> {c3: r21, c8: 0, c4: 0, c5: 0, c0: 0, c1: 0, c6: 0, c7: r21, c2: -r21} 
[c0, c1, c2 + c7, c3 + c7, c4, c5, c6, c8] --> {c3: r22, c8: 0, c4: 0, c5: 0, c0: 0, c1: 0, c6: 0, c7: -r22, c2: r22} 
[c0, c1, c2 - c6, c3, c4 + c6, c5, c7, c8] --> {c3: 0, c8: 0, c4: -r23, c5: 0, c0: 0, c1: 0, c6: r23, c7: 0, c2: r23} 
[c0, c1, c2 + c6, c3, c4 + c6, c5, c7, c8] --> {c3: 0, c8: 0, c4: r24, c5: 0, c0: 0, c1: 0, c6: -r24, c7: 0, c2: r24} 
```

We testing out procedure we compute the ideal of the Roman surface parametrized by `g`.

```python
eqg = sage_ideal( [y[i] - g[i] for i in range( 4 )] ).elimination_ideal( x ).gens()
print( eqg )
print( str( eqg.subs( {y[0]:1} ) ).replace( 'y1', 'x' ).replace( 'y2', 'y' ).replace( 'y3', 'z' ) )
```

Output:
```    
    [y1^2*y2^2 - y0*y1*y2*y3 + y1^2*y3^2 + y2^2*y3^2]
    [x^2*y^2 + x^2*z^2 + y^2*z^2 - x*y*z] 
```
 
For each of the 24 solutions in `sol_lst` we recover the projective automorphism `U`.  
Each of the 24 symmetries of the [Roman surface](https://en.wikipedia.org/wiki/Roman_surface) corresponds to the symmetries of a tetrahedron.

```python
for sol in sol_lst:
    
    # compute the projective isomorphism in terms of parametrized matrix U
    Ef = sage_matrix( sage_QQ, list( Mf ) + list( Kf.T ) )
    Egr = sage_matrix( list( Mgr.subs( sol ) ) + list( Kf.T ) )
    UpI = Egr * ~Ef
    assert ( UpI.submatrix( 4, 4 ) - sage_identity_matrix( 2 ) ).is_zero()
    U = UpI.submatrix( 0, 0, 4, 4 )
    U = U / sage_gcd( U.list() )
    assert U.dimensions() == ( 4, 4 )
    
    # verify whether U*f is a parametrization for Y for all (c0,...,c7)
    Uf = list( U * sage_vector( f ) )
    eqg_sub = [ eq.subs( {y[i]:Uf[i] for i in range( 4 )} ) for eq in eqg ]
    assert eqg_sub == [0]

    # output U with corresponding solution
    print( 'U =', list( U ), ', sol =', sol )
```

Output:
```
U = [(1, 0, 0, 0), (0,  0,  0,  1), (0,  1,  0,  0), (0,  0,  1,  0)] , sol = {c3: 0, c8: 0, c4: 0, c5: r1, c0: 0, c1: r1, c6: r1, c7: 0, c2: 0} 
U = [(1, 0, 0, 0), (0,  0,  0, -1), (0, -1,  0,  0), (0,  0,  1,  0)] , sol = {c3: 0, c8: 0, c4: 0, c5: r2, c0: 0, c1: -r2, c6: r2, c7: 0, c2: 0} 
U = [(1, 0, 0, 0), (0,  0,  0, -1), (0,  1,  0,  0), (0,  0, -1,  0)] , sol = {c3: 0, c8: 0, c4: 0, c5: r3, c0: 0, c1: -r3, c6: -r3, c7: 0, c2: 0} 
U = [(1, 0, 0, 0), (0,  0,  0,  1), (0, -1,  0,  0), (0,  0, -1,  0)] , sol = {c3: 0, c8: 0, c4: 0, c5: r4, c0: 0, c1: r4, c6: -r4, c7: 0, c2: 0} 
U = [(1, 0, 0, 0), (0,  0,  1,  0), (0,  1,  0,  0), (0,  0,  0,  1)] , sol = {c3: 0, c8: 0, c4: 0, c5: r5, c0: r5, c1: 0, c6: 0, c7: r5, c2: 0} 
U = [(1, 0, 0, 0), (0,  0, -1,  0), (0, -1,  0,  0), (0,  0,  0,  1)] , sol = {c3: 0, c8: 0, c4: 0, c5: r6, c0: -r6, c1: 0, c6: 0, c7: r6, c2: 0} 
U = [(1, 0, 0, 0), (0,  0,  1,  0), (0, -1,  0,  0), (0,  0,  0, -1)] , sol = {c3: 0, c8: 0, c4: 0, c5: r7, c0: r7, c1: 0, c6: 0, c7: -r7, c2: 0} 
U = [(1, 0, 0, 0), (0,  0, -1,  0), (0,  1,  0,  0), (0,  0,  0, -1)] , sol = {c3: 0, c8: 0, c4: 0, c5: r8, c0: -r8, c1: 0, c6: 0, c7: -r8, c2: 0} 
U = [(1, 0, 0, 0), (0,  1,  0,  0), (0,  0,  1,  0), (0,  0,  0,  1)] , sol = {c3: 0, c8: r9, c4: r9, c5: 0, c0: r9, c1: 0, c6: 0, c7: 0, c2: 0} 
U = [(1, 0, 0, 0), (0, -1,  0,  0), (0,  0, -1,  0), (0,  0,  0,  1)] , sol = {c3: 0, c8: r10, c4: r10, c5: 0, c0: -r10, c1: 0, c6: 0, c7: 0, c2: 0} 
U = [(1, 0, 0, 0), (0,  1,  0,  0), (0,  0, -1,  0), (0,  0,  0, -1)] , sol = {c3: 0, c8: r11, c4: -r11, c5: 0, c0: -r11, c1: 0, c6: 0, c7: 0, c2: 0} 
U = [(1, 0, 0, 0), (0, -1,  0,  0), (0,  0,  1,  0), (0,  0,  0, -1)] , sol = {c3: 0, c8: r12, c4: -r12, c5: 0, c0: r12, c1: 0, c6: 0, c7: 0, c2: 0} 
U = [(1, 0, 0, 0), (0,  0,  0,  1), (0,  0,  1,  0), (0,  1,  0,  0)] , sol = {c3: 0, c8: 0, c4: r13, c5: 0, c0: 0, c1: 0, c6: r13, c7: 0, c2: r13} 
U = [(1, 0, 0, 0), (0,  0,  0, -1), (0,  0, -1,  0), (0,  1,  0,  0)] , sol = {c3: 0, c8: 0, c4: -r14, c5: 0, c0: 0, c1: 0, c6: -r14, c7: 0, c2: r14} 
U = [(1, 0, 0, 0), (0,  1,  0,  0), (0,  0,  0,  1), (0,  0,  1,  0)] , sol = {c3: r15, c8: r15, c4: 0, c5: 0, c0: 0, c1: r15, c6: 0, c7: 0, c2: 0} 
U = [(1, 0, 0, 0), (0, -1,  0,  0), (0,  0,  0,  1), (0,  0, -1,  0)] , sol = {c3: -r16, c8: r16, c4: 0, c5: 0, c0: 0, c1: r16, c6: 0, c7: 0, c2: 0} 
U = [(1, 0, 0, 0), (0, -1,  0,  0), (0,  0,  0, -1), (0,  0,  1,  0)] , sol = {c3: r17, c8: r17, c4: 0, c5: 0, c0: 0, c1: -r17, c6: 0, c7: 0, c2: 0} 
U = [(1, 0, 0, 0), (0,  1,  0,  0), (0,  0,  0, -1), (0,  0, -1,  0)] , sol = {c3: -r18, c8: r18, c4: 0, c5: 0, c0: 0, c1: -r18, c6: 0, c7: 0, c2: 0} 
U = [(1, 0, 0, 0), (0,  0,  1,  0), (0,  0,  0,  1), (0,  1,  0,  0)] , sol = {c3: r19, c8: 0, c4: 0, c5: 0, c0: 0, c1: 0, c6: 0, c7: r19, c2: r19} 
U = [(1, 0, 0, 0), (0,  0, -1,  0), (0,  0,  0,  1), (0, -1,  0,  0)] , sol = {c3: r20, c8: 0, c4: 0, c5: 0, c0: 0, c1: 0, c6: 0, c7: -r20, c2: -r20} 
U = [(1, 0, 0, 0), (0,  0, -1,  0), (0,  0,  0, -1), (0,  1,  0,  0)] , sol = {c3: r21, c8: 0, c4: 0, c5: 0, c0: 0, c1: 0, c6: 0, c7: r21, c2: -r21} 
U = [(1, 0, 0, 0), (0,  0,  1,  0), (0,  0,  0, -1), (0, -1,  0,  0)] , sol = {c3: r22, c8: 0, c4: 0, c5: 0, c0: 0, c1: 0, c6: 0, c7: -r22, c2: r22} 
U = [(1, 0, 0, 0), (0,  0,  0, -1), (0,  0,  1,  0), (0, -1,  0,  0)] , sol = {c3: 0, c8: 0, c4: -r23, c5: 0, c0: 0, c1: 0, c6: r23, c7: 0, c2: r23} 
U = [(1, 0, 0, 0), (0,  0,  0,  1), (0,  0, -1,  0), (0, -1,  0,  0)] , sol = {c3: 0, c8: 0, c4: r24, c5: 0, c0: 0, c1: 0, c6: -r24, c7: 0, c2: r24} 
```




### Example 2: Projective equivalence of surfaces parametrized by bihomogeneous polynomials

We test projective equivalence of two surfaces that are represented by a parametrization with bihomogeneous polynomials.


```python
# we explicitly import the required modules
#
from surface_equivalence.class_se_ring import SERing
from surface_equivalence.find_equivalences import iso_P1xP1
from surface_equivalence.class_se_tools import SETools

SETools.filter( [] ) # replace [] with None for debug info
os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/maple/link/bin' # we need Maple

# create parameters for random parametrizations
#
d1, d2 = SERing.random_elt( [2, 3] ), SERing.random_elt( [2, 3] ) # bidegree
mon_lst = SERing.get_mon_P1xP1( d1, d2 ) # basis of monomial of bidegree (d1,d2)
mat1 = SERing.random_matrix_QQ( 4, len(mon_lst) )
while True:
    mati = SERing.random_inv_matrix_QQ( 4 )
    matL = SERing.random_inv_matrix_QQ( 2 )
    matR = SERing.random_inv_matrix_QQ( 2 )
    if mati[0, 0] != 0 and matL[0, 0] != 0 and matR[0, 0] != 0:
        break
        
flip = SERing.random_elt( [True, False] ) # True if element in Aut(P1xP1) flips the P1-factors
print 'flip     =', flip 
print 'mati     =', list(mati) 


# create random parametrizations
#
p1_lst = list( mat1 * sage_vector( mon_lst ) )
p2_lst = SERing.compose_aut_P1P1( p1_lst, matL, matR, flip )
p2_lst = list( mati * sage_vector( p2_lst ) )
print 'bidegree =', d1,d2
print 'p1_lst   =', p1_lst 
print 'p2_lst   =', p2_lst 

mato = iso_P1xP1( p1_lst, p2_lst )
print 'mato     =', list(mato)

a00 = mati[0, 0]
b00 = mato[0, 0]

( a00 / b00 ) * mato == mati

```
Output:

    flip     = True
    mati     = [(-5/11, 0, 1, 3/10), (4, -9, 0, 0), (3/5, -4/3, -5/76, -33/7), (1, 2, 1/4, 1/2)]
    bidegree = 3 3
    p1_lst   = [13/5*s^3*u^3 + s^2*t*u^3 - 3*s*t^2*u^3 - s^3*u^2*v + s^2*t*u^2*v - s*t^2*u^2*v + 1/4*t^3*u^2*v + 3/4*s^2*t*u*v^2 - s*t^2*u*v^2 - t^3*u*v^2 + 1/3*s^3*v^3 + 1/8*s^2*t*v^3 - 1/111*s*t^2*v^3 - 1/3*t^3*v^3, -2/3*s^3*u^3 + 1/5*s^2*t*u^3 + s*t^2*u^3 + 12/5*t^3*u^3 - s^3*u^2*v + s*t^2*u^2*v - 6*s^3*u*v^2 - 1/13*s^2*t*u*v^2 - s*t^2*u*v^2 + 6*t^3*u*v^2 + 5*s^3*v^3 - 17*s^2*t*v^3 - 2*s*t^2*v^3 - t^3*v^3, -s^3*u^3 + s^2*t*u^3 + s*t^2*u^3 - 4*t^3*u^3 + 1/11*s^3*u^2*v - 4*s^2*t*u^2*v + 2*s*t^2*u^2*v + 1/3*t^3*u^2*v - 9/40*s^3*u*v^2 - 1/3*s^2*t*u*v^2 - 1/4*s*t^2*u*v^2 + 1/3*t^3*u*v^2 + 5/2*s^3*v^3 + s^2*t*v^3 - 1/2*t^3*v^3, 1/2*s^2*t*u^3 - 1/3*t^3*u^3 + 2*s^3*u^2*v - 1/2*s^2*t*u^2*v - s*t^2*u^2*v + 3*t^3*u^2*v - s^3*u*v^2 + s^2*t*u*v^2 - 2*s*t^2*u*v^2 + 131/3*s^3*v^3 + 2*s^2*t*v^3 - t^3*v^3]
    p2_lst   = [31976963/4111884*s^3*u^3 + 442566851/678460860*s^2*t*u^3 + 79327517/678460860*s*t^2*u^3 - 9955163/1221229548*t^3*u^3 - 12122399551/25128180*s^3*u^2*v - 7546163291/226153620*s^2*t*u^2*v - 817668211/135692172*s*t^2*u^2*v - 15425701/26433540*t^3*u^2*v + 266767797911/25128180*s^3*u*v^2 + 7628052919/10769220*s^2*t*u*v^2 + 201784817/2937060*s*t^2*u*v^2 + 1018351253/96922980*t^3*u*v^2 - 129246969755/1675212*s^3*v^3 - 128508257761/25128180*s^2*t*v^3 - 126597503/512820*s*t^2*v^3 - 13032476179/226153620*t^3*v^3, -13723953589/133636230*s^3*u^3 - 745505287/44545410*s^2*t*u^3 - 317687029/400908690*s*t^2*u^3 - 188294741/3608178210*t^3*u^3 + 3113754401/685314*s^3*u^2*v + 120941227/228438*s^2*t*u^2*v + 35407201/2055942*s*t^2*u^2*v + 26283841/18503478*t^3*u^2*v - 887011361503/14848470*s^3*u*v^2 - 149819337727/14848470*s^2*t*u*v^2 - 12663809183/44545410*s*t^2*u*v^2 - 28416470207/400908690*t^3*u*v^2 + 596677426187/4949490*s^3*v^3 + 37903909949/707070*s^2*t*v^3 + 7888481689/4949490*s*t^2*v^3 + 58937728363/133636230*t^3*v^3, -302689209253889/2469597530400*s^3*u^3 - 25876498952281/2469597530400*s^2*t*u^3 - 494351919163/2469597530400*s*t^2*u^3 - 64757491969/7408792591200*t^3*u^3 + 2951850694437049/401045752800*s^3*u^2*v + 647334286605043/1203137258400*s^2*t*u^2*v + 92119037981627/3609411775200*s*t^2*u^2*v + 1334344009783/1546890760800*t^3*u^2*v - 264432209582610019/1737864928800*s^3*u*v^2 - 7961576227612279/744799255200*s^2*t*u*v^2 - 1026657339268511/2234397765600*s*t^2*u*v^2 - 32348959205153/1737864928800*t^3*u*v^2 + 1803918982248651413/1737864928800*s^3*v^3 + 18005314299774091/248266418400*s^2*t*v^3 + 15362016458673413/5213594786400*s*t^2*v^3 + 492966380121493/5213594786400*t^3*v^3, 824101924631/23519976480*s^3*u^3 + 47464870363/10079989920*s^2*t*u^3 + 1185896993/10079989920*s*t^2*u^3 + 286968677/18143981856*t^3*u^3 - 1055071129103/603076320*s^3*u^2*v - 317744383909/1809228960*s^2*t*u^2*v - 18877533469/5427686880*s*t^2*u^2*v + 43558783/2326151520*t^3*u^2*v + 8500878427533/290370080*s^3*u*v^2 + 25481060398727/7839992160*s^2*t*u*v^2 + 1116863242031/23519976480*s*t^2*u*v^2 - 31079181383/14111985888*t^3*u*v^2 - 123113598534193/871110240*s^3*v^3 - 49221364180763/2613330720*s^2*t*v^3 - 30031112093/124444320*s*t^2*v^3 + 55539169097/1567998432*t^3*v^3]
    mato     = [(-8990228480/101871, 0, 1798045696/9261, 899022848/15435), (7192182784/9261, -1798045696/1029, 0, 0), (1798045696/15435, -7192182784/27783, -118292480/9261, -19778502656/21609), (1798045696/9261, 3596091392/9261, 449511424/9261, 899022848/9261)]
    True

    