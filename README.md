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


### Example 1: Projective automorphisms of Roman surfaces (B1)

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

We testing the outcome below we first compute the ideal of the Roman surface parametrized by `g`.

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
 
For each of the 24 solutions in `sol_lst` obtained, we recover the corresponding projective automorphism `U`.
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




### Example 2: Projective isomorphisms between Veronese-Segre surfaces (B1)

We compute the projective isomorphisms between two Veronese-Segre surfaces.
A Veronese-Segre surface is a projective surface that can be parametrized by a map 
whose components consists of bihomogeneous polynomials. Thus the domain of this map
is the fiber product of the projective line with itself.

We start by importing the required libraries and initialize parameters.

```python
from surface_equivalence.class_se_ring import ring
from surface_equivalence.class_se_ring import SERing
from surface_equivalence.sage_interface import sage_gcd
from surface_equivalence.sage_interface import sage_matrix
from surface_equivalence.sage_interface import sage_identity_matrix
from surface_equivalence.sage_interface import sage_vector
from surface_equivalence.sage_interface import sage_ideal
from surface_equivalence.sage_interface import sage_maple
from surface_equivalence.sage_interface import sage_QQ
from surface_equivalence.sage_interface import sage_SR
from surface_equivalence.sage_interface import sage_solve
os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/maple/link/bin' # edit PATH
```

Initialize the parameters for the parametrizations `f` and `g`.
See [this file](https://github.com/niels-lubbes/surface_equivalence/blob/master/surface_equivalence/src/surface_equivalence/__main__.py)
for code for obtaining random parameters.

```python
# bidegree of f
d1, d2 = 2, 2

# coefficient matrix of f
matf = sage_matrix( sage_QQ, ring( '[(3/5, 87, -1/2, 1/4, 0, 1/8, 16/5, 0, 0), (0, 11, -1/9, 44, 0, 0, -1, 0, 1/3), (1/12, -4, 0, 2, -1, 1/4, 0, -3, -1/9), (1/3, 1/4, 1, -4, 1/2, -1, -6, 2, 22)]' ) )

# projective automorphism P^4--->P^4
matU = sage_matrix( sage_QQ, ring( '[( 1, 0, -1, -3 ), ( -1 / 18, -3, 0, 15 / 2 ), ( -1, 0, 13, 1 ), ( -3, -1 / 11, -2, 3 / 23 )]' ) )

# projective automorphism s: P^1xP^1--> P^1xP^1.
L, R  = ring( '[-1, 0, -1/87, -1/2] ' ), ring( '[2, 1, 0, 2]' )
y = [ring( 'y' + str( i ) ) for i in range( 4 )]
s = {y[0]:L[0] * y[0] + L[1] * y[1], y[1]:L[2] * y[0] + L[3] * y[1], y[2]:R[0] * y[2] + R[1] * y[3], y[3]:R[2] * y[2] + R[3] * y[3]}
```

Construct `f` and `g` from the parameters.

```python
f = list( matf * sage_vector( SERing.get_mon_P1xP1( d1, d2 ) ) )
g = list( matU * sage_vector( [ comp.subs( s ) for comp in f ] ) )
assert set( SERing.get_bidegree( f ) ) == set( SERing.get_bidegree( g ) )
print( 'f =', f )
print( 'g =', g )
```
Output:
```
f = [3/5*y0^2*y2^2 + 1/4*y0*y1*y2^2 + 16/5*y1^2*y2^2 + 87*y0^2*y2*y3 - 1/2*y0^2*y3^2 + 1/8*y0*y1*y3^2, 44*y0*y1*y2^2 - y1^2*y2^2 + 11*y0^2*y2*y3 - 1/9*y0^2*y3^2 + 1/3*y1^2*y3^2, 1/12*y0^2*y2^2 + 2*y0*y1*y2^2 - 4*y0^2*y2*y3 - y0*y1*y2*y3 - 3*y1^2*y2*y3 + 1/4*y0*y1*y3^2 - 1/9*y1^2*y3^2, 1/3*y0^2*y2^2 - 4*y0*y1*y2^2 - 6*y1^2*y2^2 + 1/4*y0^2*y2*y3 + 1/2*y0*y1*y2*y3 + 2*y1^2*y2*y3 + y0^2*y3^2 - y0*y1*y3^2 + 22*y1^2*y3^2] 
g = [-54908/37845*y0^2*y2^2 + 18683/870*y0*y1*y2^2 + 106/5*y1^2*y2^2 + 13606207/37845*y0^2*y2*y3 + 17693/870*y0*y1*y2*y3 + 91/5*y1^2*y2*y3 + 56616167/340605*y0^2*y3^2 + 235537/31320*y0*y1*y3^2 - 2794/45*y1^2*y3^2, 1631813/681210*y0^2*y2^2 - 5104643/15660*y0*y1*y2^2 - 1898/45*y1^2*y2^2 - 48113021/340605*y0^2*y2*y3 - 4976393/15660*y0*y1*y2*y3 - 1223/45*y1^2*y2*y3 - 109107127/2724840*y0^2*y3^2 - 5316293/62640*y0*y1*y3^2 + 7243/45*y1^2*y3^2, 161288/37845*y0^2*y2^2 + 37477/870*y0*y1*y2^2 - 46/5*y1^2*y2^2 - 2318353/4205*y0^2*y2*y3 + 4749/290*y0*y1*y2*y3 - 231/5*y1^2*y2*y3 - 92175587/340605*y0^2*y3^2 + 81863/31320*y0*y1*y3^2 - 11/45*y1^2*y3^2, -77790143/9574785*y0^2*y2^2 - 4185757/220110*y0*y1*y2^2 - 13019/1265*y1^2*y2^2 - 9803581988/9574785*y0^2*y2*y3 - 3213247/220110*y0*y1*y2*y3 - 5099/1265*y1^2*y2*y3 - 173517820907/344692260*y0^2*y3^2 - 1521101/344520*y0*y1*y3^2 + 164809/45540*y1^2*y3^2] 

```

We will now recover the projective automorphism defined by the matrix `matU` from only `f` and `g`.

```python
# Superset of compatible reparametrizations P^1xP^1--->P^1xP^1 consists of two families r0 and r1 parametrized by c.
# Notice that r1 flips the factors of P^1xP^1.
y = [ring( 'y' + str( i ) ) for i in range( 4 )]
c = [ring( 'c' + str( i ) ) for i in range( 8 )]
r0 = {y[0]:c[0] * y[0] + c[1] * y[1], y[1]:c[2] * y[0] + c[3] * y[1], y[2]:c[4] * y[2] + c[5] * y[3], y[3]:c[6] * y[2] + c[7] * y[3]}
r1 = {y[2]:c[0] * y[0] + c[1] * y[1], y[3]:c[2] * y[0] + c[3] * y[1], y[0]:c[4] * y[2] + c[5] * y[3], y[1]:c[6] * y[2] + c[7] * y[3]}
```

We check whether we can find a solution for `c` with compatible reparametrization `r0`.

```python
# First try r0
r=r0

# compute kernel and coefficient matrix of f
Mf = SERing.get_matrix_P1xP1( f )
Kf = Mf.right_kernel_matrix().T
assert ( Mf * Kf ).is_zero()

# compute the coefficient matrix of g composed with r
gr = [ comp.subs( r ) for comp in g ]
Mgr = SERing.get_matrix_P1xP1( gr )
assert sage_gcd( gr ) == 1

# compute c such that Mgr*Kf==0
ec_lst = ( Mgr * Kf ).list() + ring( '[(c0*c3-c1*c2)*t-1, (c4*c7-c5*c6)*t-1]' )   
sage_maple.eval( 'with(Groebner);' )
sage_maple.eval( 'gb := Basis( ' + str( ec_lst ) + ', plex(' + str( c )[1:-1] + ', t) );' )
gb_lst = ring( sage_maple.eval( 'lprint(gb);' ) )
assert gb_lst != [1]
pc_lst = sage_ideal( gb_lst ).elimination_ideal( ring( 't' ) ).primary_decomposition()
assert len(pc_lst)==1
sol_lst = ring( sage_solve( [sage_SR( comp ) for comp in pc_lst[0].gens()], [sage_SR( comp ) for comp in c], solution_dict = True ) )
sol_lst = [ sol for sol in sol_lst if sol.values() != 8 * [0] ]
assert len(sol_lst)==2
sol = sol_lst[0] # both solutions are equivalent
print( 'sol_lst =', sol_lst )
print( 'sol     =', sol )
```
Output:
```
sol_lst = [{c3: 2*r3, c4: sqrt(2)*r3, c5: -1/2*sqrt(2)*r3, c0: r3, c1: 0, c6: 0, c7: sqrt(2)*r3, c2: -2/87*r3}, {c3: 2*r4, c4: -sqrt(2)*r4, c5: 1/2*sqrt(2)*r4, c0: r4, c1: 0, c6: 0, c7: -sqrt(2)*r4, c2: -2/87*r4}]
sol     = {c3: 2*r3, c4: sqrt(2)*r3, c5: -1/2*sqrt(2)*r3, c0: r3, c1: 0, c6: 0, c7: sqrt(2)*r3, c2: -2/87*r3}
```

We now compute from the solution `sol` for `c` the matrix `U`, which defines a projective isomorphisms.



```python
#.subs({sage_SR('r1'):1,sage_SR('r2'):1})
Ef = sage_matrix( sage_QQ, list( Mf ) + list( Kf.T ) )
Egr = sage_matrix( list( Mgr.subs( sol ) ) + list( Kf.T ) )
UpI = Egr * ~Ef
U = UpI.submatrix( 0, 0, 4, 4 )
assert ( UpI.submatrix( 4, 4 ) - sage_identity_matrix( len( SERing.get_mon_P1xP1( d1, d2 ) ) - 4 ) ).is_zero()
assert U.dimensions() == ( 4, 4 )

# output U with corresponding solution
print( 'U =', list( U ), ', sol =', sol )
assert ( matU[0, 0] / U[0, 0] ) * U == matU
```

We repeat the same procedure for the case `r=r1`, but find that there are no solutions in this case.

```python
r = r1
Mf = SERing.get_matrix_P1xP1( f )
Kf = Mf.right_kernel_matrix().T
gr = [ comp.subs( r ) for comp in g ]
Mgr = SERing.get_matrix_P1xP1( gr )
ec_lst = ( Mgr * Kf ).list() + ring( '[(c0*c3-c1*c2)*t-1, (c4*c7-c5*c6)*t-1]' )   
sage_maple.eval( 'with(Groebner);' )
sage_maple.eval( 'gb := Basis( ' + str( ec_lst ) + ', plex(' + str( c )[1:-1] + ', t) );' )
gb_lst = ring( sage_maple.eval( 'lprint(gb);' ) )
assert gb_lst == [1]
```






    
