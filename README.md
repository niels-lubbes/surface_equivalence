# Surface Equivalence 


## Introduction

Surface-Equivalence is a library for computing projective isomorphisms between surfaces, if such isomorphisms exists.
The underlying theory for these for these algorithm was developed by 
[Bert Jüttler](http://www.ag.jku.at/), 
[Niels Lubbes](https://nielslubbes.com) and 
[Josef Schicho](https://www3.risc.jku.at/people/jschicho/).
We refer to [arxiv:](https://arxiv.org/abs/) for more information.

This library depends on [SageMath](https://SageMath.org) libraries. 
From some parts we use functionality of [Maple](https://www.maplesoft.com) 
and [Mathematica](https://www.wolfram.com/mathematica/).


## Installation

* Install Sage from [SageMath](https://SageMath.org).
We assume that `sage` is accessible from your commandline interface.

* Install [Maple](https://www.maplesoft.com).
We assume for some examples that `maple` is accessible from your commandline interface.

* Install [Mathematica](https://www.wolfram.com/mathematica/).
We assume for some examples that `mathematica` is accessible from your commandline interface.


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

For pasting code into a Sage terminal type `%paste`.
See [this file](https://github.com/niels-lubbes/surface_equivalence/blob/master/surface_equivalence/src/surface_equivalence/__main__.py) 
for more example usecases. 
See the [source code](https://github.com/niels-lubbes/surface_equivalence/blob/master/surface_equivalence/src/surface_equivalence)
the io-specification of each function.
The [test functions](https://github.com/niels-lubbes/surface_equivalence/blob/master/surface_equivalence/src/tests)
might be informative for how to call each function.


### Example 1: Projective automorphisms of a Roman surface (case B1)

For an explanation of the underlying theory behind the code we refer to
[(Example 3, arxiv????)](https://arxiv.org/abs/).

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
We denote the parametrizations of the two [Roman surfaces](https://en.wikipedia.org/wiki/Roman_surface) by f and g.

```python
f = ring( '[x0^2+x1^2+x2^2,x0*x1,x0*x2,x1*x2]' )
g = f

```
Our goal is to recover matrices U that define the projective automorphisms of the Roman surface.
For this purpose we first compute coefficient matrices.

```python
# compatible reparametrizations are linear and indexed by c
r = {x[0]:c[0] * y[0] + c[1] * y[1] + c[2] * y[2], x[1]:c[3] * y[0] + c[4] * y[1] + c[5] * y[2], x[2]:c[6] * y[0] + c[7] * y[1] + c[8] * y[2]}

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

We compute all solutions for c such that `Mgr*Kf==0`. This takes a few seconds.

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

We testing the outcome below we first compute the ideal of the Roman surface parametrized by g.

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
For each of the 24 solutions in `sol_lst` obtained, we recover the corresponding projective automorphism U.
Each of the 24 symmetries of the [Roman surface](https://en.wikipedia.org/wiki/Roman_surface) corresponds to the symmetries of a tetrahedron.

```python
# type "%paste" for pasting indented code into a Python or Sage terminal
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

### Example 2: Projective isomorphisms between Veronese-Segre surfaces (case B1)

We compute the projective isomorphisms between two Veronese-Segre surfaces.
A Veronese-Segre surface is a projective surface that can be parametrized by a map 
whose components consists of bihomogeneous polynomials. Thus the domain of this map
is the fiber product of the projective line with itself.

We start by importing the required libraries.

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
Initialize the parameters for the parametrizations f and g.
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
Construct f and g from the parameters.

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

#### Compute projective isomorphisms

We will now recover the projective automorphism defined by the matrix `matU` from only f and g.

```python
# Superset of compatible reparametrizations P^1xP^1--->P^1xP^1 consists of 
# two families r0 and r1 parametrized by c.
# Notice that r1 flips the factors of P^1xP^1.
y = [ring( 'y' + str( i ) ) for i in range( 4 )]
c = [ring( 'c' + str( i ) ) for i in range( 8 )]
r0 = {y[0]:c[0] * y[0] + c[1] * y[1], y[1]:c[2] * y[0] + c[3] * y[1], y[2]:c[4] * y[2] + c[5] * y[3], y[3]:c[6] * y[2] + c[7] * y[3]}
r1 = {y[2]:c[0] * y[0] + c[1] * y[1], y[3]:c[2] * y[0] + c[3] * y[1], y[0]:c[4] * y[2] + c[5] * y[3], y[1]:c[6] * y[2] + c[7] * y[3]}

```
We check whether we can find a solution for c with compatible reparametrization r0.

```python
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

We now compute from the solution sol for c the matrix U, which defines a projective isomorphisms.

```python
Ef = sage_matrix( sage_QQ, list( Mf ) + list( Kf.T ) )
Egr = sage_matrix( list( Mgr.subs( sol ) ) + list( Kf.T ) )
UpI = Egr * ~Ef
U = UpI.submatrix( 0, 0, 4, 4 )
assert ( UpI.submatrix( 4, 4 ) - sage_identity_matrix( len( SERing.get_mon_P1xP1( d1, d2 ) ) - 4 ) ).is_zero()
assert U.dimensions() == ( 4, 4 )
print( 'U =', list( U ), ', sol =', sol )

# we verify that we indeed recovered the projective isomorphism
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

### Example 4: Projective isomorphisms between surfaces (case B2)

We compute the projective isomorphisms between two projective surfaces
that are adjoint to double ruled quadrics.
For an explanation of the code we refer to [(Example 14, arxiv????)](https://arxiv.org/abs/).

We start by importing the required libraries and declaring parameters.

```python
from surface_equivalence.class_se_ring import ring
from surface_equivalence.class_se_ring import SERing
from surface_equivalence.sage_interface import sage_gcd
from surface_equivalence.sage_interface import sage_matrix
from surface_equivalence.sage_interface import sage_identity_matrix
from surface_equivalence.sage_interface import sage_vector
from surface_equivalence.sage_interface import sage_ideal
from surface_equivalence.sage_interface import sage_QQ
from surface_equivalence.sage_interface import sage_SR
from surface_equivalence.sage_interface import sage_solve
from surface_equivalence.sage_interface import sage_diff
from linear_series.class_poly_ring import PolyRing
from linear_series.class_base_points import BasePointTree
from linear_series.class_linear_series import LinearSeries

x = [ring( 'x' + str( i ) ) for i in range( 3 )]
y = [ring( 'y' + str( i ) ) for i in range( 4 )]
z = [ring( 'z' + str( i ) ) for i in range( 4 )]
c = [ring( 'c' + str( i ) ) for i in range( 8 )]

```

We initialize the parametric maps f and g.

```python
f = ring( '[x0^6*x1^2,x0*x1^5*x2^2,x1^3*x2^5,x0^5*x2^3+x0^5*x2^3+x0^5*x1*x2^2]' )
g = ring( '[y0^3*y1^2*y2^5,y1^5*y2^3*y3^2,y0^2*y1^3*y3^5,y0^5*y2^2*y3^3+y0^4*y1*y2^3*y3^2]' )
g = [g[0], g[1] + g[0], g[2], g[3] + g[2]]
assert sage_gcd( f ) == 1
assert sage_gcd( g ) == 1

```

#### Basepoint analysis for the maps f and g

We do a a basepoint analysis for f and g.

```python
bf = LinearSeries( SERing.conv( f ), PolyRing( 'x,y,z', True ) ).get_bp_tree()
bg = LinearSeries( SERing.conv( g ), PolyRing( 'x,y,v,w', True ) ).get_bp_tree()
print( 'bf =' + str(bf) )
print( 'bg =' + str(bg) )

```
Output:

```
bf = 
{ 4, <<x^2*z^6, x^5*y^2*z, x^3*y^5, x*y^2*z^5 + 2*y^3*z^5>>, QQ[x, y, z] }
chart=z, depth=0, mult=2, sol=(0, 0), { 4, <<x^2, x^5*y^2, x^3*y^5, x*y^2 + 2*y^3>>, QQ[x, y] }
    chart=t, depth=1, mult=1, sol=(0, 0), { 4, <<x^2, x^5*y^5, x^3*y^6, x*y + 2*y>>, QQ[x, y] }
        chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<x, x^9*y^5, x^8*y^6, x*y + 2*y>>, QQ[x, y] }
chart=x, depth=0, mult=3, sol=(0, 0), { 4, <<z^6, y^2*z, y^5, 2*y^3*z^5 + y^2*z^5>>, QQ[y, z] }
    chart=t, depth=1, mult=2, sol=(0, 0), { 4, <<z^3, y^2, y^5*z^2, 2*y^3*z^5 + y^2*z^4>>, QQ[y, z] }
        chart=t, depth=2, mult=1, sol=(0, 0), { 4, <<z, y^2, y^5*z^5, 2*y^3*z^6 + y^2*z^4>>, QQ[y, z] }
            chart=s, depth=3, mult=1, sol=(0, 0), { 4, <<z, y, y^9*z^5, 2*y^8*z^6 + y^5*z^4>>, QQ[y, z] }
    chart=s, depth=1, mult=1, sol=(0, 0), { 4, <<y^3*z^6, z, y^2, 2*y^5*z^5 + y^4*z^5>>, QQ[y, z] }
        chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<y^8*z^6, z, y, 2*y^9*z^5 + y^8*z^5>>, QQ[y, z] }
chart=y, depth=0, mult=3, sol=(0, 0), { 4, <<x^2*z^6, x^5*z, x^3, x*z^5 + 2*z^5>>, QQ[x, z] }
    chart=t, depth=1, mult=2, sol=(0, 0), { 4, <<x^2*z^5, x^5*z^3, x^3, x*z^3 + 2*z^2>>, QQ[x, z] }
        chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<x^5*z^5, x^6*z^3, x, x^2*z^3 + 2*z^2>>, QQ[x, z] }
            chart=t, depth=3, mult=1, sol=(0, 0), { 4, <<x^5*z^9, x^6*z^8, x, x^2*z^4 + 2*z>>, QQ[x, z] } 
bg = 
{ 4, <<x^3*y^2*v^5, x^3*y^2*v^5 + y^5*v^3*w^2, x^2*y^3*w^5, x^4*y*v^3*w^2 + x^5*v^2*w^3 + x^2*y^3*w^5>>, QQ[x, y, v, w] }
chart=xv, depth=0, mult=2, sol=(0, 0), { 4, <<y^2, y^5*w^2 + y^2, y^3*w^5, y^3*w^5 + y*w^2 + w^3>>, QQ[y, w] }
    chart=t, depth=1, mult=1, sol=(0, 0), { 4, <<y^2, y^5*w^5 + y^2, y^3*w^6, y^3*w^6 + y*w + w>>, QQ[y, w] }
        chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<y, y^9*w^5 + y, y^8*w^6, y^8*w^6 + y*w + w>>, QQ[y, w] }
chart=xw, depth=0, mult=2, sol=(0, 0), { 4, <<y^2*v^5, y^5*v^3 + y^2*v^5, y^3, y*v^3 + y^3 + v^2>>, QQ[y, v] }
    chart=s, depth=1, mult=1, sol=(0, 0), { 4, <<y^5*v^5, y^5*v^5 + y^6*v^3, y, y^2*v^3 + v^2 + y>>, QQ[y, v] }
        chart=t, depth=2, mult=1, sol=(0, 0), { 4, <<y^5*v^9, y^6*v^8 + y^5*v^9, y, y^2*v^4 + y + v>>, QQ[y, v] }
chart=yv, depth=0, mult=2, sol=(0, 0), { 4, <<x^3, x^3 + w^2, x^2*w^5, x^5*w^3 + x^2*w^5 + x^4*w^2>>, QQ[x, w] }
    chart=s, depth=1, mult=1, sol=(0, 0), { 4, <<x, w^2 + x, x^5*w^5, x^5*w^5 + x^6*w^3 + x^4*w^2>>, QQ[x, w] }
        chart=t, depth=2, mult=1, sol=(0, 0), { 4, <<x, x + w, x^5*w^9, x^6*w^8 + x^5*w^9 + x^4*w^5>>, QQ[x, w] }
chart=yw, depth=0, mult=2, sol=(0, 0), { 4, <<x^3*v^5, x^3*v^5 + v^3, x^2, x^5*v^2 + x^4*v^3 + x^2>>, QQ[x, v] }
    chart=t, depth=1, mult=1, sol=(0, 0), { 4, <<x^3*v^6, x^3*v^6 + v, x^2, x^5*v^5 + x^4*v^5 + x^2>>, QQ[x, v] }
        chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<x^8*v^6, x^8*v^6 + v, x, x^9*v^5 + x^8*v^5 + x>>, QQ[x, v] }
```

#### Computing implicit equations for testing purposes

We compute the implicit equations of the images of the maps f and g for testing purposes.
This is not required for the algorithm.

```python
eqf = sage_ideal( [z[i] - f[i] for i in range( 4 )] ).elimination_ideal( x ).gens()
assert len( eqf ) == 1
assert eqf[0].degree() == 26
eqg = sage_ideal( [z[i] - g[i] for i in range( 4 )] ).elimination_ideal( y ).gens()
assert len( eqg ) == 1
assert eqg[0].degree() == 26
print( 'eqf =', eqf )
print( 'eqg =', eqg )

```
Output:

```
eqf = [z0^12*z1^6*z2^8 + 8192*z0^13*z2^13 - 53248*z0^12*z1*z2^12*z3 + 159744*z0^11*z1^2*z2^11*z3^2 - 292864*z0^10*z1^3*z2^10*z3^3 + 366080*z0^9*z1^4*z2^9*z3^4 - 329472*z0^8*z1^5*z2^8*z3^5 + 219648*z0^7*z1^6*z2^7*z3^6 - 109824*z0^6*z1^7*z2^6*z3^7 + 41184*z0^5*z1^8*z2^5*z3^8 - 11440*z0^4*z1^9*z2^4*z3^9 + 2288*z0^3*z1^10*z2^3*z3^10 - 312*z0^2*z1^11*z2^2*z3^11 + 26*z0*z1^12*z2*z3^12 - z1^13*z3^13] 
eqg = [z0^18*z2^8 - 6*z0^17*z1*z2^8 + 15*z0^16*z1^2*z2^8 - 20*z0^15*z1^3*z2^8 + 15*z0^14*z1^4*z2^8 - 6*z0^13*z1^5*z2^8 + z0^12*z1^6*z2^8 + z1^13*z2^13 + 13*z0*z1^12*z2^12*z3 - 13*z1^13*z2^12*z3 + 78*z0^2*z1^11*z2^11*z3^2 - 156*z0*z1^12*z2^11*z3^2 + 78*z1^13*z2^11*z3^2 + 286*z0^3*z1^10*z2^10*z3^3 - 858*z0^2*z1^11*z2^10*z3^3 + 858*z0*z1^12*z2^10*z3^3 - 286*z1^13*z2^10*z3^3 + 715*z0^4*z1^9*z2^9*z3^4 - 2860*z0^3*z1^10*z2^9*z3^4 + 4290*z0^2*z1^11*z2^9*z3^4 - 2860*z0*z1^12*z2^9*z3^4 + 715*z1^13*z2^9*z3^4 + 1287*z0^5*z1^8*z2^8*z3^5 - 6435*z0^4*z1^9*z2^8*z3^5 + 12870*z0^3*z1^10*z2^8*z3^5 - 12870*z0^2*z1^11*z2^8*z3^5 + 6435*z0*z1^12*z2^8*z3^5 - 1287*z1^13*z2^8*z3^5 + 1716*z0^6*z1^7*z2^7*z3^6 - 10296*z0^5*z1^8*z2^7*z3^6 + 25740*z0^4*z1^9*z2^7*z3^6 - 34320*z0^3*z1^10*z2^7*z3^6 + 25740*z0^2*z1^11*z2^7*z3^6 - 10296*z0*z1^12*z2^7*z3^6 + 1716*z1^13*z2^7*z3^6 + 1716*z0^7*z1^6*z2^6*z3^7 - 12012*z0^6*z1^7*z2^6*z3^7 + 36036*z0^5*z1^8*z2^6*z3^7 - 60060*z0^4*z1^9*z2^6*z3^7 + 60060*z0^3*z1^10*z2^6*z3^7 - 36036*z0^2*z1^11*z2^6*z3^7 + 12012*z0*z1^12*z2^6*z3^7 - 1716*z1^13*z2^6*z3^7 + 1287*z0^8*z1^5*z2^5*z3^8 - 10296*z0^7*z1^6*z2^5*z3^8 + 36036*z0^6*z1^7*z2^5*z3^8 - 72072*z0^5*z1^8*z2^5*z3^8 + 90090*z0^4*z1^9*z2^5*z3^8 - 72072*z0^3*z1^10*z2^5*z3^8 + 36036*z0^2*z1^11*z2^5*z3^8 - 10296*z0*z1^12*z2^5*z3^8 + 1287*z1^13*z2^5*z3^8 + 715*z0^9*z1^4*z2^4*z3^9 - 6435*z0^8*z1^5*z2^4*z3^9 + 25740*z0^7*z1^6*z2^4*z3^9 - 60060*z0^6*z1^7*z2^4*z3^9 + 90090*z0^5*z1^8*z2^4*z3^9 - 90090*z0^4*z1^9*z2^4*z3^9 + 60060*z0^3*z1^10*z2^4*z3^9 - 25740*z0^2*z1^11*z2^4*z3^9 + 6435*z0*z1^12*z2^4*z3^9 - 715*z1^13*z2^4*z3^9 + 286*z0^10*z1^3*z2^3*z3^10 - 2860*z0^9*z1^4*z2^3*z3^10 + 12870*z0^8*z1^5*z2^3*z3^10 - 34320*z0^7*z1^6*z2^3*z3^10 + 60060*z0^6*z1^7*z2^3*z3^10 - 72072*z0^5*z1^8*z2^3*z3^10 + 60060*z0^4*z1^9*z2^3*z3^10 - 34320*z0^3*z1^10*z2^3*z3^10 + 12870*z0^2*z1^11*z2^3*z3^10 - 2860*z0*z1^12*z2^3*z3^10 + 286*z1^13*z2^3*z3^10 + 78*z0^11*z1^2*z2^2*z3^11 - 858*z0^10*z1^3*z2^2*z3^11 + 4290*z0^9*z1^4*z2^2*z3^11 - 12870*z0^8*z1^5*z2^2*z3^11 + 25740*z0^7*z1^6*z2^2*z3^11 - 36036*z0^6*z1^7*z2^2*z3^11 + 36036*z0^5*z1^8*z2^2*z3^11 - 25740*z0^4*z1^9*z2^2*z3^11 + 12870*z0^3*z1^10*z2^2*z3^11 - 4290*z0^2*z1^11*z2^2*z3^11 + 858*z0*z1^12*z2^2*z3^11 - 78*z1^13*z2^2*z3^11 + 13*z0^12*z1*z2*z3^12 - 156*z0^11*z1^2*z2*z3^12 + 858*z0^10*z1^3*z2*z3^12 - 2860*z0^9*z1^4*z2*z3^12 + 6435*z0^8*z1^5*z2*z3^12 - 10296*z0^7*z1^6*z2*z3^12 + 12012*z0^6*z1^7*z2*z3^12 - 10296*z0^5*z1^8*z2*z3^12 + 6435*z0^4*z1^9*z2*z3^12 - 2860*z0^3*z1^10*z2*z3^12 + 858*z0^2*z1^11*z2*z3^12 - 156*z0*z1^12*z2*z3^12 + 13*z1^13*z2*z3^12 + z0^13*z3^13 - 13*z0^12*z1*z3^13 + 78*z0^11*z1^2*z3^13 - 286*z0^10*z1^3*z3^13 + 715*z0^9*z1^4*z3^13 - 1287*z0^8*z1^5*z3^13 + 1716*z0^7*z1^6*z3^13 - 1716*z0^6*z1^7*z3^13 + 1287*z0^5*z1^8*z3^13 - 715*z0^4*z1^9*z3^13 + 286*z0^3*z1^10*z3^13 - 78*z0^2*z1^11*z3^13 + 13*z0*z1^12*z3^13 - z1^13*z3^13] 

```
#### Computing reductions of f and g

We first compute the r0 reductions of the rational maps f and g.  

```python
hf = LinearSeries.get( [8], bf )
hg = LinearSeries.get( [5, 5], bg )
assert len( hf.pol_lst ) == len( hg.pol_lst ) == 16
print( len( hf.pol_lst ), SERing.conv( hf.pol_lst ) )
print( len( hg.pol_lst ), SERing.conv( hg.pol_lst ) )

```
Output:
```
usecase_B2(422): 16 [x0*x1^5*x2^2, x0*x1^4*x2^3, x0^2*x1^4*x2^2, x1^3*x2^5, x0*x1^3*x2^4, x0^2*x1^3*x2^3, x0^3*x1^3*x2^2, x0^4*x1^3*x2, x0^2*x1^2*x2^4, x0^3*x1^2*x2^3, x0^4*x1^2*x2^2, x0^5*x1^2*x2, x0^6*x1^2, x0^4*x1*x2^3, x0^5*x1*x2^2, x0^5*x2^3] 
usecase_B2(423): 16 [y0^5*y2^2*y3^3, y0^4*y1*y2^3*y3^2, y0^4*y1*y2^2*y3^3, y0^3*y1^2*y2^5, y0^3*y1^2*y2^4*y3, y0^3*y1^2*y2^3*y3^2, y0^3*y1^2*y2^2*y3^3, y0^3*y1^2*y2*y3^4, y0^2*y1^3*y2^4*y3, y0^2*y1^3*y2^3*y3^2, y0^2*y1^3*y2^2*y3^3, y0^2*y1^3*y2*y3^4, y0^2*y1^3*y3^5, y0*y1^4*y2^3*y3^2, y0*y1^4*y2^2*y3^3, y1^5*y2^3*y3^2] 
```

Next we compute the r1 reductions of the above output maps.  

```python
hft = BasePointTree()
hft.add( 'x', ( 0, 0 ) , 2 ).add( 't', ( 0, 0 ), 1 )
hft.add( 'y', ( 0, 0 ) , 2 ).add( 't', ( 0, 0 ), 1 )
hft.add( 'z', ( 0, 0 ) , 1 )
hf = LinearSeries.get( [5], hft )
hgt = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
hgt.add( 'xv', ( 0, 0 ) , 1 )
hgt.add( 'xw', ( 0, 0 ) , 1 )
hgt.add( 'yv', ( 0, 0 ) , 1 )
hgt.add( 'yw', ( 0, 0 ) , 1 )
hg = LinearSeries.get( [3, 3], hgt )
SETools.p( len( hf.pol_lst ), SERing.conv( hf.pol_lst ) )
SETools.p( len( hg.pol_lst ), SERing.conv( hg.pol_lst ) )
assert len( hf.pol_lst ) == len( hg.pol_lst ) == 12

```
Output:
```
12 [x1^3*x2^2, x0*x1^3*x2, x1^2*x2^3, x0*x1^2*x2^2, x0^2*x1^2*x2, x0^3*x1^2, x0*x1*x2^3, x0^2*x1*x2^2, x0^3*x1*x2, x0^4*x1, x0^3*x2^2, x0^4*x2] 
12 [y0^3*y2^2*y3, y0^3*y2*y3^2, y0^2*y1*y2^3, y0^2*y1*y2^2*y3, y0^2*y1*y2*y3^2, y0^2*y1*y3^3, y0*y1^2*y2^3, y0*y1^2*y2^2*y3, y0*y1^2*y2*y3^2, y0*y1^2*y3^3, y1^3*y2^2*y3, y1^3*y2*y3^2] 
```
We again compute the r1 reductions of the newly obtained rational maps.  

```python
hft = BasePointTree()
hft.add( 'x', ( 0, 0 ) , 1 )
hft.add( 'y', ( 0, 0 ) , 1 )
hf = LinearSeries.get( [2], hft )
hg = LinearSeries( ['x', 'y', 'v', 'w' ], PolyRing( 'x,y,v,w' ) )
SETools.p( len( hf.pol_lst ), SERing.conv( hf.pol_lst ) )
SETools.p( len( hg.pol_lst ), SERing.conv( hg.pol_lst ) )
assert len( hf.pol_lst ) == len( hg.pol_lst ) == 4

```
Output:
```
4 [x1*x2, x0*x1, x0*x2, x0^2] 
4 [y0, y1, y2, y3] 
```
    
The above output maps are the (r0,r1,r2)-reductions of the rational maps f and g.    
    
#### Declare method for enforcing basepoints

We declare a method which takes as input a map gr with parameters c
and outputs conditions on c such that gr has the same basepoints as f (see `bf`):

```python
# use "%paste" to paste this code in a Sage or Python terminal.
def usecase_B2_helper_bp( gr ):
    
    eqn_lst = []

    #
    #     chart=z, depth=0, mult=2, sol=(0, 0), { 4, <<x^2, x^5*y^2, x^3*y^5, x*y^2 + 2*y^3>>, QQ[x, y] }
    #         chart=t, depth=1, mult=1, sol=(0, 0), { 4, <<x^2, x^5*y^5, x^3*y^6, x*y + 2*y>>, QQ[x, y] }
    #             chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<x, x^9*y^5, x^8*y^6, x*y + 2*y>>, QQ[x, y] }
    #
    gr0 = [ comp.subs( {x[0]:1} ) for comp in gr ]
    gr0t = [ comp.subs( {x[1]:x[1] * x[2]} ).quo_rem( x[2] ** 2 )[0] for comp in gr0 ]
    gr0ts = [ comp.subs( {x[2]:x[1] * x[2]} ).quo_rem( x[1] )[0] for comp in gr0t ]
    for a, b in [( 0, 0 ), ( 0, 1 ), ( 1, 0 )]:
        eqn_lst += [ sage_diff( comp, x[1], a, x[2], b ).subs( {x[1]:0, x[2]:0} ) for comp in gr0 ]
    eqn_lst += [ comp.subs( {x[1]:0, x[2]:0} ) for comp in gr0t ]
    eqn_lst += [ comp.subs( {x[1]:0, x[2]:0} ) for comp in gr0ts ]

    #
    #     chart=x, depth=0, mult=3, sol=(0, 0), { 4, <<z^6, y^2*z, y^5, 2*y^3*z^5 + y^2*z^5>>, QQ[y, z] }
    #         chart=t, depth=1, mult=2, sol=(0, 0), { 4, <<z^3, y^2, y^5*z^2, 2*y^3*z^5 + y^2*z^4>>, QQ[y, z] }
    #             chart=t, depth=2, mult=1, sol=(0, 0), { 4, <<z, y^2, y^5*z^5, 2*y^3*z^6 + y^2*z^4>>, QQ[y, z] }
    #                 chart=s, depth=3, mult=1, sol=(0, 0), { 4, <<z, y, y^9*z^5, 2*y^8*z^6 + y^5*z^4>>, QQ[y, z] }
    #         chart=s, depth=1, mult=1, sol=(0, 0), { 4, <<y^3*z^6, z, y^2, 2*y^5*z^5 + y^4*z^5>>, QQ[y, z] }
    #             chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<y^8*z^6, z, y, 2*y^9*z^5 + y^8*z^5>>, QQ[y, z] }
    #
    gr1 = [ comp.subs( {x[1]:1} ) for comp in gr ]
    gr1t = [ comp.subs( {x[0]:x[0] * x[2]} ).quo_rem( x[2] ** 3 )[0] for comp in gr1 ]
    gr1tt = [ comp.subs( {x[0]:x[0] * x[2]} ).quo_rem( x[2] ** 2 )[0] for comp in gr1t ]
    gr1tts = [ comp.subs( {x[2]:x[0] * x[2]} ).quo_rem( x[0] ** 1 )[0] for comp in gr1tt ]
    gr1s = [ comp.subs( {x[2]:x[0] * x[2]} ).quo_rem( x[0] ** 3 )[0] for comp in gr1 ]
    gr1ss = [ comp.subs( {x[2]:x[0] * x[2]} ).quo_rem( x[0] ** 1 )[0] for comp in gr1s ]
    for a, b in [( 0, 0 ), ( 0, 1 ), ( 1, 0 ), ( 2, 0 ), ( 1, 1 ), ( 0, 2 )]:
        eqn_lst += [ sage_diff( comp, x[0], a, x[2], b ).subs( {x[0]:0, x[2]:0} ) for comp in gr1 ]
    for a, b in [( 0, 0 ), ( 0, 1 ), ( 1, 0 )]:
        eqn_lst += [ sage_diff( comp, x[0], a, x[2], b ).subs( {x[0]:0, x[2]:0} ) for comp in gr1t ]
    eqn_lst += [ comp.subs( {x[0]:0, x[2]:0} ) for comp in gr1tt ]
    eqn_lst += [ comp.subs( {x[0]:0, x[2]:0} ) for comp in gr1tts ]
    eqn_lst += [ comp.subs( {x[0]:0, x[2]:0} ) for comp in gr1s ]
    eqn_lst += [ comp.subs( {x[0]:0, x[2]:0} ) for comp in gr1ss ]

    #
    #     chart=y, depth=0, mult=3, sol=(0, 0), { 4, <<x^2*z^6, x^5*z, x^3, x*z^5 + 2*z^5>>, QQ[x, z] }
    #         chart=t, depth=1, mult=2, sol=(0, 0), { 4, <<x^2*z^5, x^5*z^3, x^3, x*z^3 + 2*z^2>>, QQ[x, z] }
    #             chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<x^5*z^5, x^6*z^3, x, x^2*z^3 + 2*z^2>>, QQ[x, z] }
    #                 chart=t, depth=3, mult=1, sol=(0, 0), { 4, <<x^5*z^9, x^6*z^8, x, x^2*z^4 + 2*z>>, QQ[x, z] }
    #
    gr2 = [ comp.subs( {x[2]:1} ) for comp in gr ]
    gr2t = [ comp.subs( {x[0]:x[0] * x[1]} ).quo_rem( x[1] ** 3 )[0] for comp in gr2 ]
    gr2ts = [ comp.subs( {x[1]:x[0] * x[1]} ).quo_rem( x[0] ** 2 )[0] for comp in gr2t ]
    gr2tst = [ comp.subs( {x[0]:x[0] * x[1]} ).quo_rem( x[1] ** 1 )[0] for comp in gr2ts ]
    for a, b in [( 0, 0 ), ( 0, 1 ), ( 1, 0 ), ( 2, 0 ), ( 1, 1 ), ( 0, 2 )]:
        eqn_lst += [ sage_diff( comp, x[0], a, x[1], b ).subs( {x[0]:0, x[1]:0} ) for comp in gr2 ]
    for a, b in [( 0, 0 ), ( 0, 1 ), ( 1, 0 )]:
        eqn_lst += [ sage_diff( comp, x[0], a, x[1], b ).subs( {x[0]:0, x[1]:0} ) for comp in gr2t ]
    eqn_lst += [ comp.subs( {x[0]:0, x[1]:0} ) for comp in gr2ts ]
    eqn_lst += [ comp.subs( {x[0]:0, x[1]:0} ) for comp in gr2tst ]

    eqn_lst = sorted( list( set( eqn_lst ) ) )

    assert 'x0' not in str( eqn_lst )
    assert 'x1' not in str( eqn_lst )
    assert 'x2' not in str( eqn_lst )

    return eqn_lst
```
The superset of the compatible reparametrizations consist of two components r0 and r1.

```python
# we compute maps to P1xP1 from two pencils
PolyRing.reset_base_field()
bpt = BasePointTree()
bpt.add( 'y', ( 0, 0 ) , 1 )
pen1 = SERing.conv( LinearSeries.get( [1], bpt ).pol_lst )
assert set(pen1) == set(ring('[x1, x0]'))
assert set( [x[0], x[1]] ) == set( pen1 )
# thus the first pencil defines a map pen1: (x0:x1:x2) |--> (x0:x1)
bpt = BasePointTree()
bpt.add( 'x', ( 0, 0 ) , 1 )
pen2 = SERing.conv( LinearSeries.get( [1], bpt ).pol_lst )
assert set(pen2) == set(ring('[x2, x0]'))
assert set( [x[0], x[2]] ) == set( pen2 )
# thus the second pencil defines a map pen2: (x0:x1:x2) |--> (x0:x2)
# We find that
#     pen1 x pen2: P2-->P1xP1, (x0:x1:x2) |--> (x0:x1;x0:x2) and
#     pen2 x pen1: P2-->P1xP1, (x0:x1:x2) |--> (x0:x2;x0:x1)
# We obtain the following compatible reparametrizations
# by composing the maps pen1 x pen2 and pen2 x pen1
# with a parametrized map in the identity component of Aut(P1xP1).
r0 = {y[0]:c[0] * x[0] + c[1] * x[1], y[1]:c[2] * x[0] + c[3] * x[1], y[2]:c[4] * x[0] + c[5] * x[2], y[3]:c[6] * x[0] + c[7] * x[2]}
r1 = {y[0]:c[0] * x[0] + c[1] * x[2], y[1]:c[2] * x[0] + c[3] * x[2], y[2]:c[4] * x[0] + c[5] * x[1], y[3]:c[6] * x[0] + c[7] * x[1]}

```
We compute the coefficient matrix Mf for f and its kernel Kf.

```python
Mf = SERing.get_matrix_P2( f )
Kf = Mf.right_kernel_matrix().T
assert ( Mf * Kf ).is_zero()
assert Mf.dimensions() == (4, 45)
assert Kf.dimensions() == (45, 41)

```

#### Compatible reparametrization r0

We first consider the reparametrization r0 and we find the coefficient matrix Mgr00.

```python
# compose g with reparametrization r0
gcd0 = sage_gcd( [ comp.subs( r0 ) for comp in g ] )
assert gcd0 == 1
gr0 = [ comp.subs( r0 ) / gcd0 for comp in g ]
print( 'gr0 =', len( gr0 ), gcd0, gr0 )
assert SERing.get_degree( gr0 ) == 10
assert SERing.get_degree( f ) == 8

# find conditions on c so that gr0 has the same basepoints as f
eqn0_lst = usecase_B2_helper_bp( gr0 )
eqn0_lst += [ ring( '(c0*c3-c1*c2)*(c4*c7-c5*c6)*t-1' ) ]
prime0_lst = sage_ideal( eqn0_lst ).elimination_ideal( ring( 't' ) ).primary_decomposition()
print( 'eqn0_lst =', len( eqn0_lst ), eqn0_lst )
for prime0 in prime0_lst: print( '\t', prime0.gens() )
sol00 = {c[1]:0, c[2]:0, c[5]:0, c[6]:0, c[0]:1, c[4]:1}  # notice that wlog c0=c4=1
sol01 = {c[0]:0, c[3]:0, c[4]:0, c[7]:0, c[1]:1, c[5]:1}  # notice that wlog c1=c5=1
assert len( prime0_lst ) == 2
assert set( [gen.subs( sol00 ) for gen in prime0_lst[0].gens()] ) == set( [0] )
assert set( [gen.subs( sol01 ) for gen in prime0_lst[1].gens()] ) == set( [0] )

# sol00: notice that c3!=0 and c7!=0
gcd00 = sage_gcd( [ comp.subs( sol00 ) for comp in gr0] )
assert gcd00 == x[0] * x[0]
gr00 = [ comp.subs( sol00 ) / gcd00 for comp in gr0]
print( 'gr00 =', len( gr00 ), gcd00, gr00 )
assert SERing.get_degree( gr00 ) == 8
Mgr00 = SERing.get_matrix_P2( gr00 )
assert Mgr00.dimensions() == ( 4, 45 )
# find conditions for c so that Mgr00 has the same kernel as the matrix of f
p00_lst = sage_ideal( ( Mgr00 * Kf ).list() + [ring( 'c3*c7*t-1' )] ).elimination_ideal( ring( 't' ) ).primary_decomposition()
assert [p00.gens() for p00 in p00_lst] == [[2 * c[3] - c[7]]]
Mgr00 = Mgr00.subs( {c[7]:2 * c[3]} )
print( 'Mgr00 =', Mgr00.dimensions(), list( Mgr00 ) )
# found a solution: Mgr00

# sol01: notice that c2!=0 and c6!=0
gcd01 = sage_gcd( [ comp.subs( sol01 ) for comp in gr0] )
assert gcd01 == x[0] * x[0]
gr01 = [ comp.subs( sol01 ) / gcd01 for comp in gr0]
print( 'gr01 =', len( gr01 ), gcd01, gr01 )
assert SERing.get_degree( gr01 ) == 8
assert [] == sage_ideal( ( SERing.get_matrix_P2( gr01 ) * Kf ).list() + [ring( 'c2*c6*t-1' )] ).elimination_ideal( ring( 't' ) ).primary_decomposition()
# --> no solution

```
Output:
```
gr0 = 4 1 [c0^3*c2^2*c4^5*x0^10 + 3*c0^2*c1*c2^2*c4^5*x0^9*x1 + 2*c0^3*c2*c3*c4^5*x0^9*x1 + 3*c0*c1^2*c2^2*c4^5*x0^8*x1^2 + 6*c0^2*c1*c2*c3*c4^5*x0^8*x1^2 + c0^3*c3^2*c4^5*x0^8*x1^2 + c1^3*c2^2*c4^5*x0^7*x1^3 + 6*c0*c1^2*c2*c3*c4^5*x0^7*x1^3 + 3*c0^2*c1*c3^2*c4^5*x0^7*x1^3 + 2*c1^3*c2*c3*c4^5*x0^6*x1^4 + 3*c0*c1^2*c3^2*c4^5*x0^6*x1^4 + c1^3*c3^2*c4^5*x0^5*x1^5 + 5*c0^3*c2^2*c4^4*c5*x0^9*x2 + 15*c0^2*c1*c2^2*c4^4*c5*x0^8*x1*x2 + 10*c0^3*c2*c3*c4^4*c5*x0^8*x1*x2 + 15*c0*c1^2*c2^2*c4^4*c5*x0^7*x1^2*x2 + 30*c0^2*c1*c2*c3*c4^4*c5*x0^7*x1^2*x2 + 5*c0^3*c3^2*c4^4*c5*x0^7*x1^2*x2 + 5*c1^3*c2^2*c4^4*c5*x0^6*x1^3*x2 + 30*c0*c1^2*c2*c3*c4^4*c5*x0^6*x1^3*x2 + 15*c0^2*c1*c3^2*c4^4*c5*x0^6*x1^3*x2 + 10*c1^3*c2*c3*c4^4*c5*x0^5*x1^4*x2 + 15*c0*c1^2*c3^2*c4^4*c5*x0^5*x1^4*x2 + 5*c1^3*c3^2*c4^4*c5*x0^4*x1^5*x2 + 10*c0^3*c2^2*c4^3*c5^2*x0^8*x2^2 + 30*c0^2*c1*c2^2*c4^3*c5^2*x0^7*x1*x2^2 + 20*c0^3*c2*c3*c4^3*c5^2*x0^7*x1*x2^2 + 30*c0*c1^2*c2^2*c4^3*c5^2*x0^6*x1^2*x2^2 + 60*c0^2*c1*c2*c3*c4^3*c5^2*x0^6*x1^2*x2^2 + 10*c0^3*c3^2*c4^3*c5^2*x0^6*x1^2*x2^2 + 10*c1^3*c2^2*c4^3*c5^2*x0^5*x1^3*x2^2 + 60*c0*c1^2*c2*c3*c4^3*c5^2*x0^5*x1^3*x2^2 + 30*c0^2*c1*c3^2*c4^3*c5^2*x0^5*x1^3*x2^2 + 20*c1^3*c2*c3*c4^3*c5^2*x0^4*x1^4*x2^2 + 30*c0*c1^2*c3^2*c4^3*c5^2*x0^4*x1^4*x2^2 + 10*c1^3*c3^2*c4^3*c5^2*x0^3*x1^5*x2^2 + 10*c0^3*c2^2*c4^2*c5^3*x0^7*x2^3 + 30*c0^2*c1*c2^2*c4^2*c5^3*x0^6*x1*x2^3 + 20*c0^3*c2*c3*c4^2*c5^3*x0^6*x1*x2^3 + 30*c0*c1^2*c2^2*c4^2*c5^3*x0^5*x1^2*x2^3 + 60*c0^2*c1*c2*c3*c4^2*c5^3*x0^5*x1^2*x2^3 + 10*c0^3*c3^2*c4^2*c5^3*x0^5*x1^2*x2^3 + 10*c1^3*c2^2*c4^2*c5^3*x0^4*x1^3*x2^3 + 60*c0*c1^2*c2*c3*c4^2*c5^3*x0^4*x1^3*x2^3 + 30*c0^2*c1*c3^2*c4^2*c5^3*x0^4*x1^3*x2^3 + 20*c1^3*c2*c3*c4^2*c5^3*x0^3*x1^4*x2^3 + 30*c0*c1^2*c3^2*c4^2*c5^3*x0^3*x1^4*x2^3 + 10*c1^3*c3^2*c4^2*c5^3*x0^2*x1^5*x2^3 + 5*c0^3*c2^2*c4*c5^4*x0^6*x2^4 + 15*c0^2*c1*c2^2*c4*c5^4*x0^5*x1*x2^4 + 10*c0^3*c2*c3*c4*c5^4*x0^5*x1*x2^4 + 15*c0*c1^2*c2^2*c4*c5^4*x0^4*x1^2*x2^4 + 30*c0^2*c1*c2*c3*c4*c5^4*x0^4*x1^2*x2^4 + 5*c0^3*c3^2*c4*c5^4*x0^4*x1^2*x2^4 + 5*c1^3*c2^2*c4*c5^4*x0^3*x1^3*x2^4 + 30*c0*c1^2*c2*c3*c4*c5^4*x0^3*x1^3*x2^4 + 15*c0^2*c1*c3^2*c4*c5^4*x0^3*x1^3*x2^4 + 10*c1^3*c2*c3*c4*c5^4*x0^2*x1^4*x2^4 + 15*c0*c1^2*c3^2*c4*c5^4*x0^2*x1^4*x2^4 + 5*c1^3*c3^2*c4*c5^4*x0*x1^5*x2^4 + c0^3*c2^2*c5^5*x0^5*x2^5 + 3*c0^2*c1*c2^2*c5^5*x0^4*x1*x2^5 + 2*c0^3*c2*c3*c5^5*x0^4*x1*x2^5 + 3*c0*c1^2*c2^2*c5^5*x0^3*x1^2*x2^5 + 6*c0^2*c1*c2*c3*c5^5*x0^3*x1^2*x2^5 + c0^3*c3^2*c5^5*x0^3*x1^2*x2^5 + c1^3*c2^2*c5^5*x0^2*x1^3*x2^5 + 6*c0*c1^2*c2*c3*c5^5*x0^2*x1^3*x2^5 + 3*c0^2*c1*c3^2*c5^5*x0^2*x1^3*x2^5 + 2*c1^3*c2*c3*c5^5*x0*x1^4*x2^5 + 3*c0*c1^2*c3^2*c5^5*x0*x1^4*x2^5 + c1^3*c3^2*c5^5*x1^5*x2^5, c0^3*c2^2*c4^5*x0^10 + c2^5*c4^3*c6^2*x0^10 + 3*c0^2*c1*c2^2*c4^5*x0^9*x1 + 2*c0^3*c2*c3*c4^5*x0^9*x1 + 5*c2^4*c3*c4^3*c6^2*x0^9*x1 + 3*c0*c1^2*c2^2*c4^5*x0^8*x1^2 + 6*c0^2*c1*c2*c3*c4^5*x0^8*x1^2 + c0^3*c3^2*c4^5*x0^8*x1^2 + 10*c2^3*c3^2*c4^3*c6^2*x0^8*x1^2 + c1^3*c2^2*c4^5*x0^7*x1^3 + 6*c0*c1^2*c2*c3*c4^5*x0^7*x1^3 + 3*c0^2*c1*c3^2*c4^5*x0^7*x1^3 + 10*c2^2*c3^3*c4^3*c6^2*x0^7*x1^3 + 2*c1^3*c2*c3*c4^5*x0^6*x1^4 + 3*c0*c1^2*c3^2*c4^5*x0^6*x1^4 + 5*c2*c3^4*c4^3*c6^2*x0^6*x1^4 + c1^3*c3^2*c4^5*x0^5*x1^5 + c3^5*c4^3*c6^2*x0^5*x1^5 + 5*c0^3*c2^2*c4^4*c5*x0^9*x2 + 3*c2^5*c4^2*c5*c6^2*x0^9*x2 + 2*c2^5*c4^3*c6*c7*x0^9*x2 + 15*c0^2*c1*c2^2*c4^4*c5*x0^8*x1*x2 + 10*c0^3*c2*c3*c4^4*c5*x0^8*x1*x2 + 15*c2^4*c3*c4^2*c5*c6^2*x0^8*x1*x2 + 10*c2^4*c3*c4^3*c6*c7*x0^8*x1*x2 + 15*c0*c1^2*c2^2*c4^4*c5*x0^7*x1^2*x2 + 30*c0^2*c1*c2*c3*c4^4*c5*x0^7*x1^2*x2 + 5*c0^3*c3^2*c4^4*c5*x0^7*x1^2*x2 + 30*c2^3*c3^2*c4^2*c5*c6^2*x0^7*x1^2*x2 + 20*c2^3*c3^2*c4^3*c6*c7*x0^7*x1^2*x2 + 5*c1^3*c2^2*c4^4*c5*x0^6*x1^3*x2 + 30*c0*c1^2*c2*c3*c4^4*c5*x0^6*x1^3*x2 + 15*c0^2*c1*c3^2*c4^4*c5*x0^6*x1^3*x2 + 30*c2^2*c3^3*c4^2*c5*c6^2*x0^6*x1^3*x2 + 20*c2^2*c3^3*c4^3*c6*c7*x0^6*x1^3*x2 + 10*c1^3*c2*c3*c4^4*c5*x0^5*x1^4*x2 + 15*c0*c1^2*c3^2*c4^4*c5*x0^5*x1^4*x2 + 15*c2*c3^4*c4^2*c5*c6^2*x0^5*x1^4*x2 + 10*c2*c3^4*c4^3*c6*c7*x0^5*x1^4*x2 + 5*c1^3*c3^2*c4^4*c5*x0^4*x1^5*x2 + 3*c3^5*c4^2*c5*c6^2*x0^4*x1^5*x2 + 2*c3^5*c4^3*c6*c7*x0^4*x1^5*x2 + 10*c0^3*c2^2*c4^3*c5^2*x0^8*x2^2 + 3*c2^5*c4*c5^2*c6^2*x0^8*x2^2 + 6*c2^5*c4^2*c5*c6*c7*x0^8*x2^2 + c2^5*c4^3*c7^2*x0^8*x2^2 + 30*c0^2*c1*c2^2*c4^3*c5^2*x0^7*x1*x2^2 + 20*c0^3*c2*c3*c4^3*c5^2*x0^7*x1*x2^2 + 15*c2^4*c3*c4*c5^2*c6^2*x0^7*x1*x2^2 + 30*c2^4*c3*c4^2*c5*c6*c7*x0^7*x1*x2^2 + 5*c2^4*c3*c4^3*c7^2*x0^7*x1*x2^2 + 30*c0*c1^2*c2^2*c4^3*c5^2*x0^6*x1^2*x2^2 + 60*c0^2*c1*c2*c3*c4^3*c5^2*x0^6*x1^2*x2^2 + 10*c0^3*c3^2*c4^3*c5^2*x0^6*x1^2*x2^2 + 30*c2^3*c3^2*c4*c5^2*c6^2*x0^6*x1^2*x2^2 + 60*c2^3*c3^2*c4^2*c5*c6*c7*x0^6*x1^2*x2^2 + 10*c2^3*c3^2*c4^3*c7^2*x0^6*x1^2*x2^2 + 10*c1^3*c2^2*c4^3*c5^2*x0^5*x1^3*x2^2 + 60*c0*c1^2*c2*c3*c4^3*c5^2*x0^5*x1^3*x2^2 + 30*c0^2*c1*c3^2*c4^3*c5^2*x0^5*x1^3*x2^2 + 30*c2^2*c3^3*c4*c5^2*c6^2*x0^5*x1^3*x2^2 + 60*c2^2*c3^3*c4^2*c5*c6*c7*x0^5*x1^3*x2^2 + 10*c2^2*c3^3*c4^3*c7^2*x0^5*x1^3*x2^2 + 20*c1^3*c2*c3*c4^3*c5^2*x0^4*x1^4*x2^2 + 30*c0*c1^2*c3^2*c4^3*c5^2*x0^4*x1^4*x2^2 + 15*c2*c3^4*c4*c5^2*c6^2*x0^4*x1^4*x2^2 + 30*c2*c3^4*c4^2*c5*c6*c7*x0^4*x1^4*x2^2 + 5*c2*c3^4*c4^3*c7^2*x0^4*x1^4*x2^2 + 10*c1^3*c3^2*c4^3*c5^2*x0^3*x1^5*x2^2 + 3*c3^5*c4*c5^2*c6^2*x0^3*x1^5*x2^2 + 6*c3^5*c4^2*c5*c6*c7*x0^3*x1^5*x2^2 + c3^5*c4^3*c7^2*x0^3*x1^5*x2^2 + 10*c0^3*c2^2*c4^2*c5^3*x0^7*x2^3 + c2^5*c5^3*c6^2*x0^7*x2^3 + 6*c2^5*c4*c5^2*c6*c7*x0^7*x2^3 + 3*c2^5*c4^2*c5*c7^2*x0^7*x2^3 + 30*c0^2*c1*c2^2*c4^2*c5^3*x0^6*x1*x2^3 + 20*c0^3*c2*c3*c4^2*c5^3*x0^6*x1*x2^3 + 5*c2^4*c3*c5^3*c6^2*x0^6*x1*x2^3 + 30*c2^4*c3*c4*c5^2*c6*c7*x0^6*x1*x2^3 + 15*c2^4*c3*c4^2*c5*c7^2*x0^6*x1*x2^3 + 30*c0*c1^2*c2^2*c4^2*c5^3*x0^5*x1^2*x2^3 + 60*c0^2*c1*c2*c3*c4^2*c5^3*x0^5*x1^2*x2^3 + 10*c0^3*c3^2*c4^2*c5^3*x0^5*x1^2*x2^3 + 10*c2^3*c3^2*c5^3*c6^2*x0^5*x1^2*x2^3 + 60*c2^3*c3^2*c4*c5^2*c6*c7*x0^5*x1^2*x2^3 + 30*c2^3*c3^2*c4^2*c5*c7^2*x0^5*x1^2*x2^3 + 10*c1^3*c2^2*c4^2*c5^3*x0^4*x1^3*x2^3 + 60*c0*c1^2*c2*c3*c4^2*c5^3*x0^4*x1^3*x2^3 + 30*c0^2*c1*c3^2*c4^2*c5^3*x0^4*x1^3*x2^3 + 10*c2^2*c3^3*c5^3*c6^2*x0^4*x1^3*x2^3 + 60*c2^2*c3^3*c4*c5^2*c6*c7*x0^4*x1^3*x2^3 + 30*c2^2*c3^3*c4^2*c5*c7^2*x0^4*x1^3*x2^3 + 20*c1^3*c2*c3*c4^2*c5^3*x0^3*x1^4*x2^3 + 30*c0*c1^2*c3^2*c4^2*c5^3*x0^3*x1^4*x2^3 + 5*c2*c3^4*c5^3*c6^2*x0^3*x1^4*x2^3 + 30*c2*c3^4*c4*c5^2*c6*c7*x0^3*x1^4*x2^3 + 15*c2*c3^4*c4^2*c5*c7^2*x0^3*x1^4*x2^3 + 10*c1^3*c3^2*c4^2*c5^3*x0^2*x1^5*x2^3 + c3^5*c5^3*c6^2*x0^2*x1^5*x2^3 + 6*c3^5*c4*c5^2*c6*c7*x0^2*x1^5*x2^3 + 3*c3^5*c4^2*c5*c7^2*x0^2*x1^5*x2^3 + 5*c0^3*c2^2*c4*c5^4*x0^6*x2^4 + 2*c2^5*c5^3*c6*c7*x0^6*x2^4 + 3*c2^5*c4*c5^2*c7^2*x0^6*x2^4 + 15*c0^2*c1*c2^2*c4*c5^4*x0^5*x1*x2^4 + 10*c0^3*c2*c3*c4*c5^4*x0^5*x1*x2^4 + 10*c2^4*c3*c5^3*c6*c7*x0^5*x1*x2^4 + 15*c2^4*c3*c4*c5^2*c7^2*x0^5*x1*x2^4 + 15*c0*c1^2*c2^2*c4*c5^4*x0^4*x1^2*x2^4 + 30*c0^2*c1*c2*c3*c4*c5^4*x0^4*x1^2*x2^4 + 5*c0^3*c3^2*c4*c5^4*x0^4*x1^2*x2^4 + 20*c2^3*c3^2*c5^3*c6*c7*x0^4*x1^2*x2^4 + 30*c2^3*c3^2*c4*c5^2*c7^2*x0^4*x1^2*x2^4 + 5*c1^3*c2^2*c4*c5^4*x0^3*x1^3*x2^4 + 30*c0*c1^2*c2*c3*c4*c5^4*x0^3*x1^3*x2^4 + 15*c0^2*c1*c3^2*c4*c5^4*x0^3*x1^3*x2^4 + 20*c2^2*c3^3*c5^3*c6*c7*x0^3*x1^3*x2^4 + 30*c2^2*c3^3*c4*c5^2*c7^2*x0^3*x1^3*x2^4 + 10*c1^3*c2*c3*c4*c5^4*x0^2*x1^4*x2^4 + 15*c0*c1^2*c3^2*c4*c5^4*x0^2*x1^4*x2^4 + 10*c2*c3^4*c5^3*c6*c7*x0^2*x1^4*x2^4 + 15*c2*c3^4*c4*c5^2*c7^2*x0^2*x1^4*x2^4 + 5*c1^3*c3^2*c4*c5^4*x0*x1^5*x2^4 + 2*c3^5*c5^3*c6*c7*x0*x1^5*x2^4 + 3*c3^5*c4*c5^2*c7^2*x0*x1^5*x2^4 + c0^3*c2^2*c5^5*x0^5*x2^5 + c2^5*c5^3*c7^2*x0^5*x2^5 + 3*c0^2*c1*c2^2*c5^5*x0^4*x1*x2^5 + 2*c0^3*c2*c3*c5^5*x0^4*x1*x2^5 + 5*c2^4*c3*c5^3*c7^2*x0^4*x1*x2^5 + 3*c0*c1^2*c2^2*c5^5*x0^3*x1^2*x2^5 + 6*c0^2*c1*c2*c3*c5^5*x0^3*x1^2*x2^5 + c0^3*c3^2*c5^5*x0^3*x1^2*x2^5 + 10*c2^3*c3^2*c5^3*c7^2*x0^3*x1^2*x2^5 + c1^3*c2^2*c5^5*x0^2*x1^3*x2^5 + 6*c0*c1^2*c2*c3*c5^5*x0^2*x1^3*x2^5 + 3*c0^2*c1*c3^2*c5^5*x0^2*x1^3*x2^5 + 10*c2^2*c3^3*c5^3*c7^2*x0^2*x1^3*x2^5 + 2*c1^3*c2*c3*c5^5*x0*x1^4*x2^5 + 3*c0*c1^2*c3^2*c5^5*x0*x1^4*x2^5 + 5*c2*c3^4*c5^3*c7^2*x0*x1^4*x2^5 + c1^3*c3^2*c5^5*x1^5*x2^5 + c3^5*c5^3*c7^2*x1^5*x2^5, c0^2*c2^3*c6^5*x0^10 + 2*c0*c1*c2^3*c6^5*x0^9*x1 + 3*c0^2*c2^2*c3*c6^5*x0^9*x1 + c1^2*c2^3*c6^5*x0^8*x1^2 + 6*c0*c1*c2^2*c3*c6^5*x0^8*x1^2 + 3*c0^2*c2*c3^2*c6^5*x0^8*x1^2 + 3*c1^2*c2^2*c3*c6^5*x0^7*x1^3 + 6*c0*c1*c2*c3^2*c6^5*x0^7*x1^3 + c0^2*c3^3*c6^5*x0^7*x1^3 + 3*c1^2*c2*c3^2*c6^5*x0^6*x1^4 + 2*c0*c1*c3^3*c6^5*x0^6*x1^4 + c1^2*c3^3*c6^5*x0^5*x1^5 + 5*c0^2*c2^3*c6^4*c7*x0^9*x2 + 10*c0*c1*c2^3*c6^4*c7*x0^8*x1*x2 + 15*c0^2*c2^2*c3*c6^4*c7*x0^8*x1*x2 + 5*c1^2*c2^3*c6^4*c7*x0^7*x1^2*x2 + 30*c0*c1*c2^2*c3*c6^4*c7*x0^7*x1^2*x2 + 15*c0^2*c2*c3^2*c6^4*c7*x0^7*x1^2*x2 + 15*c1^2*c2^2*c3*c6^4*c7*x0^6*x1^3*x2 + 30*c0*c1*c2*c3^2*c6^4*c7*x0^6*x1^3*x2 + 5*c0^2*c3^3*c6^4*c7*x0^6*x1^3*x2 + 15*c1^2*c2*c3^2*c6^4*c7*x0^5*x1^4*x2 + 10*c0*c1*c3^3*c6^4*c7*x0^5*x1^4*x2 + 5*c1^2*c3^3*c6^4*c7*x0^4*x1^5*x2 + 10*c0^2*c2^3*c6^3*c7^2*x0^8*x2^2 + 20*c0*c1*c2^3*c6^3*c7^2*x0^7*x1*x2^2 + 30*c0^2*c2^2*c3*c6^3*c7^2*x0^7*x1*x2^2 + 10*c1^2*c2^3*c6^3*c7^2*x0^6*x1^2*x2^2 + 60*c0*c1*c2^2*c3*c6^3*c7^2*x0^6*x1^2*x2^2 + 30*c0^2*c2*c3^2*c6^3*c7^2*x0^6*x1^2*x2^2 + 30*c1^2*c2^2*c3*c6^3*c7^2*x0^5*x1^3*x2^2 + 60*c0*c1*c2*c3^2*c6^3*c7^2*x0^5*x1^3*x2^2 + 10*c0^2*c3^3*c6^3*c7^2*x0^5*x1^3*x2^2 + 30*c1^2*c2*c3^2*c6^3*c7^2*x0^4*x1^4*x2^2 + 20*c0*c1*c3^3*c6^3*c7^2*x0^4*x1^4*x2^2 + 10*c1^2*c3^3*c6^3*c7^2*x0^3*x1^5*x2^2 + 10*c0^2*c2^3*c6^2*c7^3*x0^7*x2^3 + 20*c0*c1*c2^3*c6^2*c7^3*x0^6*x1*x2^3 + 30*c0^2*c2^2*c3*c6^2*c7^3*x0^6*x1*x2^3 + 10*c1^2*c2^3*c6^2*c7^3*x0^5*x1^2*x2^3 + 60*c0*c1*c2^2*c3*c6^2*c7^3*x0^5*x1^2*x2^3 + 30*c0^2*c2*c3^2*c6^2*c7^3*x0^5*x1^2*x2^3 + 30*c1^2*c2^2*c3*c6^2*c7^3*x0^4*x1^3*x2^3 + 60*c0*c1*c2*c3^2*c6^2*c7^3*x0^4*x1^3*x2^3 + 10*c0^2*c3^3*c6^2*c7^3*x0^4*x1^3*x2^3 + 30*c1^2*c2*c3^2*c6^2*c7^3*x0^3*x1^4*x2^3 + 20*c0*c1*c3^3*c6^2*c7^3*x0^3*x1^4*x2^3 + 10*c1^2*c3^3*c6^2*c7^3*x0^2*x1^5*x2^3 + 5*c0^2*c2^3*c6*c7^4*x0^6*x2^4 + 10*c0*c1*c2^3*c6*c7^4*x0^5*x1*x2^4 + 15*c0^2*c2^2*c3*c6*c7^4*x0^5*x1*x2^4 + 5*c1^2*c2^3*c6*c7^4*x0^4*x1^2*x2^4 + 30*c0*c1*c2^2*c3*c6*c7^4*x0^4*x1^2*x2^4 + 15*c0^2*c2*c3^2*c6*c7^4*x0^4*x1^2*x2^4 + 15*c1^2*c2^2*c3*c6*c7^4*x0^3*x1^3*x2^4 + 30*c0*c1*c2*c3^2*c6*c7^4*x0^3*x1^3*x2^4 + 5*c0^2*c3^3*c6*c7^4*x0^3*x1^3*x2^4 + 15*c1^2*c2*c3^2*c6*c7^4*x0^2*x1^4*x2^4 + 10*c0*c1*c3^3*c6*c7^4*x0^2*x1^4*x2^4 + 5*c1^2*c3^3*c6*c7^4*x0*x1^5*x2^4 + c0^2*c2^3*c7^5*x0^5*x2^5 + 2*c0*c1*c2^3*c7^5*x0^4*x1*x2^5 + 3*c0^2*c2^2*c3*c7^5*x0^4*x1*x2^5 + c1^2*c2^3*c7^5*x0^3*x1^2*x2^5 + 6*c0*c1*c2^2*c3*c7^5*x0^3*x1^2*x2^5 + 3*c0^2*c2*c3^2*c7^5*x0^3*x1^2*x2^5 + 3*c1^2*c2^2*c3*c7^5*x0^2*x1^3*x2^5 + 6*c0*c1*c2*c3^2*c7^5*x0^2*x1^3*x2^5 + c0^2*c3^3*c7^5*x0^2*x1^3*x2^5 + 3*c1^2*c2*c3^2*c7^5*x0*x1^4*x2^5 + 2*c0*c1*c3^3*c7^5*x0*x1^4*x2^5 + c1^2*c3^3*c7^5*x1^5*x2^5, c0^4*c2*c4^3*c6^2*x0^10 + c0^5*c4^2*c6^3*x0^10 + c0^2*c2^3*c6^5*x0^10 + 4*c0^3*c1*c2*c4^3*c6^2*x0^9*x1 + c0^4*c3*c4^3*c6^2*x0^9*x1 + 5*c0^4*c1*c4^2*c6^3*x0^9*x1 + 2*c0*c1*c2^3*c6^5*x0^9*x1 + 3*c0^2*c2^2*c3*c6^5*x0^9*x1 + 6*c0^2*c1^2*c2*c4^3*c6^2*x0^8*x1^2 + 4*c0^3*c1*c3*c4^3*c6^2*x0^8*x1^2 + 10*c0^3*c1^2*c4^2*c6^3*x0^8*x1^2 + c1^2*c2^3*c6^5*x0^8*x1^2 + 6*c0*c1*c2^2*c3*c6^5*x0^8*x1^2 + 3*c0^2*c2*c3^2*c6^5*x0^8*x1^2 + 4*c0*c1^3*c2*c4^3*c6^2*x0^7*x1^3 + 6*c0^2*c1^2*c3*c4^3*c6^2*x0^7*x1^3 + 10*c0^2*c1^3*c4^2*c6^3*x0^7*x1^3 + 3*c1^2*c2^2*c3*c6^5*x0^7*x1^3 + 6*c0*c1*c2*c3^2*c6^5*x0^7*x1^3 + c0^2*c3^3*c6^5*x0^7*x1^3 + c1^4*c2*c4^3*c6^2*x0^6*x1^4 + 4*c0*c1^3*c3*c4^3*c6^2*x0^6*x1^4 + 5*c0*c1^4*c4^2*c6^3*x0^6*x1^4 + 3*c1^2*c2*c3^2*c6^5*x0^6*x1^4 + 2*c0*c1*c3^3*c6^5*x0^6*x1^4 + c1^4*c3*c4^3*c6^2*x0^5*x1^5 + c1^5*c4^2*c6^3*x0^5*x1^5 + c1^2*c3^3*c6^5*x0^5*x1^5 + 3*c0^4*c2*c4^2*c5*c6^2*x0^9*x2 + 2*c0^5*c4*c5*c6^3*x0^9*x2 + 2*c0^4*c2*c4^3*c6*c7*x0^9*x2 + 3*c0^5*c4^2*c6^2*c7*x0^9*x2 + 5*c0^2*c2^3*c6^4*c7*x0^9*x2 + 12*c0^3*c1*c2*c4^2*c5*c6^2*x0^8*x1*x2 + 3*c0^4*c3*c4^2*c5*c6^2*x0^8*x1*x2 + 10*c0^4*c1*c4*c5*c6^3*x0^8*x1*x2 + 8*c0^3*c1*c2*c4^3*c6*c7*x0^8*x1*x2 + 2*c0^4*c3*c4^3*c6*c7*x0^8*x1*x2 + 15*c0^4*c1*c4^2*c6^2*c7*x0^8*x1*x2 + 10*c0*c1*c2^3*c6^4*c7*x0^8*x1*x2 + 15*c0^2*c2^2*c3*c6^4*c7*x0^8*x1*x2 + 18*c0^2*c1^2*c2*c4^2*c5*c6^2*x0^7*x1^2*x2 + 12*c0^3*c1*c3*c4^2*c5*c6^2*x0^7*x1^2*x2 + 20*c0^3*c1^2*c4*c5*c6^3*x0^7*x1^2*x2 + 12*c0^2*c1^2*c2*c4^3*c6*c7*x0^7*x1^2*x2 + 8*c0^3*c1*c3*c4^3*c6*c7*x0^7*x1^2*x2 + 30*c0^3*c1^2*c4^2*c6^2*c7*x0^7*x1^2*x2 + 5*c1^2*c2^3*c6^4*c7*x0^7*x1^2*x2 + 30*c0*c1*c2^2*c3*c6^4*c7*x0^7*x1^2*x2 + 15*c0^2*c2*c3^2*c6^4*c7*x0^7*x1^2*x2 + 12*c0*c1^3*c2*c4^2*c5*c6^2*x0^6*x1^3*x2 + 18*c0^2*c1^2*c3*c4^2*c5*c6^2*x0^6*x1^3*x2 + 20*c0^2*c1^3*c4*c5*c6^3*x0^6*x1^3*x2 + 8*c0*c1^3*c2*c4^3*c6*c7*x0^6*x1^3*x2 + 12*c0^2*c1^2*c3*c4^3*c6*c7*x0^6*x1^3*x2 + 30*c0^2*c1^3*c4^2*c6^2*c7*x0^6*x1^3*x2 + 15*c1^2*c2^2*c3*c6^4*c7*x0^6*x1^3*x2 + 30*c0*c1*c2*c3^2*c6^4*c7*x0^6*x1^3*x2 + 5*c0^2*c3^3*c6^4*c7*x0^6*x1^3*x2 + 3*c1^4*c2*c4^2*c5*c6^2*x0^5*x1^4*x2 + 12*c0*c1^3*c3*c4^2*c5*c6^2*x0^5*x1^4*x2 + 10*c0*c1^4*c4*c5*c6^3*x0^5*x1^4*x2 + 2*c1^4*c2*c4^3*c6*c7*x0^5*x1^4*x2 + 8*c0*c1^3*c3*c4^3*c6*c7*x0^5*x1^4*x2 + 15*c0*c1^4*c4^2*c6^2*c7*x0^5*x1^4*x2 + 15*c1^2*c2*c3^2*c6^4*c7*x0^5*x1^4*x2 + 10*c0*c1*c3^3*c6^4*c7*x0^5*x1^4*x2 + 3*c1^4*c3*c4^2*c5*c6^2*x0^4*x1^5*x2 + 2*c1^5*c4*c5*c6^3*x0^4*x1^5*x2 + 2*c1^4*c3*c4^3*c6*c7*x0^4*x1^5*x2 + 3*c1^5*c4^2*c6^2*c7*x0^4*x1^5*x2 + 5*c1^2*c3^3*c6^4*c7*x0^4*x1^5*x2 + 3*c0^4*c2*c4*c5^2*c6^2*x0^8*x2^2 + c0^5*c5^2*c6^3*x0^8*x2^2 + 6*c0^4*c2*c4^2*c5*c6*c7*x0^8*x2^2 + 6*c0^5*c4*c5*c6^2*c7*x0^8*x2^2 + c0^4*c2*c4^3*c7^2*x0^8*x2^2 + 3*c0^5*c4^2*c6*c7^2*x0^8*x2^2 + 10*c0^2*c2^3*c6^3*c7^2*x0^8*x2^2 + 12*c0^3*c1*c2*c4*c5^2*c6^2*x0^7*x1*x2^2 + 3*c0^4*c3*c4*c5^2*c6^2*x0^7*x1*x2^2 + 5*c0^4*c1*c5^2*c6^3*x0^7*x1*x2^2 + 24*c0^3*c1*c2*c4^2*c5*c6*c7*x0^7*x1*x2^2 + 6*c0^4*c3*c4^2*c5*c6*c7*x0^7*x1*x2^2 + 30*c0^4*c1*c4*c5*c6^2*c7*x0^7*x1*x2^2 + 4*c0^3*c1*c2*c4^3*c7^2*x0^7*x1*x2^2 + c0^4*c3*c4^3*c7^2*x0^7*x1*x2^2 + 15*c0^4*c1*c4^2*c6*c7^2*x0^7*x1*x2^2 + 20*c0*c1*c2^3*c6^3*c7^2*x0^7*x1*x2^2 + 30*c0^2*c2^2*c3*c6^3*c7^2*x0^7*x1*x2^2 + 18*c0^2*c1^2*c2*c4*c5^2*c6^2*x0^6*x1^2*x2^2 + 12*c0^3*c1*c3*c4*c5^2*c6^2*x0^6*x1^2*x2^2 + 10*c0^3*c1^2*c5^2*c6^3*x0^6*x1^2*x2^2 + 36*c0^2*c1^2*c2*c4^2*c5*c6*c7*x0^6*x1^2*x2^2 + 24*c0^3*c1*c3*c4^2*c5*c6*c7*x0^6*x1^2*x2^2 + 60*c0^3*c1^2*c4*c5*c6^2*c7*x0^6*x1^2*x2^2 + 6*c0^2*c1^2*c2*c4^3*c7^2*x0^6*x1^2*x2^2 + 4*c0^3*c1*c3*c4^3*c7^2*x0^6*x1^2*x2^2 + 30*c0^3*c1^2*c4^2*c6*c7^2*x0^6*x1^2*x2^2 + 10*c1^2*c2^3*c6^3*c7^2*x0^6*x1^2*x2^2 + 60*c0*c1*c2^2*c3*c6^3*c7^2*x0^6*x1^2*x2^2 + 30*c0^2*c2*c3^2*c6^3*c7^2*x0^6*x1^2*x2^2 + 12*c0*c1^3*c2*c4*c5^2*c6^2*x0^5*x1^3*x2^2 + 18*c0^2*c1^2*c3*c4*c5^2*c6^2*x0^5*x1^3*x2^2 + 10*c0^2*c1^3*c5^2*c6^3*x0^5*x1^3*x2^2 + 24*c0*c1^3*c2*c4^2*c5*c6*c7*x0^5*x1^3*x2^2 + 36*c0^2*c1^2*c3*c4^2*c5*c6*c7*x0^5*x1^3*x2^2 + 60*c0^2*c1^3*c4*c5*c6^2*c7*x0^5*x1^3*x2^2 + 4*c0*c1^3*c2*c4^3*c7^2*x0^5*x1^3*x2^2 + 6*c0^2*c1^2*c3*c4^3*c7^2*x0^5*x1^3*x2^2 + 30*c0^2*c1^3*c4^2*c6*c7^2*x0^5*x1^3*x2^2 + 30*c1^2*c2^2*c3*c6^3*c7^2*x0^5*x1^3*x2^2 + 60*c0*c1*c2*c3^2*c6^3*c7^2*x0^5*x1^3*x2^2 + 10*c0^2*c3^3*c6^3*c7^2*x0^5*x1^3*x2^2 + 3*c1^4*c2*c4*c5^2*c6^2*x0^4*x1^4*x2^2 + 12*c0*c1^3*c3*c4*c5^2*c6^2*x0^4*x1^4*x2^2 + 5*c0*c1^4*c5^2*c6^3*x0^4*x1^4*x2^2 + 6*c1^4*c2*c4^2*c5*c6*c7*x0^4*x1^4*x2^2 + 24*c0*c1^3*c3*c4^2*c5*c6*c7*x0^4*x1^4*x2^2 + 30*c0*c1^4*c4*c5*c6^2*c7*x0^4*x1^4*x2^2 + c1^4*c2*c4^3*c7^2*x0^4*x1^4*x2^2 + 4*c0*c1^3*c3*c4^3*c7^2*x0^4*x1^4*x2^2 + 15*c0*c1^4*c4^2*c6*c7^2*x0^4*x1^4*x2^2 + 30*c1^2*c2*c3^2*c6^3*c7^2*x0^4*x1^4*x2^2 + 20*c0*c1*c3^3*c6^3*c7^2*x0^4*x1^4*x2^2 + 3*c1^4*c3*c4*c5^2*c6^2*x0^3*x1^5*x2^2 + c1^5*c5^2*c6^3*x0^3*x1^5*x2^2 + 6*c1^4*c3*c4^2*c5*c6*c7*x0^3*x1^5*x2^2 + 6*c1^5*c4*c5*c6^2*c7*x0^3*x1^5*x2^2 + c1^4*c3*c4^3*c7^2*x0^3*x1^5*x2^2 + 3*c1^5*c4^2*c6*c7^2*x0^3*x1^5*x2^2 + 10*c1^2*c3^3*c6^3*c7^2*x0^3*x1^5*x2^2 + c0^4*c2*c5^3*c6^2*x0^7*x2^3 + 6*c0^4*c2*c4*c5^2*c6*c7*x0^7*x2^3 + 3*c0^5*c5^2*c6^2*c7*x0^7*x2^3 + 3*c0^4*c2*c4^2*c5*c7^2*x0^7*x2^3 + 6*c0^5*c4*c5*c6*c7^2*x0^7*x2^3 + c0^5*c4^2*c7^3*x0^7*x2^3 + 10*c0^2*c2^3*c6^2*c7^3*x0^7*x2^3 + 4*c0^3*c1*c2*c5^3*c6^2*x0^6*x1*x2^3 + c0^4*c3*c5^3*c6^2*x0^6*x1*x2^3 + 24*c0^3*c1*c2*c4*c5^2*c6*c7*x0^6*x1*x2^3 + 6*c0^4*c3*c4*c5^2*c6*c7*x0^6*x1*x2^3 + 15*c0^4*c1*c5^2*c6^2*c7*x0^6*x1*x2^3 + 12*c0^3*c1*c2*c4^2*c5*c7^2*x0^6*x1*x2^3 + 3*c0^4*c3*c4^2*c5*c7^2*x0^6*x1*x2^3 + 30*c0^4*c1*c4*c5*c6*c7^2*x0^6*x1*x2^3 + 5*c0^4*c1*c4^2*c7^3*x0^6*x1*x2^3 + 20*c0*c1*c2^3*c6^2*c7^3*x0^6*x1*x2^3 + 30*c0^2*c2^2*c3*c6^2*c7^3*x0^6*x1*x2^3 + 6*c0^2*c1^2*c2*c5^3*c6^2*x0^5*x1^2*x2^3 + 4*c0^3*c1*c3*c5^3*c6^2*x0^5*x1^2*x2^3 + 36*c0^2*c1^2*c2*c4*c5^2*c6*c7*x0^5*x1^2*x2^3 + 24*c0^3*c1*c3*c4*c5^2*c6*c7*x0^5*x1^2*x2^3 + 30*c0^3*c1^2*c5^2*c6^2*c7*x0^5*x1^2*x2^3 + 18*c0^2*c1^2*c2*c4^2*c5*c7^2*x0^5*x1^2*x2^3 + 12*c0^3*c1*c3*c4^2*c5*c7^2*x0^5*x1^2*x2^3 + 60*c0^3*c1^2*c4*c5*c6*c7^2*x0^5*x1^2*x2^3 + 10*c0^3*c1^2*c4^2*c7^3*x0^5*x1^2*x2^3 + 10*c1^2*c2^3*c6^2*c7^3*x0^5*x1^2*x2^3 + 60*c0*c1*c2^2*c3*c6^2*c7^3*x0^5*x1^2*x2^3 + 30*c0^2*c2*c3^2*c6^2*c7^3*x0^5*x1^2*x2^3 + 4*c0*c1^3*c2*c5^3*c6^2*x0^4*x1^3*x2^3 + 6*c0^2*c1^2*c3*c5^3*c6^2*x0^4*x1^3*x2^3 + 24*c0*c1^3*c2*c4*c5^2*c6*c7*x0^4*x1^3*x2^3 + 36*c0^2*c1^2*c3*c4*c5^2*c6*c7*x0^4*x1^3*x2^3 + 30*c0^2*c1^3*c5^2*c6^2*c7*x0^4*x1^3*x2^3 + 12*c0*c1^3*c2*c4^2*c5*c7^2*x0^4*x1^3*x2^3 + 18*c0^2*c1^2*c3*c4^2*c5*c7^2*x0^4*x1^3*x2^3 + 60*c0^2*c1^3*c4*c5*c6*c7^2*x0^4*x1^3*x2^3 + 10*c0^2*c1^3*c4^2*c7^3*x0^4*x1^3*x2^3 + 30*c1^2*c2^2*c3*c6^2*c7^3*x0^4*x1^3*x2^3 + 60*c0*c1*c2*c3^2*c6^2*c7^3*x0^4*x1^3*x2^3 + 10*c0^2*c3^3*c6^2*c7^3*x0^4*x1^3*x2^3 + c1^4*c2*c5^3*c6^2*x0^3*x1^4*x2^3 + 4*c0*c1^3*c3*c5^3*c6^2*x0^3*x1^4*x2^3 + 6*c1^4*c2*c4*c5^2*c6*c7*x0^3*x1^4*x2^3 + 24*c0*c1^3*c3*c4*c5^2*c6*c7*x0^3*x1^4*x2^3 + 15*c0*c1^4*c5^2*c6^2*c7*x0^3*x1^4*x2^3 + 3*c1^4*c2*c4^2*c5*c7^2*x0^3*x1^4*x2^3 + 12*c0*c1^3*c3*c4^2*c5*c7^2*x0^3*x1^4*x2^3 + 30*c0*c1^4*c4*c5*c6*c7^2*x0^3*x1^4*x2^3 + 5*c0*c1^4*c4^2*c7^3*x0^3*x1^4*x2^3 + 30*c1^2*c2*c3^2*c6^2*c7^3*x0^3*x1^4*x2^3 + 20*c0*c1*c3^3*c6^2*c7^3*x0^3*x1^4*x2^3 + c1^4*c3*c5^3*c6^2*x0^2*x1^5*x2^3 + 6*c1^4*c3*c4*c5^2*c6*c7*x0^2*x1^5*x2^3 + 3*c1^5*c5^2*c6^2*c7*x0^2*x1^5*x2^3 + 3*c1^4*c3*c4^2*c5*c7^2*x0^2*x1^5*x2^3 + 6*c1^5*c4*c5*c6*c7^2*x0^2*x1^5*x2^3 + c1^5*c4^2*c7^3*x0^2*x1^5*x2^3 + 10*c1^2*c3^3*c6^2*c7^3*x0^2*x1^5*x2^3 + 2*c0^4*c2*c5^3*c6*c7*x0^6*x2^4 + 3*c0^4*c2*c4*c5^2*c7^2*x0^6*x2^4 + 3*c0^5*c5^2*c6*c7^2*x0^6*x2^4 + 2*c0^5*c4*c5*c7^3*x0^6*x2^4 + 5*c0^2*c2^3*c6*c7^4*x0^6*x2^4 + 8*c0^3*c1*c2*c5^3*c6*c7*x0^5*x1*x2^4 + 2*c0^4*c3*c5^3*c6*c7*x0^5*x1*x2^4 + 12*c0^3*c1*c2*c4*c5^2*c7^2*x0^5*x1*x2^4 + 3*c0^4*c3*c4*c5^2*c7^2*x0^5*x1*x2^4 + 15*c0^4*c1*c5^2*c6*c7^2*x0^5*x1*x2^4 + 10*c0^4*c1*c4*c5*c7^3*x0^5*x1*x2^4 + 10*c0*c1*c2^3*c6*c7^4*x0^5*x1*x2^4 + 15*c0^2*c2^2*c3*c6*c7^4*x0^5*x1*x2^4 + 12*c0^2*c1^2*c2*c5^3*c6*c7*x0^4*x1^2*x2^4 + 8*c0^3*c1*c3*c5^3*c6*c7*x0^4*x1^2*x2^4 + 18*c0^2*c1^2*c2*c4*c5^2*c7^2*x0^4*x1^2*x2^4 + 12*c0^3*c1*c3*c4*c5^2*c7^2*x0^4*x1^2*x2^4 + 30*c0^3*c1^2*c5^2*c6*c7^2*x0^4*x1^2*x2^4 + 20*c0^3*c1^2*c4*c5*c7^3*x0^4*x1^2*x2^4 + 5*c1^2*c2^3*c6*c7^4*x0^4*x1^2*x2^4 + 30*c0*c1*c2^2*c3*c6*c7^4*x0^4*x1^2*x2^4 + 15*c0^2*c2*c3^2*c6*c7^4*x0^4*x1^2*x2^4 + 8*c0*c1^3*c2*c5^3*c6*c7*x0^3*x1^3*x2^4 + 12*c0^2*c1^2*c3*c5^3*c6*c7*x0^3*x1^3*x2^4 + 12*c0*c1^3*c2*c4*c5^2*c7^2*x0^3*x1^3*x2^4 + 18*c0^2*c1^2*c3*c4*c5^2*c7^2*x0^3*x1^3*x2^4 + 30*c0^2*c1^3*c5^2*c6*c7^2*x0^3*x1^3*x2^4 + 20*c0^2*c1^3*c4*c5*c7^3*x0^3*x1^3*x2^4 + 15*c1^2*c2^2*c3*c6*c7^4*x0^3*x1^3*x2^4 + 30*c0*c1*c2*c3^2*c6*c7^4*x0^3*x1^3*x2^4 + 5*c0^2*c3^3*c6*c7^4*x0^3*x1^3*x2^4 + 2*c1^4*c2*c5^3*c6*c7*x0^2*x1^4*x2^4 + 8*c0*c1^3*c3*c5^3*c6*c7*x0^2*x1^4*x2^4 + 3*c1^4*c2*c4*c5^2*c7^2*x0^2*x1^4*x2^4 + 12*c0*c1^3*c3*c4*c5^2*c7^2*x0^2*x1^4*x2^4 + 15*c0*c1^4*c5^2*c6*c7^2*x0^2*x1^4*x2^4 + 10*c0*c1^4*c4*c5*c7^3*x0^2*x1^4*x2^4 + 15*c1^2*c2*c3^2*c6*c7^4*x0^2*x1^4*x2^4 + 10*c0*c1*c3^3*c6*c7^4*x0^2*x1^4*x2^4 + 2*c1^4*c3*c5^3*c6*c7*x0*x1^5*x2^4 + 3*c1^4*c3*c4*c5^2*c7^2*x0*x1^5*x2^4 + 3*c1^5*c5^2*c6*c7^2*x0*x1^5*x2^4 + 2*c1^5*c4*c5*c7^3*x0*x1^5*x2^4 + 5*c1^2*c3^3*c6*c7^4*x0*x1^5*x2^4 + c0^4*c2*c5^3*c7^2*x0^5*x2^5 + c0^5*c5^2*c7^3*x0^5*x2^5 + c0^2*c2^3*c7^5*x0^5*x2^5 + 4*c0^3*c1*c2*c5^3*c7^2*x0^4*x1*x2^5 + c0^4*c3*c5^3*c7^2*x0^4*x1*x2^5 + 5*c0^4*c1*c5^2*c7^3*x0^4*x1*x2^5 + 2*c0*c1*c2^3*c7^5*x0^4*x1*x2^5 + 3*c0^2*c2^2*c3*c7^5*x0^4*x1*x2^5 + 6*c0^2*c1^2*c2*c5^3*c7^2*x0^3*x1^2*x2^5 + 4*c0^3*c1*c3*c5^3*c7^2*x0^3*x1^2*x2^5 + 10*c0^3*c1^2*c5^2*c7^3*x0^3*x1^2*x2^5 + c1^2*c2^3*c7^5*x0^3*x1^2*x2^5 + 6*c0*c1*c2^2*c3*c7^5*x0^3*x1^2*x2^5 + 3*c0^2*c2*c3^2*c7^5*x0^3*x1^2*x2^5 + 4*c0*c1^3*c2*c5^3*c7^2*x0^2*x1^3*x2^5 + 6*c0^2*c1^2*c3*c5^3*c7^2*x0^2*x1^3*x2^5 + 10*c0^2*c1^3*c5^2*c7^3*x0^2*x1^3*x2^5 + 3*c1^2*c2^2*c3*c7^5*x0^2*x1^3*x2^5 + 6*c0*c1*c2*c3^2*c7^5*x0^2*x1^3*x2^5 + c0^2*c3^3*c7^5*x0^2*x1^3*x2^5 + c1^4*c2*c5^3*c7^2*x0*x1^4*x2^5 + 4*c0*c1^3*c3*c5^3*c7^2*x0*x1^4*x2^5 + 5*c0*c1^4*c5^2*c7^3*x0*x1^4*x2^5 + 3*c1^2*c2*c3^2*c7^5*x0*x1^4*x2^5 + 2*c0*c1*c3^3*c7^5*x0*x1^4*x2^5 + c1^4*c3*c5^3*c7^2*x1^5*x2^5 + c1^5*c5^2*c7^3*x1^5*x2^5 + c1^2*c3^3*c7^5*x1^5*x2^5] 
eqn0_lst = 26 [0, c1^2*c3^3*c7^5, 10*c0^2*c2^3*c6^3*c7^2, c1^4*c3*c5^3*c7^2 + c1^5*c5^2*c7^3 + c1^2*c3^3*c7^5, 10*c0*c1*c2^3*c6^4*c7 + 15*c0^2*c2^2*c3*c6^4*c7, 5*c0^2*c2^3*c6^4*c7, 2*c0*c1*c2^3*c6^5 + 3*c0^2*c2^2*c3*c6^5, c0^2*c2^3*c6^5, 3*c0^4*c2*c4*c5^2*c6^2 + c0^5*c5^2*c6^3 + 6*c0^4*c2*c4^2*c5*c6*c7 + 6*c0^5*c4*c5*c6^2*c7 + c0^4*c2*c4^3*c7^2 + 3*c0^5*c4^2*c6*c7^2 + 10*c0^2*c2^3*c6^3*c7^2, 12*c0^3*c1*c2*c4^2*c5*c6^2 + 3*c0^4*c3*c4^2*c5*c6^2 + 10*c0^4*c1*c4*c5*c6^3 + 8*c0^3*c1*c2*c4^3*c6*c7 + 2*c0^4*c3*c4^3*c6*c7 + 15*c0^4*c1*c4^2*c6^2*c7 + 10*c0*c1*c2^3*c6^4*c7 + 15*c0^2*c2^2*c3*c6^4*c7, 3*c0^4*c2*c4^2*c5*c6^2 + 2*c0^5*c4*c5*c6^3 + 2*c0^4*c2*c4^3*c6*c7 + 3*c0^5*c4^2*c6^2*c7 + 5*c0^2*c2^3*c6^4*c7, 4*c0^3*c1*c2*c4^3*c6^2 + c0^4*c3*c4^3*c6^2 + 5*c0^4*c1*c4^2*c6^3 + 2*c0*c1*c2^3*c6^5 + 3*c0^2*c2^2*c3*c6^5, c0^4*c2*c4^3*c6^2 + c0^5*c4^2*c6^3 + c0^2*c2^3*c6^5, c1^3*c3^2*c5^5, c1^3*c3^2*c5^5 + c3^5*c5^3*c7^2, 10*c0^3*c2^2*c4^3*c5^2, 10*c0^3*c2^2*c4^3*c5^2 + 3*c2^5*c4*c5^2*c6^2 + 6*c2^5*c4^2*c5*c6*c7 + c2^5*c4^3*c7^2, 15*c0^2*c1*c2^2*c4^4*c5 + 10*c0^3*c2*c3*c4^4*c5, 15*c0^2*c1*c2^2*c4^4*c5 + 10*c0^3*c2*c3*c4^4*c5 + 15*c2^4*c3*c4^2*c5*c6^2 + 10*c2^4*c3*c4^3*c6*c7, 5*c0^3*c2^2*c4^4*c5, 5*c0^3*c2^2*c4^4*c5 + 3*c2^5*c4^2*c5*c6^2 + 2*c2^5*c4^3*c6*c7, 3*c0^2*c1*c2^2*c4^5 + 2*c0^3*c2*c3*c4^5, 3*c0^2*c1*c2^2*c4^5 + 2*c0^3*c2*c3*c4^5 + 5*c2^4*c3*c4^3*c6^2, c0^3*c2^2*c4^5, c0^3*c2^2*c4^5 + c2^5*c4^3*c6^2, t*c1*c2*c5*c6 - t*c0*c3*c5*c6 - t*c1*c2*c4*c7 + t*c0*c3*c4*c7 - 1] 
	 [c6, c2, c1^2, c5^3] 
	 [c4, c0, c3^2, c7^4, c3*c7^3, c3*c5*c7^2 + c1*c7^3] 
gr00 = 4 x0^2 [c3^2*x0^6*x1^2, c3^5*c7^2*x0*x1^5*x2^2 + c3^2*x0^6*x1^2, c3^3*c7^5*x1^3*x2^5, c3^3*c7^5*x1^3*x2^5 + c3*c7^2*x0^5*x1*x2^2 + c7^3*x0^5*x2^3] 
Mgr00 = (4, 45) [(0, 0, 0, c3^2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, c3^2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4*c3^7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32*c3^8, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 4*c3^3, 8*c3^3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32*c3^8, 0, 0, 0)] 
```

#### Compatible reparmatrization r1

We do the same for the reparametrization r1, but here we do not find a solution.

```python
# compose g with reparametrization r1
gcd1 = sage_gcd( [ comp.subs( r1 ) for comp in g ] )
assert gcd1 == 1
gr1 = [ comp.subs( r1 ) / gcd1 for comp in g ]
print( 'gr1 =', gcd1, gr1 )
assert SERing.get_degree( gr1 ) == 10
assert SERing.get_degree( f ) == 8

# find conditions on c so that gr1 has the same basepoints as f
eqn1_lst = usecase_B2_helper_bp( gr1 )
eqn1_lst += [ ring( '(c0*c3-c1*c2)*(c4*c7-c5*c6)*t-1' ) ]
print( 'eqn1_lst =', len( eqn1_lst ), eqn1_lst )
prime1_lst = sage_ideal( eqn1_lst ).elimination_ideal( ring( 't' ) ).primary_decomposition()
for prime1 in prime1_lst: print( '\t', prime1.gens() )
sol10 = {c[0]:0, c[3]:0, c[5]:0, c[6]:0, c[1]:1, c[4]:1}  # notice that wlog c1=c4=1
sol11 = {c[1]:0, c[2]:0, c[4]:0, c[7]:0, c[0]:1, c[5]:1}  # notice that wlog c0=c5=1
assert len( prime1_lst ) == 2
assert set( [gen.subs( sol10 ) for gen in prime1_lst[0].gens()] ) == set( [0] )
assert set( [gen.subs( sol11 ) for gen in prime1_lst[1].gens()] ) == set( [0] )

# sol10: notice that c2!=0 and c7!=0
gcd10 = sage_gcd( [ comp.subs( sol10 ) for comp in gr1] )
assert gcd10 == x[0] * x[0]
gr10 = [ comp.subs( sol10 ) / gcd10 for comp in gr1]
print( 'gr10 =', len( gr10 ), gcd10, gr10 )
assert SERing.get_degree( gr10 ) == 8
assert [] == sage_ideal( ( SERing.get_matrix_P2( gr10 ) * Kf ).list() + [ring( 'c2*c7*t-1' )] ).elimination_ideal( ring( 't' ) ).primary_decomposition()
# --> no solution

# sol11: notice that c3!=0 and c6!=0
gcd11 = sage_gcd( [ comp.subs( sol11 ) for comp in gr1] )
assert gcd11 == x[0] * x[0]
gr11 = [ comp.subs( sol11 ) / gcd11 for comp in gr1]
print( 'gr11 =', len( gr11 ), gcd11, gr11 )
assert SERing.get_degree( gr11 ) == 8
assert [] == sage_ideal( ( SERing.get_matrix_P2( gr11 ) * Kf ).list() + [ring( 'c3*c6*t-1' )] ).elimination_ideal( ring( 't' ) ).primary_decomposition()
# --> no solution

```
Output:
```
gr1 = 1 [c0^3*c2^2*c4^5*x0^10 + 5*c0^3*c2^2*c4^4*c5*x0^9*x1 + 10*c0^3*c2^2*c4^3*c5^2*x0^8*x1^2 + 10*c0^3*c2^2*c4^2*c5^3*x0^7*x1^3 + 5*c0^3*c2^2*c4*c5^4*x0^6*x1^4 + c0^3*c2^2*c5^5*x0^5*x1^5 + 3*c0^2*c1*c2^2*c4^5*x0^9*x2 + 2*c0^3*c2*c3*c4^5*x0^9*x2 + 15*c0^2*c1*c2^2*c4^4*c5*x0^8*x1*x2 + 10*c0^3*c2*c3*c4^4*c5*x0^8*x1*x2 + 30*c0^2*c1*c2^2*c4^3*c5^2*x0^7*x1^2*x2 + 20*c0^3*c2*c3*c4^3*c5^2*x0^7*x1^2*x2 + 30*c0^2*c1*c2^2*c4^2*c5^3*x0^6*x1^3*x2 + 20*c0^3*c2*c3*c4^2*c5^3*x0^6*x1^3*x2 + 15*c0^2*c1*c2^2*c4*c5^4*x0^5*x1^4*x2 + 10*c0^3*c2*c3*c4*c5^4*x0^5*x1^4*x2 + 3*c0^2*c1*c2^2*c5^5*x0^4*x1^5*x2 + 2*c0^3*c2*c3*c5^5*x0^4*x1^5*x2 + 3*c0*c1^2*c2^2*c4^5*x0^8*x2^2 + 6*c0^2*c1*c2*c3*c4^5*x0^8*x2^2 + c0^3*c3^2*c4^5*x0^8*x2^2 + 15*c0*c1^2*c2^2*c4^4*c5*x0^7*x1*x2^2 + 30*c0^2*c1*c2*c3*c4^4*c5*x0^7*x1*x2^2 + 5*c0^3*c3^2*c4^4*c5*x0^7*x1*x2^2 + 30*c0*c1^2*c2^2*c4^3*c5^2*x0^6*x1^2*x2^2 + 60*c0^2*c1*c2*c3*c4^3*c5^2*x0^6*x1^2*x2^2 + 10*c0^3*c3^2*c4^3*c5^2*x0^6*x1^2*x2^2 + 30*c0*c1^2*c2^2*c4^2*c5^3*x0^5*x1^3*x2^2 + 60*c0^2*c1*c2*c3*c4^2*c5^3*x0^5*x1^3*x2^2 + 10*c0^3*c3^2*c4^2*c5^3*x0^5*x1^3*x2^2 + 15*c0*c1^2*c2^2*c4*c5^4*x0^4*x1^4*x2^2 + 30*c0^2*c1*c2*c3*c4*c5^4*x0^4*x1^4*x2^2 + 5*c0^3*c3^2*c4*c5^4*x0^4*x1^4*x2^2 + 3*c0*c1^2*c2^2*c5^5*x0^3*x1^5*x2^2 + 6*c0^2*c1*c2*c3*c5^5*x0^3*x1^5*x2^2 + c0^3*c3^2*c5^5*x0^3*x1^5*x2^2 + c1^3*c2^2*c4^5*x0^7*x2^3 + 6*c0*c1^2*c2*c3*c4^5*x0^7*x2^3 + 3*c0^2*c1*c3^2*c4^5*x0^7*x2^3 + 5*c1^3*c2^2*c4^4*c5*x0^6*x1*x2^3 + 30*c0*c1^2*c2*c3*c4^4*c5*x0^6*x1*x2^3 + 15*c0^2*c1*c3^2*c4^4*c5*x0^6*x1*x2^3 + 10*c1^3*c2^2*c4^3*c5^2*x0^5*x1^2*x2^3 + 60*c0*c1^2*c2*c3*c4^3*c5^2*x0^5*x1^2*x2^3 + 30*c0^2*c1*c3^2*c4^3*c5^2*x0^5*x1^2*x2^3 + 10*c1^3*c2^2*c4^2*c5^3*x0^4*x1^3*x2^3 + 60*c0*c1^2*c2*c3*c4^2*c5^3*x0^4*x1^3*x2^3 + 30*c0^2*c1*c3^2*c4^2*c5^3*x0^4*x1^3*x2^3 + 5*c1^3*c2^2*c4*c5^4*x0^3*x1^4*x2^3 + 30*c0*c1^2*c2*c3*c4*c5^4*x0^3*x1^4*x2^3 + 15*c0^2*c1*c3^2*c4*c5^4*x0^3*x1^4*x2^3 + c1^3*c2^2*c5^5*x0^2*x1^5*x2^3 + 6*c0*c1^2*c2*c3*c5^5*x0^2*x1^5*x2^3 + 3*c0^2*c1*c3^2*c5^5*x0^2*x1^5*x2^3 + 2*c1^3*c2*c3*c4^5*x0^6*x2^4 + 3*c0*c1^2*c3^2*c4^5*x0^6*x2^4 + 10*c1^3*c2*c3*c4^4*c5*x0^5*x1*x2^4 + 15*c0*c1^2*c3^2*c4^4*c5*x0^5*x1*x2^4 + 20*c1^3*c2*c3*c4^3*c5^2*x0^4*x1^2*x2^4 + 30*c0*c1^2*c3^2*c4^3*c5^2*x0^4*x1^2*x2^4 + 20*c1^3*c2*c3*c4^2*c5^3*x0^3*x1^3*x2^4 + 30*c0*c1^2*c3^2*c4^2*c5^3*x0^3*x1^3*x2^4 + 10*c1^3*c2*c3*c4*c5^4*x0^2*x1^4*x2^4 + 15*c0*c1^2*c3^2*c4*c5^4*x0^2*x1^4*x2^4 + 2*c1^3*c2*c3*c5^5*x0*x1^5*x2^4 + 3*c0*c1^2*c3^2*c5^5*x0*x1^5*x2^4 + c1^3*c3^2*c4^5*x0^5*x2^5 + 5*c1^3*c3^2*c4^4*c5*x0^4*x1*x2^5 + 10*c1^3*c3^2*c4^3*c5^2*x0^3*x1^2*x2^5 + 10*c1^3*c3^2*c4^2*c5^3*x0^2*x1^3*x2^5 + 5*c1^3*c3^2*c4*c5^4*x0*x1^4*x2^5 + c1^3*c3^2*c5^5*x1^5*x2^5, c0^3*c2^2*c4^5*x0^10 + c2^5*c4^3*c6^2*x0^10 + 5*c0^3*c2^2*c4^4*c5*x0^9*x1 + 3*c2^5*c4^2*c5*c6^2*x0^9*x1 + 2*c2^5*c4^3*c6*c7*x0^9*x1 + 10*c0^3*c2^2*c4^3*c5^2*x0^8*x1^2 + 3*c2^5*c4*c5^2*c6^2*x0^8*x1^2 + 6*c2^5*c4^2*c5*c6*c7*x0^8*x1^2 + c2^5*c4^3*c7^2*x0^8*x1^2 + 10*c0^3*c2^2*c4^2*c5^3*x0^7*x1^3 + c2^5*c5^3*c6^2*x0^7*x1^3 + 6*c2^5*c4*c5^2*c6*c7*x0^7*x1^3 + 3*c2^5*c4^2*c5*c7^2*x0^7*x1^3 + 5*c0^3*c2^2*c4*c5^4*x0^6*x1^4 + 2*c2^5*c5^3*c6*c7*x0^6*x1^4 + 3*c2^5*c4*c5^2*c7^2*x0^6*x1^4 + c0^3*c2^2*c5^5*x0^5*x1^5 + c2^5*c5^3*c7^2*x0^5*x1^5 + 3*c0^2*c1*c2^2*c4^5*x0^9*x2 + 2*c0^3*c2*c3*c4^5*x0^9*x2 + 5*c2^4*c3*c4^3*c6^2*x0^9*x2 + 15*c0^2*c1*c2^2*c4^4*c5*x0^8*x1*x2 + 10*c0^3*c2*c3*c4^4*c5*x0^8*x1*x2 + 15*c2^4*c3*c4^2*c5*c6^2*x0^8*x1*x2 + 10*c2^4*c3*c4^3*c6*c7*x0^8*x1*x2 + 30*c0^2*c1*c2^2*c4^3*c5^2*x0^7*x1^2*x2 + 20*c0^3*c2*c3*c4^3*c5^2*x0^7*x1^2*x2 + 15*c2^4*c3*c4*c5^2*c6^2*x0^7*x1^2*x2 + 30*c2^4*c3*c4^2*c5*c6*c7*x0^7*x1^2*x2 + 5*c2^4*c3*c4^3*c7^2*x0^7*x1^2*x2 + 30*c0^2*c1*c2^2*c4^2*c5^3*x0^6*x1^3*x2 + 20*c0^3*c2*c3*c4^2*c5^3*x0^6*x1^3*x2 + 5*c2^4*c3*c5^3*c6^2*x0^6*x1^3*x2 + 30*c2^4*c3*c4*c5^2*c6*c7*x0^6*x1^3*x2 + 15*c2^4*c3*c4^2*c5*c7^2*x0^6*x1^3*x2 + 15*c0^2*c1*c2^2*c4*c5^4*x0^5*x1^4*x2 + 10*c0^3*c2*c3*c4*c5^4*x0^5*x1^4*x2 + 10*c2^4*c3*c5^3*c6*c7*x0^5*x1^4*x2 + 15*c2^4*c3*c4*c5^2*c7^2*x0^5*x1^4*x2 + 3*c0^2*c1*c2^2*c5^5*x0^4*x1^5*x2 + 2*c0^3*c2*c3*c5^5*x0^4*x1^5*x2 + 5*c2^4*c3*c5^3*c7^2*x0^4*x1^5*x2 + 3*c0*c1^2*c2^2*c4^5*x0^8*x2^2 + 6*c0^2*c1*c2*c3*c4^5*x0^8*x2^2 + c0^3*c3^2*c4^5*x0^8*x2^2 + 10*c2^3*c3^2*c4^3*c6^2*x0^8*x2^2 + 15*c0*c1^2*c2^2*c4^4*c5*x0^7*x1*x2^2 + 30*c0^2*c1*c2*c3*c4^4*c5*x0^7*x1*x2^2 + 5*c0^3*c3^2*c4^4*c5*x0^7*x1*x2^2 + 30*c2^3*c3^2*c4^2*c5*c6^2*x0^7*x1*x2^2 + 20*c2^3*c3^2*c4^3*c6*c7*x0^7*x1*x2^2 + 30*c0*c1^2*c2^2*c4^3*c5^2*x0^6*x1^2*x2^2 + 60*c0^2*c1*c2*c3*c4^3*c5^2*x0^6*x1^2*x2^2 + 10*c0^3*c3^2*c4^3*c5^2*x0^6*x1^2*x2^2 + 30*c2^3*c3^2*c4*c5^2*c6^2*x0^6*x1^2*x2^2 + 60*c2^3*c3^2*c4^2*c5*c6*c7*x0^6*x1^2*x2^2 + 10*c2^3*c3^2*c4^3*c7^2*x0^6*x1^2*x2^2 + 30*c0*c1^2*c2^2*c4^2*c5^3*x0^5*x1^3*x2^2 + 60*c0^2*c1*c2*c3*c4^2*c5^3*x0^5*x1^3*x2^2 + 10*c0^3*c3^2*c4^2*c5^3*x0^5*x1^3*x2^2 + 10*c2^3*c3^2*c5^3*c6^2*x0^5*x1^3*x2^2 + 60*c2^3*c3^2*c4*c5^2*c6*c7*x0^5*x1^3*x2^2 + 30*c2^3*c3^2*c4^2*c5*c7^2*x0^5*x1^3*x2^2 + 15*c0*c1^2*c2^2*c4*c5^4*x0^4*x1^4*x2^2 + 30*c0^2*c1*c2*c3*c4*c5^4*x0^4*x1^4*x2^2 + 5*c0^3*c3^2*c4*c5^4*x0^4*x1^4*x2^2 + 20*c2^3*c3^2*c5^3*c6*c7*x0^4*x1^4*x2^2 + 30*c2^3*c3^2*c4*c5^2*c7^2*x0^4*x1^4*x2^2 + 3*c0*c1^2*c2^2*c5^5*x0^3*x1^5*x2^2 + 6*c0^2*c1*c2*c3*c5^5*x0^3*x1^5*x2^2 + c0^3*c3^2*c5^5*x0^3*x1^5*x2^2 + 10*c2^3*c3^2*c5^3*c7^2*x0^3*x1^5*x2^2 + c1^3*c2^2*c4^5*x0^7*x2^3 + 6*c0*c1^2*c2*c3*c4^5*x0^7*x2^3 + 3*c0^2*c1*c3^2*c4^5*x0^7*x2^3 + 10*c2^2*c3^3*c4^3*c6^2*x0^7*x2^3 + 5*c1^3*c2^2*c4^4*c5*x0^6*x1*x2^3 + 30*c0*c1^2*c2*c3*c4^4*c5*x0^6*x1*x2^3 + 15*c0^2*c1*c3^2*c4^4*c5*x0^6*x1*x2^3 + 30*c2^2*c3^3*c4^2*c5*c6^2*x0^6*x1*x2^3 + 20*c2^2*c3^3*c4^3*c6*c7*x0^6*x1*x2^3 + 10*c1^3*c2^2*c4^3*c5^2*x0^5*x1^2*x2^3 + 60*c0*c1^2*c2*c3*c4^3*c5^2*x0^5*x1^2*x2^3 + 30*c0^2*c1*c3^2*c4^3*c5^2*x0^5*x1^2*x2^3 + 30*c2^2*c3^3*c4*c5^2*c6^2*x0^5*x1^2*x2^3 + 60*c2^2*c3^3*c4^2*c5*c6*c7*x0^5*x1^2*x2^3 + 10*c2^2*c3^3*c4^3*c7^2*x0^5*x1^2*x2^3 + 10*c1^3*c2^2*c4^2*c5^3*x0^4*x1^3*x2^3 + 60*c0*c1^2*c2*c3*c4^2*c5^3*x0^4*x1^3*x2^3 + 30*c0^2*c1*c3^2*c4^2*c5^3*x0^4*x1^3*x2^3 + 10*c2^2*c3^3*c5^3*c6^2*x0^4*x1^3*x2^3 + 60*c2^2*c3^3*c4*c5^2*c6*c7*x0^4*x1^3*x2^3 + 30*c2^2*c3^3*c4^2*c5*c7^2*x0^4*x1^3*x2^3 + 5*c1^3*c2^2*c4*c5^4*x0^3*x1^4*x2^3 + 30*c0*c1^2*c2*c3*c4*c5^4*x0^3*x1^4*x2^3 + 15*c0^2*c1*c3^2*c4*c5^4*x0^3*x1^4*x2^3 + 20*c2^2*c3^3*c5^3*c6*c7*x0^3*x1^4*x2^3 + 30*c2^2*c3^3*c4*c5^2*c7^2*x0^3*x1^4*x2^3 + c1^3*c2^2*c5^5*x0^2*x1^5*x2^3 + 6*c0*c1^2*c2*c3*c5^5*x0^2*x1^5*x2^3 + 3*c0^2*c1*c3^2*c5^5*x0^2*x1^5*x2^3 + 10*c2^2*c3^3*c5^3*c7^2*x0^2*x1^5*x2^3 + 2*c1^3*c2*c3*c4^5*x0^6*x2^4 + 3*c0*c1^2*c3^2*c4^5*x0^6*x2^4 + 5*c2*c3^4*c4^3*c6^2*x0^6*x2^4 + 10*c1^3*c2*c3*c4^4*c5*x0^5*x1*x2^4 + 15*c0*c1^2*c3^2*c4^4*c5*x0^5*x1*x2^4 + 15*c2*c3^4*c4^2*c5*c6^2*x0^5*x1*x2^4 + 10*c2*c3^4*c4^3*c6*c7*x0^5*x1*x2^4 + 20*c1^3*c2*c3*c4^3*c5^2*x0^4*x1^2*x2^4 + 30*c0*c1^2*c3^2*c4^3*c5^2*x0^4*x1^2*x2^4 + 15*c2*c3^4*c4*c5^2*c6^2*x0^4*x1^2*x2^4 + 30*c2*c3^4*c4^2*c5*c6*c7*x0^4*x1^2*x2^4 + 5*c2*c3^4*c4^3*c7^2*x0^4*x1^2*x2^4 + 20*c1^3*c2*c3*c4^2*c5^3*x0^3*x1^3*x2^4 + 30*c0*c1^2*c3^2*c4^2*c5^3*x0^3*x1^3*x2^4 + 5*c2*c3^4*c5^3*c6^2*x0^3*x1^3*x2^4 + 30*c2*c3^4*c4*c5^2*c6*c7*x0^3*x1^3*x2^4 + 15*c2*c3^4*c4^2*c5*c7^2*x0^3*x1^3*x2^4 + 10*c1^3*c2*c3*c4*c5^4*x0^2*x1^4*x2^4 + 15*c0*c1^2*c3^2*c4*c5^4*x0^2*x1^4*x2^4 + 10*c2*c3^4*c5^3*c6*c7*x0^2*x1^4*x2^4 + 15*c2*c3^4*c4*c5^2*c7^2*x0^2*x1^4*x2^4 + 2*c1^3*c2*c3*c5^5*x0*x1^5*x2^4 + 3*c0*c1^2*c3^2*c5^5*x0*x1^5*x2^4 + 5*c2*c3^4*c5^3*c7^2*x0*x1^5*x2^4 + c1^3*c3^2*c4^5*x0^5*x2^5 + c3^5*c4^3*c6^2*x0^5*x2^5 + 5*c1^3*c3^2*c4^4*c5*x0^4*x1*x2^5 + 3*c3^5*c4^2*c5*c6^2*x0^4*x1*x2^5 + 2*c3^5*c4^3*c6*c7*x0^4*x1*x2^5 + 10*c1^3*c3^2*c4^3*c5^2*x0^3*x1^2*x2^5 + 3*c3^5*c4*c5^2*c6^2*x0^3*x1^2*x2^5 + 6*c3^5*c4^2*c5*c6*c7*x0^3*x1^2*x2^5 + c3^5*c4^3*c7^2*x0^3*x1^2*x2^5 + 10*c1^3*c3^2*c4^2*c5^3*x0^2*x1^3*x2^5 + c3^5*c5^3*c6^2*x0^2*x1^3*x2^5 + 6*c3^5*c4*c5^2*c6*c7*x0^2*x1^3*x2^5 + 3*c3^5*c4^2*c5*c7^2*x0^2*x1^3*x2^5 + 5*c1^3*c3^2*c4*c5^4*x0*x1^4*x2^5 + 2*c3^5*c5^3*c6*c7*x0*x1^4*x2^5 + 3*c3^5*c4*c5^2*c7^2*x0*x1^4*x2^5 + c1^3*c3^2*c5^5*x1^5*x2^5 + c3^5*c5^3*c7^2*x1^5*x2^5, c0^2*c2^3*c6^5*x0^10 + 5*c0^2*c2^3*c6^4*c7*x0^9*x1 + 10*c0^2*c2^3*c6^3*c7^2*x0^8*x1^2 + 10*c0^2*c2^3*c6^2*c7^3*x0^7*x1^3 + 5*c0^2*c2^3*c6*c7^4*x0^6*x1^4 + c0^2*c2^3*c7^5*x0^5*x1^5 + 2*c0*c1*c2^3*c6^5*x0^9*x2 + 3*c0^2*c2^2*c3*c6^5*x0^9*x2 + 10*c0*c1*c2^3*c6^4*c7*x0^8*x1*x2 + 15*c0^2*c2^2*c3*c6^4*c7*x0^8*x1*x2 + 20*c0*c1*c2^3*c6^3*c7^2*x0^7*x1^2*x2 + 30*c0^2*c2^2*c3*c6^3*c7^2*x0^7*x1^2*x2 + 20*c0*c1*c2^3*c6^2*c7^3*x0^6*x1^3*x2 + 30*c0^2*c2^2*c3*c6^2*c7^3*x0^6*x1^3*x2 + 10*c0*c1*c2^3*c6*c7^4*x0^5*x1^4*x2 + 15*c0^2*c2^2*c3*c6*c7^4*x0^5*x1^4*x2 + 2*c0*c1*c2^3*c7^5*x0^4*x1^5*x2 + 3*c0^2*c2^2*c3*c7^5*x0^4*x1^5*x2 + c1^2*c2^3*c6^5*x0^8*x2^2 + 6*c0*c1*c2^2*c3*c6^5*x0^8*x2^2 + 3*c0^2*c2*c3^2*c6^5*x0^8*x2^2 + 5*c1^2*c2^3*c6^4*c7*x0^7*x1*x2^2 + 30*c0*c1*c2^2*c3*c6^4*c7*x0^7*x1*x2^2 + 15*c0^2*c2*c3^2*c6^4*c7*x0^7*x1*x2^2 + 10*c1^2*c2^3*c6^3*c7^2*x0^6*x1^2*x2^2 + 60*c0*c1*c2^2*c3*c6^3*c7^2*x0^6*x1^2*x2^2 + 30*c0^2*c2*c3^2*c6^3*c7^2*x0^6*x1^2*x2^2 + 10*c1^2*c2^3*c6^2*c7^3*x0^5*x1^3*x2^2 + 60*c0*c1*c2^2*c3*c6^2*c7^3*x0^5*x1^3*x2^2 + 30*c0^2*c2*c3^2*c6^2*c7^3*x0^5*x1^3*x2^2 + 5*c1^2*c2^3*c6*c7^4*x0^4*x1^4*x2^2 + 30*c0*c1*c2^2*c3*c6*c7^4*x0^4*x1^4*x2^2 + 15*c0^2*c2*c3^2*c6*c7^4*x0^4*x1^4*x2^2 + c1^2*c2^3*c7^5*x0^3*x1^5*x2^2 + 6*c0*c1*c2^2*c3*c7^5*x0^3*x1^5*x2^2 + 3*c0^2*c2*c3^2*c7^5*x0^3*x1^5*x2^2 + 3*c1^2*c2^2*c3*c6^5*x0^7*x2^3 + 6*c0*c1*c2*c3^2*c6^5*x0^7*x2^3 + c0^2*c3^3*c6^5*x0^7*x2^3 + 15*c1^2*c2^2*c3*c6^4*c7*x0^6*x1*x2^3 + 30*c0*c1*c2*c3^2*c6^4*c7*x0^6*x1*x2^3 + 5*c0^2*c3^3*c6^4*c7*x0^6*x1*x2^3 + 30*c1^2*c2^2*c3*c6^3*c7^2*x0^5*x1^2*x2^3 + 60*c0*c1*c2*c3^2*c6^3*c7^2*x0^5*x1^2*x2^3 + 10*c0^2*c3^3*c6^3*c7^2*x0^5*x1^2*x2^3 + 30*c1^2*c2^2*c3*c6^2*c7^3*x0^4*x1^3*x2^3 + 60*c0*c1*c2*c3^2*c6^2*c7^3*x0^4*x1^3*x2^3 + 10*c0^2*c3^3*c6^2*c7^3*x0^4*x1^3*x2^3 + 15*c1^2*c2^2*c3*c6*c7^4*x0^3*x1^4*x2^3 + 30*c0*c1*c2*c3^2*c6*c7^4*x0^3*x1^4*x2^3 + 5*c0^2*c3^3*c6*c7^4*x0^3*x1^4*x2^3 + 3*c1^2*c2^2*c3*c7^5*x0^2*x1^5*x2^3 + 6*c0*c1*c2*c3^2*c7^5*x0^2*x1^5*x2^3 + c0^2*c3^3*c7^5*x0^2*x1^5*x2^3 + 3*c1^2*c2*c3^2*c6^5*x0^6*x2^4 + 2*c0*c1*c3^3*c6^5*x0^6*x2^4 + 15*c1^2*c2*c3^2*c6^4*c7*x0^5*x1*x2^4 + 10*c0*c1*c3^3*c6^4*c7*x0^5*x1*x2^4 + 30*c1^2*c2*c3^2*c6^3*c7^2*x0^4*x1^2*x2^4 + 20*c0*c1*c3^3*c6^3*c7^2*x0^4*x1^2*x2^4 + 30*c1^2*c2*c3^2*c6^2*c7^3*x0^3*x1^3*x2^4 + 20*c0*c1*c3^3*c6^2*c7^3*x0^3*x1^3*x2^4 + 15*c1^2*c2*c3^2*c6*c7^4*x0^2*x1^4*x2^4 + 10*c0*c1*c3^3*c6*c7^4*x0^2*x1^4*x2^4 + 3*c1^2*c2*c3^2*c7^5*x0*x1^5*x2^4 + 2*c0*c1*c3^3*c7^5*x0*x1^5*x2^4 + c1^2*c3^3*c6^5*x0^5*x2^5 + 5*c1^2*c3^3*c6^4*c7*x0^4*x1*x2^5 + 10*c1^2*c3^3*c6^3*c7^2*x0^3*x1^2*x2^5 + 10*c1^2*c3^3*c6^2*c7^3*x0^2*x1^3*x2^5 + 5*c1^2*c3^3*c6*c7^4*x0*x1^4*x2^5 + c1^2*c3^3*c7^5*x1^5*x2^5, c0^4*c2*c4^3*c6^2*x0^10 + c0^5*c4^2*c6^3*x0^10 + c0^2*c2^3*c6^5*x0^10 + 3*c0^4*c2*c4^2*c5*c6^2*x0^9*x1 + 2*c0^5*c4*c5*c6^3*x0^9*x1 + 2*c0^4*c2*c4^3*c6*c7*x0^9*x1 + 3*c0^5*c4^2*c6^2*c7*x0^9*x1 + 5*c0^2*c2^3*c6^4*c7*x0^9*x1 + 3*c0^4*c2*c4*c5^2*c6^2*x0^8*x1^2 + c0^5*c5^2*c6^3*x0^8*x1^2 + 6*c0^4*c2*c4^2*c5*c6*c7*x0^8*x1^2 + 6*c0^5*c4*c5*c6^2*c7*x0^8*x1^2 + c0^4*c2*c4^3*c7^2*x0^8*x1^2 + 3*c0^5*c4^2*c6*c7^2*x0^8*x1^2 + 10*c0^2*c2^3*c6^3*c7^2*x0^8*x1^2 + c0^4*c2*c5^3*c6^2*x0^7*x1^3 + 6*c0^4*c2*c4*c5^2*c6*c7*x0^7*x1^3 + 3*c0^5*c5^2*c6^2*c7*x0^7*x1^3 + 3*c0^4*c2*c4^2*c5*c7^2*x0^7*x1^3 + 6*c0^5*c4*c5*c6*c7^2*x0^7*x1^3 + c0^5*c4^2*c7^3*x0^7*x1^3 + 10*c0^2*c2^3*c6^2*c7^3*x0^7*x1^3 + 2*c0^4*c2*c5^3*c6*c7*x0^6*x1^4 + 3*c0^4*c2*c4*c5^2*c7^2*x0^6*x1^4 + 3*c0^5*c5^2*c6*c7^2*x0^6*x1^4 + 2*c0^5*c4*c5*c7^3*x0^6*x1^4 + 5*c0^2*c2^3*c6*c7^4*x0^6*x1^4 + c0^4*c2*c5^3*c7^2*x0^5*x1^5 + c0^5*c5^2*c7^3*x0^5*x1^5 + c0^2*c2^3*c7^5*x0^5*x1^5 + 4*c0^3*c1*c2*c4^3*c6^2*x0^9*x2 + c0^4*c3*c4^3*c6^2*x0^9*x2 + 5*c0^4*c1*c4^2*c6^3*x0^9*x2 + 2*c0*c1*c2^3*c6^5*x0^9*x2 + 3*c0^2*c2^2*c3*c6^5*x0^9*x2 + 12*c0^3*c1*c2*c4^2*c5*c6^2*x0^8*x1*x2 + 3*c0^4*c3*c4^2*c5*c6^2*x0^8*x1*x2 + 10*c0^4*c1*c4*c5*c6^3*x0^8*x1*x2 + 8*c0^3*c1*c2*c4^3*c6*c7*x0^8*x1*x2 + 2*c0^4*c3*c4^3*c6*c7*x0^8*x1*x2 + 15*c0^4*c1*c4^2*c6^2*c7*x0^8*x1*x2 + 10*c0*c1*c2^3*c6^4*c7*x0^8*x1*x2 + 15*c0^2*c2^2*c3*c6^4*c7*x0^8*x1*x2 + 12*c0^3*c1*c2*c4*c5^2*c6^2*x0^7*x1^2*x2 + 3*c0^4*c3*c4*c5^2*c6^2*x0^7*x1^2*x2 + 5*c0^4*c1*c5^2*c6^3*x0^7*x1^2*x2 + 24*c0^3*c1*c2*c4^2*c5*c6*c7*x0^7*x1^2*x2 + 6*c0^4*c3*c4^2*c5*c6*c7*x0^7*x1^2*x2 + 30*c0^4*c1*c4*c5*c6^2*c7*x0^7*x1^2*x2 + 4*c0^3*c1*c2*c4^3*c7^2*x0^7*x1^2*x2 + c0^4*c3*c4^3*c7^2*x0^7*x1^2*x2 + 15*c0^4*c1*c4^2*c6*c7^2*x0^7*x1^2*x2 + 20*c0*c1*c2^3*c6^3*c7^2*x0^7*x1^2*x2 + 30*c0^2*c2^2*c3*c6^3*c7^2*x0^7*x1^2*x2 + 4*c0^3*c1*c2*c5^3*c6^2*x0^6*x1^3*x2 + c0^4*c3*c5^3*c6^2*x0^6*x1^3*x2 + 24*c0^3*c1*c2*c4*c5^2*c6*c7*x0^6*x1^3*x2 + 6*c0^4*c3*c4*c5^2*c6*c7*x0^6*x1^3*x2 + 15*c0^4*c1*c5^2*c6^2*c7*x0^6*x1^3*x2 + 12*c0^3*c1*c2*c4^2*c5*c7^2*x0^6*x1^3*x2 + 3*c0^4*c3*c4^2*c5*c7^2*x0^6*x1^3*x2 + 30*c0^4*c1*c4*c5*c6*c7^2*x0^6*x1^3*x2 + 5*c0^4*c1*c4^2*c7^3*x0^6*x1^3*x2 + 20*c0*c1*c2^3*c6^2*c7^3*x0^6*x1^3*x2 + 30*c0^2*c2^2*c3*c6^2*c7^3*x0^6*x1^3*x2 + 8*c0^3*c1*c2*c5^3*c6*c7*x0^5*x1^4*x2 + 2*c0^4*c3*c5^3*c6*c7*x0^5*x1^4*x2 + 12*c0^3*c1*c2*c4*c5^2*c7^2*x0^5*x1^4*x2 + 3*c0^4*c3*c4*c5^2*c7^2*x0^5*x1^4*x2 + 15*c0^4*c1*c5^2*c6*c7^2*x0^5*x1^4*x2 + 10*c0^4*c1*c4*c5*c7^3*x0^5*x1^4*x2 + 10*c0*c1*c2^3*c6*c7^4*x0^5*x1^4*x2 + 15*c0^2*c2^2*c3*c6*c7^4*x0^5*x1^4*x2 + 4*c0^3*c1*c2*c5^3*c7^2*x0^4*x1^5*x2 + c0^4*c3*c5^3*c7^2*x0^4*x1^5*x2 + 5*c0^4*c1*c5^2*c7^3*x0^4*x1^5*x2 + 2*c0*c1*c2^3*c7^5*x0^4*x1^5*x2 + 3*c0^2*c2^2*c3*c7^5*x0^4*x1^5*x2 + 6*c0^2*c1^2*c2*c4^3*c6^2*x0^8*x2^2 + 4*c0^3*c1*c3*c4^3*c6^2*x0^8*x2^2 + 10*c0^3*c1^2*c4^2*c6^3*x0^8*x2^2 + c1^2*c2^3*c6^5*x0^8*x2^2 + 6*c0*c1*c2^2*c3*c6^5*x0^8*x2^2 + 3*c0^2*c2*c3^2*c6^5*x0^8*x2^2 + 18*c0^2*c1^2*c2*c4^2*c5*c6^2*x0^7*x1*x2^2 + 12*c0^3*c1*c3*c4^2*c5*c6^2*x0^7*x1*x2^2 + 20*c0^3*c1^2*c4*c5*c6^3*x0^7*x1*x2^2 + 12*c0^2*c1^2*c2*c4^3*c6*c7*x0^7*x1*x2^2 + 8*c0^3*c1*c3*c4^3*c6*c7*x0^7*x1*x2^2 + 30*c0^3*c1^2*c4^2*c6^2*c7*x0^7*x1*x2^2 + 5*c1^2*c2^3*c6^4*c7*x0^7*x1*x2^2 + 30*c0*c1*c2^2*c3*c6^4*c7*x0^7*x1*x2^2 + 15*c0^2*c2*c3^2*c6^4*c7*x0^7*x1*x2^2 + 18*c0^2*c1^2*c2*c4*c5^2*c6^2*x0^6*x1^2*x2^2 + 12*c0^3*c1*c3*c4*c5^2*c6^2*x0^6*x1^2*x2^2 + 10*c0^3*c1^2*c5^2*c6^3*x0^6*x1^2*x2^2 + 36*c0^2*c1^2*c2*c4^2*c5*c6*c7*x0^6*x1^2*x2^2 + 24*c0^3*c1*c3*c4^2*c5*c6*c7*x0^6*x1^2*x2^2 + 60*c0^3*c1^2*c4*c5*c6^2*c7*x0^6*x1^2*x2^2 + 6*c0^2*c1^2*c2*c4^3*c7^2*x0^6*x1^2*x2^2 + 4*c0^3*c1*c3*c4^3*c7^2*x0^6*x1^2*x2^2 + 30*c0^3*c1^2*c4^2*c6*c7^2*x0^6*x1^2*x2^2 + 10*c1^2*c2^3*c6^3*c7^2*x0^6*x1^2*x2^2 + 60*c0*c1*c2^2*c3*c6^3*c7^2*x0^6*x1^2*x2^2 + 30*c0^2*c2*c3^2*c6^3*c7^2*x0^6*x1^2*x2^2 + 6*c0^2*c1^2*c2*c5^3*c6^2*x0^5*x1^3*x2^2 + 4*c0^3*c1*c3*c5^3*c6^2*x0^5*x1^3*x2^2 + 36*c0^2*c1^2*c2*c4*c5^2*c6*c7*x0^5*x1^3*x2^2 + 24*c0^3*c1*c3*c4*c5^2*c6*c7*x0^5*x1^3*x2^2 + 30*c0^3*c1^2*c5^2*c6^2*c7*x0^5*x1^3*x2^2 + 18*c0^2*c1^2*c2*c4^2*c5*c7^2*x0^5*x1^3*x2^2 + 12*c0^3*c1*c3*c4^2*c5*c7^2*x0^5*x1^3*x2^2 + 60*c0^3*c1^2*c4*c5*c6*c7^2*x0^5*x1^3*x2^2 + 10*c0^3*c1^2*c4^2*c7^3*x0^5*x1^3*x2^2 + 10*c1^2*c2^3*c6^2*c7^3*x0^5*x1^3*x2^2 + 60*c0*c1*c2^2*c3*c6^2*c7^3*x0^5*x1^3*x2^2 + 30*c0^2*c2*c3^2*c6^2*c7^3*x0^5*x1^3*x2^2 + 12*c0^2*c1^2*c2*c5^3*c6*c7*x0^4*x1^4*x2^2 + 8*c0^3*c1*c3*c5^3*c6*c7*x0^4*x1^4*x2^2 + 18*c0^2*c1^2*c2*c4*c5^2*c7^2*x0^4*x1^4*x2^2 + 12*c0^3*c1*c3*c4*c5^2*c7^2*x0^4*x1^4*x2^2 + 30*c0^3*c1^2*c5^2*c6*c7^2*x0^4*x1^4*x2^2 + 20*c0^3*c1^2*c4*c5*c7^3*x0^4*x1^4*x2^2 + 5*c1^2*c2^3*c6*c7^4*x0^4*x1^4*x2^2 + 30*c0*c1*c2^2*c3*c6*c7^4*x0^4*x1^4*x2^2 + 15*c0^2*c2*c3^2*c6*c7^4*x0^4*x1^4*x2^2 + 6*c0^2*c1^2*c2*c5^3*c7^2*x0^3*x1^5*x2^2 + 4*c0^3*c1*c3*c5^3*c7^2*x0^3*x1^5*x2^2 + 10*c0^3*c1^2*c5^2*c7^3*x0^3*x1^5*x2^2 + c1^2*c2^3*c7^5*x0^3*x1^5*x2^2 + 6*c0*c1*c2^2*c3*c7^5*x0^3*x1^5*x2^2 + 3*c0^2*c2*c3^2*c7^5*x0^3*x1^5*x2^2 + 4*c0*c1^3*c2*c4^3*c6^2*x0^7*x2^3 + 6*c0^2*c1^2*c3*c4^3*c6^2*x0^7*x2^3 + 10*c0^2*c1^3*c4^2*c6^3*x0^7*x2^3 + 3*c1^2*c2^2*c3*c6^5*x0^7*x2^3 + 6*c0*c1*c2*c3^2*c6^5*x0^7*x2^3 + c0^2*c3^3*c6^5*x0^7*x2^3 + 12*c0*c1^3*c2*c4^2*c5*c6^2*x0^6*x1*x2^3 + 18*c0^2*c1^2*c3*c4^2*c5*c6^2*x0^6*x1*x2^3 + 20*c0^2*c1^3*c4*c5*c6^3*x0^6*x1*x2^3 + 8*c0*c1^3*c2*c4^3*c6*c7*x0^6*x1*x2^3 + 12*c0^2*c1^2*c3*c4^3*c6*c7*x0^6*x1*x2^3 + 30*c0^2*c1^3*c4^2*c6^2*c7*x0^6*x1*x2^3 + 15*c1^2*c2^2*c3*c6^4*c7*x0^6*x1*x2^3 + 30*c0*c1*c2*c3^2*c6^4*c7*x0^6*x1*x2^3 + 5*c0^2*c3^3*c6^4*c7*x0^6*x1*x2^3 + 12*c0*c1^3*c2*c4*c5^2*c6^2*x0^5*x1^2*x2^3 + 18*c0^2*c1^2*c3*c4*c5^2*c6^2*x0^5*x1^2*x2^3 + 10*c0^2*c1^3*c5^2*c6^3*x0^5*x1^2*x2^3 + 24*c0*c1^3*c2*c4^2*c5*c6*c7*x0^5*x1^2*x2^3 + 36*c0^2*c1^2*c3*c4^2*c5*c6*c7*x0^5*x1^2*x2^3 + 60*c0^2*c1^3*c4*c5*c6^2*c7*x0^5*x1^2*x2^3 + 4*c0*c1^3*c2*c4^3*c7^2*x0^5*x1^2*x2^3 + 6*c0^2*c1^2*c3*c4^3*c7^2*x0^5*x1^2*x2^3 + 30*c0^2*c1^3*c4^2*c6*c7^2*x0^5*x1^2*x2^3 + 30*c1^2*c2^2*c3*c6^3*c7^2*x0^5*x1^2*x2^3 + 60*c0*c1*c2*c3^2*c6^3*c7^2*x0^5*x1^2*x2^3 + 10*c0^2*c3^3*c6^3*c7^2*x0^5*x1^2*x2^3 + 4*c0*c1^3*c2*c5^3*c6^2*x0^4*x1^3*x2^3 + 6*c0^2*c1^2*c3*c5^3*c6^2*x0^4*x1^3*x2^3 + 24*c0*c1^3*c2*c4*c5^2*c6*c7*x0^4*x1^3*x2^3 + 36*c0^2*c1^2*c3*c4*c5^2*c6*c7*x0^4*x1^3*x2^3 + 30*c0^2*c1^3*c5^2*c6^2*c7*x0^4*x1^3*x2^3 + 12*c0*c1^3*c2*c4^2*c5*c7^2*x0^4*x1^3*x2^3 + 18*c0^2*c1^2*c3*c4^2*c5*c7^2*x0^4*x1^3*x2^3 + 60*c0^2*c1^3*c4*c5*c6*c7^2*x0^4*x1^3*x2^3 + 10*c0^2*c1^3*c4^2*c7^3*x0^4*x1^3*x2^3 + 30*c1^2*c2^2*c3*c6^2*c7^3*x0^4*x1^3*x2^3 + 60*c0*c1*c2*c3^2*c6^2*c7^3*x0^4*x1^3*x2^3 + 10*c0^2*c3^3*c6^2*c7^3*x0^4*x1^3*x2^3 + 8*c0*c1^3*c2*c5^3*c6*c7*x0^3*x1^4*x2^3 + 12*c0^2*c1^2*c3*c5^3*c6*c7*x0^3*x1^4*x2^3 + 12*c0*c1^3*c2*c4*c5^2*c7^2*x0^3*x1^4*x2^3 + 18*c0^2*c1^2*c3*c4*c5^2*c7^2*x0^3*x1^4*x2^3 + 30*c0^2*c1^3*c5^2*c6*c7^2*x0^3*x1^4*x2^3 + 20*c0^2*c1^3*c4*c5*c7^3*x0^3*x1^4*x2^3 + 15*c1^2*c2^2*c3*c6*c7^4*x0^3*x1^4*x2^3 + 30*c0*c1*c2*c3^2*c6*c7^4*x0^3*x1^4*x2^3 + 5*c0^2*c3^3*c6*c7^4*x0^3*x1^4*x2^3 + 4*c0*c1^3*c2*c5^3*c7^2*x0^2*x1^5*x2^3 + 6*c0^2*c1^2*c3*c5^3*c7^2*x0^2*x1^5*x2^3 + 10*c0^2*c1^3*c5^2*c7^3*x0^2*x1^5*x2^3 + 3*c1^2*c2^2*c3*c7^5*x0^2*x1^5*x2^3 + 6*c0*c1*c2*c3^2*c7^5*x0^2*x1^5*x2^3 + c0^2*c3^3*c7^5*x0^2*x1^5*x2^3 + c1^4*c2*c4^3*c6^2*x0^6*x2^4 + 4*c0*c1^3*c3*c4^3*c6^2*x0^6*x2^4 + 5*c0*c1^4*c4^2*c6^3*x0^6*x2^4 + 3*c1^2*c2*c3^2*c6^5*x0^6*x2^4 + 2*c0*c1*c3^3*c6^5*x0^6*x2^4 + 3*c1^4*c2*c4^2*c5*c6^2*x0^5*x1*x2^4 + 12*c0*c1^3*c3*c4^2*c5*c6^2*x0^5*x1*x2^4 + 10*c0*c1^4*c4*c5*c6^3*x0^5*x1*x2^4 + 2*c1^4*c2*c4^3*c6*c7*x0^5*x1*x2^4 + 8*c0*c1^3*c3*c4^3*c6*c7*x0^5*x1*x2^4 + 15*c0*c1^4*c4^2*c6^2*c7*x0^5*x1*x2^4 + 15*c1^2*c2*c3^2*c6^4*c7*x0^5*x1*x2^4 + 10*c0*c1*c3^3*c6^4*c7*x0^5*x1*x2^4 + 3*c1^4*c2*c4*c5^2*c6^2*x0^4*x1^2*x2^4 + 12*c0*c1^3*c3*c4*c5^2*c6^2*x0^4*x1^2*x2^4 + 5*c0*c1^4*c5^2*c6^3*x0^4*x1^2*x2^4 + 6*c1^4*c2*c4^2*c5*c6*c7*x0^4*x1^2*x2^4 + 24*c0*c1^3*c3*c4^2*c5*c6*c7*x0^4*x1^2*x2^4 + 30*c0*c1^4*c4*c5*c6^2*c7*x0^4*x1^2*x2^4 + c1^4*c2*c4^3*c7^2*x0^4*x1^2*x2^4 + 4*c0*c1^3*c3*c4^3*c7^2*x0^4*x1^2*x2^4 + 15*c0*c1^4*c4^2*c6*c7^2*x0^4*x1^2*x2^4 + 30*c1^2*c2*c3^2*c6^3*c7^2*x0^4*x1^2*x2^4 + 20*c0*c1*c3^3*c6^3*c7^2*x0^4*x1^2*x2^4 + c1^4*c2*c5^3*c6^2*x0^3*x1^3*x2^4 + 4*c0*c1^3*c3*c5^3*c6^2*x0^3*x1^3*x2^4 + 6*c1^4*c2*c4*c5^2*c6*c7*x0^3*x1^3*x2^4 + 24*c0*c1^3*c3*c4*c5^2*c6*c7*x0^3*x1^3*x2^4 + 15*c0*c1^4*c5^2*c6^2*c7*x0^3*x1^3*x2^4 + 3*c1^4*c2*c4^2*c5*c7^2*x0^3*x1^3*x2^4 + 12*c0*c1^3*c3*c4^2*c5*c7^2*x0^3*x1^3*x2^4 + 30*c0*c1^4*c4*c5*c6*c7^2*x0^3*x1^3*x2^4 + 5*c0*c1^4*c4^2*c7^3*x0^3*x1^3*x2^4 + 30*c1^2*c2*c3^2*c6^2*c7^3*x0^3*x1^3*x2^4 + 20*c0*c1*c3^3*c6^2*c7^3*x0^3*x1^3*x2^4 + 2*c1^4*c2*c5^3*c6*c7*x0^2*x1^4*x2^4 + 8*c0*c1^3*c3*c5^3*c6*c7*x0^2*x1^4*x2^4 + 3*c1^4*c2*c4*c5^2*c7^2*x0^2*x1^4*x2^4 + 12*c0*c1^3*c3*c4*c5^2*c7^2*x0^2*x1^4*x2^4 + 15*c0*c1^4*c5^2*c6*c7^2*x0^2*x1^4*x2^4 + 10*c0*c1^4*c4*c5*c7^3*x0^2*x1^4*x2^4 + 15*c1^2*c2*c3^2*c6*c7^4*x0^2*x1^4*x2^4 + 10*c0*c1*c3^3*c6*c7^4*x0^2*x1^4*x2^4 + c1^4*c2*c5^3*c7^2*x0*x1^5*x2^4 + 4*c0*c1^3*c3*c5^3*c7^2*x0*x1^5*x2^4 + 5*c0*c1^4*c5^2*c7^3*x0*x1^5*x2^4 + 3*c1^2*c2*c3^2*c7^5*x0*x1^5*x2^4 + 2*c0*c1*c3^3*c7^5*x0*x1^5*x2^4 + c1^4*c3*c4^3*c6^2*x0^5*x2^5 + c1^5*c4^2*c6^3*x0^5*x2^5 + c1^2*c3^3*c6^5*x0^5*x2^5 + 3*c1^4*c3*c4^2*c5*c6^2*x0^4*x1*x2^5 + 2*c1^5*c4*c5*c6^3*x0^4*x1*x2^5 + 2*c1^4*c3*c4^3*c6*c7*x0^4*x1*x2^5 + 3*c1^5*c4^2*c6^2*c7*x0^4*x1*x2^5 + 5*c1^2*c3^3*c6^4*c7*x0^4*x1*x2^5 + 3*c1^4*c3*c4*c5^2*c6^2*x0^3*x1^2*x2^5 + c1^5*c5^2*c6^3*x0^3*x1^2*x2^5 + 6*c1^4*c3*c4^2*c5*c6*c7*x0^3*x1^2*x2^5 + 6*c1^5*c4*c5*c6^2*c7*x0^3*x1^2*x2^5 + c1^4*c3*c4^3*c7^2*x0^3*x1^2*x2^5 + 3*c1^5*c4^2*c6*c7^2*x0^3*x1^2*x2^5 + 10*c1^2*c3^3*c6^3*c7^2*x0^3*x1^2*x2^5 + c1^4*c3*c5^3*c6^2*x0^2*x1^3*x2^5 + 6*c1^4*c3*c4*c5^2*c6*c7*x0^2*x1^3*x2^5 + 3*c1^5*c5^2*c6^2*c7*x0^2*x1^3*x2^5 + 3*c1^4*c3*c4^2*c5*c7^2*x0^2*x1^3*x2^5 + 6*c1^5*c4*c5*c6*c7^2*x0^2*x1^3*x2^5 + c1^5*c4^2*c7^3*x0^2*x1^3*x2^5 + 10*c1^2*c3^3*c6^2*c7^3*x0^2*x1^3*x2^5 + 2*c1^4*c3*c5^3*c6*c7*x0*x1^4*x2^5 + 3*c1^4*c3*c4*c5^2*c7^2*x0*x1^4*x2^5 + 3*c1^5*c5^2*c6*c7^2*x0*x1^4*x2^5 + 2*c1^5*c4*c5*c7^3*x0*x1^4*x2^5 + 5*c1^2*c3^3*c6*c7^4*x0*x1^4*x2^5 + c1^4*c3*c5^3*c7^2*x1^5*x2^5 + c1^5*c5^2*c7^3*x1^5*x2^5 + c1^2*c3^3*c7^5*x1^5*x2^5] 
eqn1_lst = 26 [0, c1^2*c3^3*c7^5, c1^4*c3*c5^3*c7^2 + c1^5*c5^2*c7^3 + c1^2*c3^3*c7^5, 10*c0*c1*c2^3*c6^4*c7 + 15*c0^2*c2^2*c3*c6^4*c7, 5*c0^2*c2^3*c6^4*c7, c1^2*c2^3*c6^5 + 6*c0*c1*c2^2*c3*c6^5 + 3*c0^2*c2*c3^2*c6^5, 2*c0*c1*c2^3*c6^5 + 3*c0^2*c2^2*c3*c6^5, c0^2*c2^3*c6^5, 12*c0^3*c1*c2*c4^2*c5*c6^2 + 3*c0^4*c3*c4^2*c5*c6^2 + 10*c0^4*c1*c4*c5*c6^3 + 8*c0^3*c1*c2*c4^3*c6*c7 + 2*c0^4*c3*c4^3*c6*c7 + 15*c0^4*c1*c4^2*c6^2*c7 + 10*c0*c1*c2^3*c6^4*c7 + 15*c0^2*c2^2*c3*c6^4*c7, 3*c0^4*c2*c4^2*c5*c6^2 + 2*c0^5*c4*c5*c6^3 + 2*c0^4*c2*c4^3*c6*c7 + 3*c0^5*c4^2*c6^2*c7 + 5*c0^2*c2^3*c6^4*c7, 6*c0^2*c1^2*c2*c4^3*c6^2 + 4*c0^3*c1*c3*c4^3*c6^2 + 10*c0^3*c1^2*c4^2*c6^3 + c1^2*c2^3*c6^5 + 6*c0*c1*c2^2*c3*c6^5 + 3*c0^2*c2*c3^2*c6^5, 4*c0^3*c1*c2*c4^3*c6^2 + c0^4*c3*c4^3*c6^2 + 5*c0^4*c1*c4^2*c6^3 + 2*c0*c1*c2^3*c6^5 + 3*c0^2*c2^2*c3*c6^5, c0^4*c2*c4^3*c6^2 + c0^5*c4^2*c6^3 + c0^2*c2^3*c6^5, c1^3*c3^2*c5^5, c1^3*c3^2*c5^5 + c3^5*c5^3*c7^2, 15*c0^2*c1*c2^2*c4^4*c5 + 10*c0^3*c2*c3*c4^4*c5, 15*c0^2*c1*c2^2*c4^4*c5 + 10*c0^3*c2*c3*c4^4*c5 + 15*c2^4*c3*c4^2*c5*c6^2 + 10*c2^4*c3*c4^3*c6*c7, 5*c0^3*c2^2*c4^4*c5, 5*c0^3*c2^2*c4^4*c5 + 3*c2^5*c4^2*c5*c6^2 + 2*c2^5*c4^3*c6*c7, 3*c0*c1^2*c2^2*c4^5 + 6*c0^2*c1*c2*c3*c4^5 + c0^3*c3^2*c4^5, 3*c0*c1^2*c2^2*c4^5 + 6*c0^2*c1*c2*c3*c4^5 + c0^3*c3^2*c4^5 + 10*c2^3*c3^2*c4^3*c6^2, 3*c0^2*c1*c2^2*c4^5 + 2*c0^3*c2*c3*c4^5, 3*c0^2*c1*c2^2*c4^5 + 2*c0^3*c2*c3*c4^5 + 5*c2^4*c3*c4^3*c6^2, c0^3*c2^2*c4^5, c0^3*c2^2*c4^5 + c2^5*c4^3*c6^2, t*c1*c2*c5*c6 - t*c0*c3*c5*c6 - t*c1*c2*c4*c7 + t*c0*c3*c4*c7 - 1] 
	 [c6, c0, c5^2, c3^3] 
	 [c4, c2, c7^2, c1^3] 
gr10 = 4 x0^2 [c2^2*x0^5*x2^3, c2^5*c7^2*x0^6*x1^2 + c2^2*x0^5*x2^3, c2^3*c7^5*x0*x1^5*x2^2, c2^3*c7^5*x0*x1^5*x2^2 + c2*c7^2*x0^2*x1^2*x2^4 + c7^3*x1^3*x2^5] 
gr11 = 4 x0^2 [c3^2*x0*x1^5*x2^2, c3^5*c6^2*x1^3*x2^5 + c3^2*x0*x1^5*x2^2, c3^3*c6^5*x0^5*x2^3, c3^3*c6^5*x0^5*x2^3 + c6^3*x0^6*x1^2 + c3*c6^2*x0^4*x1^3*x2] 
```

#### Computing projective isomorphisms

We recover from Mgr00 the projective isomorphisms in terms of a parametrized matrix U.

```python
Mgr = Mgr00
Ef = sage_matrix( sage_QQ, list( Mf ) + list( Kf.T ) )
Egr = sage_matrix( list( Mgr ) + list( Kf.T ) )
UpI = Egr * ~Ef
assert ( UpI.submatrix( 4, 4 ) - sage_identity_matrix( 41 ) ).is_zero()
U = UpI.submatrix( 0, 0, 4, 4 )
U = U / sage_gcd( U.list() )
print( 'U =\n' + str( U ) )

# verify whether U*f is a parametrization for X for all (c0,...,c7)
Uf = list( U * sage_vector( f ) )
eqg_sub = [ eq.subs( {z[i]:Uf[i] for i in range( 4 )} ) for eq in eqg ]
assert eqg_sub == [0]

```
Output:

```
U = 
[      1       0       0       0]
[      1  4*c3^5       0       0]
[      0       0 32*c3^6       0]
[      0       0 32*c3^6    4*c3] 

```

### Example 5: Projective automorphisms of rational normal scroll (case B4)
    
We compute the projective automorphism of the 
rational normal scrolls that is parametrized 
the birational map f.
Further explanation of this example can be found 
in the accompanying [arxiv article (not yet online)](https://arxiv.org).

We start by importing the required libraries.

```python
from surface_equivalence.class_se_ring import ring
from surface_equivalence.class_se_ring import SERing
from surface_equivalence.sage_interface import sage_gcd
from surface_equivalence.sage_interface import sage_matrix
from surface_equivalence.sage_interface import sage_identity_matrix
from surface_equivalence.sage_interface import sage_vector
from surface_equivalence.sage_interface import sage_ideal
from surface_equivalence.sage_interface import sage_QQ
from surface_equivalence.sage_interface import sage_SR
from surface_equivalence.sage_interface import sage_solve
from linear_series.class_poly_ring import PolyRing
from linear_series.class_base_points import BasePointTree
from linear_series.class_linear_series import LinearSeries

os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/mathematica/link/bin' # edit path

```
We let f to be the parametrization of a rational normal scroll and 
we assume that `f==g`.

```python
f = ring('[-x0*x1^2 + x1^3, x1^2*x2, x1*x2^2, x0*x1*x2, -x0*x2^2 + x2^3]')
g = f

```
We start by doing a basepoint analysis for f.

```python
print( LinearSeries( SERing.conv( f ), PolyRing( 'x,y,z', True ) ).get_bp_tree() )
```
Output:
```
{ 5, <<x^3 - x^2*z, x^2*y, x*y^2, x*y*z, y^3 - y^2*z>>, QQ[x, y, z] }
chart=z, depth=0, mult=1, sol=(0, 1), { 5, <<x^3 - x^2, x^2*y, x*y^2, x*y, y^3 - y^2>>, QQ[x, y] }
chart=z, depth=0, mult=1, sol=(1, 0), { 5, <<x^3 - x^2, x^2*y, x*y^2, x*y, y^3 - y^2>>, QQ[x, y] }
chart=z, depth=0, mult=2, sol=(0, 0), { 5, <<x^3 - x^2, x^2*y, x*y^2, x*y, y^3 - y^2>>, QQ[x, y] }
```
We compute the generators of bigraded Cox rings associated to f.
    
```python
PolyRing.reset_base_field()
p1 = ( 0, 0 ); p2 = ( 1, 0 ); p3 = ( 0, 1 )

# e0-e1
bpt = BasePointTree(); bpt.add( 'z', p1, 1 )
f0p1 = SERing.conv( LinearSeries.get( [1], bpt ).pol_lst )

# e0-e2-e3
bpt = BasePointTree(); bpt.add( 'z', p2, 1 ); bpt.add( 'z', p3, 1 )
f1m2 = SERing.conv( LinearSeries.get( [1], bpt ).pol_lst )

# 2e0-e1-e2-e3
bpt = BasePointTree(); bpt.add( 'z', p1, 1 ); bpt.add( 'z', p2, 1 ); bpt.add( 'z', p3, 1 )
f1m1 = SERing.conv( LinearSeries.get( [2], bpt ).pol_lst )

# 3e0-2e1-e2-e3
bpt = BasePointTree(); bpt.add( 'z', p1, 2 ); bpt.add( 'z', p2, 1 ); bpt.add( 'z', p3, 1 )
f1m0 = SERing.conv( LinearSeries.get( [3], bpt ).pol_lst )

print( 'f0p1 =', len( f0p1 ), f0p1 )
print( 'f1m2 =', len( f1m2 ), f1m2 )
print( 'f1m1 =', len( f1m1 ), f1m1 )
print( 'f1m0 =', len( f1m0 ), f1m0 )

```
Output:
```
f0p1 = 2 [x1, x2] 
f1m2 = 1 [-x0 + x1 + x2] 
f1m1 = 3 [-x0*x1 + x1^2, x1*x2, -x0*x2 + x2^2] 
f1m0 = 5 [-x0*x1^2 + x1^3, x1^2*x2, x1*x2^2, x0*x1*x2, -x0*x2^2 + x2^3]
```

From the output we obtain the following generators for the Cox ring associated to f:
u0, u1, u2 and u3 with weights (0,1), (0,1), (1,-2) and (1,-1), respectively. 
We compute the list M1m0 of monomials of weight (1,0) in these generators.

```python
U = [ring( 'x1' ), ring( 'x2' ), ring( 'x1+x2-x0' ), ring( 'x1*x2' )]
u = ring( 'u0,u1,u2,u3' )
w_lst = [( 0, 1 ), ( 0, 1 ), ( 1, -2 ), ( 1, -1 )]
M1m0 = SERing.get_wmon_lst( u, w_lst, 1, 0 )
print( 'M1m0 =', M1m0 )

```
Output:
```
M1m0 = [u1*u3, u0*u3, u1^2*u2, u0*u1*u2, u0^2*u2]
```

Let F be the map whose components are defined by M1m0.
We compose F with a parametrized map

`P: (u0,u1,u2,u3)|-->(c0*u0+c1*u1, c2*u0+c3*u1, c4*u2, c5*u3+c6*u0*u2+c7*u1*u2)`.
 
We obtain the composition PoF after we replace u0, u1, u2 and u3 
with x1, x2, x1+x2-x0 and x1*x2, respectively.


```python
F = [ comp.subs( {u[i]:U[i] for i in range( 4 )} ) for comp in M1m0 ]
z = [ring( 'z' + str( i ) ) for i in range( 6 )]
c = [ring( 'c' + str( i ) ) for i in range( 8 )]
dctZ = { u[i]:z[i] for i in range( 4 )}
dctP = { z[0]:c[0] * u[0] + c[1] * u[1], z[1]:c[2] * u[0] + c[3] * u[1], z[2]:c[4] * u[2], z[3]:c[5] * u[3] + c[6] * u[0] * u[2] + c[7] * u[1] * u[2]}
PoF = [ comp.subs( dctZ ).subs( dctP ) for comp in M1m0 ]
PoF = [ comp.subs( {u[i]:U[i] for i in range( 4 )} ) for comp in PoF]
PoF = [ comp / sage_gcd( PoF ) for comp in PoF]
print( 'PoF =', PoF )

```
Output:
```
PoF = [-c2*c6*x0*x1^2 + c2*c6*x1^3 - c3*c6*x0*x1*x2 - c2*c7*x0*x1*x2 + c2*c5*x1^2*x2 + c2*c6*x1^2*x2 + c3*c6*x1^2*x2 + c2*c7*x1^2*x2 - c3*c7*x0*x2^2 + c3*c5*x1*x2^2 + c3*c6*x1*x2^2 + c2*c7*x1*x2^2 + c3*c7*x1*x2^2 + c3*c7*x2^3, -c0*c6*x0*x1^2 + c0*c6*x1^3 - c1*c6*x0*x1*x2 - c0*c7*x0*x1*x2 + c0*c5*x1^2*x2 + c0*c6*x1^2*x2 + c1*c6*x1^2*x2 + c0*c7*x1^2*x2 - c1*c7*x0*x2^2 + c1*c5*x1*x2^2 + c1*c6*x1*x2^2 + c0*c7*x1*x2^2 + c1*c7*x1*x2^2 + c1*c7*x2^3, -c2^2*c4*x0*x1^2 + c2^2*c4*x1^3 - 2*c2*c3*c4*x0*x1*x2 + c2^2*c4*x1^2*x2 + 2*c2*c3*c4*x1^2*x2 - c3^2*c4*x0*x2^2 + 2*c2*c3*c4*x1*x2^2 + c3^2*c4*x1*x2^2 + c3^2*c4*x2^3, -c0*c2*c4*x0*x1^2 + c0*c2*c4*x1^3 - c1*c2*c4*x0*x1*x2 - c0*c3*c4*x0*x1*x2 + c0*c2*c4*x1^2*x2 + c1*c2*c4*x1^2*x2 + c0*c3*c4*x1^2*x2 - c1*c3*c4*x0*x2^2 + c1*c2*c4*x1*x2^2 + c0*c3*c4*x1*x2^2 + c1*c3*c4*x1*x2^2 + c1*c3*c4*x2^3, -c0^2*c4*x0*x1^2 + c0^2*c4*x1^3 - 2*c0*c1*c4*x0*x1*x2 + c0^2*c4*x1^2*x2 + 2*c0*c1*c4*x1^2*x2 - c1^2*c4*x0*x2^2 + 2*c0*c1*c4*x1*x2^2 + c1^2*c4*x1*x2^2 + c1^2*c4*x2^3] 
```
We now recover the matrix M that defines a projective automorphism induced by P.

```python
# type "%paste" to paste indented code in a Sage or Python terminal
M = []
for pol in [ comp.subs( dctZ ).subs( dctP ) for comp in M1m0]:
    row = []
    for mon in M1m0:
        row += [ pol.coefficient( mon ) ]
    M += [row]
M = sage_matrix( M )
MF = list( M * sage_vector( F ) )
assert MF == PoF
print( 'M =\n' + str( M ) )

```
Output:
```
M =
[              c3*c5               c2*c5               c3*c7       c3*c6 + c2*c7               c2*c6]
[              c1*c5               c0*c5               c1*c7       c1*c6 + c0*c7               c0*c6]
[                  0                   0             c3^2*c4          2*c2*c3*c4             c2^2*c4]
[                  0                   0            c1*c3*c4 c1*c2*c4 + c0*c3*c4            c0*c2*c4]
[                  0                   0             c1^2*c4          2*c0*c1*c4             c0^2*c4]
```
We compute the inverse Q of the map F
and compute the composition QoPoF in order to obtain a compatible reparametrization from 
the projective plane to the projective plane.

```python
t = ring( 't' )
x = ring( 'x0,x1,x2' )
id = [ F[i] * z[0] - z[i] * F[0] for i in range( 5 ) ] + [t * F[0] - 1]
I1 = sage_ideal( id ).elimination_ideal( [t, x[2] ] ).gens()
I2 = sage_ideal( id ).elimination_ideal( [t, x[1] ] ).gens()
I1 = [ elt for elt in I1 if elt.degree( x[0] ) == 1 and elt.degree( x[1] ) == 1 ][0]
I2 = [ elt for elt in I2 if elt.degree( x[0] ) == 1 and elt.degree( x[2] ) == 1 ][0]
Q0 = I1.coefficient( x[1] )
Q1 = -I1.coefficient( x[0] )
Q2 = I2.coefficient( x[0] )
Q = [Q0, Q1, Q2]
QoF = [comp.subs( {z[i]:F[i] for i in range( 5 )} ) for comp in Q]
QoF = [comp / sage_gcd( QoF ) for comp in QoF]
assert QoF == [x[0], x[1], x[2]]
QoPoF = [ comp.subs( {z[i]:PoF[i] for i in range( 5 )} ) for comp in Q]
QoPoF = [ comp / sage_gcd( QoPoF ) for comp in QoPoF ]
print( 'Q =', Q )
print( 'QoPoF =', len(QoPoF), QoPoF )

```
Output:

```
Q = [z0 + z1 - z3, z1, z0]
QoPoF = 3 [c0*c2*c4*x0*x1^2 - c0*c2*c4*x1^3 + c1*c2*c4*x0*x1*x2 + c0*c3*c4*x0*x1*x2 - c0*c2*c4*x1^2*x2 - c1*c2*c4*x1^2*x2 - c0*c3*c4*x1^2*x2 + c1*c3*c4*x0*x2^2 - c1*c2*c4*x1*x2^2 - c0*c3*c4*x1*x2^2 - c1*c3*c4*x1*x2^2 - c1*c3*c4*x2^3 - c0*c6*x0*x1^2 - c2*c6*x0*x1^2 + c0*c6*x1^3 + c2*c6*x1^3 - c1*c6*x0*x1*x2 - c3*c6*x0*x1*x2 - c0*c7*x0*x1*x2 - c2*c7*x0*x1*x2 + c0*c5*x1^2*x2 + c2*c5*x1^2*x2 + c0*c6*x1^2*x2 + c1*c6*x1^2*x2 + c2*c6*x1^2*x2 + c3*c6*x1^2*x2 + c0*c7*x1^2*x2 + c2*c7*x1^2*x2 - c1*c7*x0*x2^2 - c3*c7*x0*x2^2 + c1*c5*x1*x2^2 + c3*c5*x1*x2^2 + c1*c6*x1*x2^2 + c3*c6*x1*x2^2 + c0*c7*x1*x2^2 + c1*c7*x1*x2^2 + c2*c7*x1*x2^2 + c3*c7*x1*x2^2 + c1*c7*x2^3 + c3*c7*x2^3, -c0*c6*x0*x1^2 + c0*c6*x1^3 - c1*c6*x0*x1*x2 - c0*c7*x0*x1*x2 + c0*c5*x1^2*x2 + c0*c6*x1^2*x2 + c1*c6*x1^2*x2 + c0*c7*x1^2*x2 - c1*c7*x0*x2^2 + c1*c5*x1*x2^2 + c1*c6*x1*x2^2 + c0*c7*x1*x2^2 + c1*c7*x1*x2^2 + c1*c7*x2^3, -c2*c6*x0*x1^2 + c2*c6*x1^3 - c3*c6*x0*x1*x2 - c2*c7*x0*x1*x2 + c2*c5*x1^2*x2 + c2*c6*x1^2*x2 + c3*c6*x1^2*x2 + c2*c7*x1^2*x2 - c3*c7*x0*x2^2 + c3*c5*x1*x2^2 + c3*c6*x1*x2^2 + c2*c7*x1*x2^2 + c3*c7*x1*x2^2 + c3*c7*x2^3] 

```

From the compatible reparametrizations QoPoF we compute the
projective automorphisms U in terms of a 5x5 matrix parametrized by c.


```python
assert f == f1m0
gr = [ comp.subs( {x[i]:QoPoF[i] for i in range( 3 )} ) for comp in f]
gcd_gr = sage_gcd( gr )
gr = [ comp / gcd_gr for comp in gr ]
Mf = SERing.get_matrix_P2( f )
Mgr = SERing.get_matrix_P2( gr )
Kf = Mf.right_kernel_matrix().T
assert ( Mf * Kf ).is_zero()
assert ( Mgr * Kf ).is_zero()
Ef = sage_matrix( sage_QQ, list( Mf ) + list( Kf.T ) )
Egr = sage_matrix( list( Mgr ) + list( Kf.T ) )
UpI = Egr * ~Ef
assert ( UpI.submatrix( 5, 5 ) - sage_identity_matrix( 5 ) ).is_zero()
U = UpI.submatrix( 0, 0, 5, 5 )

# verify whether U*f is a parametrization for X for all (c0,...,c7)
Uf = list( U * sage_vector( f ) )
eqX = sage_ideal( [ z[i] - f[i] for i in range( 5 )] ).elimination_ideal( [x[0], x[1], x[2]] ).gens()
eqXs = [ eq.subs( {z[i]:Uf[i] for i in range( 5 )} ) for eq in eqX ]
assert eqXs == [0, 0, 0]

print( 'U = '+str(U) )

```
Output:
```
U =
[                                                                                c0^2*c4 - c0*c6                                            c0^2*c4 + 2*c0*c1*c4 - c0*c5 - c0*c6 - c1*c6 - c0*c7                                            2*c0*c1*c4 + c1^2*c4 - c1*c5 - c1*c6 - c0*c7 - c1*c7                                                                     -2*c0*c1*c4 + c1*c6 + c0*c7                                                                                 c1^2*c4 - c1*c7]
[                                                                                          c0*c6                                                                   c0*c5 + c0*c6 + c1*c6 + c0*c7                                                                   c1*c5 + c1*c6 + c0*c7 + c1*c7                                                                                  -c1*c6 - c0*c7                                                                                           c1*c7]
[                                                                                          c2*c6                                                                   c2*c5 + c2*c6 + c3*c6 + c2*c7                                                                   c3*c5 + c3*c6 + c2*c7 + c3*c7                                                                                  -c3*c6 - c2*c7                                                                                           c3*c7]
[                                                                      -c0*c2*c4 + c0*c6 + c2*c6 -c0*c2*c4 - c1*c2*c4 - c0*c3*c4 + c0*c5 + c2*c5 + c0*c6 + c1*c6 + c2*c6 + c3*c6 + c0*c7 + c2*c7 -c1*c2*c4 - c0*c3*c4 - c1*c3*c4 + c1*c5 + c3*c5 + c1*c6 + c3*c6 + c0*c7 + c1*c7 + c2*c7 + c3*c7                                             c1*c2*c4 + c0*c3*c4 - c1*c6 - c3*c6 - c0*c7 - c2*c7                                                                       -c1*c3*c4 + c1*c7 + c3*c7]
[                                                                                c2^2*c4 - c2*c6                                            c2^2*c4 + 2*c2*c3*c4 - c2*c5 - c2*c6 - c3*c6 - c2*c7                                            2*c2*c3*c4 + c3^2*c4 - c3*c5 - c3*c6 - c2*c7 - c3*c7                                                                     -2*c2*c3*c4 + c3*c6 + c2*c7                                                                                 c3^2*c4 - c3*c7] 

```
    
### Example 6: Projective isomorphisms between conic bundles (case B5)


We compute the projective isomorphisms between surfaces that are conic bundles.
Further explanation of this example can be found 
in the accompanying [arxiv article (not yet online)](https://arxiv.org).

We start by importing the required libraries.

```python
from surface_equivalence.class_se_ring import ring
from surface_equivalence.class_se_ring import SERing
from surface_equivalence.sage_interface import sage_gcd
from surface_equivalence.sage_interface import sage_matrix
from surface_equivalence.sage_interface import sage_identity_matrix
from surface_equivalence.sage_interface import sage_vector
from surface_equivalence.sage_interface import sage_ideal
from surface_equivalence.sage_interface import sage_QQ
from surface_equivalence.sage_interface import sage_SR
from surface_equivalence.sage_interface import sage_solve
from surface_equivalence.sage_interface import sage_Compositions
from linear_series.class_poly_ring import PolyRing
from linear_series.class_base_points import BasePointTree
from linear_series.class_linear_series import LinearSeries

os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/mathematica/link/bin' # edit path

```
We compute the projective isomorphism between two 
conic bundles X and Y that are parametrized by the 
birational maps ff and gg, respectively.
The domain of ff is the projective plane in coordinates (x0:x1:x2)
and the domain of gg is the fiber product of the projective line with itself with 
coordinates (y0:y1;y2:y3).

```python
ff = ring("[-x0^2*x1^2 + x1^4, x1^3*x2, -x0^2*x1^2 + x0*x1^3, x1^2*x2^2, x0*x1^2*x2, x1*x2^3, x0*x1*x2^2, x0^2*x1*x2, -x0^2*x2^2 + x2^4, -x0^2*x2^2 + x0*x2^3]")
gg = ring("[y0^2*y1^2*y2*y3, y0*y1^3*y3^2, y0^2*y1^2*y2^2, y1^4*y2*y3 - y1^4*y3^2, y0*y1^3*y2*y3, y0*y1^3*y2^2, y1^4*y2^2 - y1^4*y3^2, y0^3*y1*y2*y3 - y0^2*y1^2*y3^2, y0^3*y1*y2^2, y0^4*y2^2 - y0^2*y1^2*y3^2]")

```

#### Graded coordinate ring associated to ff

```python
print( LinearSeries( SERing.conv( ff ), PolyRing( 'x,y,z' ) ).get_bp_tree() )

```
Output:
```
{ 10, <<x^4 - x^2*z^2, x^3*y, x^3*z - x^2*z^2, x^2*y^2, x^2*y*z, x*y^3, x*y^2*z, x*y*z^2, y^4 - y^2*z^2, y^3*z - y^2*z^2>>, QQ[x, y, z] }
chart=z, depth=0, mult=1, sol=(0, 1), { 10, <<x^4 - x^2, x^3*y, x^3 - x^2, x^2*y^2, x^2*y, x*y^3, x*y^2, x*y, y^4 - y^2, y^3 - y^2>>, QQ[x, y] }
chart=z, depth=0, mult=1, sol=(1, 0), { 10, <<x^4 - x^2, x^3*y, x^3 - x^2, x^2*y^2, x^2*y, x*y^3, x*y^2, x*y, y^4 - y^2, y^3 - y^2>>, QQ[x, y] }
chart=z, depth=0, mult=2, sol=(0, 0), { 10, <<x^4 - x^2, x^3*y, x^3 - x^2, x^2*y^2, x^2*y, x*y^3, x*y^2, x*y, y^4 - y^2, y^3 - y^2>>, QQ[x, y] }
```
We construct linear series associated to ff in order to determine
the generators of a graded coordinate ring associated to ff.

```python
PolyRing.reset_base_field()

# basepoints in chart x0!=0;
p1 = ( 0, 0 ); p2 = ( 0, 1 ); p3 = ( 1, 0 )

# 0f+p = e0-e1
bp_tree = BasePointTree(); bp_tree.add( 'z', p1, 1 )
f0p1 = SERing.conv( LinearSeries.get( [1], bp_tree ).pol_lst )

# 1f-3p = e0+e1-e2-e3
bp_tree = BasePointTree(); bp_tree.add( 'z', p2, 1 ); bp_tree.add( 'z', p3, 1 )
f1m3 = SERing.conv( LinearSeries.get( [1], bp_tree ).pol_lst )

# 1f-2p = 2e0-e2-e3
bp_tree = BasePointTree(); bp_tree.add( 'z', p2, 1 ); bp_tree.add( 'z', p3, 1 )
f1m2 = SERing.conv( LinearSeries.get( [2], bp_tree ).pol_lst )

# 1f-1p = 3e0-e1-e2-e3
bp_tree = BasePointTree(); bp_tree.add( 'z', p1, 1 ); bp_tree.add( 'z', p2, 1 ); bp_tree.add( 'z', p3, 1 )
f1m1 = SERing.conv( LinearSeries.get( [3], bp_tree ).pol_lst )

# 1f-0p = 4e0-2e1-e2-e3
bp_tree = BasePointTree(); bp_tree.add( 'z', p1, 2 ); bp_tree.add( 'z', p2, 1 ); bp_tree.add( 'z', p3, 1 )
f1m0 = SERing.conv( LinearSeries.get( [4], bp_tree ).pol_lst )

# 2f-4p = 4e0-2e2-2e3
bp_tree = BasePointTree(); bp_tree.add( 'z', p2, 2 ); bp_tree.add( 'z', p3, 2 )
f2m4 = SERing.conv( LinearSeries.get( [4], bp_tree ).pol_lst )

print( 'f0p1 =', len( f0p1 ), f0p1 )
print( 'f1m3 =', len( f1m3 ), f1m3 )
print( 'f1m2 =', len( f1m2 ), f1m2 )
print( 'f1m1 =', len( f1m1 ), f1m1 )
print( 'f1m0 =', len( f1m0 ), f1m0 )
print( 'f2m4 =', len( f2m4 ), f2m4 )

```
Output:
```
f0p1 = 2 [x1, x2] 
f1m3 = 1 [-x0 + x1 + x2] 
f1m2 = 4 [-x0^2 + x1^2 + x0*x2, x1*x2, -x0^2 + x0*x1 + x0*x2, -x0*x2 + x2^2] 
f1m1 = 7 [-x0^2*x1 + x1^3, x1^2*x2, -x0^2*x1 + x0*x1^2, x1*x2^2, x0*x1*x2, -x0^2*x2 + x2^3, -x0^2*x2 + x0*x2^2] 
f1m0 = 10 [-x0^2*x1^2 + x1^4, x1^3*x2, -x0^2*x1^2 + x0*x1^3, x1^2*x2^2, x0*x1^2*x2, x1*x2^3, x0*x1*x2^2, x0^2*x1*x2, -x0^2*x2^2 + x2^4, -x0^2*x2^2 + x0*x2^3] 
f2m4 = 9 [3*x0^4 - 4*x0^3*x1 + x1^4 - 4*x0^3*x2 + 4*x0^2*x1*x2 - x0^2*x2^2 + 2*x0*x2^3, -x0^3*x2 + x1^3*x2 + 2*x0^2*x2^2 - x0*x2^3, 2*x0^4 - 3*x0^3*x1 + x0*x1^3 - 3*x0^3*x2 + 3*x0^2*x1*x2 + x0*x2^3, x1^2*x2^2, -x0^3*x2 + x0*x1^2*x2 + 2*x0^2*x2^2 - x0*x2^3, x0^4 - 2*x0^3*x1 + x0^2*x1^2 - 2*x0^3*x2 + 2*x0^2*x1*x2 + x0^2*x2^2, x0^3*x2 - x0^2*x1*x2 - 2*x0^2*x2^2 + x0*x2^3 + x1*x2^3, x0^3*x2 - x0^2*x1*x2 - 2*x0^2*x2^2 + x0*x1*x2^2 + x0*x2^3, x0^2*x2^2 - 2*x0*x2^3 + x2^4] 
```

By inspection we recover the generators of graded ring associated to ff.

```python
U = U0, U1, U2, U3, U4 = ring( 'x1' ), ring( 'x2' ), ring( 'x1+x2-x0' ), ring( 'x1*x2' ), ring( '(x1+x2-x0)^2' )

```

We compute smallest d such that there is a relation between the generators of weight (2,d). 
The idea is to compare the number of monomials of given bi-weight 
with the dimension predicted by the Riemann-Roch formula.

```python
u = u0, u1, u2, u3, u4 = ring( 'u0,u1,u2,u3,u4' )
for d in reversed( [-i for i in range( 8 )] ):
    w_lst = [( 0, 1 ), ( 0, 1 ), ( 1, -3 ), ( 1, -2 ), ( 1, -2 )]
    print( '\tweight=', ( 2, d ), ',\t#monomials =', len( SERing.get_wmon_lst( u, w_lst, 2, d ) ), ',\tRiemann-Roch =', 29 + 5 * d )

```
Output:
```
weight= (2, -7) ,	#monomials = 0  ,	Riemann-Roch = -6 
weight= (2, -6) ,	#monomials = 1  ,	Riemann-Roch = -1 
weight= (2, -5) ,	#monomials = 4  ,	Riemann-Roch = 4 
weight= (2, -4) ,	#monomials = 10 ,	Riemann-Roch = 9 
weight= (2, -3) ,	#monomials = 16 ,	Riemann-Roch = 14
weight= (2, -2) ,	#monomials = 22 ,	Riemann-Roch = 19
weight= (2, -1) ,	#monomials = 28 ,	Riemann-Roch = 24
weight= (2,  0) ,	#monomials = 34 ,	Riemann-Roch = 29

```
We conclude from the above output that d is equal to -4
and we proceed by computing the polynomial relation between the monomials of weight (2,-4).

```python
T2m4 = ring( '[u3^2,u3*u4,u4^2,u0*u2*u3,u0*u2*u4,u1*u2*u3,u1*u2*u4,u0^2*u2^2,u0*u1*u2^2,u1^2*u2^2]' )
T1m0 = ring( '[u1^2*u4,u1^2*u3,u1^3*u2,u0*u1*u4,u0*u1*u3,u0*u1^2*u2,u0^2*u4,u0^2*u3,u0^2*u1*u2,u0^3*u2]' )
a = a0, a1, a2, a3, a4, a5, a6, a7, a8, a9 = [elt.subs( {u[i]:U[i] for i in range( 5 )} ) for elt in T2m4 ]
mata = sage_matrix( sage_QQ, SERing.get_matrix_P2( a ) )
kera = mata.transpose().right_kernel().matrix()
assert kera * sage_vector( a ) == sage_vector( [0] )
assert str(kera) == '[ 0  1  0  0  0  0  0  0 -1  0]'
assert a1 - a8 == 0

```
Thus the relation between the monomials a0,...,a9 of weight (2,-4) in 
the graded coordinate ring associated to ff is `a1 - a8==0`.

#### Graded coordinate ring associated to gg

We now proceed to do the same procedure for the map gg as we did for ff.
We start with a basepoint analysis.

```python
print( LinearSeries( SERing.conv( gg ), PolyRing( 'x,y,v,w' ) ).get_bp_tree() )

```
Output:
```
{ 10, <<x^2*y^2*v*w, x*y^3*w^2, x^2*y^2*v^2, y^4*v*w - y^4*w^2, x*y^3*v*w, x*y^3*v^2, y^4*v^2 - y^4*w^2, x^3*y*v*w - x^2*y^2*w^2, x^3*y*v^2, x^4*v^2 - x^2*y^2*w^2>>, QQ[x, y, v, w] }
chart=xw, depth=0, mult=2, sol=(0, 0), { 10, <<y^2*v, y^3, y^2*v^2, y^4*v - y^4, y^3*v, y^3*v^2, y^4*v^2 - y^4, -y^2 + y*v, y*v^2, -y^2 + v^2>>, QQ[y, v] }
    chart=t, depth=1, mult=1, sol=(1, 0), { 10, <<y^2*v, y^3*v, y^2*v^2, y^4*v^3 - y^4*v^2, y^3*v^2, y^3*v^3, y^4*v^4 - y^4*v^2, -y^2 + y, y*v, -y^2 + 1>>, QQ[y, v] }
chart=yv, depth=0, mult=1, sol=(0, 1), { 10, <<x^2*w, x*w^2, x^2, -w^2 + w, x*w, x, -w^2 + 1, x^3*w - x^2*w^2, x^3, x^4 - x^2*w^2>>, QQ[x, w] } 
```
We construct linear series in order to determine
the generators of a graded coordinate ring associated to gg.

```python
tree_211 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
tree_211.add( 'xw', ( 0, 0 ), 2 ).add( 't', ( 1, 0 ), 1 )
tree_211.add( 'yv', ( 0, 1 ), 1 )

tree_422 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
tree_422.add( 'xw', ( 0, 0 ), 4 ).add( 't', ( 1, 0 ), 2 )
tree_422.add( 'yv', ( 0, 1 ), 2 )

# 1g+0q = 4l0+2l1-2e1-e2-e3
g1m0 = SERing.conv( LinearSeries.get( [4, 2], tree_211 ).pol_lst )
# 1g-3q = (l0+l1-e1-e2-e3) + (b-e1)
g1m3 = SERing.conv( LinearSeries.get( [1, 2], tree_211 ).pol_lst )
# 1g-2q = 2l0+2l1-2e1-e2-e3
g1m2 = SERing.conv( LinearSeries.get( [2, 2], tree_211 ).pol_lst )
# 1g-1q = 3l0+2l1-2e1-e2-e3
g1m1 = SERing.conv( LinearSeries.get( [3, 2], tree_211 ).pol_lst )
# 2g-4q = 4l0+4l1-4e1-2e2-2e3
g2m4 = SERing.conv( LinearSeries.get( [4, 4], tree_422 ).pol_lst )

print( 'g1m0 =', len( g1m0 ), g1m0 )
print( 'g1m3 =', len( g1m3 ), g1m3 )
print( 'g1m2 =', len( g1m2 ), g1m2 )
print( 'g1m1 =', len( g1m1 ), g1m1 )
print( 'g2m4 =', len( g2m4 ), g2m4 )

```
Output:
```
g1m0 = 10 [y0^4*y2^2 - y0^2*y1^2*y3^2, y0^3*y1*y2^2, y0^3*y1*y2*y3 - y0^2*y1^2*y3^2, y0^2*y1^2*y2^2, y0^2*y1^2*y2*y3, y0*y1^3*y2^2, y0*y1^3*y2*y3, y0*y1^3*y3^2, y1^4*y2^2 - y1^4*y3^2, y1^4*y2*y3 - y1^4*y3^2] 
g1m3 = 1 [y0*y2^2 + y1*y2^2 - y1*y2*y3] 
g1m2 = 4 [y0^2*y2^2 + y1^2*y2*y3 - y1^2*y3^2, y0*y1*y2^2, y0*y1*y2*y3 + y1^2*y2*y3 - y1^2*y3^2, y1^2*y2^2 - y1^2*y2*y3] 
g1m1 = 7 [y0^3*y2^2 - y0*y1^2*y3^2, y0^2*y1*y2^2, y0^2*y1*y2*y3 - y0*y1^2*y3^2, y0*y1^2*y2^2, y0*y1^2*y2*y3, y1^3*y2^2 - y1^3*y3^2, y1^3*y2*y3 - y1^3*y3^2] 
g2m4 = 9 [y0^4*y2^4 + 2*y1^4*y2^3*y3 + 4*y0*y1^3*y2^2*y3^2 - y1^4*y2^2*y3^2 - 4*y0*y1^3*y2*y3^3 - 4*y1^4*y2*y3^3 + 3*y1^4*y3^4, y0^3*y1*y2^4 - y1^4*y2^3*y3 + 2*y1^4*y2^2*y3^2 - y1^4*y2*y3^3, y0^3*y1*y2^3*y3 + y1^4*y2^3*y3 + 3*y0*y1^3*y2^2*y3^2 - 3*y0*y1^3*y2*y3^3 - 3*y1^4*y2*y3^3 + 2*y1^4*y3^4, y0^2*y1^2*y2^4, y0^2*y1^2*y2^3*y3 - y1^4*y2^3*y3 + 2*y1^4*y2^2*y3^2 - y1^4*y2*y3^3, y0^2*y1^2*y2^2*y3^2 + 2*y0*y1^3*y2^2*y3^2 + y1^4*y2^2*y3^2 - 2*y0*y1^3*y2*y3^3 - 2*y1^4*y2*y3^3 + y1^4*y3^4, y0*y1^3*y2^4 + y1^4*y2^3*y3 - y0*y1^3*y2^2*y3^2 - 2*y1^4*y2^2*y3^2 + y1^4*y2*y3^3, y0*y1^3*y2^3*y3 + y1^4*y2^3*y3 - y0*y1^3*y2^2*y3^2 - 2*y1^4*y2^2*y3^2 + y1^4*y2*y3^3, y1^4*y2^4 - 2*y1^4*y2^3*y3 + y1^4*y2^2*y3^2] 
```
By inspection we recover the generators of graded ring associated to gg.

```python
V = V0, V1, V2, V3, V4 = ring( 'y0' ), ring( 'y1' ), ring( 'y0*y2^2+y1*y2^2-y1*y2*y3' ), ring( 'y0*y1*y2^2' ), ring( 'y0^2*y2^2+y1^2*y2^2-y1^2*y3^2' )

```
Since ff and gg parametrize conic bundles that are projectively isomorphic
we know that that there exists a relation between the monomials of bidegree (2,-4)
in the coordinate associated to gg.

```python
b = b0, b1, b2, b3, b4, b5, b6, b7, b8, b9 = [elt.subs( {u[i]:V[i] for i in range( 5 )} ) for elt in T2m4 ]
matb = sage_matrix( sage_QQ, SERing.get_matrix_P1xP1( b ) )
kerb = matb.transpose().right_kernel().matrix()
assert kerb * sage_vector( b ) == sage_vector( [0] )
assert str(kerb) == '[  1 1/2   0  -1   0  -1   0   0 1/2   0]'
assert 2 * b0 + b1 - 2 * b3 - 2 * b5 + b8 == 0

```

#### Computing projective isomorphisms

Let F and G be the maps whose components form a basis for the 
monomials of weight (1,0).

```python
F = [ elt.subs( {u[i]:U[i] for i in range( 5 )} ) for elt in T1m0]
G = [ elt.subs( {u[i]:V[i] for i in range( 5 )} ) for elt in T1m0]

```
We compute the inverse Q of the map G.

```python
y = ring( 'y0,y1,y2,y3' )
z = ring( 'z0,z1,z2,z3,z4,z5,z6,z7,z8,z9' )
t = ring( 't' )
id = [ G[i] * z[0] - z[i] * G[0] for i in range( 10 ) ] + [t * G[0] - 1]
I01 = sage_ideal( id ).elimination_ideal( [t, y[2], y[3] ] ).gens()
I23 = sage_ideal( id ).elimination_ideal( [t, y[0], y[1] ] ).gens()
I01 = [ elt for elt in I01 if elt.degree( y[0] ) == 1 and elt.degree( y[1] ) == 1 ][0]
I23 = [ elt for elt in I23 if elt.degree( y[2] ) == 1 and elt.degree( y[3] ) == 1 ][0]
Q0 = I01.coefficient( y[1] )
Q1 = -I01.coefficient( y[0] )
Q2 = I23.coefficient( y[3] )
Q3 = -I23.coefficient( y[2] )
Q = [Q0, Q1, Q2, Q3]
assert Q == ring("[-z9, -z8, -z8, -z6 - 2*z7 + z8 + z9]")

# check whether Q is the inverse of G
QoG = [q.subs( { z[i]:G[i] for i in range( 10 ) } ) for q in Q ]
gcd01 = sage_gcd( QoG[0], QoG[1] )
gcd23 = sage_gcd( QoG[2], QoG[3] )
QoG = [QoG[0] / gcd01, QoG[1] / gcd01, QoG[2] / gcd23, QoG[3] / gcd23]
assert QoG == [y[0], y[1], y[2], y[3]]

```
We compute the composition QoPoF of F with P and Q.

```python
c = c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12 = [ring( 'c' + str( i ) ) for i in range( 13 )]
dctZ = { u[i]:z[i] for i in range( 5 )}
dctP = { z[0]:c0 * u0 + c1 * u1, z[1]:c2 * u0 + c3 * u1, z[2]:c4 * u2, z[3]:c5 * u3 + c6 * u4 + c7 * u0 * u2 + c8 * u1 * u2, z[4]:c9 * u3 + c10 * u4 + c11 * u0 * u2 + c12 * u1 * u2 }
PoF = [ comp.subs( dctZ ).subs( dctP ) for comp in T1m0]
PoF = [ comp.subs( {u[i]:U[i] for i in range( 5 )} ) for comp in PoF]
PoF = [ comp / sage_gcd( PoF ) for comp in PoF]
QoPoF = [ comp.subs( {z[i]:PoF[i] for i in range( 10 )} ) for comp in Q]
gcd01 = sage_gcd( QoPoF[0], QoPoF[1] )
gcd23 = sage_gcd( QoPoF[2], QoPoF[3] )
QoPoF = [QoPoF[0] / gcd01, QoPoF[1] / gcd01, QoPoF[2] / gcd23, QoPoF[3] / gcd23]
print( 'QoPoF =', len( QoPoF ), QoPoF )

```
Output:
```
QoPoF = 4 [-c0*x1 - c1*x2, -c2*x1 - c3*x2, c2*c4*x0*x1 - c2*c4*x1^2 + c3*c4*x0*x2 - c2*c4*x1*x2 - c3*c4*x1*x2 - c3*c4*x2^2, -c0*c4*x0*x1 - c2*c4*x0*x1 + c0*c4*x1^2 + c2*c4*x1^2 - c1*c4*x0*x2 - c3*c4*x0*x2 + c0*c4*x1*x2 + c1*c4*x1*x2 + c2*c4*x1*x2 + c3*c4*x1*x2 + c1*c4*x2^2 + c3*c4*x2^2 - 2*c6*x0^2 - c10*x0^2 + 4*c6*x0*x1 + 2*c7*x0*x1 + 2*c10*x0*x1 + c11*x0*x1 - 2*c6*x1^2 - 2*c7*x1^2 - c10*x1^2 - c11*x1^2 + 4*c6*x0*x2 + 2*c8*x0*x2 + 2*c10*x0*x2 + c12*x0*x2 - 2*c5*x1*x2 - 4*c6*x1*x2 - 2*c7*x1*x2 - 2*c8*x1*x2 - c9*x1*x2 - 2*c10*x1*x2 - c11*x1*x2 - c12*x1*x2 - 2*c6*x2^2 - 2*c8*x2^2 - c10*x2^2 - c12*x2^2] 
```

We create a list rel_lst of equations in c.

```python
# use %paste to paste indented code into a Sage or Python terminal
b = T2m4
rel_g4m2 = 2 * b[0] + b[1] - 2 * b[3] - 2 * b[5] + b[8]
rel_g4m2 = rel_g4m2.subs( dctZ ).subs( dctP ).subs( {u[i]:U[i] for i in range( 5 )} )
rel_lst = []; x = ring( '[x0,x1,x2]' )
for exp in sage_Compositions( 4 + 3, length = 3 ): 
    rel_lst += [rel_g4m2.coefficient( {x[i]:exp[i] - 1 for i in range( 3 )} )]
rel_lst += [ ( c0 * c3 - c1 * c2 ) * c4 * ( c5 * c10 - c9 * c6 ) * ring('t') - 1 ]

```
We solve the equations rel_lst for c by computing a primary decomposition.

```python
prime_lst = sage_ideal( rel_lst ).elimination_ideal( t ).primary_decomposition()
assert len(prime_lst) == 4
prime0, prime1, prime2, prime3 = [prime.gens() for prime in prime_lst]
sol0 = sage_solve( [sage_SR( gen ) for gen in prime0], [sage_SR( elt ) for elt in c], solution_dict = True )
sol1 = sage_solve( [sage_SR( gen ) for gen in prime1], [sage_SR( elt ) for elt in c], solution_dict = True )
sol2 = sage_solve( [sage_SR( gen ) for gen in prime2], [sage_SR( elt ) for elt in c], solution_dict = True )
sol3 = sage_solve( [sage_SR( gen ) for gen in prime3], [sage_SR( elt ) for elt in c], solution_dict = True )
print('prime0 =', prime0)
print('sol0   =', sol0)
print('prime1 =', prime1)
print('sol1   =', sol1)
print('prime2 =', prime2)
print('sol2   =', sol2)
print('prime3 =', prime3)
print('sol3   =', sol3)

```
Output:
```
prime0 = [c8, c7, 2*c6 + c10, c5, c2, c1, c3*c11 - c0*c12, 2*c9*c10 - c11*c12, 2*c3*c4 - c12, 2*c0*c4 - c11] 
sol0   = {c10: -2*r2, c12: -4*r2*r3/r4, c6: r2, c8: 0, c2: 0, c4: r1, c0: 1/2*r4/r1, c11: r4, c7: 0, c9: r3, c3: -2*r2*r3/(r1*r4), c5: 0, c1: 0} 
prime1 = [c8, c7, 2*c6 + c10, c5, c3, c0, c1*c11 - c2*c12, 2*c9*c10 - c11*c12, 2*c2*c4 - c11, 2*c1*c4 - c12] 
sol1   = {c10: -2*r6, c12: r8, c6: r6, c8: 0, c2: -2*r6*r7/(r5*r8), c4: r5, c0: 0, c11: -4*r6*r7/r8, c7: 0, c9: r7, c3: 0, c5: 0, c1: 1/2*r8/r5} 
prime2 = [c8, c7, c6, 2*c5 + c9, c3, c0, c1*c11 - c2*c12, 2*c9*c10 - c11*c12, 2*c2*c4 - c11, 2*c1*c4 - c12] 
sol2   = {c10: r11, c12: r12, c6: 0, c8: 0, c2: -2*r10*r11/(r12*r9), c4: r9, c0: 0, c11: -4*r10*r11/r12, c7: 0, c9: -2*r10, c3: 0, c5: r10, c1: 1/2*r12/r9} 
prime3 = [c8, c7, c6, 2*c5 + c9, c2, c1, c3*c11 - c0*c12, 2*c9*c10 - c11*c12, 2*c3*c4 - c12, 2*c0*c4 - c11] 
sol3   = {c10: r15, c12: -4*r14*r15/r16, c6: 0, c8: 0, c2: 0, c4: r13, c0: 1/2*r16/r13, c11: r16, c7: 0, c9: -2*r14, c3: -2*r14*r15/(r13*r16), c5: r14, c1: 0} 
```

After normalization we may assume that some of the nonzero c is equal to 1. 
We add the corresponding equations in order to simplify the solutions.

```python
prime0 += [c0 - 1, c4 - 1]
prime1 += [c1 - 1, c4 - 1]
prime2 += [c1 - 1, c4 - 1]
prime3 += [c0 - 1, c4 - 1]
sol0 = sage_solve( [sage_SR( gen ) for gen in prime0], [sage_SR( elt ) for elt in c], solution_dict = True )
sol1 = sage_solve( [sage_SR( gen ) for gen in prime1], [sage_SR( elt ) for elt in c], solution_dict = True )
sol2 = sage_solve( [sage_SR( gen ) for gen in prime2], [sage_SR( elt ) for elt in c], solution_dict = True )
sol3 = sage_solve( [sage_SR( gen ) for gen in prime3], [sage_SR( elt ) for elt in c], solution_dict = True )
print('sol0   =', sol0)
print('sol1   =', sol1)
print('sol2   =', sol2)
print('sol3   =', sol3)

```
Output
```
sol0 = {c10: -2*r17, c12: -2*r17*r18, c6: r17, c8: 0, c2: 0, c4: 1, c0: 1, c11: 2, c7: 0, c9: r18, c3: -r17*r18, c5: 0, c1: 0} 
sol1 = {c10: -2*r19, c12: 2, c6: r19, c8: 0, c2: -r19*r20, c4: 1, c0: 0, c11: -2*r19*r20, c7: 0, c9: r20, c3: 0, c5: 0, c1: 1} 
sol2 = {c10: r22, c12: 2, c6: 0, c8: 0, c2: -r21*r22, c4: 1, c0: 0, c11: -2*r21*r22, c7: 0, c9: -2*r21, c3: 0, c5: r21, c1: 1} 
sol3 = {c10: r24, c12: -2*r23*r24, c6: 0, c8: 0, c2: 0, c4: 1, c0: 1, c11: 2, c7: 0, c9: -2*r23, c3: -r23*r24, c5: r23, c1: 0} 
```
We set declare the above solutions manualy:

```python    
r0, r1 = ring( 'r0,r1' )
sol0 = {c0:1, c1:0, c2:0     , c3:-r0*r1, c4:1, c5:0 , c6:r0, c7:0, c8:0, c9:r1   , c10:-2*r0, c11:2       , c12:-2*r0*r1}
sol1 = {c0:0, c1:1, c2:-r0*r1, c3:0     , c4:1, c5:0 , c6:r0, c7:0, c8:0, c9:r1   , c10:-2*r0, c11:-2*r0*r1, c12:2       }
sol2 = {c0:0, c1:1, c2:-r0*r1, c3:0     , c4:1, c5:r0, c6:0 , c7:0, c8:0, c9:-2*r0, c10:r1   , c11:-2*r0*r1, c12:2       }
sol3 = {c0:1, c1:0, c2:0     , c3:-r0*r1, c4:1, c5:r0, c6:0 , c7:0, c8:0, c9:-2*r0, c10:r1   , c11:2       , c12:-2*r0*r1}
sol_lst = [sol0, sol1, sol2, sol3]

```
For each of the four solutions we compose the corresponding
compatible reparametrization with gg. Notice that after 
we substitute a solution the components have a non-trivial 
greatest common divisor.

```python
y = ring( '[y0,y1,y2,y3]' )
gr_lst = []
for sol in sol_lst:
    gr = [ comp.subs( {y[i]:QoPoF[i] for i in range( 4 )} ).subs( sol ) for comp in gg ]
    gcd_gr = sage_gcd( gr )
    gr_lst += [[ comp / gcd_gr for comp in gr ]]
    print( gcd_gr )    

```
Output:
```
r0^2*r1^3*x1*x2^3 
r0^2*r1^3*x1^3*x2 
r0^2*r1^3*x0^2*x1^2 - 2*r0^2*r1^3*x0*x1^3 + r0^2*r1^3*x1^4 - 2*r0^2*r1^3*x0*x1^2*x2 + 2*r0^2*r1^3*x1^3*x2 + r0^2*r1^3*x1^2*x2^2 
r0^2*r1^3*x0^2*x2^2 - 2*r0^2*r1^3*x0*x1*x2^2 + r0^2*r1^3*x1^2*x2^2 - 2*r0^2*r1^3*x0*x2^3 + 2*r0^2*r1^3*x1*x2^3 + r0^2*r1^3*x2^4 
```

In order to test our result we compute the implicit equations for the image Y of gg. 
However, this is not required for computing the projective isomorphisms between X and Y.

```python
z = ring( 'z0,z1,z2,z3,z4,z5,z6,z7,z8,z9' )
y = ring( 'y0,y1,y2,y3' )
igg = SERing.R.ideal( [ z[i] - gg[i] for i in range( 10 )  ] ).elimination_ideal( [y[i] for i in range( 4 )] )
print( 'igg =', list(igg.gens()) )

```
Output:
```
igg = [z1*z8 - z8^2 + z2*z9, z6*z7 + z0*z8 - z2*z8 - z3*z8 + z5*z8 - z8^2 + z2*z9 - z4*z9 + z5*z9 - z6*z9, z5*z7 + z2*z8 - z4*z8 - z5*z9, z4*z7 + z0*z8 - z8^2 + z2*z9 - z4*z9, z3*z7 + z0*z8 - z2*z8 + z4*z8 - z8^2 + z2*z9 - z3*z9 - z4*z9 + z5*z9, z2*z7 - z0*z8 + z8^2 - z2*z9, z1*z7 - z7*z8 + z0*z9 - z1*z9, z0*z7 + z7*z8 - z0*z9, z4*z5 - z5^2 - z0*z6 + z2*z6 - z4*z8 + z5*z8 + z3*z9 - z6*z9, z2*z5 - z2*z8 - z6*z8 + z5*z9, z1*z5 - z5^2 + z2*z6, z0*z5 - z2*z8 - z3*z8 + z5*z9, z4^2 - z5^2 + z2*z6, z3*z4 + z3*z5 - z4*z6, z2*z4 - z2*z8 - z3*z8 + z5*z9, z1*z4 - z5^2 + z2*z6 - z4*z8 + z5*z8 + z3*z9 - z6*z9, z0*z4 - z2*z8 + z5*z9, z2*z3 - z0*z6 - z4*z8 + z5*z8 + z3*z9 - z6*z9, z1*z3 - z3*z5 - z1*z6 + z4*z6, z0*z3 + z4*z8 - z5*z8 - z3*z9 + z6*z9, z2^2 - z5*z8, z1*z2 - z2*z8 + z5*z9, z0*z2 - z4*z8, z1^2 - z5^2 + z2*z6 + z5*z8 - z8^2 + z2*z9 - z6*z9, z0*z1 - z0*z8 + z4*z9, z0^2 - z8^2 + z2*z9, z5^2*z8 - z2*z6*z8 - z5*z8^2 + z2*z8*z9 + z6*z8*z9 - z5*z9^2, z3*z5*z8 - z4*z6*z8 - z3*z8^2 + z6*z8^2 + z0*z6*z9 - z2*z6*z9 + z4*z8*z9 - z5*z8*z9 - z3*z9^2 + z6*z9^2, z3^2*z8 + 2*z0*z6*z8 - z2*z6*z8 + 2*z4*z8^2 - 2*z5*z8^2 - 2*z3*z5*z9 + z5*z6*z9 - 2*z3*z8*z9 + 2*z6*z8*z9, 2*z3*z5^2 - z5^2*z6 - 2*z0*z6^2 + z2*z6^2 - 4*z4*z6*z8 + 3*z5*z6*z8 - 2*z3*z8^2 + 2*z6*z8^2 - z3^2*z9 + 2*z0*z6*z9 - 2*z2*z6*z9 + 4*z3*z6*z9 - 3*z6^2*z9 + 2*z4*z8*z9 - 2*z5*z8*z9 - 2*z3*z9^2 + 2*z6*z9^2] 
```

For each gr in gr_lst we recover 
the corresponding projective isomorphism between X and Y
in terms of a parametrized matrix U.
 
```python
# compute the coefficient matrix of ff and its kernel
mff = SERing.get_matrix_P2( ff )
kff = mff.right_kernel_matrix().T
assert ( mff * kff ).is_zero()

# Compute isomorphism U for each gr in gr_lst
for gr in gr_lst:

    # compute the coefficient matrix of gr
    mgr = SERing.get_matrix_P2( gr )
    mgk = mgr * kff
    assert mgk.is_zero()  # because the surfaces X and Y in P^9 are already linearly normal
    
    # compute isomorphism U
    Ef = sage_matrix( sage_QQ, mff.rows() + kff.T.rows() )
    Egr = sage_matrix( mgr.rows() + kff.T.rows() )
    UpI = Egr * Ef.inverse()
    assert ( UpI.submatrix( 10, 10 ) - sage_identity_matrix( 5 ) ).is_zero()
    U = UpI.submatrix( 0, 0, 10, 10 )
    print( 'U =', U.dimensions(), list( U ) )
    
    # check U by using the equations of Y
    Uff = list( U * sage_vector( ff ) )
    iggs = igg.subs( {z[i]:Uff[i] for i in range( 10 )} )
    assert iggs.is_zero()

```
Output:
```
U = (10, 10) [(-r0, r0^2*r1 - r0*r1 - 2*r0, 2*r0, 2*r0^2*r1 - r0*r1 - r0, -2*r0^2*r1 + r0*r1 + 2*r0, r0^2*r1, -2*r0^2*r1, r0^2*r1, 0, 0), (-r0, 2*r0^2*r1 - 2*r0*r1 - 2*r0, 2*r0, -r0^3*r1^2 + 2*r0^2*r1^2 + 4*r0^2*r1 - r0*r1^2 - 2*r0*r1 - r0, -4*r0^2*r1 + 2*r0*r1 + 2*r0, -2*r0^3*r1^2 + 2*r0^2*r1^2 + 2*r0^2*r1, 2*r0^3*r1^2 - 2*r0^2*r1^2 - 4*r0^2*r1, 2*r0^2*r1, -r0^3*r1^2, 2*r0^3*r1^2), (0, r0^2*r1, 0, 2*r0^2*r1, -2*r0^2*r1, r0^2*r1, -2*r0^2*r1, r0^2*r1, 0, 0), (0, -r0^2*r1, 0, r0^3*r1^2 - 2*r0^2*r1^2 - 2*r0^2*r1, 2*r0^2*r1, r0^3*r1^3 + 2*r0^3*r1^2 - r0^2*r1^3 - 2*r0^2*r1^2 - r0^2*r1, -2*r0^3*r1^2 + 2*r0^2*r1^2 + 2*r0^2*r1, -r0^2*r1, r0^3*r1^3 + r0^3*r1^2, -r0^3*r1^3 - 2*r0^3*r1^2), (0, r0^2*r1, 0, -r0^3*r1^2 + r0^2*r1^2 + 2*r0^2*r1, -2*r0^2*r1, -2*r0^3*r1^2 + r0^2*r1^2 + r0^2*r1, 2*r0^3*r1^2 - r0^2*r1^2 - 2*r0^2*r1, r0^2*r1, -r0^3*r1^2, 2*r0^3*r1^2), (0, 0, 0, -r0^3*r1^2, 0, -2*r0^3*r1^2, 2*r0^3*r1^2, 0, -r0^3*r1^2, 2*r0^3*r1^2), (0, -r0^2*r1, 0, 2*r0^3*r1^2 - 2*r0^2*r1^2 - 2*r0^2*r1, 2*r0^2*r1, 2*r0^3*r1^3 + 4*r0^3*r1^2 - r0^2*r1^3 - 2*r0^2*r1^2 - r0^2*r1, -4*r0^3*r1^2 + 2*r0^2*r1^2 + 2*r0^2*r1, -r0^2*r1, 2*r0^3*r1^3 + 2*r0^3*r1^2, -2*r0^3*r1^3 - 4*r0^3*r1^2), (r0 - 1, -r0^2*r1 + 2*r0*r1 + 2*r0 - r1 - 1, -2*r0 + 1, -2*r0^2*r1 + 2*r0*r1 + r0, 2*r0^2*r1 - 2*r0*r1 - 2*r0, -r0^2*r1, 2*r0^2*r1, -r0^2*r1, 0, 0), (-r0, -2*r0, 2*r0, -r0, 2*r0, 0, 0, 0, 0, 0), (2*r0 - 2, -r0^2*r1 + 2*r0*r1 + 4*r0 - r1 - 2, -4*r0 + 2, -2*r0^2*r1 + 2*r0*r1 + 2*r0, 2*r0^2*r1 - 2*r0*r1 - 4*r0, -r0^2*r1, 2*r0^2*r1, -r0^2*r1, 0, 0)] 
U = (10, 10) [(0, r0^2*r1, 0, 2*r0^2*r1 - r0*r1 - r0, -2*r0^2*r1, r0^2*r1 - r0*r1 - 2*r0, -2*r0^2*r1 + r0*r1 + 2*r0, r0^2*r1, -r0, 2*r0), (-r0^3*r1^2, -2*r0^3*r1^2 + 2*r0^2*r1^2 + 2*r0^2*r1, 2*r0^3*r1^2, -r0^3*r1^2 + 2*r0^2*r1^2 + 4*r0^2*r1 - r0*r1^2 - 2*r0*r1 - r0, 2*r0^3*r1^2 - 2*r0^2*r1^2 - 4*r0^2*r1, 2*r0^2*r1 - 2*r0*r1 - 2*r0, -4*r0^2*r1 + 2*r0*r1 + 2*r0, 2*r0^2*r1, -r0, 2*r0), (0, r0^2*r1, 0, 2*r0^2*r1, -2*r0^2*r1, r0^2*r1, -2*r0^2*r1, r0^2*r1, 0, 0), (r0^3*r1^3 + r0^3*r1^2, r0^3*r1^3 + 2*r0^3*r1^2 - r0^2*r1^3 - 2*r0^2*r1^2 - r0^2*r1, -r0^3*r1^3 - 2*r0^3*r1^2, r0^3*r1^2 - 2*r0^2*r1^2 - 2*r0^2*r1, -2*r0^3*r1^2 + 2*r0^2*r1^2 + 2*r0^2*r1, -r0^2*r1, 2*r0^2*r1, -r0^2*r1, 0, 0), (-r0^3*r1^2, -2*r0^3*r1^2 + r0^2*r1^2 + r0^2*r1, 2*r0^3*r1^2, -r0^3*r1^2 + r0^2*r1^2 + 2*r0^2*r1, 2*r0^3*r1^2 - r0^2*r1^2 - 2*r0^2*r1, r0^2*r1, -2*r0^2*r1, r0^2*r1, 0, 0), (-r0^3*r1^2, -2*r0^3*r1^2, 2*r0^3*r1^2, -r0^3*r1^2, 2*r0^3*r1^2, 0, 0, 0, 0, 0), (2*r0^3*r1^3 + 2*r0^3*r1^2, 2*r0^3*r1^3 + 4*r0^3*r1^2 - r0^2*r1^3 - 2*r0^2*r1^2 - r0^2*r1, -2*r0^3*r1^3 - 4*r0^3*r1^2, 2*r0^3*r1^2 - 2*r0^2*r1^2 - 2*r0^2*r1, -4*r0^3*r1^2 + 2*r0^2*r1^2 + 2*r0^2*r1, -r0^2*r1, 2*r0^2*r1, -r0^2*r1, 0, 0), (0, -r0^2*r1, 0, -2*r0^2*r1 + 2*r0*r1 + r0, 2*r0^2*r1, -r0^2*r1 + 2*r0*r1 + 2*r0 - r1 - 1, 2*r0^2*r1 - 2*r0*r1 - 2*r0, -r0^2*r1, r0 - 1, -2*r0 + 1), (0, 0, 0, -r0, 0, -2*r0, 2*r0, 0, -r0, 2*r0), (0, -r0^2*r1, 0, -2*r0^2*r1 + 2*r0*r1 + 2*r0, 2*r0^2*r1, -r0^2*r1 + 2*r0*r1 + 4*r0 - r1 - 2, 2*r0^2*r1 - 2*r0*r1 - 4*r0, -r0^2*r1, 2*r0 - 2, -4*r0 + 2)] 
U = (10, 10) [(0, 0, 0, r0^2*r1 - r0*r1, 0, -r0*r1 - r0, r0*r1, 0, 0, 0), (0, -r0^3*r1^2 + 2*r0^2*r1^2 - r0*r1^2, 0, 2*r0^2*r1^2 + 2*r0^2*r1 - 2*r0*r1^2 - 2*r0*r1, -2*r0^2*r1^2 + 2*r0*r1^2, -r0*r1^2 - 2*r0*r1 - r0, 2*r0*r1^2 + 2*r0*r1, -r0*r1^2, 0, 0), (0, 0, 0, r0^2*r1, 0, 0, 0, 0, 0, 0), (r0^3*r1^3 - r0^2*r1^3, r0^3*r1^3 + r0^3*r1^2 - 2*r0^2*r1^3 - 2*r0^2*r1^2, -r0^3*r1^3 + 2*r0^2*r1^3, -r0^2*r1^3 - 2*r0^2*r1^2 - r0^2*r1, 2*r0^2*r1^3 + 2*r0^2*r1^2, 0, 0, 0, 0, 0), (0, -r0^3*r1^2 + r0^2*r1^2, 0, r0^2*r1^2 + r0^2*r1, -r0^2*r1^2, 0, 0, 0, 0, 0), (0, -r0^3*r1^2, 0, 0, 0, 0, 0, 0, 0, 0), (2*r0^3*r1^3 - r0^2*r1^3, 2*r0^3*r1^3 + 2*r0^3*r1^2 - 2*r0^2*r1^3 - 2*r0^2*r1^2, -2*r0^3*r1^3 + 2*r0^2*r1^3, -r0^2*r1^3 - 2*r0^2*r1^2 - r0^2*r1, 2*r0^2*r1^3 + 2*r0^2*r1^2, 0, 0, 0, 0, 0), (0, 0, 0, -r0^2*r1 + 2*r0*r1 - r1, 0, 2*r0*r1 + r0 - 2*r1 - 1, -2*r0*r1 + 2*r1, 0, -r1 - 1, 2*r1 + 1), (0, 0, 0, 0, 0, -r0, 0, 0, 0, 0), (0, 0, 0, -r0^2*r1 + 2*r0*r1 - r1, 0, 2*r0*r1 + 2*r0 - 2*r1 - 2, -2*r0*r1 + 2*r1, 0, -r1 - 2, 2*r1 + 2)] 
U = (10, 10) [(0, -r0*r1 - r0, 0, r0^2*r1 - r0*r1, r0*r1, 0, 0, 0, 0, 0), (0, -r0*r1^2 - 2*r0*r1 - r0, 0, 2*r0^2*r1^2 + 2*r0^2*r1 - 2*r0*r1^2 - 2*r0*r1, 2*r0*r1^2 + 2*r0*r1, -r0^3*r1^2 + 2*r0^2*r1^2 - r0*r1^2, -2*r0^2*r1^2 + 2*r0*r1^2, -r0*r1^2, 0, 0), (0, 0, 0, r0^2*r1, 0, 0, 0, 0, 0, 0), (0, 0, 0, -r0^2*r1^3 - 2*r0^2*r1^2 - r0^2*r1, 0, r0^3*r1^3 + r0^3*r1^2 - 2*r0^2*r1^3 - 2*r0^2*r1^2, 2*r0^2*r1^3 + 2*r0^2*r1^2, 0, r0^3*r1^3 - r0^2*r1^3, -r0^3*r1^3 + 2*r0^2*r1^3), (0, 0, 0, r0^2*r1^2 + r0^2*r1, 0, -r0^3*r1^2 + r0^2*r1^2, -r0^2*r1^2, 0, 0, 0), (0, 0, 0, 0, 0, -r0^3*r1^2, 0, 0, 0, 0), (0, 0, 0, -r0^2*r1^3 - 2*r0^2*r1^2 - r0^2*r1, 0, 2*r0^3*r1^3 + 2*r0^3*r1^2 - 2*r0^2*r1^3 - 2*r0^2*r1^2, 2*r0^2*r1^3 + 2*r0^2*r1^2, 0, 2*r0^3*r1^3 - r0^2*r1^3, -2*r0^3*r1^3 + 2*r0^2*r1^3), (-r1 - 1, 2*r0*r1 + r0 - 2*r1 - 1, 2*r1 + 1, -r0^2*r1 + 2*r0*r1 - r1, -2*r0*r1 + 2*r1, 0, 0, 0, 0, 0), (0, -r0, 0, 0, 0, 0, 0, 0, 0, 0), (-r1 - 2, 2*r0*r1 + 2*r0 - 2*r1 - 2, 2*r1 + 2, -r0^2*r1 + 2*r0*r1 - r1, -2*r0*r1 + 2*r1, 0, 0, 0, 0, 0)] 
```
