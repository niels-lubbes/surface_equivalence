# Orbital 


## Introduction

Surface-Equivalence is a library for computing projective equivalences between surfaces, if such equivalence exists.

This library depends on [SageMath](https://SageMath.org) libraries and uses [Maple](https://www.maplesoft.com) for Groebner basis computations.

## Installation

* Install Sage from [SageMath](https://SageMath.org).
We assume that `sage` is accessible from your commandline interface.

* Install [Maple](https://www.maplesoft.com).
We assume that `maple` is accessible from your commandline interface.

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


### Example 1: Projective equivalence of surfaces parametrized by bihomogeneous polynomials

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

    