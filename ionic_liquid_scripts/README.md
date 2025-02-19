Scripts to generate a set of ionic liquids

The files this notebook generates are:
 
 - il0 through il9 have a series of modified imidazole with 10
   different functionalized imidazole along with bistriflimide (TFSI)
   as the anion (300 of each)

 - a functionalized imidazole plus TFSI, in combinations with various
   concentrations of tetraglyme.
   - il_neutral10: 10% tetraglyme
   - il_neutral20: 20% tetraglyme
   - il_neutral30: 30% tetraglyme
   - il_neutral40 - 40% tetraglyme
   - il_neutral50 - 50% tetraglyme

 - des12-16 are a series of deep eutectic solvent simulations, with
   choline chloride as the minor component, with ethylene glycol as
   the major component.
 
    - des12: 250/500
    - des13: 200/600
    - des14: 175/700
    - des15: 150/750
    - des16: 125/750
    
 - Note that to get packmol to run in a reasonable amount of time,
this scripts use a really big box, so there are some enormous gaps in
many of the structures that would need quite a lot of NPT
equilibration.
 
 - Also, the placement isn’t smart enough to put the anions and
cations close to each other, and instead the positions are randomized;
the files have been energy minimized, but there will likely be quite a
bit of heat/box shrinkage as anions and cations get closer.
 
Note: it's known that OpenFF currently gives very high viscosities for
ionic liquds (there is no charge scaling), so be aware.