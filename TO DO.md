TO DO

- resolve density function better
- record the number of collisions an electron has undergone
- make an electron transport code by
    - calculating all excitation rates
    - callback at 10km intervalls to save results etc into histogram to save fluxes
- more energies 12 16 kev


- markov chain

- cdf of radial distribution for 1/e calculation, look into kernel density estimation

- height distirbution of spwaning  height of escaping electorns


__________________________________
Discussions with Tima (Mai 2026)

- check propagation without mangetic field to check
- Tima is working in space domain, not time, i.e. sampled number of mean free paths is calculated using 
    r = sum(dz / sum(sigma n_i) coss(theta))

- check in boris mover, that the condition of distance between 2 points << mfp is satisfied
- check step size for small pitch angle in boris mover!
- precalcute and tabulate phase functions, cross sections etc
- in ionizing collisions: sample secondary energy => momentun and energy conservation should give deflection angle
- inelastic collisions: shifted in energy. Björn has other opinions, see itakawa paper on measurements
    f(E) = 0.43 E (ln(0.1 E))^2
- initial conditions: insert electrons from point source with isotropic pitch angle distribution instead of common gyrocenter?
- isotropic spherical and isotropic in a plane is not the same distribution? look up!
- check mirror point and free path: much shorter delta t should result in longer path length at the mirror point? could be bc tima is working in the space domaing ( see point 2), here in the time domain it should be less of a problem. still, check!!
- energy cost per ionization:
    - plot cost of every interaction (ie by cross section): cost should be altitude dependent, and vary quite a bit by species.
    - only when summed up over all species, it should be roughly constant.
    - it should be closer to 32 eV if all cross sections are used for total cross section (also double ionization).
- 200m steps in atmosphere, interpolate linearly
- to resolve collisions quicker, first choose gas, then class of scattering (elastic, inelastic or ionizing) and then the specific process => less cross sections to evaluation means large speed ups.
- precalculate random numbers
- tabulate and linearly inperpolate phase angles
- statistics in cross sections can be dramatically improved:
    - trace electrons => electrons flux will still be standart statistics
    - at every collision, record all cross sections as weights => increases statistics
    - still need to select one particular scattering process to continue the electron path
- phase functions have three terms, can be used to quickly probe pahse function of O2, N2
- O phase functions can be probed by 2 random numbers (see notes)
- check ionization corss section and secondary electron distribution ( Green at al 1972, in zotero)

verification methods:
- mirroring altitude
- compare backscattered energy with Reimei satellite or vision rocket (if difference is less than factor 2, very good)
- compare with timas profiles (see email)
- compare vertical, non-converging magnetic field with dipole field for istropic flux
    => difference shoudl be around factor of 1.3-1.4


