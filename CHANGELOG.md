# PynAbi Changelog

## 0.1.2

 * Feature: `SCFMixing` and `SCFDirectMinimization` require `Tolerance` to be specified
 * Feature: when `NonSelfConsistentCal` is not explicitly set, PynAbi will check the existence of `Tolerance` (Abinit implicitly assumes a SCF calculation)
 * Documentation: for `Smearing` and `Metal`
 * Fix: `typat`, `natom`, `nband` spelled correctly
 * Fix: compatibility check also on base dataset
 * Fix: `SymmetricGrid` now correctly writes `kptopt` and `shiftk` 

## 0.1.1

 * Fix: `CanDelay` in same dataset of one of its `Delayed` properties must delay that property
 * Feature: you can now append stampable to groups of datsets to cut down repeated code

## 0.1.0

 * __Breaking__: refactor of `kspace` -> there are now three specialized (mutually exclusive) classes `ManualGrid`, `SymmetricGrid`, `AutomaticGrid`, `Path`
 * Refactor using `IndexedWithDefault` to cut down boilerplate code
 * Refactor of the delaying code: `DelayedInfo` has been added and a delayed value now can be of type different from input
 * `SCFDirectMinimization` does not take steps as argument (already taken by `MaxSteps`)
 * Reject delayed value only if already defined in same dataset

## 0.0.3 and previous

 * File handling thorugh `AbIn` and `AbOut`
 * Atom basis and lattices with many helper functions
 * Definition of the reciprocal space
 * Occupation of bands and spin types
 * Basic relaxation: implemention of `ionmov`, `occopt` with their asscoated variables
 * Units of energy and length of Abinit
 * (non) SCF calculation
 * Partial check for imcompatible definitions and values