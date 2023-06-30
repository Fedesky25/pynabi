# PynAbi Changelog

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