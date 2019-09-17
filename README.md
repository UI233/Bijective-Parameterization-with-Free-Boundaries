# Bijective-Parameterization-with-Free-Boundaries

## Introduction

This project is a rough implementation of paper [bijective](http://faculty.cs.tamu.edu/schaefer/research/bijective.pdf) used for mesh parameterization.

This project can be successfully complied under G++ 7.3.0.

## Usage

*./bijective filename [optimization error]*

The default error is 0.005.

## Implementation

The initial parameterization is obtained by Tutte's method with circular boudary.
And the energy is optimized by L-BFGS.

## Requirement

The mesh to be parameterized should be equiped with a valid boundary and it should have no hole.
