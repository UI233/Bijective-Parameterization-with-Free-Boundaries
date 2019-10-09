# Bijective Parameterization

The goal of this repository is to replicate the algorithm mentioned in [bijective](http://faculty.cs.tamu.edu/schaefer/research/bijective.pdf) parameterization. Making use of this project, you can parameterize a 3D mesh with specified boundary into a planar mesh with free boundary.

## Introduction

### What is Parameterization?

Parameterization is a map which maps the surface of a 3D object into the 2D plane. Formally speaking, it is a function f: M -> R<sup>2</sup> where M is a surface in R<sup>3</sup>.

A typical example of parameterization is the map of the *Earth* which projects the surface of sphere (our planet) into a 2D plane(a sheet of paper). As you can see in this example, an application of parameterization is for applying plausible texture to the surface of a 3D object.

### Brief Introduction To Bijective Parameterization

For every triangle in mesh, parameterization can be viewed as a linear transformation locally, which maps from R<sup>2</sup> to R<sup>2</sup>, so local 2 singular values are well-defined. Making use of these 2 singular values, the distortion energy, which measures how largely the triangle is distorted by the parameterization, can be defined over every single triangle (eg. MIPS energy &sigma;<sub>1</sub>/&sigma;<sub>2</sub> + &sigma;<sub>2</sub>/&sigma;<sub>1</sub>). Plus, to avoid boundary edges from degeneration, s specified energy is proposed over every boundary edge in this paper. In order to an obtain nice parameterization, the sum of these two kinds of energy should be minimal.

In my implementation, the energy is optimized using L-BFGS algorithm.

## Prerequisites

This project was developed in  Ubuntu 18.04 using *OpenMesh* and *Eigen* library so you can directly build this project using *CMake*(version >= 3.0) in Linux. 

In other platforms with compiler supporting C++14(or higher), you have to replace the *OpenMesh* library in my repository with corresponding version for your environment and modify the path of *OpenMesh* library in CMakeLists.txt to build this project.

## Constrains About Mesh

The mesh to be parameterized should have boundary, which is better homeomorphic to a circle.
Moreover the orientation of every triangle needs to be counterclockwise or the program may get stuck during optimization step.

## Getting Started

### How To Build

Once these prerequisites are satisfied, you can build this project directly via CMake.

![Build Project](https://github.com/UI233/imageTemp/blob/master/build.gif)

### Usage

`Bijective FILE [ERROR]`

eg. `Bijective cow.obj 0.05`

ERROR term describes the terminating condition in the stage of optimization .

Default error is 0.05.

## Demo

<table border="1">
    <tr>
        <td>
        Raw
        </td>
        <td>
            Parameterization
        </td>
    </tr>
    <tr>
        <td>
            <img src = "https://github.com/UI233/imageTemp/blob/master/cow.jpg">
        </td>
        <td>
            <img src ="https://github.com/UI233/imageTemp/blob/master/param.jpg">
        </td>
    </tr>
</table>

## References

\[1\] [Bijective-Parameterization-with-Free-Boundaries](http://faculty.cs.tamu.edu/schaefer/research/bijective.pdf)

\[2\] [L-BFGS](https://en.wikipedia.org/wiki/Limited-memory_BFGS)