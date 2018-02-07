# MeshStitching

This applications transfer the salient geometric features from a target mesh to a source mesh while trying to preserve the input triangulation.

This application was used to solve misregistration and adapt to  geometric and topological changes in:

https://dl.acm.org/citation.cfm?id=3073679 

To compile in linux, install PNG and Embree on your machine. Then, run the provided Makefile.

In Windows, move the \include and \lib folders from 4Windows.zip to the project root directory. Build the project from MeshStitching.sln, and copy the content of \dll to the \x64\Release.

Run the application from command line as:

Stitching --source Source.ply --target Target.ply

Add the parameter --output Sticthed.ply to run the application in default configuration and skip the user interface.

For an explanatory video on the usage of the interface please visit,

http://www.cs.jhu.edu/~fpradan1/code/