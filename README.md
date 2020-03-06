# MPI_PhotoFilter

App for applying filters on PNM and PGM images through a distributed system, using OpenMPI, in C. The app can apply one or more of the following filters: smooth, blur, sharpen, mean, emboss. The purpose of using MPI is to divide the image matrix between a different number of processes, to achieve high scalability in order to reduce the processing time. For the color, PNM images, a derived MPI struct type was constructed.

An experiment could be run using the following format:

mpirun -np N ./filter image_in.pnm image_out.pnm filter1 filter2 ... filterX
