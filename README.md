# Image-resize-using-Bicubic-Interpolation

Description:

The bicubic interpolation algorithm resizes images in PNG format, allowing both reduction and enlargement of the image. The algorithm operates on the original resolution of the image and can adjust its size by a user-specified factor.

Need for Parallelization:

Parallelization could be beneficial due to the significant weight of the imageInterpolate function in the execution time. This function represents approximately 80% of the total execution time and is called three times for each color channel.

Implementation Descriptions for Parallel Variants:

OpenMP:

Run: ./resize Example.png <factor>

Iterations of the for loops in the imageInterpolate function are distributed among multiple threads, collapsing them into a single loop.
The "arr" matrix is private since it is the local matrix of each process.
After testing various scheduling types, it was concluded that static scheduling is the most suitable as iterations are balanced.
Pthread:

Run: ./resize Example.png <factor> <num_threads>

For simplicity, the imageInterpolate function is not called three times in the main function. Instead, the processing of the three color channels is performed in a function specific to the threads.
pthread_struct structure: This structure contains parameters needed for the function executed by each thread, including thread ID, total available cores, image dimensions, scaling factor, mutexes for accessing matrices calculating pixel neighbors, matrices for color channels, and matrices resulting from interpolation.
Each thread receives a portion of the image and performs bicubic interpolation for each color channel. Mutexes are used to avoid data access contention.
MPI:

Run: mpirun -np <num_proc> ./resize Example.png <factor>

The implementation idea starts with the MPI_Scatterv function, a version of MPI_Scatter that allows the distribution of variable-sized data chunks among processes.
The master process calculates the size of each color channel matrix piece for each process, along with offsets indicating where each process's processing will start. It broadcasts this data.
Each process calls imageInterpolate to process its matrix piece, sending the processed matrix back to the master process, which assembles the results.
MPI_Pthread:

Run: mpirun -np <num_proc> ./resize Example.png <factor> <num_threads>

Image portions for each color channel are distributed among processes using MPI, and each thread processes a part of the assigned piece.
MPI_OpenMP:

Run: mpirun -np <num_proc> ./resize Example.png <factor>

Image portions for each color channel are distributed among processes using MPI, and iterations in the imageInterpolate function are collapsed and evenly distributed among threads.

Best speedup : Pthread and OpenMP




