#include <iostream>
#include <cmath>
#include <fstream> //for testing with excel
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;


int main(int argc, char* argv[])
{

    ofstream testingdoc;
    testingdoc.open("testingdoc.csv");
    if (argc < 2)
    {
        printf("Enter value for grid_size a parameter at run time.\n");
        return 0;
    }
    int fineness = atoi(argv[1]); //measure of how coarse the grid is.
    double h = 1.0/fineness;
    double epsillon = 1;
    double w = 1.87;
    double errorTolerance = 1E-8;
    double old;
    int converged = 0;
    //initialise and obtain all parameters.
    int a, b, x1, x2, y1, y2, xDesired, yDesired;
    double q1, q2, inputx1, inputx2, inputy1, inputy2,inputxDesired,inputyDesired, timeStart, timeFinish;

    //INPUT INITIAL PARAMETERS//
    a = 4;
    b = 6;
    q1 = 1;
    inputx1 = 3.1415926535897932;
    inputy1 = 2.7182818284590452;
    q2 = -2;
    inputx2 = 1.6180339887487848;
    inputy2 = 4.6692016091029906;
    inputxDesired = 3.0;
    inputyDesired = 3.0;

    //PARAMETER MODIFICATION FOR FUTURE USAGES, DESCRETATION, ETC..//
    a = fineness * a;
    b = fineness * b;
    inputxDesired = (fineness * inputxDesired);
    xDesired = (int)inputxDesired;
    inputyDesired = (fineness * inputyDesired);
    yDesired = (int)inputyDesired;
    inputx1 = (fineness * inputx1);
    x1 = (int)inputx1;
    inputy1 = (fineness * inputy1);
    y1 = (int)inputy1;
    inputx2 = (fineness * inputx2);
    x2 = (int)inputx2;
    inputy2 = (fineness * inputy2);
    y2 = (int)inputy2;


    MPI_Init(NULL,NULL);
    timeStart = MPI_Wtime();
    //Create data for the grid of potentials, grid of densities and a checker pattern for future decomposition.
    double potential;
    double* density = new double[(b+2)*(a+2)];
    double* potentialGrid = new double[(b+2)*(a+2)];
    int* checkerPattern = new int[(b+2)*(a+2)];
    //allocate initial values to above arrays.
    for(int i = 0; i<((b+2)*(a+2)); i++)
    {
        potentialGrid[i] = 0.0;
        density[i] = 0.0;
        checkerPattern[i] = 0;
    }
    //Fill in checkerPattern array for decomposition of problem//
    bool oddNumber = true;
    for(int y = 1; y<b+1; y++)
    {
        for(int x = 1; x<a+1; x++)
        {
            if (oddNumber == true)
            {
                checkerPattern[(y)*(a+2)+x] = 1;
                oddNumber = false;
            }
            else
            {
                oddNumber = true;
            }
        }
    }
    //Obtain initial parameters//
    int world_size,world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    //Error if bad number of processors//
    if ((b+2)*fineness%world_size != 0 && world_rank == 0)
    {
        cout << "Uneven distribution of processors, " << world_size <<" for number of grid points, " << (a+2)*(b+2)*fineness << ", change number of processors to a factor of this number." << endl;
        MPI_Finalize();
        return 0;
    }
    //Assign initial densities//
    density[(b-y1+1)*(a+2)+(x1+1)] = q1/pow(h,2);  //note: goes from 0-4 for a,b=5
    density[(b-y2+1)*(a+2)+(x2+1)] = q2/pow(h,2);
    //Set up work per processor//
    int blocksize = ((a+2)*(b+2))/world_size;
    int* blockLength = new int[world_size];
    int* blockDisplacement = new int[world_size];
    //Create arrays for holding information on displacements and block lengths for MPI function usages.//
    for (int i = 0; i < world_size; i++)
    {
        blockLength[i] = blocksize;
        blockDisplacement[i] = blocksize*i;
    }

    bool conclude = false;
    while(conclude == false)
    {
        //decompose grid problem into alternating spaces to avoid overlap in parallel processes.
        //first part
        for(int y = b*(world_rank*1.0/world_size); y<b*((world_rank+1.0)/world_size); y++) //Distribute work evenly over the y-axis//
        {
            for(int x = 1; x<a+1; x++)
            {
                if (checkerPattern[(y+1)*(a+2)+x] == 1)
                {
                    //Implementation of algorithm from question//
                    //Chose to maintain (a+2) and (y+1) to clarify the differences between the original problem and the 0 surrounding buffer I introduced.
                    potential = 0.25*(potentialGrid[(y+1)*(a+2)+(x-1)]+potentialGrid[(y+1)*(a+2)+(x+1)]+potentialGrid[(y+1-1)*(a+2)+(x)]+potentialGrid[(y+1+1)*(a+2)+(x)]+(pow(h,2)*density[(y+1)*(a+2)+(x)])/epsillon);
                    old = potentialGrid[(y+1)*(a+2)+x];
                    potentialGrid[(y+1)*(a+2)+x] = potentialGrid[(y+1)*(a+2)+x]+w*(potential - potentialGrid[(y+1)*(a+2)+x]);
                    //Convergence test
                    if (abs(old - potentialGrid[(y+1)*(a+2)+x]) < errorTolerance)
                    {
                        converged += 1;
                    }
                    else
                    {
                        converged = 0;
                    }
                }
            }
        }
        MPI_Allgatherv(&potentialGrid[(blocksize)*world_rank],blocksize,MPI_DOUBLE,&potentialGrid[0],blockLength,blockDisplacement,MPI_DOUBLE,MPI_COMM_WORLD);
        //Update all ranks grid arrays with new fully updates values

        //Second part of decomposition, same as previously but for alternating positions.
        for(int y = b*(world_rank*1.0/world_size); y<b*((world_rank+1.0)/world_size); y++)
        {
            for(int x = 1; x<a+1; x++)
            {
                if (checkerPattern[(y+1)*(a+2)+x] == 0)
                {
                    potential = 0.25*(potentialGrid[(y+1)*(a+2)+(x-1)]+potentialGrid[(y+1)*(a+2)+(x+1)]+potentialGrid[(y+1-1)*(a+2)+(x)]+potentialGrid[(y+1+1)*(a+2)+(x)]+(pow(h,2)*density[(y+1)*(a+2)+(x)])/epsillon);
                    old = potentialGrid[(y+1)*(a+2)+x];
                    potentialGrid[(y+1)*(a+2)+x] = potentialGrid[(y+1)*(a+2)+x]+w*(potential - potentialGrid[(y+1)*(a+2)+x]);
                    if (abs(old - potentialGrid[(y+1)*(a+2)+x]) < errorTolerance)
                    {
                        converged += 1;
                    }
                    else
                    {
                        converged = 0;
                    }
                }
            }
        }
        //reduce convergence parameters to the minimal one and check if converged or not.
        MPI_Allgatherv(&potentialGrid[(blocksize)*world_rank],blocksize,MPI_DOUBLE,&potentialGrid[0],blockLength,blockDisplacement,MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Allreduce(&converged,&converged,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

        if (converged > ((b+2)*(a+2)))
        {
            conclude = true;
        }
        if (world_rank == 0)
        {
        //cout << potentialGrid[(b-yDesired+1)*(a+2)+xDesired+1] << endl;
        }
    }


    //print matrices to excel document for inerpretartion/analysis.
    if (world_rank == 0)
    {
        for(int y = 1; y<(b+1); y++)
        {
            for(int x = 1; x<(a+1); x++)
            {
                testingdoc << potentialGrid[y*(a+2)+x] << ",";
            }
            testingdoc << "\n";
        }
        testingdoc.close();
    }
    timeFinish = MPI_Wtime();
    if (world_rank == 0)
    {
        cout << endl << "THE VALUE AT (" << inputxDesired/fineness << "," << inputyDesired/fineness <<") IS: " << potentialGrid[(b-yDesired+1)*(a+2)+xDesired+1] << endl;
        cout << "Time taken = " << timeFinish-timeStart << "s"<< endl;
    }

    MPI_Finalize();
    return 0;
}
