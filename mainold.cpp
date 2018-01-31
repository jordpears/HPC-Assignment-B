#include <iostream>
#include <cmath>
#include <fstream> //for testing with excel
#include <mpi.h>
#include <stdio.h>

//appears to work if you look at the one below it but uncertain.

using namespace std;


int main(int argc, char** argv)
{

    ofstream testingdoc;
    testingdoc.open("testingdoc.csv");
    int fineness = 120; //measure of how coarse the grid is.

    double h = 1.0/fineness;
    double epsillon = 1;
    double w = 1.0;
    double errorTolerance = 1E-7;
    double old;
    int converged = 0;
    //initialise and obtain all parameters.
    int a, b, x1, x2, y1, y2, xdes, ydes;
    double q1, q2, tempx1, tempx2, tempy1, tempy2,tempxdes,tempydes;

    a = 4;
    b = 6;
    q1 = 1;
    tempx1 = 3.1415926535897932;
    tempy1 = 2.7182818284590452;
    q2 = -2;
    tempx2 = 1.6180339887487848;
    tempy2 = 4.6692016091029906; //the problem

    tempxdes = 3.0;
    tempydes = 3.0;

    a = fineness * a;
    b = fineness * b;
    tempxdes = (fineness * tempxdes);
    xdes = (int)tempxdes;
    tempydes = (fineness * tempydes);
    ydes = (int)tempydes;
    tempx1 = (fineness * tempx1);
    x1 = (int)tempx1;
    tempy1 = (fineness * tempy1);
    y1 = (int)tempy1;
    tempx2 = (fineness * tempx2);
    x2 = (int)tempx2;
    tempy2 = (fineness * tempy2);
    y2 = (int)tempy2;


    MPI_Init(NULL,NULL);

    double potential;
    double* density = new double[(b+2)*(a+2)];
    double* grid = new double[(b+2)*(a+2)];
    int* checkerPattern = new int[(b+2)*(a+2)];

    for(int i = 0; i<((b+2)*(a+2)); i++)
    {
        grid[i] = 0.0;
        density[i] = 0.0;
        checkerPattern[i] = 0;
    }

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

    int world_size,world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if ((b+2)*fineness%world_size != 0 && world_rank == 0)
    {
        cout << "Uneven distribution of processors, " << world_size <<" for number of grid points, " << (a+2)*(b+2)*fineness << ", change number of processors to a factor of this number." << endl;
        MPI_Finalize();
        return 0;
    }

    density[(b-y1)*(a+2)+(x1+1)] = q1/pow(h,2);  //note goes from 0-4 for a,b=5
    density[(b-y2)*(a+2)+(x2+1)] = q2/pow(h,2);

    int blocksize = ((a+2)*(b+2))/world_size;
    int* blen = new int[world_size];
    int* displ = new int[world_size];

    for (int i = 0; i < world_size; i++)
    {
        blen[i] = blocksize;
        displ[i] = blocksize*i;
    }

    bool conclude = false;


    while(conclude == false)
    {

        //black part#

        /*

        if (world_rank == 0)
        {
        cout << "=================beg" << endl;
        {
            for(int y = 1; y<(b+1); y++)
            {
                for(int x = 1; x<(a+1); x++)
                {
                    cout << grid[y*(a+2)+x] << ",";
                }
                cout << endl;
            }
        }
        cout << "==========" << endl;
        }
        */


        for(int y = b*(world_rank*1.0/world_size); y<b*((world_rank+1.0)/world_size); y++)//for(int y = 1; y<(b+1)/world_size; y++)
        {
            for(int x = 1; x<a+1; x++)
            {
                if (checkerPattern[(y+1)*(a+2)+x] == 1)
                {
                    potential = 0.25*(grid[(y+1)*(a+2)+(x-1)]+grid[(y+1)*(a+2)+(x+1)]+grid[(y+1-1)*(a+2)+(x)]+grid[(y+1+1)*(a+2)+(x)]+(pow(h,2)*density[(y+1)*(a+2)+(x)])/epsillon);
                    old = grid[(y+1)*(a+2)+x];
                    grid[(y+1)*(a+2)+x] = grid[(y+1)*(a+2)+x]+w*(potential - grid[(y+1)*(a+2)+x]);
                    if (abs(old - grid[(y+1)*(a+2)+x]) < errorTolerance && old != 0)
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
        MPI_Allgatherv(&grid[(blocksize)*world_rank],blocksize,MPI_DOUBLE,&grid[0],blen,displ,MPI_DOUBLE,MPI_COMM_WORLD);

        //red part

        for(int y = b*(world_rank*1.0/world_size); y<b*((world_rank+1.0)/world_size); y++)//for(int y = 1; y<(b+1)/world_size; y++)
        {
            for(int x = 1; x<a+1; x++)
            {
                if (checkerPattern[(y+1)*(a+2)+x] == 0)
                {
                    potential = 0.25*(grid[(y+1)*(a+2)+(x-1)]+grid[(y+1)*(a+2)+(x+1)]+grid[(y+1-1)*(a+2)+(x)]+grid[(y+1+1)*(a+2)+(x)]+(pow(h,2)*density[(y+1)*(a+2)+(x)])/epsillon);
                    old = grid[(y+1)*(a+2)+x];
                    grid[(y+1)*(a+2)+x] = grid[(y+1)*(a+2)+x]+w*(potential - grid[(y+1)*(a+2)+x]);
                    if (abs(old - grid[(y+1)*(a+2)+x]) < errorTolerance && old != 0)
                    {
                        converged += 1;
                    }
                    else
                    {
                        converged = 0;
                    }
                }
            }


        MPI_Allgatherv(&grid[(blocksize)*world_rank],blocksize,MPI_DOUBLE,&grid[0],blen,displ,MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Allreduce(&converged,&converged,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

        if (converged > ((b+2)*(a+2)))
        {
            conclude = true;
        }
    }


    /*
    // serial algorithm
    if (world_rank == 0)
    {
        while(conclude == false)
        {
            for(int y = 1; y<b+1; y++)
            {
                for(int x = 1; x<a+1; x++)
                {
                    potential = 0.25*(grid[(y)*(a+2)+(x-1)]+grid[(y)*(a+2)+(x+1)]+grid[(y-1)*(a+2)+(x)]+grid[(y+1)*(a+2)+(x)]+(pow(h,2)*density[(y)*(a+2)+(x)])/epsillon);
                    old = grid[(y)*(a+2)+x];
                    grid[(y)*(a+2)+x] = grid[(y)*(a+2)+x]+w*(potential - grid[(y)*(a+2)+x]);
                    if (abs(old - grid[(y)*(a+2)+x]) < errorTolerance)
                    {
                        converged += 1;
                    }
                    else
                    {
                        converged = 0;
                    }
                }
            }
            if (converged > ((b+2)*(a+2)))
            {
                conclude = true;
            }
        }
    }
    */

//print matrices for inerpretartion/analysis.

    if (world_rank == 0)
    {
        for(int y = 1; y<(b+1); y++)
        {
            for(int x = 1; x<(a+1); x++)
            {
                //cout << grid[y*(a+2)+x] << ",";
                testingdoc << grid[y*(a+2)+x] << ",";
            }
            testingdoc << "\n";
            //cout << endl;
        }
        testingdoc.close();
    }


    //cout << world_rank << endl;
    if (world_rank == 0)
    {
        cout << endl << "THE VALUE AT (" << tempxdes/fineness << "," << tempydes/fineness <<") IS: " << grid[ydes*(a+2)+xdes] << endl;
    }

    MPI_Finalize();


    return 0;
}
