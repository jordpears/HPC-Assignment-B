#include <iostream>
#include <cmath>
#include <fstream> //for testing with excel

//appears to work if you look at the one below it but uncertain.

using namespace std;

int main()
{
    ofstream testingdoc;
    testingdoc.open("testingdoc.csv");
    int fineness = 100; //measure of how coarse the grid is.
    double h = 1.0/fineness;
    double epsillon = 1;
    double w = 1.99;
    double errorTolerance = 1E-7;
    double old;
    int converged = 0;
    double potential;
    //initialise and obtain all parameters.
    int a, b, x1, x2, y1, y2, xdes, ydes;
    double q1, q2, tempx1, tempx2, tempy1, tempy2,tempxdes,tempydes;
//    cout << "Input a: ";
//    cin >> a;
//    cout << "Input b: ";
//    cin >> b;
//    cout << "Input q1: ";
//    cin >> q1;
//    cout << "Input x1: ";
//    cin >> tempx1;
//    cout << "Input y1: ";
//    cin >> tempy1;
//    cout << "Input q2: ";
//    cin >> q2;
//    cout << "Input x2: ";
//    cin >> tempx2;
//    cout << "Input y2: ";
//    cin >> tempy2;

    a = 4;b=9;q1=1;tempx1=0;tempy1=0;q2=-1;tempx2=3;tempy2=8; //testing
    a = 4; b = 6; q1 = 1; tempx1 = 3.1415926535897932; tempy1 = 2.7182818284590452; q2 = -2; tempx2 = 1.6180339887487848; tempy2 = 4.6692016091029906; //the problem

    tempxdes = 3.0;
    tempydes = 3.0;

    a = fineness * a;
    b = fineness * b;
    tempxdes = (fineness * tempxdes)+0.5;
    xdes = (int)tempxdes;
    tempydes = (fineness * tempydes)+0.5;
    ydes = (int)tempydes;
    tempx1 = (fineness * tempx1)+0.5;
    x1 = (int)tempx1;
    tempy1 = (fineness * tempy1)+0.5;
    y1 = (int)tempy1;
    tempx2 = (fineness * tempx2)+0.5;
    x2 = (int)tempx2;
    tempy2 = (fineness * tempy2)+0.5;
    y2 = (int)tempy2;
    double density[(b+2)*(a+2)];
    double grid[(b+2)*(a+2)];

    for(int i = 0; i<((b+2)*(a+2)); i++)
    {
        grid[i] = 0.0;
        density[i] = 0.0;
    }

    density[(y1+1)*(a+2)+(x1+1)] = q1/pow(h,2);  //note goes from 0-4 for a,b=5
    density[(y2+1)*(a+2)+(x2+1)] = q2/pow(h,2);

    while(converged < (a+2)*(b+2))
    {
        for(int y = 1; y<b+1; y++)
        {
            for(int x = 1; x<a+1; x++)
            {
                potential = 0.25*(grid[(y)*(a+2)+(x-1)]+grid[(y)*(a+2)+(x+1)]+grid[(y-1)*(a+2)+(x)]+grid[(y+1)*(a+2)+(x)]+(pow(h,2)*density[(y)*(a+2)+(x)])/epsillon);
                old = grid[(y)*(a+2)+x];
                grid[(y)*(a+2)+x] = grid[(y)*(a+2)+x]+w*(potential - grid[(y)*(a+2)+x]);
                if (old - grid[(y)*(a+2)+x] < errorTolerance)
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

    //print matrices for inerpretartion/analysis.
    for(int y = 1; y<(b+1); y++)
    {
        for(int x = 1; x<(a+1); x++)
        {
            cout << density[y*(a+2)+x] << ",";
            testingdoc << grid[y*(a+2)+x] << ",";
        }
        testingdoc << "\n";
        cout << endl;
    }
    testingdoc.close();

    cout << endl << "THE VALUE AT" << tempxdes << "," << tempydes <<"IS: " << grid[ydes*(a+2)+xdes];

    return 0;
}
