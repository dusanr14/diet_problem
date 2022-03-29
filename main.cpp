#define SC_INCLUDE_FX
#include <iostream>
#include <fstream>
#include <systemc>
#include <deque>
#include <vector>
#include <cmath>


#define NUMOFVAR 50
#define NUMOFSLACK 50
#define ROWSIZE (NUMOFSLACK+1)
#define COLSIZE (NUMOFSLACK+NUMOFVAR+1)

#define WIDTH 36
#define FIXED_POINT 16
using namespace std;

typedef sc_dt::sc_fixed_fast<WIDTH,FIXED_POINT> num_t;
typedef sc_dt::sc_fix_fast num_t_matrix[ROWSIZE][COLSIZE];


//typedef std::vector<num_t> array_temp;
//typedef std::vector<double> orig_array_t;


//for the matrix
//
//typedef std::vector< vector<num_t> > mat_t;
//typedef std::vector< vector<float> > orig_mat_t;

num_t wv_fixed[ROWSIZE][COLSIZE];
//od_profesora//num_t* wv_fixed;
num_t pivot_fixed;

//mat_t fixed_matrix;
//second parameter changer=d, d = matrix member (not array), i push back in matrix whole array
void copy2fix(num_t wv_fixed[ROWSIZE][COLSIZE],const float wv[ROWSIZE][COLSIZE], int W, int F)
{
	//cout<<"Pocinje copy2fix"<<endl;
	for (int j=0; j<ROWSIZE; j++)
	{
		for(int i=0; i<COLSIZE; i++)
		{
		num_t d;
		d = wv[j][i];
		//cout<<wv[j][i]<<" "<<d<<" ";
		if (d.overflow_flag())
			std::cout << "Overflow in conversion.\n";
		
		//PROBLEM!!!
		wv_fixed[j][i] = d;
		//od_profesora//wv_fixed[i*ROWSIZE+j] = d;
		//cout<<wv_fixed[j][i]<<endl;
		//cout<<wv_fixed[j][i]<<" "<<wv[j][i]<< endl;
		}

	}
	//while(1);
}

bool passCheck(const num_t wv_fixed[ROWSIZE][COLSIZE], const float wv[ROWSIZE][COLSIZE],
			    double delta)
{
	
	
	for (size_t j = 0; j != ROWSIZE; ++j)
	{
		for (size_t i = 0; i != COLSIZE; ++i)
		{
		if (std::abs(wv_fixed[j][i] - wv[j][i]) > delta)
			return false;
		}
	}
	return true;
}

bool checkOptimality(float wv[ROWSIZE][COLSIZE])
{
    for(int i=0;i<COLSIZE-1;i++)
    {
        if(wv[ROWSIZE-1][i]<0)//min> max<
            return false;
    }
    return true;
}
bool isUnbounded(float wv[ROWSIZE][COLSIZE],int pivotCol)
{
    for(int j=0;j<ROWSIZE-1;j++)
    {
        if(wv[j][pivotCol]>0)
            return false;
    }
    return true;
}
void print(float wv[ROWSIZE][COLSIZE])
{
    for(int j=0;j<ROWSIZE;j++)
        {
            for(int i=0;i<COLSIZE;i++)
            {
                cout<<wv[j][i]<<" ";
            }
            cout<<endl;
        }
        cout<<endl<<endl<<endl;
}
void makeMatrix(float wv[ROWSIZE][COLSIZE])
{

	fstream myFile;
    myFile.open("baza.txt",ios::in); //otvaram fajl u read modu
	if(myFile.is_open())
    {
        for(int j = 0; j < ROWSIZE; j++)
        {
            for(int i = 0; i< NUMOFVAR; i++)
            {
              myFile >> wv[j][i];
            }
        }
		for(int j = 0;j< NUMOFSLACK;j++)
		{
			myFile >> wv[j][COLSIZE-1];
		}
    }
    myFile.close();

}
int findPivotCol(float wv[ROWSIZE][COLSIZE])
{
     float minnegval=wv[ROWSIZE-1][0];
       int loc=0;
        for(int i=1;i<COLSIZE-1;i++)
        {
            if(wv[ROWSIZE-1][i]<minnegval)
            {
                minnegval=wv[ROWSIZE-1][i];
                loc=i;
            }
        }
        return loc;
}

int findPivotRow(float wv[ROWSIZE][COLSIZE],int pivotCol)
{
    float rat[ROWSIZE-1];
    for(int j=0;j<ROWSIZE-1;j++)
        {
            if(wv[j][pivotCol]>0)
            {
                rat[j]=wv[j][COLSIZE-1]/wv[j][pivotCol];
            }
            else
            {
                rat[j]=0;
            }
        }

        float minpozval=99999999;
        int loc=0;
        for(int j=0;j<ROWSIZE-1;j++)
        {
            if(rat[j]>0)
            {
                if(rat[j]<minpozval)
                {
                    minpozval=rat[j];
                    loc=j;
                }
            }
        }
        return loc;
}
void doPivoting(num_t wv[ROWSIZE][COLSIZE],int pivotRow,int pivotCol,num_t pivot)
{
    num_t newRow[COLSIZE];
    num_t pivotColVal[ROWSIZE];
    num_t nr;
    num_t pcv;
    for(int i=0;i<COLSIZE;i++)
        {
            nr = wv[pivotRow][i]/pivot;
            newRow[i]= nr;
        }

        for(int j=0;j<ROWSIZE;j++)
        {
            pcv = wv[j][pivotCol];
            pivotColVal[j]=pcv;
        }

        for(int j=0;j<ROWSIZE;j++)
        {
            if(j==pivotRow)
            {
                for(int i=0;i<COLSIZE;i++)
                {
                    wv[j][i]=wv[j][i]/pivot;
                }
            }
            else
            {
                for(int i=0;i<COLSIZE;i++)
                {
                    wv[j][i]=wv[j][i]-newRow[i]*pivotColVal[j];
                }
            }
        }
}
void solutions(float wv[ROWSIZE][COLSIZE])
{
    for(int i=0;i<NUMOFVAR; i++)  //every basic column has the values, get it form B array
     {
        int count0 = 0;
        int index = 0;
        for(int j=0; j<ROWSIZE-1; j++)
        {
            if(wv[j][i]==0.0)
            {
                count0 = count0+1;
            }
            else if(wv[j][i]==1)
            {
                index = j;
            }


        }

        if(count0 == ROWSIZE - 2 )
        {
            cout<<"variable"<<i+1<<": "<<wv[index][COLSIZE-1]<<endl;  //every basic column has the values, get it form B array
        }
        else
        {
            cout<<"variable"<<i+1<<": "<<0<<endl;
        }
    }

    cout<<""<<endl;
    cout<<endl<<"Optimal solution is "<<wv[ROWSIZE-1][COLSIZE-1]<<endl;
}
void simplexCalculate(float wv[ROWSIZE][COLSIZE])
{

    //float minnegval;
    //float minpozval;
    //int loc;
    int pivotRow;
    int pivotCol;
    bool unbounded=false;
    float pivot;

    //float solVar[NUMOFVAR];

    while(!checkOptimality(wv))
    {
        pivotCol=findPivotCol(wv);

        if(isUnbounded(wv,pivotCol))
        {
            unbounded=true;
            break;
        }


        pivotRow=findPivotRow(wv,pivotCol);

        pivot=wv[pivotRow][pivotCol];
        
	// ulaze u funkciju promenim da su fixed
	
	copy2fix(wv_fixed,wv, WIDTH,FIXED_POINT);
	pivot_fixed = pivot;
        doPivoting(wv_fixed,pivotRow,pivotCol,pivot_fixed);
        
        for(int j=0;j<ROWSIZE; j++)
	{
		for(int i =0;i<COLSIZE;i++)
		{
			wv[j][i]=wv_fixed[j][i];
		}
	}
	
	
    }
    //Ispisivanje rezultata
    if(unbounded)
    {
        cout<<"Unbounded"<<endl;
    }
    else
    {
        //print(wv);

        solutions(wv);

    }
}

int sc_main(int argc, char*argv[])
{
 
    float wv[ROWSIZE][COLSIZE];
	for(int j=0;j<ROWSIZE; j++)
	{
		for(int i =0;i<COLSIZE;i++)
		{
			wv[j][i]=0;
		}
	}

	makeMatrix(wv);
	for(int j=0;j<ROWSIZE-1;j++)
	{
		{
			wv[j][NUMOFVAR+j]=1;
		}
	}


        simplexCalculate(wv);
    return 0;
}

//napravim python skriptu koja kompajlira, menja w i h, i onda vidim kad dobijem zadovoljavajuci rezultat
