#include <iostream>
#include <fstream>
#define NUMOFVAR 3
#define NUMOFSLACK 3
#define ROWSIZE (NUMOFSLACK+1)
#define COLSIZE (NUMOFSLACK+NUMOFVAR+1)

using namespace std;
bool checkOptimality(float wv[ROWSIZE][COLSIZE])
{
    for(int i=0;i<COLSIZE-1;i++)
    {
        if(wv[ROWSIZE-1][i]<0)
            return false;
    }
    return true;
}
bool isUnbounded(float wv[ROWSIZE][COLSIZE],int pivotCol)
{
    for(int j=0;j<ROWSIZE-1;j++)
    {
        if(wv[j][pivotCol]>=0)
            return false;
    }
    return true;
}
int main()
{
    fstream myFile;
    myFile.open("baza.txt",ios::in); //otvaram fajl u read modu
    float wv[ROWSIZE][COLSIZE];

    if(myFile.is_open())
    {


        for(int j = 0; j < ROWSIZE; j++)
        {
            for(int i = 0; i< COLSIZE; i++)
            {
              myFile >> wv[j][i];
            }
        }
    }

    /* ispis matrice, kao neka provera
    for(int j=0;j<ROWSIZE;j++)
        {
            for(int i=0;i<COLSIZE;i++)
            {
                cout<<wv[j][i]<<" ";
            }
            cout<<endl;
        }
        cout<<endl<<endl<<endl;
    */

    float rat[ROWSIZE-1];
    float minnegval;
    float minpozval;
    int loc;
    int pivotRow;
    int pivotCol;
    bool unbounded=false;
    float pivot;
    float newRow[COLSIZE];
    float pivotColVal[ROWSIZE];
    //float solVar[NUMOFVAR];

    while(!checkOptimality(wv))
    {
        minnegval=wv[ROWSIZE-1][0];
        loc=0;
        for(int i=1;i<COLSIZE-1;i++)
        {
            if(wv[ROWSIZE-1][i]<minnegval)
            {
                minnegval=wv[ROWSIZE-1][i];
                loc=i;
            }
        }
        pivotCol=loc;

        if(isUnbounded(wv,pivotCol))
        {
            unbounded=true;
            break;
        }

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

        minpozval=rat[0];
        loc=0;
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
        pivotRow=loc;

        pivot=wv[pivotRow][pivotCol];

        for(int i=0;i<COLSIZE;i++)
        {
            newRow[i]=wv[pivotRow][i]/pivot;
        }

        for(int j=0;j<ROWSIZE;j++)
        {
            pivotColVal[j]=wv[j][pivotCol];
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
    //Ispisivanje rezultata
    if(unbounded)
    {
        cout<<"Unbounded"<<endl;
    }
    else
    {
        for(int j=0;j<ROWSIZE;j++)
        {
            for(int i=0;i<COLSIZE;i++)
            {
                cout<<wv[j][i]<<" ";
            }
            cout<<endl;
        }

        cout<<endl<<"Optimal solution is "<<wv[ROWSIZE-1][COLSIZE-1]<<endl;
    }

    return 0;
}
