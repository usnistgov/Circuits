
// BP2022.cc
// Rene Peralta
// June 3, 2022
// This is a faster implementation of the BP algorithm.
// The function bool reachable(unsigned long long T, int K, int S) which 
// returns true if T is the sum of K elements among Base[S..BaseSize-1],
// is implemented via 
//   i) splitting the bases into two groups A_Base and B_Base
//  ii) computing all i element sums from L1;
// iii) for each such sum S, see if T+S is reachable as a sum of
//      K-i elements in L2;

/* The idea of the heuristic is to
  - keep a "base" of available signals (initially the
    base is just the set of variables x0, ..., xn);
  - for each required matrix output (I called them
    "Targets") keep a "distance" from the base to the
    output. e.g. Dist[3]+1 is the smallest number of base
    elements that I need to sum in order to obtain the
    third row of the matrix;
  - greedily pick a new basis element by adding two existing
    basis elements;

  The current criteria for picking the new basis element is
   - if a target is the sum of two basis elements, pick those
   - otherwise, pick the one that minimizes the sum of new distances
   - resolve ties by maximizing the euclidean norm of the 
   vector of new distances; 
*/

#include <math.h>
#include <ctype.h>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <unordered_set>
#include <set>

using namespace std;

const int DistanceBound=8;
const int MaxBaseSize=200;

string VarName[MaxBaseSize];
int NumInputs;
int NumTargets;
int ProgramSize;
unsigned long long Target[MaxBaseSize];
int Dist[MaxBaseSize]; //distance from current base to Target[i]
bool M[MaxBaseSize][MaxBaseSize];
string TargetName[MaxBaseSize];
unsigned long long Base[MaxBaseSize];
int BaseSize;
int TargetsFound;
int MaxDistance;

//array of indices into Base[MaxBaseSize];
int IndexArray[MaxBaseSize]; //an array of indices into Base[MaxBaseSize];

set <unsigned long long> A_BaseSums[DistanceBound];
unordered_set <unsigned long long> B_BaseSums[DistanceBound];

void InitBase();
void ReadTargetMatrix(); //reads matrix M[][], inits Dist[], inits Target[],

bool is_target(unsigned long long x);
int which_target(unsigned long long x);
bool is_base(unsigned long long x);
int Distance(int u); //calculates the distance from the base to Target[u]
int TotalDistance(); //returns the sum of distances to targets
void PickNewBaseElement();
void SetLowSums();
void SetSums(int IndArray[], int offset, int IndArraySize, set <unsigned long long> Sums[]);
void SetSums(int IndArray[], int offset, int IndArraySize, unordered_set <unsigned long long> Sums[]);


ifstream TheMatrix;
ofstream out_file;


int
main(int argc, char *argv[])
{
 assert(argc==2);
          
  ProgramSize = 0;
  TheMatrix.open(argv[1]); 
  ReadTargetMatrix();
  InitBase(); 
  assert(BaseSize == NumInputs);
  cout << "MaxDistance is " << MaxDistance << endl;
  assert(MaxDistance < DistanceBound);
  
  SetSums(IndexArray,0, BaseSize/2,A_BaseSums);
  SetSums(IndexArray,(BaseSize/2) , BaseSize - BaseSize/2,B_BaseSums);
  
  while (TargetsFound < NumTargets) 
  {
    PickNewBaseElement();
  }
}//main

void InitBase()
{
  TargetsFound = 0;
  Base[0] = 1;
  for (int i = 1; i < NumInputs; i++) Base[i] = 2*Base[i-1];
  BaseSize = NumInputs; //initial base is just the xi's
  for (int i = 0; i < NumTargets; i++) 
    if (Dist[i] == 0) 
    {
      TargetsFound++;
      cout << "output " << TargetName[i] << " = " ;
      int safecount = 0;
      for (int u = 0; u < NumInputs; u++) 
      {
        if (M[i][u] == 1) 
        {
          cout << VarName[u] ;
          safecount++;
          assert(safecount < 2);
        }
      }
      cout << endl;
    }
}

int TotalDistance() //returns the sum of distances to targets
{
  int D = 0;
  int t;
  for (int i = 0; i < NumTargets; i++) 
  {
    t = Distance(i);
    D = D + t;
  }
  return D;
}

unsigned long long NewBase; //global variable containing a candidate new base

void PickNewBaseElement()
{
  int MinDistance;
  unsigned long long TheBest;
  int ThisDist;
  int ThisNorm, OldNorm;
  int besti,bestj, d;
  bool easytarget;

  MinDistance = BaseSize*NumTargets; //i.e. something big
  OldNorm = 0; //i.e. something small
  //try all pairs of bases
  for (int i = 0; i < BaseSize - 1; i++)
  {
    for (int j = i+1; j < BaseSize; j++)
    {
      NewBase = Base[i] ^ Base[j];
      //if NewBase is not new continue
      if (is_base(NewBase)) continue;
      //if NewBase is target then choose it
      easytarget = false;
      if (is_target(NewBase))
      {
        easytarget = true;
        besti = i;
        bestj = j;
        TheBest = NewBase;
        break;
      }
      ThisDist = TotalDistance();
      if (ThisDist <= MinDistance)
      {
        //calculate Norm
        ThisNorm = 0;
        for (int b = 0; b < NumTargets; b++)
        {
          d = Distance(b); 
          ThisNorm = ThisNorm + d*d;
        }
        //resolve tie in favor of largest norm
        if ((ThisDist < MinDistance) || (ThisNorm > OldNorm) )
        {
          besti = i;
          bestj = j;
          TheBest = NewBase;
          MinDistance = ThisDist; 
          OldNorm = ThisNorm;
        }
      }
    }
    if (easytarget) break;
  }
  //update Dist array
  NewBase = TheBest; // Distance uses NewBase
  int newD = 0; //new distance
  int newN = 0; //new norm
  for (int i = 0; i < NumTargets; i++) 
  {
    Dist[i] = Distance(i);
    newD = newD + Dist[i];
    newN = newN + Dist[i]*Dist[i];
  }
  Base[BaseSize] = TheBest;
  BaseSize++;
  // update A_BaseSums and B_BaseSums
  SetSums(IndexArray,0, BaseSize/2,A_BaseSums);
  
  SetSums(IndexArray,(BaseSize/2) , BaseSize - BaseSize/2,B_BaseSums);
  //output linear program
  ProgramSize++;
  //cout << ProgramSize << " : " ;
  //cout << "X" << BaseSize-1 << " = X" << besti << " + X" << bestj;
  cout << "X" << BaseSize-1;
  if (besti < NumInputs) cout << " = " << VarName[besti];
  else cout << " = X" << besti;
  if (bestj < NumInputs) cout << " + " << VarName[bestj];
  else cout << " + X" << bestj;
  cout << endl;
  if (is_target(TheBest)) 
    cout << "output " << TargetName[which_target(TheBest)] << " = " << "X" << BaseSize-1 << endl;
  // if a target is found update counter
  if (is_target(TheBest)) TargetsFound++;
} //PickNewBaseElement()

//reads matrix M[][], inits Dist[], inits Target[]
void ReadTargetMatrix() 
{
  //TheMatrix.open("matrix848");
  
  TheMatrix >> NumTargets;
  TheMatrix >> NumInputs;
  cout << NumTargets << " equations and ";
  cout << NumInputs << " variables" << endl;
  assert(NumInputs < 8*sizeof(unsigned long long) );
  for (int i = 0; i < NumInputs; i++) TheMatrix >> VarName[i] ;
  for (int i = 0; i < NumInputs; i++) cout << VarName[i] << " ";
  cout << endl;

  int bit;
  MaxDistance = 0;
  for (int i = 0; i < NumTargets; i++)
  //read row i
  {
    TheMatrix >> TargetName[i];
    unsigned long long PowerOfTwo  = 1;
    Dist[i] = -1; //initial distance from Target[i] is Hamming weight - 1
    for (int j = 0; j < NumInputs; j++) 
    {
      TheMatrix >> bit;
      if (bit) 
      {
        M[i][j] = true;
        Dist[i]++; 
        Target[i] = Target[i] + PowerOfTwo;
      }
      else M[i][j] = false;
      
      PowerOfTwo = PowerOfTwo * 2;
    }
    if (Dist[i] > MaxDistance) MaxDistance = Dist[i];
  }
} //ReadTargetMatrix()

bool is_target(unsigned long long x)
{
  for (int i = 0; i < NumTargets; i++)
    if (x == Target[i]) return true;
  return false;
} //is_target

int which_target(unsigned long long x)
{
  for (int i = 0; i < NumTargets; i++)
    if (x == Target[i]) return i;
  cout << "error can't find target" << endl;
  exit(0);
} //which_target

bool is_base(unsigned long long x)
{
  //sanity check, shouldn't ask if 0 is base
  if (x==0) { cout << "asking if 0 is in Base " <<endl ; exit(0); }
  
  for (int i = 0; i < BaseSize; i++) if (x == Base[i]) return true;
  return false;
} //is_base

// Distance is 1 less than the number of elements
// in the base that I need to add in order to get Target[u].
// The next function calculates the distance from the base,
// augmented by NewBase, to Target[u]. Uses the following observations:
// Adding to the base can only decrease distance. 
// Also, since NewBase is the sum of two old base 
// elements, the distance from the augmented base 
// to Target[u] can decrease at most by 1. If the
// the distance decreases, then NewBase must be one
// of the summands.
  
// returns true if T is the sum of an element in A_BaseSums[i-1]
// and an element in A_BaseSums[K-i-1] for some i in [0,K]
// note that BaseSums[i] contains the sum of i+1 bases.
// need to deal with i = 0 and i = K case separately
bool reachableNew(unsigned long long T, int K)
{
   unsigned long long tempTarget;
   assert(K <= MaxDistance + 1);
   if (K==0) return is_base(T);

   assert(K > 0);
   
   //i = 0 see if T is in B_BaseSums[K-1]
   if (B_BaseSums[K-1].find(T) != B_BaseSums[K-1].end()) return true;

   //i = K see if T is in A_BaseSums[K-1]
   if (A_BaseSums[K-1].find(T) != A_BaseSums[K-1].end()) return true;
   
   for (int a = 1; a < K; a++)
   {
     //traverse A_BaseSums[a-1]
     for (set<unsigned long long>::iterator it= A_BaseSums[a-1].begin(); 
                                         it!=A_BaseSums[a-1].end(); ++it)
     {
       tempTarget = T ^ (*it);
       if (B_BaseSums[K-a-1].find(tempTarget) != B_BaseSums[K-a-1].end()) 
       {
         return true;
       }
     }
   }
   return false;
}

int Distance(int u) 
{
  //if Target[u] is in augmented base return 0;
  if (is_base(Target[u]) || (NewBase == Target[u])) return 0;
  
  // Try all combinations of Dist[u]-1 base elements until one sums 
  // to Target[u] + NewBase. If this is true, then Target[u] is the
  // sum of Dist[u] elements in the augmented base, and therefore
  // the distance decreases by 1.
  
  if (Dist[u] == 0) return 0;
  if (reachableNew(Target[u] ^ NewBase,Dist[u]-1))
  {
    return (Dist[u]-1);
  }
  else 
  {
    return Dist[u]; //keep old distance 
  }
} //Distance(int u) 

bool FirstKtuple;

void initIndexArray(int offset, int K, int theIndexArray[])
{
  //the array will contain K non-zero elements
  FirstKtuple = true;
  for (int i = 0; i < MaxBaseSize; i++)
  {
    if (i < K) theIndexArray[i] = i+offset;
      else theIndexArray[i] = 0;
  }
}//initIndexArray()

// advances IndexArray of K elements in Base
// H bound is the highest base index allowed
bool stepIndexArray(int K, int theArray[], int Hbound)
{
  assert(K <= BaseSize/2);
  if (FirstKtuple) { FirstKtuple = false; return true;}
  //find largest index that can be increased
  int increaseMe;
  increaseMe = K-1;
  for (int i = 0; theArray[K-1-i] > (Hbound -1 - i); i++) increaseMe--; 
  if (increaseMe < 0) return false; 
  theArray[increaseMe]++;

  for (int i = increaseMe+1; i < K; i++) 
    theArray[i] = theArray[i-1] + 1;

  return true;
} //stepIndexArray(int K)

unsigned long long sumBasesInIndexArray( int K, int theArray[])
{
  unsigned long long theSum = 0;
  for (int i = 0; i < K; i++) theSum = theSum^(Base[theArray[i]]);
  return theSum;
}

void SetSums(int IndArray[], int offset, int IndArraySize, set <unsigned long long> Sums[])
{
  bool step;
  unsigned long long sum;

  int MD = 0;
  for (int i = 0; i < NumTargets; i++) if (Dist[i] > MD) MD = Dist[i];

  for (int d = 0; d <= MD; d++)
  {
    if ((d+1)  > IndArraySize ) break;
    initIndexArray(offset,d+1, IndArray);
    while (1)
    {

      step = stepIndexArray(d+1,IndArray,IndArraySize+offset-1);
      if (step) 
      {
        sum = sumBasesInIndexArray(d+1,IndArray) ;
        
        Sums[d].insert(sum);
      }
      else break;
    }
  }
}//SetSums

void SetSums(int IndArray[], int offset, int IndArraySize, 
                             unordered_set <unsigned long long> Sums[])
{
  bool step;
  unsigned long long sum;

  int MD = 0;
  for (int i = 0; i < NumTargets; i++) if (Dist[i] > MD) MD = Dist[i];
  for (int d = 0; d <= MD; d++)
  {
    if ((d+1)  > IndArraySize ) break;
    initIndexArray(offset,d+1, IndArray);
    while (1)
    {

      step = stepIndexArray(d+1,IndArray,IndArraySize+offset-1);
      if (step) 
      {
        sum = sumBasesInIndexArray(d+1,IndArray) ;
        
        Sums[d].insert(sum);
      }
      else break;
    }
  }
}//SetSums

