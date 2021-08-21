#include <Rcpp.h>
using namespace Rcpp;

int getSum(std::vector<int>& fenwick, int index);
void incrementFenwick(std::vector<int>& fenwick, int n, int index, int val);

// [[Rcpp::export]]
double TauSVIV(IntegerVector sv, IntegerVector iv) {
  
  //sv: a representation of the sorted real vector X. 
  //sv is a vector of length X, where the ith value in sv contains the index of the ith largest value in X.
  //Ties are represented by having the tied indices in sv being negative.
  //n missing values (R NaN) results in a string of n 0s at the end of sv.
  //Done by Tau.R, SV()
  
  //iv: an 'integer' representation of real vector Y: 
  //the first value in iv is the size of Y; afterwards, each value in iv corresponds to a value in Y;
  //but that the lowest value in Y becomes '1' in IV; the 2nd lowest value in Y becomes '2' in IV; and so on.
  //Missing values (NaN) are represented by '0'.
  //This allows this function to know for sure the ranking of each value in Y within Y, which is what matters in this
  //non-parametric correlation test.
  //Done by Tau.R, IV()
  
  //See Christensen (2005)
  //c: concordant pairs, d: discordant pairs
  
  long long ex = 0, ey = 0, c = 0, d = 0, n = 0;
  
  //The first value of iv signals the (1-based) index of the largest value in iv;
  //Hence, maxy is the largest value in iv.
  int maxy = iv[iv[0] - 1];
 
  //An array to keep track of how many times a value in IV has been encountered,
  //to handle duplicate Y values.
  std::vector<int> yFreq (maxy, 0);
  
  //Fenwick tree
  std::vector<int> fw (maxy + 1, 0);
  
  int nx = sv.size();
  
  //ci, di: concordant and discordant pairs found on the ith comparison
  //fy: number of times the Y[i] value has been encountered
  //rxi: the rank of X[i] in X
  int ci = 0, fy = 0, di = 0, rxi = 0, yi = 0;
  
  int i = 0;
  
  //0 values in sv represent an NaN in vector X
  while (i < nx && sv[i] != 0) {
  
    rxi = sv[i];
    if (rxi > 0) {
      //Not a tie
      yi = iv[rxi];
      if (yi != 0) {
        
        ci = getSum(fw, yi - 1); fy = yFreq[yi - 1];
        di = n - ci - fy;
        n++; c+= ci; d += di; ey += fy;
        yFreq[yi - 1] = fy + 1;
        incrementFenwick(fw, maxy + 1, yi, 1);
        
      }
      i++;
    }
    
    else {
      
      yi = iv[-rxi];
      std::vector<int> localfreq (maxy, 0);
      std::vector<int> locallist (maxy, 0);
      int nlist = 0, nstreak = 0;
      bool streak = true;
      
      if (yi != 0) {
        locallist[0] = yi;
        localfreq[yi - 1] = 1;
        nlist++; nstreak++;
        ci = getSum(fw, yi - 1); fy = yFreq[yi - 1]; di = n - ci - fy;
        c+= ci; d += di; ey += fy;
      }
      i++;
      
      while (((i < nx && sv[i] != 0) && streak)) {
        rxi = sv[i];
        if (rxi < 0) {
          streak = false;
          yi = iv[-rxi];
        }
        else{
          yi = iv[rxi];
        }
        
        if (yi != 0) {
          int localfy = localfreq[yi - 1];
          ex += nstreak - localfy;
          if (localfy == 0) {
            locallist[nlist] = yi; nlist++;
          }
          localfreq[yi - 1] = localfy + 1;
          
          ci = getSum(fw, yi - 1); fy = yFreq[yi - 1]; di = n - ci - fy;
          c+= ci; d += di; ey += fy;
          nstreak++;
        }
        
        i++;
      }
      
      n += nstreak;
      int j = 0, nlistj = 0;
      while (j < nlist) {
        nlistj = locallist[j];
        incrementFenwick(fw, maxy + 1, nlistj, localfreq[nlistj - 1]);
        yFreq[nlistj - 1] += localfreq[nlistj - 1];
        j++;
      }
      
    }
    
    
  }
  
  //Christensen (2005)'s tau-b formula
  return ((c - d)/(double) (sqrt((c + d + ex) * (c + d + ey))));
  
}

void incrementFenwick(std::vector<int>& fenwick, int n, int index, int val) {
  //Updates the cumulative sum array (as stored in the fenwick structure), corresponding to an increase in 'val' at 'index'
  //Array size = n
  while (index < n) {
    fenwick[index] += val;
    index += (index & (-index));
  }
}

int getSum(std::vector<int>& fenwick, int index) {
  //Gets the cumulative sum of an array as represented in a fenwick structure, up to index
  int sum = 0;
  
  while (index > 0) {
    sum += fenwick[index];
    index -= (index & (-index));
  }
  
  return sum;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/

