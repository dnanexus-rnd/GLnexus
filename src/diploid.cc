#include <utility>
#include <math.h>
#include <assert.h>
using namespace std;

namespace diploid {

// n_gt = (nA+1) choose 2 = (nA+1)!/2/(nA-1)! = (nA+1)(nA)/2
int genotypes(int n_allele) {
    return (n_allele+1)*n_allele/2;
}

int alleles_gt(int n_allele, int i, int j) {
    /*
      0 1 2
    0 0 1 2 
    1   3 4
    2     5

      0 1 2 3
    0 0 1 2 3
    1   4 5 6
    2     7 8
    3       9

    // t[0,0] = 0
    // t[1,1] = n_allele
    // t[2,2] = n_allele + (n_allele - 1)
    // t[3,3] = n_allele + (n_allele - 1) + (n_allele - 2)
    // t[i,i] = i*n_allele-(i-1)(i)/2
    */

    assert(i >= 0 && i < n_allele);
    assert(j >= 0 && j < n_allele);

    if (j < i) {
        swap(i, j);
    }
    return (i == 0 ? 0 : (n_allele*i-(i-1)*i/2)) + (j-i);
}

// given a genotype index, return a pair with the indices of the constituent n_allele
pair<int,int> gt_alleles(int n_allele, int gt) {
    /*
    0 0
    0 1
    0 2
    1 1
    1 2
    2 2

      0 1 2
    0 0 1 2
    1   3 4
    2     5
    */

    /*
    http://stackoverflow.com/a/243342
    WTF? really!? 
    TODO: make a lookup table for small values of n_allele
    */

    assert(gt >= 0 && gt < genotypes(n_allele));

    double m = n_allele;
    double row = (-2*m - 1 + sqrt( (4*m*(m+1) - 8*(double)gt - 7) )) / -2;
    if( row == (double)(int) row ) row -= 1;
    int irow = (int) row;
    
    int col = gt - n_allele * irow + irow*(irow+1) / 2;

    return pair<int,int>(irow, col);
}

namespace trio {
// Given genotype indices of two parents and one child, return:
//   0 if one of the child alleles is present in one parent and the other child
//   allele is present in the other parent
//   1 if one child allele is found in a parent, but the other child allele is
//   not found in the other parent
//   2 if neither child allele is found in the parents
// Note, this is not simply number of alleles not found in either parent, e.g.
// p1=0/0, p2=1/1, ch=1/1 => 1
int mendelian_inconsistencies(int n_allele, int gt_p1, int gt_p2, int gt_ch) {
    auto a_p1 = gt_alleles(n_allele, gt_p1);
    auto a_p2 = gt_alleles(n_allele, gt_p2);
    auto a_ch = gt_alleles(n_allele, gt_ch);
    auto a1 = a_ch.first, a2 = a_ch.second;

    if (a1 == a_p1.first || a1 == a_p1.second) {
        return (a2 == a_p2.first || a2 == a_p2.second) ? 0 : 1;
    } else if (a1 == a_p2.first || a1 == a_p2.second) {
        return (a2 == a_p1.first || a2 == a_p1.second) ? 0 : 1;
    } else {
        return (a2 == a_p1.first || a2 == a_p1.second || a2 == a_p2.first || a2 == a_p2.second) ? 1 : 2;
    }
}
} // namespace trio

} // namespace diploid
