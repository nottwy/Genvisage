#include <bits/stdc++.h>
#include <map>
#include <iostream>
#include <fstream>
#include <queue>
#include <string>
#include <algorithm>
using namespace std;

typedef double val_type;

struct feature{
    vector<val_type > vals;
    int correct;
    int fid;

    feature(vector<val_type > _vals, int _correct, int _fid):vals(_vals),correct(_correct),fid(_fid){}
    feature(){}
    bool operator<(const feature& a) const  //sort from high to low
    {
        return correct > a.correct;
    }
};

struct topG{
    int gid;
    int correct;
    topG(int _gid=0, int _correct =0):gid(_gid),correct(_correct){}
    bool operator<(const topG& rhs) const // sort from low to high
    {
        return correct < rhs.correct;
    }
};

struct top_1d{
    int fid; // feature id
    int correct_p;
    int correct_n;
    int correct; // number of correctly classified genes
    top_1d(int _fid=0, int _correct=0,int _correct_p=0,int _correct_n=0):fid(_fid),
           correct(_correct), correct_p(_correct_p), correct_n(_correct_n){}
    bool operator<(const top_1d& rhs) const
    {
        return correct < rhs.correct;
    }
};

struct top_2d{
    int f1;
    int f2;
    int miss;
    int miss_p;
    int miss_n;
    top_2d(int _f1=0, int _f2=0, int _miss=0, int _miss_p=0, int _miss_n=0):f1(_f1),f2(_f2),miss(_miss),miss_p(_miss_p),miss_n(_miss_n){}
    bool operator<(const top_2d& a) const
    {
        return miss > a.miss;
    }
};

struct point{
	double x;
	double y;
	point(double _x, double _y):x(_x),y(_y){}
};


int  FEATURE_CONSIDER, MATRIX_GENE_NUM,GENE_NUM, FEATURE_NUM,TOPK;
long FEATURE_CONSIDER_TOTAL;
double DELTA,REST;
bool EARLY_TERMINATE, SORT_G, SORT_F,HORIZONTAL, WEIGHTED, SAMPLING,SAMPLING_OPT,TRASNFORM,PRINT,HIST; // WEIGHTED: pos v.s. neg
string DELIMITER, MATRIXF, EXPF, EXPF2, output_dir;
int printF1, printF2, histSize, buckets;


vector<vector<val_type > > matrix;
vector<int > genes; // indicate the positive or negative of this gene
vector<feature > features; // genes are not sorted: initialization for ff; $features is not used afterwards
vector<vector<val_type > > raw; // matrix[fid][gid]
vector<feature > ff; // genes are sorted or not -> use ff in the algorithm 
vector<vector<val_type > > sample_f;
map<string, int> gene_id;  // gene name -> gid
vector<string > featureName;
std::vector<pair<int, int> > prints;
int pos_num,pos_weight;