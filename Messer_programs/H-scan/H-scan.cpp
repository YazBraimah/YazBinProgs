// g++ -O3 ./H-scan.cpp -o H-scan

#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;


class snp_data
{
public:
  int n; // number of samples
  int l; // number of snps

  vector<string> seq; // genotypes
  vector<int>    pos; // snp positions

  snp_data(ifstream& infile, int ms)
  {
    if(ms==0) // regular input file
      {
	// count sample size and initialize window

	string line;
	getline(infile,line); n = count(line.begin(),line.end(), ','); infile.seekg(0,ios::beg);
	seq.resize(n);
 
	// read datafile line by line

	while(getline(infile,line) && infile.good())
	  {
	    string token; 
	    istringstream ss(line); 

	    getline(ss,token,','); pos.push_back(atoi(token.c_str()));

	    for(int i=0; i<n; i++)
	      { 
		getline(ss,token,','); 
		seq[i].append(token);
	      }
	  }
	l = seq[0].length();
      }	

    else // MS input file
      {
	string line,sub;

	while(line.find("positions") == string::npos && !infile.eof()) { getline(infile,line); } 

	istringstream iss(line); iss >> sub;

	while(iss >> sub) 
	  {
	    int p = int(1e6*atof(sub.c_str()));
	    pos.push_back(p);
	  }

	while(getline(infile,line))
	  {
	    if(line.find("0") != string::npos || line.find("1") != string::npos)
	      {
		seq.push_back(line);
	      }
	  }
	n = seq.size();
	l = seq[0].length();
      }
  }
};
  

vector<double> calc_H(snp_data &S, int x, int g, int pair, int dist)
{
  vector<double> H; 
  H.push_back(0.0);

  // cycle over all pairs

  for(int i=0; i<S.n; i++)
    {
      for(int j=i+1; j<S.n; j++)
	{
	  double h = 0.0;
	  
	  int x_r = x; // extend from x to the right
	  
	  while((x_r+1<S.l) && 
		(S.pos[x_r+1]-S.pos[x_r]<g) &&
		(S.seq[i].at(x_r+1)!='.') &&
		(S.seq[j].at(x_r+1)!='.') &&
		((S.seq[i].at(x_r+1)==S.seq[j].at(x_r+1)) ||
		 (S.seq[i].at(x_r+1)=='N') ||
		 (S.seq[j].at(x_r+1)=='N') ) ) { x_r++; }
	  
	  int x_l = x; // extend from x to the left
	  
	  while((x_l-1>=0) &&
		(S.seq[i].at(x_l-1)==S.seq[j].at(x_l-1)) && 
		(S.pos[x_l]-S.pos[x_l-1]<g) &&
		(S.seq[i].at(x_l-1)!='.') &&
		(S.seq[j].at(x_l-1)!='.') &&
		((S.seq[i].at(x_l-1)==S.seq[j].at(x_l-1)) ||
		 (S.seq[i].at(x_l-1)=='N') ||
		 (S.seq[j].at(x_l-1)=='N') ) ) { x_l--; }
	  

	  if(x_r!=x_l)
	    {
	      switch(dist) // tract length distance estimation
		{
		case 0: // physical distance in bp (default) 
		  h = (double)(S.pos[x_r]-S.pos[x_l]-1);
		  break;
		case 1: // number of segregating sites in sample
		  h = (double)(x_r-x_l-1);
		  break;
		case 2: // product of the above
		  h = (double)(S.pos[x_r]-S.pos[x_l]-1)*(double)(x_r-x_l-1);
		  break;	
		default: exit(1);			
		}
	    }
	  
	  H[0] += h; 
	  if(pair==1) { H.push_back(h); }
	}
    }
  H[0] = 2.0*H[0]/((double)S.n*(double)(S.n-1));
  return H;
};


int main(int argc, char **argv)
{
  ifstream snp_stream; // SNP data file stream

  int g = 1e9;      // max gapsize (default: 1e9 bp)
  int s = 1;        // moving window shift (default: 1 SNP)
  int left  = 0;    // left end of analysis range (default: 0)
  int right = 1e9;  // right end of analysis range (default: 1e9)
  int ms = 0;       // ms input mode [1: ms mode, 0: regular mode (default)]
  int pair = 0;     // pairwise distances [0: do not report pairwise distances (default), 1: report pairwise distances]
  int dist = 0;     // distance calculation: [0: physical distance in bp (default), 1: number of segregating sites in sample, 2: product of both]
  
  int opt;

  while ((opt = getopt(argc, argv, "i:s:g:l:r:d:pm")) != -1) 
    {
      switch (opt) 
	{
	case 'i':
	  snp_stream.open(optarg, std::ifstream::in);
	  break;
	case 's':
	  s = atoi(optarg);
	  break;
	case 'g':
	  g = atof(optarg);
	  break;
	case 'l':
	  left = atoi(optarg);
	  break;
	case 'r':
	  right = atoi(optarg);
	  break;
	case 'd':
	  dist = atoi(optarg);
	  break;			
	case 'p':
	  pair = 1;
	  break;
	case 'm':
	  ms = 1;
	  break;
	default: exit(1);
      	}
    }

  if (!snp_stream.is_open()) { cerr << "invalid SNP file" << endl ; exit(1); }

  snp_data S(snp_stream,ms);
  
  if(pair==0) { cout << "x\tH"; cout << endl; }
  
  for(int x=0; x<S.l; x+=s)
    {
      if(S.pos[x]>=left && S.pos[x]<=right)
	{
	  vector<double> H = calc_H(S,x,g,pair,dist);

	  cout << S.pos[x];
	  for(int i=0; i<H.size(); i++) { cout << "\t"<< H[i]; } 
	  cout << endl;
	}
    }
  return 0;
}

