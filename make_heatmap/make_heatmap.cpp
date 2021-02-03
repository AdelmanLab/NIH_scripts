//Created by Adam Burkholder, National Institute of Environmental Health Sciences, 2011-2015

#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <sstream>
#include <tr1/unordered_map>
#include <errno.h>
#include <set>
#include <algorithm>
#include <getopt.h>
#include <cstring>
#include <cstdlib>
#include <cmath>

#ifndef SINGLE
#include <pthread.h>
#endif

using namespace std;
using tr1::unordered_map;

struct chr_entry {								//stores features in gene list file, bin start and end locations
	vector<string> id;
	vector<long> bins_start;
	vector<long> bins_end;
	vector<vector<pair<long,long> > > bins;
};

struct ptr_entry {								//stores pointers to single feature in gene list file
	string *id;
	long *bins_start;
	long *bins_end;
	vector<pair<long,long> > *bins;
};

struct totals_info {							//stores counts intersecting each bin for a given feature from the gene list file
	string desc;
	string chr;
	long start;
	long end;
	string strand;
	vector<vector<double> > total;
	vector<vector<long> > count;
	vector<long> length;
};

struct hit_entry {								//stores information about single lines in hit file
	string chr,strand;
	long location;
	long location2;
	double value;
};

class opt_parser {
public:
	static void usage(void) {
		cout << "Usage: make_heatmap [opts] -b c ... [output] [bin start] [size] [count]\n"
				"Usage: make_heatmap [opts] -b v ... [output] [bin count]\n"
				"Usage: make_heatmap [opts] [hits] [genes] [output] [bins]\n"
				"Usage: make_heatmap [opts] -p [plus hits] -m [minus hits] [genes] ...\n"
		        "Available Options:\n"
		        "  --help                     produce this help message\n"
				"  -p [ --plus ] arg          specify file containing plus strand hits\n"
				"  -m [ --minus ] arg         specify file containing minus strand hits\n"
#ifndef SINGLE
				"  -t [ --threads ] arg (=1)  specify number of threads to use\n"
#endif
		        "  -b [ --bintype ] arg (=c)  specify method of defining relative bin locations:\n"
		        "                               c  command line, fixed bin size\n"
		        "                               v  command line, variable bin size\n"
				"                               f  file\n"
		        "  -h [ --hittype ] arg (=G)  specify hit file type:\n"
		        "                               G  standard bedGraph (0-based start coordinate)\n"
				"                               g  bedgraph (1-based start coordinate)\n"
		        "                               b  basic bed\n"
				"                               e  extended bed\n"
		        "                               c  cppmatch\n"
		        "  -l [ --hitloc ] arg (=p)   specify location within hit to match to gene\n"
		        "                             list:\n"
				"                               p  physical start\n"
				"                               d  physical end\n"
				"                               s  genetic start, requires -p, -m, -h c, or -h e\n"
		        "                               e  genetic end, requires -p, -m, -h c, or -h e\n"
		        "                               c  center\n"
		        "  -a [ --anchor ] arg (=s)   specify location within each gene to anchor\n"
		        "                             relative bin locations:\n"
		        "                               s  genetic start, requires stranded gene list\n"
		        "                               e  genetic end, requires stranded gene list\n"
				"                               p  physical start\n"
				"                               d  physical end\n"
		        "                               u  user-specified, contained in second column\n"
		        "                                  of gene list file\n"
		        "  -v [ --binvalue ] arg (=t) specify method of determining value for each bin\n"
		        "                               t  compute total of hit values\n"
		        "                               a  compute average of hit values\n"
				"                               d  compute density of hit values\n"
				"                               c  compute mean per-nt coverage\n"
		        "  -d [ --binloc ] arg (=g)   specify method of determining location of each bin\n"
		        "                               g  interpret relative bin locations as genetic\n"
		        "                                  distance, requires stranded gene list\n"
				"                               p  interpret relative bin locations as physical\n"
				"                                  distance\n"
		        "  -s [ --strands ] arg (=b)  specify handling of strand identifiers:\n"
				"                               b  utilize hits from both strands, also for use\n"
				"                                  with data that lacks strand information\n"
				"                               s  utilize only same-strand hits, requires\n"
		        "                                  -p, -m, -h c, or -h e\n"
		        "                               o  utilize only opposite-strand hits, requires\n"
		        "                                  -p, -m, -h c, or -h e\n"
				"  --nostrand                  indicates no strand specific methods are to be\n"
				"                              used, applies -s b, -l p, -a p, and -d p\n"
				"  --nohead                    suppresses printing of header to output file\n";
		return;
	}
	int s,t,b,h,l,a,v,d,o;
	char *plushits,*minushits,*hits,*genelist,*output,*binfile;
	long start,size,count;
	opt_parser(int argc,char **args) {
		istringstream temp;
		struct option long_options[]={
				{ "help",0,NULL,'u'},
				{"strands",1,NULL,'s'},
				{"threads",1,NULL,'t'},
				{"bintype",1,NULL,'b'},
				{"hittype",1,NULL,'h'},
				{"hitloc",1,NULL,'l'},
				{"anchor",1,NULL,'a'},
				{"binvalue",1,NULL,'v'},
				{"plus",1,NULL,'p'},
				{"minus",1,NULL,'m'},
				{"binloc",1,NULL,'d'},
				{"nostrand",0,NULL,'n'},
				{"nohead",0,NULL,'o'}
		};
		s=0;
		t=1;
		b=1;
		h=4;
		l=2;
		a=0;
		v=0;
		d=0;
		o=0;
		plushits=NULL;
		minushits=NULL;
		hits=NULL;
		genelist=NULL;
		output=NULL;
		binfile=NULL;
		start=0;
		size=0;
		count=0;
		int lflag=0;
		int sflag=0;
		int opt,dummy;
		bool a_flag=0;
		while((opt=getopt_long(argc,args,"us:b:h:l:a:v:p:m:t:d:",long_options,&dummy))!=-1) {
			switch(opt) {
			case 'u':
				usage();
				exit(0);
			case 's':
				sflag=1;
				if(strcmp(optarg,"b")==0) {
					s=0;
				}
				else if(strcmp(optarg,"s")==0) {
					s=1;
				}
				else if(strcmp(optarg,"o")==0) {
					s=2;
				}
				else {
					cout << "Error: \"" << optarg << "\" is not a supported argument for the \"-s\" option\n";
					usage();
					exit(1);
				}
				break;
			case 'b':
				if(strcmp(optarg,"f")==0) {
					b=0;
				}
				else if(strcmp(optarg,"c")==0) {
					b=1;
				}
				else if(strcmp(optarg,"v")==0) {
					b=2;
				}
				else {
					cout << "Error: \"" << optarg << "\" is not a supported argument for the \"-b\" option\n";
					usage();
					exit(1);
				}
				break;
			case 'h':
				if(strcmp(optarg,"g")==0) {
					h=0;
				}
				else if(strcmp(optarg,"b")==0) {
					h=1;
				}
				else if(strcmp(optarg,"e")==0) {
					h=2;
				}
				else if(strcmp(optarg,"c")==0) {
					h=3;
				}
				else if(strcmp(optarg,"G")==0) {
					h=4;
				}
				else {
					cout << "Error \"" << optarg << "\" is not a supported argument for the \"-h\" option\n";
					usage();
					exit(1);
				}
				break;
			case 'l':
				lflag=1;
				if(strcmp(optarg,"s")==0) {
					l=0;
				}
				else if(strcmp(optarg,"e")==0) {
					l=1;
				}
				else if(strcmp(optarg,"p")==0) {
					l=2;
				}
				else if(strcmp(optarg,"d")==0) {
					l=3;
				}
				else if(strcmp(optarg,"c")==0) {
					l=4;
				}
				else {
					cout << "Error \"" << optarg << "\" is not a supported argument for the \"-l\" option\n";
					usage();
					exit(1);
				}
				break;
			case 'a':
				a_flag=1;
				if(strcmp(optarg,"s")==0) {
					a=0;
				}
				else if(strcmp(optarg,"e")==0) {
					a=1;
				}
				else if(strcmp(optarg,"p")==0) {
					a=2;
				}
				else if(strcmp(optarg,"d")==0) {
					a=3;
				}
				else if(strcmp(optarg,"u")==0) {
					a=4;
				}
				else {
					cout << "Error: \"" << optarg << "\" is not a supported argument for the \"-a\" option\n";
					usage();
					exit(1);
				}
				break;
			case 'd':
				if(strcmp(optarg,"g")==0) {
					d=0;
				}
				else if(strcmp(optarg,"p")==0) {
					d=1;
				}
				else {
					cout << "Error: \"" << optarg << "\" is not a supported argument for the \"-d\" option\n";
				}
				break;
			case 'v':
				if(strcmp(optarg,"t")==0) {
					v=0;
				}
				else if(strcmp(optarg,"a")==0) {
					v=1;
				}
				else if(strcmp(optarg,"d")==0) {
					v=2;
				}
				else if(strcmp(optarg,"c")==0) {
					v=3;
				}
				else {
					cout << "Error: \"" << optarg << "\" is not a supported argument for the \"-v\" option\n";
					usage();
					exit(1);
				}
				break;
			case 'p':
				plushits=optarg;
				break;
			case 'm':
				minushits=optarg;
				break;
			case 't':
				temp.str(optarg);
				temp >> t;
				if(temp.fail() || t<1) {
					cout << "Error: -t argument must be an integer value greater than 0\n";
					usage();
					exit(1);
				}
				break;
			case 'n':
				s=0;
				l=2;
				a=2;
				d=1;
				break;
			case 'o':
				o=1;
				break;
			case '?':
				usage();
				exit(1);
			}
		}
		try {
			istringstream temp;
			temp.exceptions(ios::failbit);
			switch(argc-optind) {
			case 0:
				if(plushits==NULL && minushits==NULL) {					//check for appropriate number of arguments given specified options
					throw 0;
				}
				else {
					throw 1;
				}
				break;
			case 1:
				if(plushits==NULL && minushits==NULL) {
					throw 1;
				}
				else {
					throw 2;
				}
				break;
			case 2:
				if(plushits==NULL && minushits==NULL) {
					throw 2;
				}
				else {
					switch(b) {
					case 0:
						throw 3;
						break;
					case 1:
						throw 4;
						break;
					case 2:
						throw 5;
						break;
					}
				}
				break;
			case 3:
				if(plushits==NULL && minushits==NULL) {
					switch(b) {
					case 0:
						throw 3;
						break;
					case 1:
						throw 4;
						break;
					case 2:
						throw 5;
						break;
					}
				}
				else {
					switch(b) {
					case 0:
						genelist=args[optind];
						output=args[optind+1];
						binfile=args[optind+2];
						break;
					case 1:
						throw 6;
						break;
					case 2:
						genelist=args[optind];
						output=args[optind+1];
						temp.str(args[optind+2]);
						temp.clear();
						temp >> count;
						break;
					}
				}
				break;
			case 4:
				if(plushits==NULL && minushits==NULL) {
					switch(b) {
					case 0:
						hits=args[optind];
						genelist=args[optind+1];
						output=args[optind+2];
						binfile=args[optind+3];
						break;
					case 1:
						throw 6;
						break;
					case 2:
						hits=args[optind];
						genelist=args[optind+1];
						output=args[optind+2];
						temp.str(args[optind+3]);
						temp.clear();
						temp >> count;
						break;
					}
				}
				else {
					switch(b) {
					case 1:
						throw 5;
						break;
					case 0:
					case 2:
						throw 7;
						break;
					}
				}
				break;
			case 5:
				if(plushits==NULL && minushits==NULL) {
					switch(b) {
					case 1:
						throw 5;
						break;
					case 0:
					case 2:
						throw 7;
						break;
					}
				}
				else {
					switch(b) {
					case 1:
						genelist=args[optind];
						output=args[optind+1];
						temp.str(args[optind+2]);
						temp >> start;
						temp.str(args[optind+3]);
						temp.clear();
						temp >> size;
						temp.str(args[optind+4]);
						temp.clear();
						temp >> count;
						break;
					case 0:
					case 2:
						throw 7;
						break;
					}
				}
				break;
			case 6:
				if(plushits==NULL && minushits==NULL) {
					switch(b) {
					case 1:
						hits=args[optind];
						genelist=args[optind+1];
						output=args[optind+2];
						temp.str(args[optind+3]);
						temp.clear();
						temp >> start;
						temp.str(args[optind+4]);
						temp.clear();
						temp >> size;
						temp.str(args[optind+5]);
						temp.clear();
						temp >> count;
						break;
					case 0:
					case 2:
						throw 7;
						break;
					}
				}
				else {
					throw 7;
				}
				break;
			default:
				throw 7;
			}
		}
		catch(ios::failure &f) {
			cout << "Error: bin start, size, and count must be integer values\n";
			usage();
			exit(1);
		}
		catch(int e) {
			switch(e) {
			case 0:
				cout << "Error: hit file name must be specified\n";
				break;
			case 1:
				cout << "Error: gene list file name must be specified\n";
				break;
			case 2:
				cout << "Error: output file name must be specified\n";
				break;
			case 3:
				cout << "Error: bin file name must be specified\n";
				break;
			case 4:
				cout << "Error: bin start location must be specified\n";
				break;
			case 5:
				cout << "Error: bin count must be specified\n";
				break;
			case 6:
				cout << "Error: bin size must be specified\n";
				break;
			case 7:
				cout << "Error: too many arguments provided given selected options\n";
				break;
			}
			usage();
			exit(1);
		}
		if((h!=2 && h!=3) && (plushits==NULL && minushits==NULL)) {
			if(s!=0) {
				cout << "Error: specified hit file type requires -p or -m for strand-specific matching\n";
				usage();
				exit(1);
			}
			if(l<2) {
				cout << "Error: specified hit file type requires -p or -m for genetic definition of\n"
						"       hit locations\n";
				usage();
				exit(1);
			}

		}
		if(a_flag==1 && b==2) {
			cout << "Warning: specified anchor position will be ignored due to use of variable\n"
					"         bin size\n";
		}
		if(lflag==1 && v==3) {
			cout << "Warning: coverage mode selected, specified hit location will be ignored\n";
		}
		if(plushits!=NULL || minushits!=NULL) {
			if(sflag==0 && lflag==0) {
				cout << "Warning: strand-specific hit file(s) provided, but default physical start\n"
						"         will be utilized as hit location and strand-specific matching\n"
						"         has not been enabled\n";
			}
			else if(lflag==0 && v!=3) {
				cout << "Warning: strand-specific hit file(s) provided, but default physical start\n"
						"         will be utilized as hit location\n";
			}
			else if(sflag==0) {
				cout << "Warning: strand-specific hit file(s) provided, but strand-specific matching\n"
						"         has not been enabled\n";
			}
		}
	}
};

class bin_parser {																						//calculates and stores relative bin start and end locations
	static bool comp_func(pair<long,long> a,pair<long,long> b) {
		return(a.first<b.first);
	}
	char delim;
	void get_delim(ifstream &ifs) {
		char c;
		while(1) {
			ifs.get(c);
			if(c=='\n') {
				delim='\n';
				ifs.seekg(0);
				break;
			}
			if(c=='\r') {
				ifs.get(c);
				if(c=='\n') {
					delim='\n';
					ifs.seekg(0);
					break;
				}
				else {
					delim='\r';
					ifs.seekg(0);
					break;
				}
			}
		}
	}
public:
	vector<pair<long,long> > bins;
	bin_parser(opt_parser &op) {
		if(op.b==0) {																					//read bin start/end locations from file if requested
			ifstream binfile(op.binfile);
			string binline;
			long start,end;
			get_delim(binfile);
			do {
				getline(binfile,binline,delim);
				istringstream binline_stream(binline);
				binline_stream >> start >> end;
				if(binline_stream.fail()) {
					cout << "Bin file contains formatting errors, exiting: " << binline << endl;
					exit(1);
				}
				bins.push_back(pair<long,long>(start,end));
			} while(!binfile.eof());
			sort(bins.begin(),bins.end(),comp_func);
			binfile.close();
		}
		else if(op.b==1) {																				//generate bin start/end locations from specified start, size, and count
			for(long i=0;i<op.count;i++) {
				bins.push_back(pair<long,long>(op.start+(op.size*i),op.start+(op.size*i)+op.size-1));
			}
		}
		else {																							//for variable bin size, just make appropriately sized vector
			bins.resize(op.count,pair<long,long>(0,0));
		}
	}
};

class file_reader {																			//opens hit file for reading
protected:
	ifstream one;
	string strand;
	char delim;
	char get_delim(ifstream &ifs) {
		char c,d;
		while(1) {
			ifs.get(c);
			if(c=='\n') {
				d='\n';
				ifs.seekg(0);
				break;
			}
			if(c=='\r') {
				ifs.get(c);
				if(c=='\n') {
					d='\n';
					ifs.seekg(0);
					break;
				}
				else {
					d='\r';
					ifs.seekg(0);
					break;
				}
			}
		}
		return(d);
	}
#ifndef SINGLE
	pthread_mutex_t linelock;
#endif
public:
	vector<string> one_hitlist;
	vector<string> two_hitlist;
	file_reader(opt_parser &op) {
#ifndef SINGLE
		pthread_mutex_init(&linelock,NULL);
#endif
		if(op.hits!=NULL) {
			char *pch;
			pch=strtok(op.hits,",");
			while(pch!=NULL) {
				one_hitlist.push_back(pch);
				pch=strtok(NULL,",");
			}
			one.open(one_hitlist[0].c_str());
			if(one.fail()) {
				cout << "Error: could not open hit file \"" << one_hitlist[0] << "\"\n";
				exit(1);
			}
			delim=get_delim(one);
		}
		else if(op.plushits!=NULL) {
			char *pch;
			pch=strtok(op.plushits,",");
			while(pch!=NULL) {
				one_hitlist.push_back(pch);
				pch=strtok(NULL,",");
			}
			pch=strtok(op.minushits,",");
			while(pch!=NULL) {
				two_hitlist.push_back(pch);
				pch=strtok(NULL,",");
			}
			one.open(one_hitlist[0].c_str());
			if(one.fail()) {
				cout << "Error: could not open hit file \"" << one_hitlist[0] << "\"\n";
				exit(1);
			}
			strand="plus";
			delim=get_delim(one);
		}
		else {
			char *pch;
			pch=strtok(op.minushits,",");
			while(pch!=NULL) {
				one_hitlist.push_back(pch);
				pch=strtok(NULL,",");
			}
			one.open(one_hitlist[0].c_str());
			if(one.fail()) {
				cout << "Error: could not open hit file \"" << one_hitlist[0] << "\"\n";
				exit(1);
			}
			strand="minus";
			delim=get_delim(one);
		}
	}
	virtual void read(string &line,string &str,int &ret,int *index) {								//for single hit file passed without -p or -m
#ifndef SINGLE
		pthread_mutex_lock(&linelock);
#endif
		do {
			getline(one,line,delim);
		} while(line=="" && one.good());
		if(!one.good() && (*index)<((int)one_hitlist.size()-1)) {
			one.close();
			(*index)++;
			one.open(one_hitlist[*index].c_str());
			if(one.fail()) {
				cout << "Error: could not open hit file \"" << one_hitlist[*index] << "\"\n";
				exit(1);
			}
			delim=get_delim(one);
			do {
				getline(one,line,delim);
			} while(line=="" && one.good());
		}
		ret=one.good();
#ifndef SINGLE
		pthread_mutex_unlock(&linelock);
#endif
		return;
	}
	virtual ~file_reader() {
		one.close();
	}
};

class file_reader_one : public file_reader {											//for single hit file passed with -p or -m
public:
	file_reader_one(opt_parser &op) : file_reader(op) { }
	void read(string &line,string &str,int &ret,int *index) {
#ifndef SINGLE
		pthread_mutex_lock(&linelock);
#endif
		do {
			getline(one,line,delim);
		} while(line=="" && one.good());
		if(!one.good() && (*index)<((int)one_hitlist.size()-1)) {
			one.close();
			(*index)++;
			one.open(one_hitlist[*index].c_str());
			if(one.fail()) {
				cout << "Error: could not open hit file \"" << one_hitlist[*index] << "\"\n";
				exit(1);
			}
			delim=get_delim(one);
			do {
				getline(one,line,delim);
			} while(line=="" && one.good());
		}
		str=strand;																		//supply strand and return value to calling function
		ret=one.good();
#ifndef SINGLE
		pthread_mutex_unlock(&linelock);
#endif
		return;
	}
};

class file_reader_two : public file_reader {											//for two hit files, passed with -p and -m
	ifstream two;
	ifstream *current;
	char delim2;
	char *current_delim;
public:
	file_reader_two(opt_parser &op) : file_reader(op) {
		two.open(two_hitlist[0].c_str());
		if(two.fail()) {
			cout << "Error: could not open hit file \"" << two_hitlist[0] << "\"\n";
			exit(1);
		}
		strand="plus";
		delim2=get_delim(two);
		current=&one;
		current_delim=&delim;
	}
	void read(string &line,string &str,int &ret,int *index) {
#ifndef SINGLE
		pthread_mutex_lock(&linelock);
#endif
		do {
			getline(*current,line,*current_delim);														//read a line at a time from the plus strand hit file, then switch to minus
		} while(line=="" && current->good());
		if(!current->good()) {
			if(current==&one) {
				current=&two;
				current_delim=&delim2;
				strand="minus";
				do {
					getline(*current,line,*current_delim);
				} while(line=="" && current->good());
			}
			else if((*index)<((int)one_hitlist.size()-1)) {
				one.close();
				two.close();
				(*index)++;
				one.open(one_hitlist[*index].c_str());
				if(one.fail()) {
					cout << "Error: could not open hit file \"" << one_hitlist[*index] << "\"\n";
					exit(1);
				}
				delim=get_delim(one);
				two.open(two_hitlist[*index].c_str());
				if(two.fail()) {
					cout << "Error: could not open hit file \"" << two_hitlist[*index] << "\"\n";
					exit(1);
				}
				delim2=get_delim(two);
				strand="plus";
				current=&one;
				current_delim=&delim;
				do {
					getline(*current,line,*current_delim);
				} while(line=="" && current->good());
			}
		}
		str=strand;
		ret=current->good();
#ifndef SINGLE
		pthread_mutex_unlock(&linelock);
#endif
		return;
	}
	~file_reader_two() {
		two.close();
	}
};

class data {																							//inherited by genelist_parser and hit_parser, stores features from gene list file, bin start/end locations, and counts per bin
protected:
	static unordered_map<string,chr_entry> db;
	static unordered_map<string,unordered_map<string,chr_entry> > db_split;
	static map<string,totals_info> table;
	static unordered_map<string,string> strand_map;
	static int max_index;
	static int index;
	static file_reader *fr;
};

unordered_map<string,chr_entry> data::db;
unordered_map<string,unordered_map<string,chr_entry> > data::db_split;
map<string,totals_info> data::table;
unordered_map<string,string> data::strand_map;
int data::index;
int data::max_index;
file_reader *data::fr;

class genelist_parser : public data {																//reads all lines from gene list file, generates specific bin start/end locations per feature
	vector<ofstream*> outfile;
	char delim;
	void get_delim(ifstream &ifs) {
		char c;
		while(1) {
			ifs.get(c);
			if(c=='\n') {
				delim='\n';
				ifs.seekg(0);
				break;
			}
			if(c=='\r') {
				ifs.get(c);
				if(c=='\n') {
					delim='\n';
					ifs.seekg(0);
					break;
				}
				else {
					delim='\r';
					ifs.seekg(0);
					break;
				}
			}
		}
	}
public:
	genelist_parser(opt_parser &op,bin_parser &bp) {
		ifstream genelist;
		string line;
		string id;
		string desc;
		string chr;
		long physical_start;
		long physical_end;
		string strand;
		long anchor_val;
		index=0;
		genelist.open(op.genelist);
		if(genelist.fail()) {
			cout << "Error: could not open gene list file \"" << op.genelist << "\"\n";
			exit(1);
		}
		get_delim(genelist);
		char *pch;
		pch=strtok(op.output,",");
		while(pch!=NULL) {
			outfile.push_back(new ofstream(pch));
			if(outfile.back()->fail()) {
				cout << "Error: could not create output file \"" << pch << "\"\n";
				exit(1);
			}
			pch=strtok(NULL,",");
		}
		if((int)outfile.size()!=max_index && outfile.size()!=1) {
			cout << "Error: when multiple output files are specified, the count must be equal to the input file count\n";
			exit(1);
		}
		if(op.s==0) {																												//for strand-independent matching
			if(op.d==0 || op.a<2) {																									//if bin distance or anchor utilize strand information
				do {
					getline(genelist,line,delim);
					istringstream temp1(line);
					temp1 >> id >> desc >> chr >> physical_start >> physical_end >> strand;
					if(strand=="+") {
						strand="plus";
					}
					if(strand=="-") {
						strand="minus";
					}
					if(temp1.fail()) {																								//skip lines with bad formatting or duplicate id's
						if(line!="") cout << "Gene list file contains bad line, skipping: " << line << endl;
					}
					else if(db.find(id)!=db.end()) {
						cout << "Gene list file contains duplicate unique identifier, skipping: " << line << endl;
					}
					else {
						db[chr].id.push_back(id);
						db[chr].bins.push_back(vector<pair<long,long> >());
						if(strand=="plus") {
							switch(op.a) {																							//determine correct field for anchor position based on strand
							case 0:
							case 2:
								anchor_val=physical_start;
								break;
							case 1:
							case 3:
								anchor_val=physical_end;
								break;
							case 4:
								istringstream temp2(desc);
								temp2 >> anchor_val;
								if(temp2.fail()) {
									cout << "Gene list file contains bad anchor value, skipping: " << line << endl;
									continue;
								}
								break;
							}
							if(op.b!=2) {																							//for fixed bin size, add relative bin locations to anchor to determine specific bin start/end
								for(vector<pair<long,long> >::iterator i=bp.bins.begin();i!=bp.bins.end();i++) {
									db[chr].bins.back().push_back(pair<long,long>(anchor_val+i->first,anchor_val+i->second));
								}
							}
							else {																														//for variable bin size
								float size=(float)(physical_end-physical_start+1)/(float)op.count;
								if(size>=1) {
									long inc[2]={(long)(floor(size)),(long)(ceil(size))};
									long bigs=(physical_end-physical_start+1)-(inc[0]*op.count);														//distribute remainder of (interval size)/(count) across by adding 1 to size of as many bins as necessary, starting with most downstream bin
									for(long i=0;i<(op.count-bigs);i++) {
										db[chr].bins.back().push_back(pair<long,long>(physical_start+(i*inc[0]),physical_start+((i+1)*inc[0])-1));
									}
									long big_start=physical_start+((op.count-bigs)*inc[0]);
									for(long i=0;i<bigs;i++) {
										db[chr].bins.back().push_back(pair<long,long>(big_start+(i*inc[1]),big_start+((i+1)*inc[1])-1));
									}
								}
								else {
									long inc[2]={(long)(floor(1/size)),(long)(ceil(1/size))};															//if bin count is larger than interval size, create "fractionally" sized bins by creating 1/size copies per position
									long extras=op.count-((physical_end-physical_start+1)*inc[0]);														//create one extra copy of as many of positions as necessary to meet requested count, starting with the most upstream positions
									for(long i=0;i<extras;i++) {
										for(long j=0;j<inc[1];j++) {
											db[chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
										}
									}
									for(long i=extras;i<(physical_end-physical_start+1);i++) {
										for(long j=0;j<inc[0];j++) {
											db[chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
										}
									}
								}
							}
						}
						else if(strand=="minus") {																													//for minus strand features, genetic start and end are swapped
							switch(op.a) {
							case 0:
							case 3:
								anchor_val=physical_end;
								break;
							case 1:
							case 2:
								anchor_val=physical_start;
								break;
							case 4:
								istringstream temp2(desc);
								temp2 >> anchor_val;
								if(temp2.fail()) {
									cout << "Gene list file contains bad anchor value, skipping: " << line << endl;
									continue;
								}
								break;
							}
							if(op.d==0) {																															//flip orientation of bins for minus strand features, if genetic bin distance is specified, but keep sorted low to high coordinate to expedite matching
								if(op.b!=2) {
									for(vector<pair<long,long> >::reverse_iterator i=bp.bins.rbegin();i<bp.bins.rend();i++) {
										db[chr].bins.back().push_back(pair<long,long>(anchor_val-i->second,anchor_val-i->first));
									}
								}
								else {
									float size=(float)(physical_end-physical_start+1)/(float)op.count;
									if(size>=1) {
										long inc[2]={(long)(floor(size)),(long)(ceil(size))};
										long bigs=(physical_end-physical_start+1)-(inc[0]*op.count);
										for(long i=0;i<bigs;i++) {
											db[chr].bins.back().push_back(pair<long,long>(physical_start+(i*inc[1]),physical_start+((i+1)*inc[1])-1));
										}
										long small_start=physical_start+(bigs*inc[1]);
										for(long i=0;i<(op.count-bigs);i++) {
											db[chr].bins.back().push_back(pair<long,long>(small_start+(i*inc[0]),small_start+((i+1)*inc[0])-1));
										}
									}
									else {
										long inc[2]={(long)(floor(1/size)),(long)(ceil(1/size))};
										long extras=op.count-((physical_end-physical_start+1)*inc[0]);
										for(long i=0;i<(physical_end-physical_start+1-extras);i++) {
											for(long j=0;j<inc[0];j++) {
												db[chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
											}
										}
										for(long i=(physical_end-physical_start+1-extras);i<(physical_end-physical_start+1);i++) {
											for(long j=0;j<inc[1];j++) {
												db[chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
											}
										}
									}
								}
							}
							else {																																//otherwise, determine bin start/end as for plus strand features
								if(op.b!=2) {
									for(vector<pair<long,long> >::iterator i=bp.bins.begin();i!=bp.bins.end();i++) {
										db[chr].bins.back().push_back(pair<long,long>(anchor_val+i->first,anchor_val+i->second));
									}
								}
								else {
									float size=(float)(physical_end-physical_start+1)/(float)op.count;
									if(size>=1) {
										long inc[2]={(long)(floor(size)),(long)(ceil(size))};
										long bigs=(physical_end-physical_start+1)-(inc[0]*op.count);
										for(long i=0;i<(op.count-bigs);i++) {
											db[chr].bins.back().push_back(pair<long,long>(physical_start+(i*inc[0]),physical_start+((i+1)*inc[0])-1));
										}
										long big_start=physical_start+((op.count-bigs)*inc[0]);
										for(long i=0;i<bigs;i++) {
											db[chr].bins.back().push_back(pair<long,long>(big_start+(i*inc[1]),big_start+((i+1)*inc[1])-1));
										}
									}
									else {
										long inc[2]={(long)(floor(1/size)),(long)(ceil(1/size))};
										long extras=op.count-((physical_end-physical_start+1)*inc[0]);
										for(long i=0;i<extras;i++) {
											for(long j=0;j<inc[1];j++) {
												db[chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
											}
										}
										for(long i=extras;i<(physical_end-physical_start+1);i++) {
											for(long j=0;j<inc[0];j++) {
												db[chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
											}
										}
									}
								}
							}
						}
						else {
							cout << "Gene list file contains bad strand identifier, skipping: " << strand << endl;
							continue;
						}
						db[chr].bins_start.push_back(db[chr].bins.back().front().first);									//store overall bin start and end locations for easy access
						db[chr].bins_end.push_back(db[chr].bins.back().back().second);
						table[id].chr=chr;																					//create entry for feature in output table
						table[id].desc=desc;
						table[id].start=physical_start;
						table[id].end=physical_end;
						table[id].strand=strand;
						for(size_t i=0;i<bp.bins.size();i++) {																//initialize total and count of intersections to 0
							for(int j=0;j<max_index;j++) {
								if(i==0) {
									table[id].total.push_back(vector<double>());
									table[id].count.push_back(vector<long>());
								}
								table[id].total[j].push_back(0);
								table[id].count[j].push_back(0);
							}
							table[id].length.push_back(db[chr].bins.back()[i].second-db[chr].bins.back()[i].first+1);		//determine all bin lengths for density calculations
						}
					}
				} while(!genelist.eof());
			}
			else {																											//for completely strand-independent operation, determine bins as with physical anchor, genetic bin distance on plus strand
				do {
					getline(genelist,line,delim);
					istringstream temp1(line);
					temp1 >> id >> desc >> chr >> physical_start >> physical_end;
					if(temp1.fail()) {
						if(line!="") cout << "Gene list file contains bad line, skipping: " << line << endl;
					}
					else if(db.find(id)!=db.end()) {
						cout << "Gene list file contains duplicate unique identifier, skipping: " << line << endl;
					}
					else {
						db[chr].id.push_back(id);
						db[chr].bins.push_back(vector<pair<long,long> >());
						switch(op.a) {
						case 2:
							anchor_val=physical_start;
							break;
						case 3:
							anchor_val=physical_end;
							break;
						case 4:
							istringstream temp2(desc);
							temp2 >> anchor_val;
							if(temp2.fail()) {
								cout << "Gene list file contains bad anchor value, skipping: " << line << endl;
								continue;
								}
							break;
						}
						if(op.b!=2) {
							for(vector<pair<long,long> >::iterator i=bp.bins.begin();i!=bp.bins.end();i++) {
								db[chr].bins.back().push_back(pair<long,long>(anchor_val+i->first,anchor_val+i->second));
							}
						}
						else {
							float size=(float)(physical_end-physical_start+1)/(float)op.count;
							if(size>=1) {
								long inc[2]={(long)(floor(size)),(long)(ceil(size))};
								long bigs=(physical_end-physical_start+1)-(inc[0]*op.count);
								for(long i=0;i<(op.count-bigs);i++) {
									db[chr].bins.back().push_back(pair<long,long>(physical_start+(i*inc[0]),physical_start+((i+1)*inc[0])-1));
								}
								long big_start=physical_start+((op.count-bigs)*inc[0]);
								for(long i=0;i<bigs;i++) {
									db[chr].bins.back().push_back(pair<long,long>(big_start+(i*inc[1]),big_start+((i+1)*inc[1])-1));
								}
							}
							else {
								long inc[2]={(long)(floor(1/size)),(long)(ceil(1/size))};
								long extras=op.count-((physical_end-physical_start+1)*inc[0]);
								for(long i=0;i<extras;i++) {
									for(long j=0;j<inc[1];j++) {
										db[chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
									}
								}
								for(long i=extras;i<(physical_end-physical_start+1);i++) {
									for(long j=0;j<inc[0];j++) {
										db[chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
									}
								}
							}
						}
						db[chr].bins_start.push_back(db[chr].bins.back().front().first);
						db[chr].bins_end.push_back(db[chr].bins.back().back().second);
						table[id].chr=chr;
						table[id].desc=desc;
						table[id].start=physical_start;
						table[id].end=physical_end;
						temp1 >> table[id].strand;
						if(table[id].strand.empty()) {
							table[id].strand="NA";
						}
						for(size_t i=0;i<bp.bins.size();i++) {
							for(int j=0;j<max_index;j++) {
								if(i==0) {
									table[id].total.push_back(vector<double>());
									table[id].count.push_back(vector<long>());
								}
								table[id].total[j].push_back(0);
								table[id].count[j].push_back(0);
							}
							table[id].length.push_back(db[chr].bins.back()[i].second-db[chr].bins.back()[i].first+1);
						}
					}
				} while(!genelist.eof());
			}
		}
		else {																														//for same- or opposite-strand matching, same as above, but separate plus and minus strand features
			do {
				getline(genelist,line,delim);
				istringstream temp1(line);
				temp1 >> id >> desc >> chr >> physical_start >> physical_end >> strand;
				if(temp1.fail()) {
					if(line!="") cout << "Gene list file contains bad line, skipping: " << line << endl;
				}
				else if(db_split.find(id)!=db_split.end()) {
					cout << "Gene list file contains duplicate unique identifier, skipping: " << line << endl;
				}
				else {
					if(strand=="+") {
						strand="plus";
					}
					if(strand=="-") {
						strand="minus";
					}
					db_split[strand][chr].id.push_back(id);
					db_split[strand][chr].bins.push_back(vector<pair<long,long> >());
					if(strand=="plus") {
						switch(op.a) {
						case 0:
						case 2:
							anchor_val=physical_start;
							break;
						case 1:
						case 3:
							anchor_val=physical_end;
							break;
						case 4:
							istringstream temp2(desc);
							temp2 >> anchor_val;
							if(temp2.fail()) {
								cout << "Gene list file contains bad anchor value, skipping: " << line << endl;
								continue;
							}
							break;
						}
						if(op.b!=2) {
							for(vector<pair<long,long> >::iterator i=bp.bins.begin();i!=bp.bins.end();i++) {
								db_split[strand][chr].bins.back().push_back(pair<long,long>(anchor_val+i->first,anchor_val+i->second));
							}
						}
						else {
							float size=(float)(physical_end-physical_start+1)/(float)op.count;
							if(size>=1) {
								long inc[2]={(long)(floor(size)),(long)(ceil(size))};
								long bigs=(physical_end-physical_start+1)-(inc[0]*op.count);
								for(long i=0;i<(op.count-bigs);i++) {
									db_split[strand][chr].bins.back().push_back(pair<long,long>(physical_start+(i*inc[0]),physical_start+((i+1)*inc[0])-1));
								}
								long big_start=physical_start+((op.count-bigs)*inc[0]);
								for(long i=0;i<bigs;i++) {
									db_split[strand][chr].bins.back().push_back(pair<long,long>(big_start+(i*inc[1]),big_start+((i+1)*inc[1])-1));
								}
							}
							else {
								long inc[2]={(long)(floor(1/size)),(long)(ceil(1/size))};
								long extras=op.count-((physical_end-physical_start+1)*inc[0]);
								for(long i=0;i<extras;i++) {
									for(long j=0;j<inc[1];j++) {
										db_split[strand][chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
									}
								}
								for(long i=extras;i<(physical_end-physical_start+1);i++) {
									for(long j=0;j<inc[0];j++) {
										db_split[strand][chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
									}
								}
							}
						}
					}
					else if(strand=="minus") {
						switch(op.a) {
						case 0:
						case 3:
							anchor_val=physical_end;
							break;
						case 1:
						case 2:
							anchor_val=physical_start;
							break;
						case 4:
							istringstream temp2(desc);
							temp2 >> anchor_val;
							if(temp2.fail()) {
								cout << "Gene list file contains bad anchor value, skipping: " << line << endl;
								continue;
							}
							break;
						}
						if(op.d==1) {
							if(op.b!=2) {
								for(vector<pair<long,long> >::iterator i=bp.bins.begin();i!=bp.bins.end();i++) {
									db_split[strand][chr].bins.back().push_back(pair<long,long>(anchor_val+i->first,anchor_val+i->second));
								}
							}
							else {
								float size=(float)(physical_end-physical_start+1)/(float)op.count;
								if(size>=1) {
									long inc[2]={(long)(floor(size)),(long)(ceil(size))};
									long bigs=(physical_end-physical_start+1)-(inc[0]*op.count);
									for(long i=0;i<(op.count-bigs);i++) {
										db_split[strand][chr].bins.back().push_back(pair<long,long>(physical_start+(i*inc[0]),physical_start+((i+1)*inc[0])-1));
									}
									long big_start=physical_start+((op.count-bigs)*inc[0]);
									for(long i=0;i<bigs;i++) {
										db_split[strand][chr].bins.back().push_back(pair<long,long>(big_start+(i*inc[1]),big_start+((i+1)*inc[1])-1));
									}
								}
								else {
									long inc[2]={(long)(floor(1/size)),(long)(ceil(1/size))};
									long extras=op.count-((physical_end-physical_start+1)*inc[0]);
									for(long i=0;i<extras;i++) {
										for(long j=0;j<inc[1];j++) {
											db_split[strand][chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
										}
									}
									for(long i=extras;i<(physical_end-physical_start+1);i++) {
										for(long j=0;j<inc[0];j++) {
											db_split[strand][chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
										}
									}
								}
							}
						}
						else {
							if(op.b!=2) {
								for(vector<pair<long,long> >::reverse_iterator i=bp.bins.rbegin();i<bp.bins.rend();i++) {
									db_split[strand][chr].bins.back().push_back(pair<long,long>(anchor_val-i->second,anchor_val-i->first));
								}
							}
							else {
								float size=(float)(physical_end-physical_start+1)/(float)op.count;
								if(size>=1) {
									long inc[2]={(long)(floor(size)),(long)(ceil(size))};
									long bigs=(physical_end-physical_start+1)-(inc[0]*op.count);
									for(long i=0;i<bigs;i++) {
										db_split[strand][chr].bins.back().push_back(pair<long,long>(physical_start+(i*inc[1]),physical_start+((i+1)*inc[1])-1));
									}
									long small_start=physical_start+(bigs*inc[1]);
									for(long i=0;i<(op.count-bigs);i++) {
										db_split[strand][chr].bins.back().push_back(pair<long,long>(small_start+(i*inc[0]),small_start+((i+1)*inc[0])-1));
									}
								}
								else {
									long inc[2]={(long)(floor(1/size)),(long)(ceil(1/size))};
									long extras=op.count-((physical_end-physical_start+1)*inc[0]);
									for(long i=0;i<(physical_end-physical_start+1-extras);i++) {
										for(long j=0;j<inc[0];j++) {
											db_split[strand][chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
										}
									}
									for(long i=(physical_end-physical_start+1-extras);i<(physical_end-physical_start+1);i++) {
										for(long j=0;j<inc[1];j++) {
											db_split[strand][chr].bins.back().push_back(pair<long,long>(physical_start+i,physical_start+i));
										}
									}
								}
							}
						}
					}
					else {
						cout << "Gene list file contains bad strand identifier, skipping: " << strand << endl;
						continue;
					}
					db_split[strand][chr].bins_start.push_back(db_split[strand][chr].bins.back().front().first);
					db_split[strand][chr].bins_end.push_back(db_split[strand][chr].bins.back().back().second);
					table[id].chr=chr;
					table[id].desc=desc;
					table[id].start=physical_start;
					table[id].end=physical_end;
					table[id].strand=strand;
					for(size_t i=0;i<bp.bins.size();i++) {
						for(int j=0;j<max_index;j++) {
							if(i==0) {
								table[id].total.push_back(vector<double>());
								table[id].count.push_back(vector<long>());
							}
							table[id].total[j].push_back(0);
							table[id].count[j].push_back(0);
						}
						table[id].length.push_back(db_split[strand][chr].bins.back()[i].second-db_split[strand][chr].bins.back()[i].first+1);
					}
				}
			} while(!genelist.eof());
			if(op.s==1) {											//lookup table to determine strand of features to check for intersections, given strand of hit file entry
				strand_map["plus"]="plus";							//for same strand matching
				strand_map["minus"]="minus";
			}
			else {													//for opposite strand
				strand_map["plus"]="minus";
				strand_map["minus"]="plus";
			}
		}
		genelist.close();
	}
	void print_header(opt_parser &op,bin_parser &bp) {				//write options specified to output file
		for(size_t k=0;k<outfile.size();k++) {
			if(op.o==0) {
				*(outfile[k]) << "Match Type: ";
				switch(op.s) {
				case 0:
					*(outfile[k]) << "both strands\n";
					break;
				case 1:
					*(outfile[k]) << "same strand\n";
					break;
				case 2:
					*(outfile[k]) << "opposite strand\n";
					break;
				}
				*(outfile[k]) << "Anchor: ";
				switch(op.a) {
				case 0:
					*(outfile[k]) << "genetic start\n";
					break;
				case 1:
					*(outfile[k]) << "genetic end\n";
					break;
				case 2:
					*(outfile[k]) << "physical start\n";
					break;
				case 3:
					*(outfile[k]) << "physical end\n";
					break;
				case 4:
					*(outfile[k]) << "user-specified\n";
					break;
				}
				*(outfile[k]) << "Hit Locations: ";
				switch(op.l) {
				case 0:
					*(outfile[k]) << "genetic start\n";
					break;
				case 1:
					*(outfile[k]) << "genetic end\n";
					break;
				case 2:
					*(outfile[k]) << "physical start\n";
					break;
				case 3:
					*(outfile[k]) << "physical end\n";
					break;
				case 4:
					*(outfile[k]) << "center\n";
					break;
				}
				*(outfile[k]) << "Bin Locations: ";
				switch(op.d) {
				case 0:
					*(outfile[k]) << "genetic\n";
					break;
				case 1:
					*(outfile[k]) << "physical\n";
					break;
				}
				*(outfile[k]) << "Bin Values: ";
				switch(op.v) {
				case 0:
					*(outfile[k]) << "total\n";
					break;
				case 1:
					*(outfile[k]) << "average\n";
					break;
				case 2:
					*(outfile[k]) << "density\n";
				}
				*(outfile[k]) << "Bin Size: ";
				switch(op.b) {
				case 0:
				case 1:
					*(outfile[k]) << "fixed\n";
					break;
				case 2:
					*(outfile[k]) << "variable\n";
					break;
				}
				*(outfile[k]) << "Gene List File: " << op.genelist << '\n';
				*(outfile[k]) << "Hit File(s):";
				if(max_index>1 && outfile.size()==1) {
					*(outfile[k]) << "\t\t\t\t\t\t";
					for(int j=0;j<max_index;j++) {
						if(op.hits!=NULL) {
							*(outfile[k]) << fr->one_hitlist[j];
						}
						else {
							if(op.plushits==NULL || op.minushits==NULL) {
								*(outfile[k]) << fr->one_hitlist[j];
							}
							else {
								*(outfile[k]) << fr->one_hitlist[j];
								*(outfile[k]) << " " << fr->two_hitlist[j];
							}
						}
						for(vector<pair<long,long> >::iterator i=bp.bins.begin();i!=bp.bins.end();i++) {
							*(outfile[k]) << "\t";
						}
					}
				}
				else {
					if(op.hits!=NULL) {
						*(outfile[k]) << " " << fr->one_hitlist[k];
					}
					else {
						if(op.plushits==NULL || op.minushits==NULL) {
							*(outfile[k]) << " " << fr->one_hitlist[k];
						}
						else {
							*(outfile[k]) << " " << fr->one_hitlist[k];
							*(outfile[k]) << " " << fr->two_hitlist[k];
						}
					}
				}
				*(outfile[k]) << '\n';
			}
			*(outfile[k]) << "Gene ID\tDescription\tChromosome\tGene Start\tGene End\tStrand";
			if(max_index>1 && outfile.size()==1) {
				for(int j=0;j<max_index;j++) {
					switch(op.b) {
					case 0:
					case 1:
						for(vector<pair<long,long> >::iterator i=bp.bins.begin();i!=bp.bins.end();i++) {
							*(outfile[k]) << '\t' << i->first << ":" << i->second;
						}
						break;
					case 2:
						for(long i=0;i<op.count;i++) {
							*(outfile[k]) << '\t' << i;
						}
						break;
					}
				}
			}
			else {
				switch(op.b) {
				case 0:
				case 1:
					for(vector<pair<long,long> >::iterator i=bp.bins.begin();i!=bp.bins.end();i++) {
						*(outfile[k]) << '\t' << i->first << ":" << i->second;
					}
					break;
				case 2:
					for(long i=0;i<op.count;i++) {
						*(outfile[k]) << '\t' << i;
					}
					break;
				}
			}
			*(outfile[k]) << "\n";
		}
		return;
	}
	void print_results(opt_parser &op) {																																	//write per-bin counts to output file
		for(size_t l=0;l<outfile.size();l++) {
			if(op.v==0) {
				for(map<string,totals_info>::iterator i=table.begin();i!=table.end();i++) {
					*(outfile[l]) << i->first << '\t' << i->second.desc << '\t' << i->second.chr << '\t' << i->second.start << '\t' << i->second.end << '\t' << i->second.strand;
					if(max_index>1 && outfile.size()==1) {
						for(int m=0;m<max_index;m++) {
							if(op.d==0 && i->second.strand=="minus") {																													//for genetic bin distance, print bins of minus strand features in reverse order
								for(vector<double>::reverse_iterator j=i->second.total[m].rbegin();j<i->second.total[m].rend();j++) {
									*(outfile[l]) << '\t' << *j;
								}
							}
							else {
								for(vector<double>::iterator j=i->second.total[m].begin();j!=i->second.total[m].end();j++) {
									*(outfile[l]) << '\t' << *j;
								}
							}
						}
					}
					else {
						if(op.d==0 && i->second.strand=="minus") {																													//for genetic bin distance, print bins of minus strand features in reverse order
							for(vector<double>::reverse_iterator j=i->second.total[l].rbegin();j<i->second.total[l].rend();j++) {
								*(outfile[l]) << '\t' << *j;
							}
						}
						else {
							for(vector<double>::iterator j=i->second.total[l].begin();j!=i->second.total[l].end();j++) {
								*(outfile[l]) << '\t' << *j;
							}
						}
					}
					*(outfile[l]) << '\n';
				}
			}
			else if(op.v==1) {																																					//if bin average is requested, divide total by number of hit file entries intersecting bin
				double average;
				for(map<string,totals_info>::iterator i=table.begin();i!=table.end();i++) {
					if(op.d==0 && i->second.strand=="minus") {
						vector<double>::reverse_iterator j=i->second.total[l].rbegin();
						vector<long>::reverse_iterator k=i->second.count[l].rbegin();
						*(outfile[l]) << i->first << '\t' << i->second.desc << '\t' << i->second.chr << '\t' << i->second.start << '\t' << i->second.end << '\t' << i->second.strand;
						if(max_index>1 && outfile.size()==1) {
							for(int m=0;m<max_index;m++) {
								for(;j<i->second.total[m].rend();j++) {
									if(*k!=0) {
										average=(*j)/(*k);
									}
									else {
										average=0;
									}
									*(outfile[l]) << '\t' << average;
									k++;
								}
								j=i->second.total[m+1].rbegin();
								k=i->second.count[m+1].rbegin();
							}
						}
						else {
							for(;j<i->second.total[l].rend();j++) {
								if(*k!=0) {
									average=(*j)/(*k);
								}
								else {
									average=0;
								}
								*(outfile[l]) << '\t' << average;
								k++;
							}
						}
					}
					else {
						vector<double>::iterator j=i->second.total[l].begin();
						vector<long>::iterator k=i->second.count[l].begin();
						*(outfile[l]) << i->first << '\t' << i->second.desc << '\t' << i->second.chr << '\t' << i->second.start << '\t' << i->second.end << '\t' << i->second.strand;
						if(max_index>1 && outfile.size()==1) {
							for(int m=0;m<max_index;m++) {
								for(;j!=i->second.total[m].end();j++) {
									if(*k!=0) {
										average=(*j)/(*k);
									}
									else {
										average=0;
									}
									*(outfile[l]) << '\t' << average;
									k++;
								}
								j=i->second.total[m+1].begin();
								k=i->second.count[m+1].begin();
							}
						}
						else {
							for(;j!=i->second.total[l].end();j++) {
								if(*k!=0) {
									average=(*j)/(*k);
								}
								else {
									average=0;
								}
								*(outfile[l]) << '\t' << average;
								k++;
							}
						}
					}
					*(outfile[l]) << '\n';
				}
			}
			else {																																								//if bin density is requested, divide total by bin size
				double density;
				for(map<string,totals_info>::iterator i=table.begin();i!=table.end();i++) {
					if(op.d==0 && i->second.strand=="minus") {
						vector<double>::reverse_iterator j=i->second.total[l].rbegin();
						vector<long>::reverse_iterator k=i->second.length.rbegin();
						*(outfile[l]) << i->first << '\t' << i->second.desc << '\t' << i->second.chr << '\t' << i->second.start << '\t' << i->second.end << '\t' << i->second.strand;
						if(max_index>1 && outfile.size()==1) {
							for(int m=0;m<max_index;m++) {
								for(;j!=i->second.total[m].rend();j++) {
									if(*k!=0) {
										density=(*j)/(*k);
									}
									else {
										density=0;
									}
									*(outfile[l]) << '\t' << density;
									k++;
								}
								j=i->second.total[m+1].rbegin();
								k=i->second.length.rbegin();
							}
						}
						else {
							for(;j!=i->second.total[l].rend();j++) {
								if(*k!=0) {
									density=(*j)/(*k);
								}
								else {
									density=0;
								}
								*(outfile[l]) << '\t' << density;
								k++;
							}
						}
					}
					else {
						vector<double>::iterator j=i->second.total[l].begin();
						vector<long>::iterator k=i->second.length.begin();
						*(outfile[l]) << i->first << '\t' << i->second.desc << '\t' << i->second.chr << '\t' << i->second.start << '\t' << i->second.end << '\t' << i->second.strand;
						if(max_index>1 && outfile.size()==1) {
							for(int m=0;m<max_index;m++) {
								for(;j!=i->second.total[m].end();j++) {
									if(*k!=0) {
										density=(*j)/(*k);
									}
									else {
										density=0;
									}
									*(outfile[l]) << '\t' << density;
									k++;
								}
								j=i->second.total[m+1].begin();
								k=i->second.length.begin();
							}
						}
						else {
							for(;j!=i->second.total[l].end();j++) {
								if(*k!=0) {
									density=(*j)/(*k);
								}
								else {
									density=0;
								}
								*(outfile[l]) << '\t' << density;
								k++;
							}
						}
					}
					*(outfile[l]) << '\n';
				}
			}
			outfile[l]->close();
		}
		return;
	}
};

class hit_parser : public data {
	static bool comp_func_ub(long a,pair<long,long> b) {
		return(a<=b.second);
	}
protected:
#ifndef SINGLE
	pthread_mutex_t tablelock;
#endif
	int s,l,v;
public:
	hit_parser(opt_parser &op) {
#ifndef SINGLE
		pthread_mutex_init(&tablelock,NULL);
#endif
		s=op.s;
		l=op.l;
		v=op.v;
		if(op.plushits==NULL && op.minushits==NULL) {
			fr=new file_reader(op);
		}
		else if(op.plushits==NULL || op.minushits==NULL) {
			fr=new file_reader_one(op);
		}
		else {
			fr=new file_reader_two(op);
			if(fr->one_hitlist.size()!=fr->two_hitlist.size()) {
				cout << "Error: the counts of plus and minus strand input files must be equal\n";
				exit(1);
			}
		}
		max_index=fr->one_hitlist.size();
	}
	virtual ~hit_parser() {
		delete fr;
	}
	virtual int update(hit_entry*)=0;
	inline void set_location(hit_entry *he,string &str,long &start,long &end) {							//determine location to be intersected with bins
		if(str.empty() || str=="plus") {																//for plus strand or strand-independent
			switch(l) {
			case 0:
			case 2:
				he->location=start;
				break;
			case 1:
			case 3:
				he->location=end;
				break;
			case 4:
				he->location=(start+end)/2;
				break;
			}
		}
		else {																							//for minus strand
			switch(l) {
			case 0:
			case 3:
				he->location=end;
				break;
			case 1:
			case 2:
				he->location=start;
				break;
			case 4:
				he->location=(start+end)/2;
				break;
			}
		}
	}
	void query(void) {																																//performs intersection of hit location and bins of all features
		ptr_entry arrays;
		hit_entry he;
		if(v!=3) {
			if(s==0) {																																	//strand-independent matching
				while(update(&he)) {
					if(db.find(he.chr)!=db.end()) {
						size_t max=db[he.chr].id.size();
						arrays.id=&db[he.chr].id[0];																									//store pointers to bin information to reduce lookups
						arrays.bins_start=&db[he.chr].bins_start[0];
						arrays.bins_end=&db[he.chr].bins_end[0];
						arrays.bins=&db[he.chr].bins[0];
						for(size_t i=0;i<max;i++) {
							if(he.location<arrays.bins_start[i] || he.location>arrays.bins_end[i]) continue;											//move to next gene list feature if hit locations falls outside of overall bin start and end
							vector<pair<long,long> >::iterator j=upper_bound(arrays.bins[i].begin(),arrays.bins[i].end(),he.location,comp_func_ub);		//find first bin with end coordinate greater than or equal to hit location
							if(he.location>=j->first) {																									//ensure hit location is also greater than or equal to bin start
#ifndef SINGLE
								pthread_mutex_lock(&tablelock);
								table[arrays.id[i]].total[index][j-arrays.bins[i].begin()]+=he.value;															//add value to bin total, increment intersection count
								table[arrays.id[i]].count[index][j-arrays.bins[i].begin()]++;
								pthread_mutex_unlock(&tablelock);
#else
								table[arrays.id[i]].total[index][j-arrays.bins[i].begin()]+=he.value;															//add value to bin total, increment intersection count
								table[arrays.id[i]].count[index][j-arrays.bins[i].begin()]++;
#endif
							}
						}
					}
				}
			}
			else {
				while(update(&he)) {																													//for strand specific matching check same or opposite strand features only
					string str=strand_map[he.strand];
					if(db_split.find(str)!=db_split.end()) {
						if(db_split[str].find(he.chr)!=db_split[str].end()) {
							size_t max=db_split[str][he.chr].id.size();
							arrays.id=&db_split[str][he.chr].id[0];
							arrays.bins_start=&db_split[str][he.chr].bins_start[0];
							arrays.bins_end=&db_split[str][he.chr].bins_end[0];
							arrays.bins=&db_split[str][he.chr].bins[0];
							for(size_t i=0;i<max;i++) {
								if(he.location<arrays.bins_start[i] || he.location>arrays.bins_end[i]) continue;
								vector<pair<long,long> >::iterator j=upper_bound(arrays.bins[i].begin(),arrays.bins[i].end(),he.location,comp_func_ub);
								if(he.location>=j->first) {
#ifndef SINGLE
									pthread_mutex_lock(&tablelock);
									table[arrays.id[i]].total[index][j-arrays.bins[i].begin()]+=he.value;
									table[arrays.id[i]].count[index][j-arrays.bins[i].begin()]++;
									pthread_mutex_unlock(&tablelock);
#else
									table[arrays.id[i]].total[index][j-arrays.bins[i].begin()]+=he.value;
									table[arrays.id[i]].count[index][j-arrays.bins[i].begin()]++;
#endif
								}
							}
						}
					}
				}
			}
		}
		else {
			long int_length;
			if(s==0) {																																	//strand-independent matching
				while(update(&he)) {
					if(db.find(he.chr)!=db.end()) {
						size_t max=db[he.chr].id.size();
						arrays.id=&db[he.chr].id[0];																									//store pointers to bin information to reduce lookups
						arrays.bins_start=&db[he.chr].bins_start[0];
						arrays.bins_end=&db[he.chr].bins_end[0];
						arrays.bins=&db[he.chr].bins[0];
						for(size_t i=0;i<max;i++) {
							if(he.location2<arrays.bins_start[i] || he.location>arrays.bins_end[i]) continue;											//move to next gene list feature if hit locations falls outside of overall bin start and end
							vector<pair<long,long> >::iterator j=upper_bound(arrays.bins[i].begin(),arrays.bins[i].end()-1,he.location2,comp_func_ub);		//find first bin with end coordinate greater than or equal to hit location
							int_length=std::min(he.location2,j->second)-std::max(he.location,j->first)+1;
							while(int_length>0) {
#ifndef SINGLE
								pthread_mutex_lock(&tablelock);
								table[arrays.id[i]].total[index][j-arrays.bins[i].begin()]+=(he.value*int_length);
								table[arrays.id[i]].count[index][j-arrays.bins[i].begin()]++;
								pthread_mutex_unlock(&tablelock);
#else
								table[arrays.id[i]].total[index][j-arrays.bins[i].begin()]+=(he.value*int_length);
								table[arrays.id[i]].count[index][j-arrays.bins[i].begin()]++;
#endif
								if(j==arrays.bins[i].begin()) break;
								j--;
								int_length=std::min(he.location2,j->second)-std::max(he.location,j->first)+1;
							}
						}
					}
				}
			}
			else {
				while(update(&he)) {
					string str=strand_map[he.strand];
					if(db_split.find(str)!=db_split.end()) {
						if(db_split[str].find(he.chr)!=db_split[str].end()) {
							size_t max=db_split[str][he.chr].id.size();
							arrays.id=&db_split[str][he.chr].id[0];
							arrays.bins_start=&db_split[str][he.chr].bins_start[0];
							arrays.bins_end=&db_split[str][he.chr].bins_end[0];
							arrays.bins=&db_split[str][he.chr].bins[0];
							for(size_t i=0;i<max;i++) {
								if(he.location2<arrays.bins_start[i] || he.location>arrays.bins_end[i]) continue;
								vector<pair<long,long> >::iterator j=upper_bound(arrays.bins[i].begin(),arrays.bins[i].end()-1,he.location2,comp_func_ub);
								int_length=std::min(he.location2,j->second)-std::max(he.location,j->first)+1;
								while(int_length>0) {
#ifndef SINGLE
									pthread_mutex_lock(&tablelock);
									table[arrays.id[i]].total[index][j-arrays.bins[i].begin()]+=(he.value*int_length);
									table[arrays.id[i]].count[index][j-arrays.bins[i].begin()]++;
									pthread_mutex_unlock(&tablelock);
#else
									table[arrays.id[i]].total[index][j-arrays.bins[i].begin()]+=(he.value*int_length);
									table[arrays.id[i]].count[index][j-arrays.bins[i].begin()]++;
#endif
									if(j==arrays.bins[i].begin()) break;
									j--;
									int_length=std::min(he.location2,j->second)-std::max(he.location,j->first)+1;
								}
							}
						}
					}
				}
			}
		}
		return;
	}
};



class hit_parser_g : public hit_parser {
public:
	hit_parser_g(opt_parser &op) : hit_parser(op) { }									//interprets bedgraph file input
	int update(hit_entry *he) {
		int ret;
		string line,str;
		long start,end;
		fr->read(line,str,ret,&index);
		istringstream linestream(line);
		linestream >> he->chr >> start >> end >> he->value;
		if(linestream.fail() && ret==1) {
			if(he->chr!="track" && line!="\r") {
				cout << "Hit file contains bad line, skipping: " << line << "\n";
			}
			return(update(he));															//loop recursively until a "good" line is discovered
		}
		he->strand=str;
		set_location(he,str,start,end);
		return(ret);
	}
};

class hit_parser_G : public hit_parser {
public:
	hit_parser_G(opt_parser &op) : hit_parser(op) { }									//interprets bedgraph file input
	int update(hit_entry *he) {
		int ret;
		string line,str;
		long start,end;
		fr->read(line,str,ret,&index);
		istringstream linestream(line);
		linestream >> he->chr >> start >> end >> he->value;
		if(linestream.fail() && ret==1) {
			if(he->chr!="track" && line!="\r") {
				cout << "Hit file contains bad line, skipping: " << line << "\n";
			}
			return(update(he));															//loop recursively until a "good" line is discovered
		}
		he->strand=str;
		start++;
		set_location(he,str,start,end);
		return(ret);
	}
};

class hit_parser_b : public hit_parser {												//interprets basic bed file input
public:
	hit_parser_b(opt_parser &op) : hit_parser(op) { }
	int update(hit_entry *he) {
		int ret;
		string line,str;
		long start,end;
		fr->read(line,str,ret,&index);
		istringstream linestream(line);
		linestream >> he->chr >> start >> end;
		if(linestream.fail() && ret==1) {
			if(he->chr!="track" && line!="\r") {
				cout << "Hit file contains bad line, skipping: " << line << "\n";
			}
			return(update(he));
		}
		start++;
		he->value=1;
		he->strand=str;
		set_location(he,str,start,end);
		return(ret);
	}
};

class hit_parser_c : public hit_parser {												//interprets cppmatch file input (unstranded)
public:
	hit_parser_c(opt_parser &op) : hit_parser(op) { }
	int update (hit_entry *he) {
		int ret;
		string line,str,temp;
		long start,end;
		fr->read(line,str,ret,&index);
		istringstream linestream(line);
		linestream >> temp >> he->value >> he->chr >> start >> end;
		if(linestream.fail() && ret==1) {
			cout << "Hit file contains bad line, skipping: " << line << "\n";
			return(update(he));
		}
		set_location(he,str,start,end);
		return(ret);
	}
};

class hit_parser_cs : public hit_parser {												//interprets cppmatch file input (stranded)
public:
	hit_parser_cs(opt_parser &op) : hit_parser(op) { }
	int update(hit_entry *he) {
		int ret;
		string line,str,temp;
		long start,end;
		fr->read(line,str,ret,&index);
		istringstream linestream(line);
		linestream >> temp >> he->value >> he->chr >> start >> end >> str;
		if(str=="+") {
			str="plus";
		}
		if(str=="-") {
			str="minus";
		}
		if(str!="plus" && str!="minus") {
			linestream.clear(linestream.failbit);
		}
		if(linestream.fail() && ret==1) {
			cout << "Hit file contains bad line, skipping: " << line << "\n";
			return(update(he));
		}
		he->strand=str;
		set_location(he,str,start,end);
		return(ret);
	}
};

class hit_parser_e : public hit_parser {												//interprets extended bed file input (unstranded)
public:
	hit_parser_e(opt_parser &op) : hit_parser(op) { }
	int update(hit_entry *he) {
		int ret;
		string line,str,temp,temp2;
		long start,end;
		fr->read(line,str,ret,&index);
		he->strand=str;
		istringstream linestream(line);
		linestream >> he->chr >> start >> end >> temp >> temp2;
		if(linestream.fail() && ret==1) {
			if(he->chr!="track" && line!="\r") {
				cout << "Hit file contains bad line, skipping: " << line << "\n";
			}
			return(update(he));
		}
		start++;
		he->value=1;
		set_location(he,str,start,end);
		return(ret);
	}
};

class hit_parser_es : public hit_parser {												//interprets extended bed file input (stranded)
public:
	hit_parser_es(opt_parser &op) : hit_parser(op) { }
	int update(hit_entry *he) {
		int ret;
		string line,str,temp,temp2;
		char prestr;
		long start,end;
		fr->read(line,str,ret,&index);
		istringstream linestream(line);
		linestream >> he->chr >> start >> end >> temp >> temp2 >> prestr;
		switch(prestr) {
		case '-':
			he->strand="minus";
			break;
		case '+':
			he->strand="plus";
			break;
		default:
			linestream.clear(linestream.failbit);
			break;
		}
		if(linestream.fail() && ret==1) {
			if(he->chr!="track" && line!="\r") {
				cout << "Hit file contains bad line, skipping: " << line << "\n";
			}
			return(update(he));
		}
		start++;
		he->value=1;
		set_location(he,str,start,end);
		return(ret);
	}
};

class hit_parser_gc : public hit_parser {
public:
	hit_parser_gc(opt_parser &op) : hit_parser(op) { }									//interprets bedgraph file input
	int update(hit_entry *he) {
		int ret;
		string line,str;
		fr->read(line,str,ret,&index);
		istringstream linestream(line);
		linestream >> he->chr >> he->location >> he->location2 >> he->value;
		if(linestream.fail() && ret==1) {
			if(he->chr!="track" && line!="\r") {
				cout << "Hit file contains bad line, skipping: " << line << "\n";
			}
			return(update(he));															//loop recursively until a "good" line is discovered
		}
		he->strand=str;
		return(ret);
	}
};

class hit_parser_Gc : public hit_parser {
public:
	hit_parser_Gc(opt_parser &op) : hit_parser(op) { }									//interprets bedgraph file input
	int update(hit_entry *he) {
		int ret;
		string line,str;
		fr->read(line,str,ret,&index);
		istringstream linestream(line);
		linestream >> he->chr >> he->location >> he->location2 >> he->value;
		if(linestream.fail() && ret==1) {
			if(he->chr!="track" && line!="\r") {
				cout << "Hit file contains bad line, skipping: " << line << "\n";
			}
			return(update(he));															//loop recursively until a "good" line is discovered
		}
		he->strand=str;
		he->location++;
		return(ret);
	}
};

class hit_parser_bc : public hit_parser {												//interprets basic bed file input
public:
	hit_parser_bc(opt_parser &op) : hit_parser(op) { }
	int update(hit_entry *he) {
		int ret;
		string line,str;
		fr->read(line,str,ret,&index);
		istringstream linestream(line);
		linestream >> he->chr >> he->location >> he->location2;
		if(linestream.fail() && ret==1) {
			if(he->chr!="track" && line!="\r") {
				cout << "Hit file contains bad line, skipping: " << line << "\n";
			}
			return(update(he));
		}
		he->location++;
		he->value=1;
		he->strand=str;
		return(ret);
	}
};

class hit_parser_cc : public hit_parser {												//interprets cppmatch file input (unstranded)
public:
	hit_parser_cc(opt_parser &op) : hit_parser(op) { }
	int update (hit_entry *he) {
		int ret;
		string line,str,temp;
		fr->read(line,str,ret,&index);
		istringstream linestream(line);
		linestream >> temp >> he->value >> he->chr >> he->location >> he->location2;
		if(linestream.fail() && ret==1) {
			cout << "Hit file contains bad line, skipping: " << line << "\n";
			return(update(he));
		}
		return(ret);
	}
};

class hit_parser_csc : public hit_parser {												//interprets cppmatch file input (stranded)
public:
	hit_parser_csc(opt_parser &op) : hit_parser(op) { }
	int update(hit_entry *he) {
		int ret;
		string line,str,temp;
		fr->read(line,str,ret,&index);
		istringstream linestream(line);
		linestream >> temp >> he->value >> he->chr >> he->location >> he->location2 >> str;
		if(str=="+") {
			str="plus";
		}
		if(str=="-") {
			str="minus";
		}
		if(str!="plus" && str!="minus") {
			linestream.clear(linestream.failbit);
		}
		if(linestream.fail() && ret==1) {
			cout << "Hit file contains bad line, skipping: " << line << "\n";
			return(update(he));
		}
		he->strand=str;
		return(ret);
	}
};

class hit_parser_ec : public hit_parser {												//interprets extended bed file input (unstranded)
public:
	hit_parser_ec(opt_parser &op) : hit_parser(op) { }
	int update(hit_entry *he) {
		int ret;
		string line,str,temp,temp2;
		fr->read(line,str,ret,&index);
		he->strand=str;
		istringstream linestream(line);
		linestream >> he->chr >> he->location >> he->location2 >> temp >> temp2;
		if(linestream.fail() && ret==1) {
			if(he->chr!="track" && line!="\r") {
				cout << "Hit file contains bad line, skipping: " << line << "\n";
			}
			return(update(he));
		}
		he->location++;
		he->value=1;
		return(ret);
	}
};

class hit_parser_esc : public hit_parser {												//interprets extended bed file input (stranded)
public:
	hit_parser_esc(opt_parser &op) : hit_parser(op) { }
	int update(hit_entry *he) {
		int ret;
		string line,str,temp,temp2;
		char prestr;
		fr->read(line,str,ret,&index);
		istringstream linestream(line);
		linestream >> he->chr >> he->location >> he->location2 >> temp >> temp2 >> prestr;
		switch(prestr) {
		case '-':
			he->strand="minus";
			break;
		case '+':
			he->strand="plus";
			break;
		default:
			linestream.clear(linestream.failbit);
			break;
		}
		if(linestream.fail() && ret==1) {
			if(he->chr!="track" && line!="\r") {
				cout << "Hit file contains bad line, skipping: " << line << "\n";
			}
			return(update(he));
		}
		he->location++;
		he->value=1;
		return(ret);
	}
};

#ifndef SINGLE
void *t_query(void *hit) {															//reads lines from the hit file(s) and performs intersections, exits when no more lines are available
	hit_parser *hp=reinterpret_cast<hit_parser*>(hit);
	hp->query();
	pthread_exit(NULL);
}
#endif

int main(int argc,char** args) {
	opt_parser op(argc,args);
	bin_parser bp(op);
	hit_parser *hp;
	switch(op.h) {																	//create appropriate object given input file type
	case 0:
		if(op.v==3) {
			hp=new hit_parser_gc(op);
		}
		else {
			hp=new hit_parser_g(op);
		}
		break;
	case 1:
		if(op.v==3) {
			hp=new hit_parser_bc(op);
		}
		else {
			hp=new hit_parser_b(op);
		}
		break;
	case 2:
		if(op.hits!=NULL && op.l<2) {
			if(op.v==3) {
				hp=new hit_parser_esc(op);
			}
			else {
				hp=new hit_parser_es(op);
			}
		}
		else {
			if(op.v==3) {
				hp=new hit_parser_ec(op);
			}
			else {
				hp=new hit_parser_e(op);
			}
		}
		break;
	case 3:
		if(op.hits!=NULL && op.l<2) {
			if(op.v==3) {
				hp=new hit_parser_csc(op);
			}
			else {
				hp=new hit_parser_cs(op);
			}
		}
		else {
			if(op.v==3) {
				hp=new hit_parser_cc(op);
			}
			else {
				hp=new hit_parser_c(op);
			}
		}
		break;
	case 4:
		if(op.v==3) {
			hp=new hit_parser_Gc(op);
		}
		else {
			hp=new hit_parser_G(op);
		}
		break;
	default:
		cout << "Error: incompatible input file format and hit location specified\n";
		opt_parser::usage();
		return(1);
	}
	genelist_parser glp(op,bp);
#ifndef SINGLE
	if(op.t>1) {																	//create specified number of threads to read hit file(s) and perform intersections
		pthread_t tid[op.t];
		int ti;
		for(ti=0;ti<op.t;ti++) {
			pthread_create(&tid[ti],NULL,t_query,reinterpret_cast<void*>(hp));
		}
		for(ti=0;ti<op.t;ti++) {													//wait until all threads exit
			pthread_join(tid[ti],NULL);
		}
	}
	else {																			//for single thread, just perform intersections
		hp->query();
	}
#else
	hp->query();
#endif
	glp.print_header(op,bp);														//print results
	glp.print_results(op);
	delete hp;
	return(0);
}
