#ifndef VERSION_H
#define VERSION_H

#include <iostream>
#include <cstring>

using namespace std;

//file version.h

class version
{
 private:

	char acc[20];
	int ver;

 public:
	version();
	friend ostream &operator<<(ostream &ous, const version &v);
	bool operator<(const version &v) const {return (strcmp(acc, v.acc)< 0);}
	bool operator==(const version &v) const {return (strcmp(acc, v.acc) == 0);}
	version& operator=(const version &v);
	void get(char *ll);  //to get from  L22785.1  GI:438913
	void quickGet(const char *ll);

	/** the version number is larger*/
	bool newer(const version &v) const {return ver > v.ver;}
};

#endif
