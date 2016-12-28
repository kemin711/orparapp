#include "version.h"
#include <cstdlib>
#include <cctype>


version::version()
	:ver(1)
{
	acc[0]='\0';
}

void version::quickGet(const char *ll) {
	//only from L22785.2  
	//if missing the version .2 then give it a default verion 1
	int i = 0;
	while (ll[i] != '.' && ll[i] != '\0') acc[i] = ll[i++];
	acc[i]='\0';
	if (ll[i] == '.') ver = atoi(ll+i+1);
}
	
void version::get(char *ll) {
	char *pp = ll+1;
	while (!isspace(*pp)) pp++;
	*pp = '\0';
	quickGet(ll);
}

version& version::operator=(const version &v)
{
	if (this != &v) {
		strcpy(acc, v.acc);
		ver = v.ver;
	}
	return *this;
}

ostream &operator<<(ostream &ous, const version &v)
{
	ous << v.acc << '.' << v.ver;
	return ous;
}

