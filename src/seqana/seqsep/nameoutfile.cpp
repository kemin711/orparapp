#include "nameoutfile.h"
#include <cstring>

/* to name output file for this and other outputfile names 
 * ipvf: input file name;  s: file for this species
 * nons: files for other species;  
 * name: name of the organism 
 * */

void nameFile(const char ipf[], string &s, string &nons, const string &name) {
	const char *p = strchr(ipf, '.');
	if (p) {
		s = string(ipf, p-ipf);
	}
	else s = ipf;
	nons = s;
	s += '.';
	s += name;
	nons += ".rest";
}

