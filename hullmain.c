#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>

/*
 * Ken Clarkson wrote this.  Copyright (c) 1995 by AT&T..
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

#define POINTSITES 1

#include "hull.h"

point	site_blocks[MAXBLOCKS];
int	num_blocks;

/* int	getopt(int, char**, char*); */
extern char	*optarg;
extern int optind;
extern int opterr;

static long num_sites;
static short vd = 0;
static int dim;

/*FILE *INFILE, *OUTFILE, *DFILE = stderr, *TFILE;*/ FILE *INFILE, *OUTFILE, *TFILE, *DFILE;

long site_numm(site p) {

	long i,j;

	if (vd && p==infinity) return -1;
	if (!p) return -2;
	for (i=0; i<num_blocks; i++)
		if ((j=p-site_blocks[i])>=0 && j < BLOCKSIZE*dim) 
			return j/dim+BLOCKSIZE*i;
	return -3;
}


site new_site (site p, long j) {

	assert(num_blocks+1<MAXBLOCKS);
	if (0==(j%BLOCKSIZE)) {
		assert(num_blocks < MAXBLOCKS);
		return(site_blocks[num_blocks++]=(site)malloc(BLOCKSIZE*site_size));
	} else
		return p+dim;
}

site read_next_site(long j){

	int i=0, k=0;
	static char buf[100], *s;

	if (j!=-1) p = new_site(p,j);
	if (j!=0) while ((s=fgets(buf,sizeof(buf),INFILE))) {
 		if (buf[0]=='%') continue;
		for (k=0; buf[k] && isspace(buf[k]); k++);
		if (buf[k]) break;
	}
	if (!s) return 0;
	if (j!=0) fprintf(TFILE, "%s", buf+k);
	while (buf[k]) {
		while (buf[k] && isspace(buf[k])) k++;
		if (buf[k] && j!=-1) {
			if (sscanf(buf+k,"%lf",p+i)==EOF) {
				fprintf(DFILE, "bad input line: %s\n", buf);
				exit(1);
			}
			p[i] = floor(mult_up*p[i]+0.5);   
			mins[i] = (mins[i]<p[i]) ? mins[i] : p[i];
			maxs[i] = (maxs[i]>p[i]) ? maxs[i] : p[i];
		}
		if (buf[k]) i++;
		while (buf[k] && !isspace(buf[k])) k++;
	}

	if (!dim) dim = i;
	if (i!=dim) {DEB(-10,inconsistent input);DEBTR(-10); exit(1);}	
	return p;
}




/*
site read_next_site(long j){

	int i;
	double pi;

	p = new_site(p,j);
	for (i=0; (i<dim) && (fscanf(INFILE,"%lf",p+i)!=EOF); i++) {
		pi = p[i] *= mult_up;
		p[i] = floor(p[i]+0.5);   
		if (abs(p[i]-pi)>1) {DEB(-3,uhoh);DEBTR(-3); exit(1);}
		mins[i] = (mins[i]<p[i]) ? mins[i] : p[i];
		maxs[i] = (maxs[i]>p[i]) ? maxs[i] : p[i];
	}	
	
	if (i==0) return NULL;
	assert(i==dim);
	return p;
}
*/

site get_site_offline(long i) {

	if (i>=num_sites) return NULL;
	else return site_blocks[i/BLOCKSIZE]+(i%BLOCKSIZE)*dim;
}


long *shufmat;
void make_shuffle(void){
	long i,t,j;
	static long mat_size = 0;

	if (mat_size<=num_sites) {
		mat_size = num_sites+1;
		shufmat = (long*)malloc(mat_size*sizeof(long));
	}
	for (i=0;i<=num_sites;i++) shufmat[i] = i;
	for (i=0;i<num_sites;i++){
		t = shufmat[i];
		shufmat[i] = shufmat[j = i + (num_sites-i)*double_rand()];
		shufmat[j] = t;
	}
}

static long (*shuf)(long);
long noshuffle(long i) {return i;}
long shufflef(long i) {return shufmat[i];}

static site (*get_site_n)(long);
site get_next_site(void) {
	static long s_num = 0;
	return (*get_site_n)((*shuf)(s_num++)); }


void errline(char *s) {fprintf(stderr, s); fprintf(stderr,"\n"); return;} void tell_options(void) {

	errline("options:");
	errline( "-m mult  multiply by mult before rounding;");
	errline( "-d compute delaunay triangulation;");
	errline( "-s seed  shuffle with srand(seed);");
	errline( "-r 	  shuffle with srand(time(0));");
	errline( "-i<name> read input from <name>;");
	errline( "-X<name> chatter to <name>;");
	errline( "-f<fmt> main output in <fmt>:");
	errline("	ps->postscript, mp->metapost, cpr->cpr format, off->OFF format, vn->vertex numbers(default)");
	errline( "-A alpha shape, find suitable alpha");
	errline( "-aa<alpha> alpha shape, alpha=<alpha>;");
	errline( "-af<fmt> output alpha shape in <fmt>;");
	errline( "-oo  opt==o (default) list of simplices;");
	errline( "-oF<name>  prefix of output files is <name>;");
	errline( "-oN  no main output;");
	errline( "-ov volumes of Voronoi facets");
	errline( "-oh incidence histograms");
}


void echo_command_line(FILE *F, int argc, char **argv) {
	fprintf(F,"%%");
	while (--argc>=0) 
		fprintf(F, "%s%s", *argv++, (argc>0) ? " " : "");
		fprintf(F,"\n");
}

char *output_forms[] = {"vn", "ps", "mp", "cpr", "off"};

out_func *out_funcs[] = {&vlist_out, &ps_out, &mp_out, &cpr_out, &off_out};


int set_out_func(char *s) {

	int i;

	for (i=0;i< sizeof(out_funcs)/(sizeof (out_func*)); i++)
		if (strcmp(s,output_forms[i])==0) return i;
	tell_options();
	return 0;
}

void make_output(simplex *root, void *(*visit_gen)(simplex*, visit_func* visit),
		visit_func* visit,
		out_func* out_funcp,
		FILE *F){

	out_funcp(0,0,F,-1);
	visit(0, out_funcp);
	visit_gen(root, visit);
	out_funcp(0,0,F,1);
	fclose(F);
}



void main(int argc, char **argv) {
    DFILE = stderr;

	long	seed = 0;
	short	shuffle = 0,
		pine = 0,
		output = 1,
		hist = 0,
		vol = 0,
		ahull = 0,
		ofn = 0,
		ifn = 0;
	int	option;
	double	alpha = 0;
	char	ofile[50] = "",
		ifile[50] = "",
		ofilepre[50] = "";
	FILE	*T;
	int	main_out_form=0, alpha_out_form=0;
	simplex *root;
	fg 	*faces_gr;

	mult_up = 1;

	while ((option = getopt(argc, argv, "i:m:rs:do:X:a:Af:")) != EOF) {
		switch (option) {
		case 'm' :
			sscanf(optarg,"%lf",&mult_up);
			DEBEXP(-4,mult_up);
			break;
		case 'r':
			shuffle = 1;
			break;
		case 's':
			seed = atol(optarg);
			shuffle = 1;
			break;
		case 'd' :
			vd = 1;
			break;
		case 'i' :
			strcpy(ifile, optarg);
			break;
		case 'X' : 
			DFILE = efopen(optarg, "w");
			break;
		case 'f' :
			main_out_form = set_out_func(optarg);
			break;
		case 'A':
			vd = ahull = 1;
			break;
		case 'a' :
			vd = ahull = 1;
			switch(optarg[0]) {
				case 'a': sscanf(optarg+1,"%lf",&alpha); break;
				case 'f':alpha_out_form=set_out_func(optarg+1);
					break;
				case '\0': break;
				default: tell_options();
			 }
			break;
		case 'o': switch (optarg[0]) {
				case 'o': output=1; break;
				case 'N': output=0; break;
				case 'v': vd = vol = 1; break;
				case 'h': hist = 1; break;
				case 'F': strcpy(ofile, optarg+1); break;
				default: errline("illegal output option");
				exit(1);
			}
			break;
		default :
			tell_options();
			exit(1);
		}
	}

	ifn = (strlen(ifile)!=0); 
	INFILE = ifn ? efopen(ifile, "r") : stdin;
	fprintf(DFILE, "reading from %s\n", ifn ? ifile : "stdin");

	ofn = (strlen(ofile)!=0);

	strcpy(ofilepre, ofn ? ofile : (ifn ? ifile : "hout") );

	if (output) {
		if (ofn && main_out_form > 0) {
			strcat(ofile, "."); 
			strcat(ofile, output_forms[main_out_form]);
		}
		OUTFILE = ofn ? efopen(ofile, "w") : stdout;
		fprintf(DFILE, "main output to %s\n", ofn ? ofile : "stdout");
	} else fprintf(DFILE, "no main output\n");

	TFILE = efopen(tmpnam(tmpfilenam), "w");

	read_next_site(-1);
/*	fprintf(DFILE,"dim=%d\n",dim);fflush(DFILE); */
	if (dim > MAXDIM) panic("dimension bound MAXDIM exceeded"); 

	point_size = site_size = sizeof(Coord)*dim;

	if (shuffle) {
		fprintf(DFILE, "reading sites...");
		for (num_sites=0; read_next_site(num_sites); num_sites++);
		fprintf(DFILE,"done; num_sites=%d\n", num_sites);fflush(DFILE);
		fprintf(DFILE,"shuffling...");
		init_rand(seed);
		make_shuffle();
		shuf = &shufflef;
		get_site_n = get_site_offline;
	} else {
		fprintf(DFILE,"not shuffling\n");
		shuf = &noshuffle;
		get_site_n = read_next_site;
	}

	fprintf(DFILE, "finding %s...",
		vd ? "Delaunay triangulation" : "convex hull");

	root = build_convex_hull(get_next_site, site_numm, dim, vd);

	fclose(TFILE);
	fprintf(DFILE, "done\n"); fflush(DFILE);

	if (output) {
		out_func* mof = out_funcs[main_out_form];
		visit_func *pr = facets_print;
		
		if (main_out_form==0) echo_command_line(OUTFILE,argc,argv);
		else if (vd) pr = ridges_print;
		
		make_output(root, visit_hull, pr, mof, OUTFILE);
	}

	if (ahull) {
		char ahname[50];
		out_func* aof = out_funcs[alpha_out_form];
		if (alpha_out_form==0)
			sprintf(ahname, "%s-alf", ofilepre);
		else	sprintf(ahname, "%s-alf.%s", ofilepre,
				output_forms[alpha_out_form]);

		T = efopen(ahname,"w");
		fprintf(DFILE, "finding alpha shape; output to %s\n", ahname);
		fflush(DFILE);
		if (alpha==0) alpha=find_alpha(root);
		DEBEXP(-10, alpha)
		if (alpha_out_form==0) echo_command_line(T,argc,argv);	
		alph_test(0,0,&alpha);
		make_output(root, visit_outside_ashape, afacets_print, aof, T);
	}


	if (vol) {
		char volnam[50];
		sprintf(volnam, "%s-vol", ofilepre);
		fprintf(DFILE, "finding volumes; output to %s\n", volnam);
		fflush(DFILE);
		find_volumes(faces_gr=build_fg(root), efopen(volnam,"w"));
	}

	if (vd && hist) {
		char hnam[50];
		sprintf(hnam, "%s-hist", ofilepre);
		fprintf(DFILE,"finding incidence histograms; output to %s\n", hnam);
		fflush(DFILE);
		if (!faces_gr) faces_gr = build_fg(root);
		print_hist_fg(root, faces_gr, efopen(hnam,"w"));
	}

	free_hull_storage();

	exit(0);
}
