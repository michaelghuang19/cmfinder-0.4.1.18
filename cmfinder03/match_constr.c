#include "squid.h"
#include <stdlib.h>
#include <string.h>
#include "match_constr.h"
#include "global.h"

#define MAXSHIFT 50

MatchPtr** ReadMatchConstr(char* matchfile, int nseq, SQINFO* sqinfo, double threshold, Range* range)
{
	FILE *fin;
	char buffer[MAXLINE]; 
	char name1[MAXLINE];
	char name2[MAXLINE];
	int  start1;
	int  start2;
	int  stop1;
	int  stop2;
	int  i, j;
	double e_val;  
	char* last_name1=NULL;  
	char* last_name2=NULL;  

	MatchPtr** mat;  

	mat = (MatchPtr**) MallocOrDie (nseq * sizeof(MatchPtr*));
	for (j = 0; j < nseq; j++){    
		mat[j] = (MatchPtr*) MallocOrDie ( nseq * sizeof(MatchPtr));
		memset(mat[j], 0, nseq*sizeof(MatchPtr));    
	}  

	/* Count the number of seqs */
	if ( (fin = fopen(matchfile, "r") ) == NULL) 
		Die("Fail to read file %s", matchfile);

	i=-1;
	j=-1;
	while(fgets(buffer,MAXLINE, fin) > 0) {      
		//Read the file
		if (strlen(buffer) < 5) continue;    
		memset(name1, 0, sizeof(char)*100);
		memset(name2, 0, sizeof(char)*100);

		sscanf(buffer, "%s :\t %d - %d\t %lf", name1, &start1, &stop1, &e_val);

		if (fgets(buffer,MAXLINE, fin) > 0) {      
			sscanf(buffer, "%s :\t %d - %d", name2, &start2, &stop2);
		}
		else break;
		if (e_val > threshold) continue;    

		//Find corresponding seq index
		if (last_name1 == NULL  || strcmp(last_name1, name1)) {      
			for(i=0; i < nseq; i++) {
				if(strncmp(sqinfo[i].name, name1,strlen(sqinfo[i].name))==0) 
					break;	
			}
			if (i== nseq) Warn("Unknow sequence name %s in file %s", name1, matchfile);
		}
		if (last_name2 == NULL  || strcmp(last_name2, name2)) {      
			for(j=0; j < nseq; j++) {
				if(strncmp(sqinfo[j].name, name2, strlen(sqinfo[j].name))==0) 
					break;	
			}      
			if (j== nseq) Warn("Unknow sequence name %s in file %s", name2, matchfile);      
		}
		if (i==nseq || j==nseq) continue;    
		if (range){
			if (range[i].start != -1 && (stop1 < range[i].start || start1 > range[i].stop)) continue;
			if (range[j].start != -1 && (stop2 < range[j].start || start2 > range[j].stop)) continue;	
		}

		//Insert to the matrix
		{
			MatchPtr temp;      
			MatchPtr prev;
			MatchPtr curr;
			temp = (MatchPtr) malloc(sizeof(struct MatchNode));
			temp->node.seq_id1 = i;
			temp->node.seq_id2 = j;
			temp->node.start1 = start1;
			temp->node.start2 = start2;
			temp->node.stop1 =  stop1;
			temp->node.stop2 =  stop2;
			temp->node.e_val =  e_val;
			temp->node.valid = 1;

			prev = NULL;
			curr = mat[i][j];      
			while(curr) {
				if (curr->node.start1 > start1) break;	
				//Conflicting anchoring points
				if (!CheckMatch(start1, stop1,  start2, stop2, curr, threshold)){
					if (temp->node.e_val < 10 * curr->node.e_val){
						curr->node.valid =0;
					}
					else if (temp->node.e_val < 10 * curr->node.e_val){
						temp->node.valid =0;
					}
					else{
						curr->node.valid =0;
						temp->node.valid =0;	    
					}	  
				}
				prev = curr;	
				curr = curr->next;
			}
			if (prev == NULL)
				mat[i][j] = temp;
			temp->next = curr;          
			printf("match %d (%d %d) - %d (%d %d)\n", i, start1,stop1, j, start2, stop2);      
		}        
	}  
	return mat;  
}

int CheckMatch(int start1, int stop1, int start2, int stop2, MatchPtr mc, double threshold)
{ 
	MatchPtr curr;
	MatchRecord* r;  
	for (curr=mc; curr!= NULL; curr = curr->next){
		r = &(curr->node);    
		if (r->e_val > threshold || !curr->node.valid) continue;    
		if (stop1 < r->start1 && stop2 < r->start2) {      
			break;      
		}
		else if (start1 > r->stop1 && start2 > r->stop2) {          
			continue;      
		}
		else if(stop2 > r->stop2 && stop1 < r->start1) {       
			return 0;      
		}
		else if(stop1 > r->stop1 && stop2 < r->start2){            
			return 0;      
		}
		else if (start1 < r->start1 && start2 > r->stop2){
			return 0;      
		}
		else if (start2 < r->start2 && start1 > r->stop1){
			return 0;      
		}    
		else {
			int ostart1, ostop1, ostart2, ostop2;
			ostart1 =  start1 > r->start1 ? start1 : r->start1;
			ostart2 =  start2 > r->start2 ? start2 : r->start2;
			ostop1  =  stop1  < r->stop1 ?  stop1  : r->stop1;
			ostop2  =  stop2  < r->stop2 ?  stop2  : r->stop2;
			if (abs(ostart1 - ostart2 - (r->start1 - r->start2)) > MAXSHIFT)
				return 0;
			if (abs(ostop1 -  ostop2 -  (r->stop1 - r->stop2)) > MAXSHIFT)
				return 0;
		}    
	}  
	return 1;        
}



double CheckConserve(int seq_id, int start, int stop, MatchPtr mc)
{ 
	MatchPtr curr = mc;
	MatchRecord* r;    
	int     mstart=0, mstop=0;
	double  max = 0;
	while(curr){
		int ostart;
		int ostop;
		double olap;
		r = &(curr->node);
		if (seq_id == r->seq_id1) {
			mstart = r->start1;
			mstop = r->stop1;
		}
		else if (seq_id == r->seq_id2) {
			mstart = r->start2;
			mstop = r->stop2;
		}
		else{
			return max;
		}
		ostart =  mstart > start ? mstart : start;
		ostop  =  mstop  < stop ? mstop : stop;
		olap= ((double)(ostop - ostart))/ (stop - start);
		if (olap > max){
			//printf("olap %.2f %d / %d\n", olap, ostop - ostart, stop - start);
			max = olap;
		}
		curr = curr->next;    
	}  
	return max;
}
