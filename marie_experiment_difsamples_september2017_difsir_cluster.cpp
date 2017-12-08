////////////////// RANGE EXPANSION MODEL  //////////////////////
////////////////// FRANCISCO ENCINAS VISO ///////////////////////////////////////////////////

//Test for the genome architecture amd the migration/colonization function of the range expansion model


#include <iostream>
#include <cstdlib>
#include <vector>
#include <time.h>
#include <math.h>  
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>
#include <iomanip> 
#include <algorithm>
#include <boost/format.hpp>
#include <time.h>
#include <sstream>




using namespace std;
using boost::format;


struct indiv{ //properties of individual plants
        vector<int> genome;
        int patch_x;
        int patch_y;                       // Patch (or site) in the landscape
        int chosen;                   //for when an adult survives to the next year; to make sure the same plant doesn't survive twice. 10(chosen, -10(not chosen) / Also to stablish which individuals survive in a local extinction
        double fitness;
        int id;
        int S_allele1(void){return genome[(genome.size()/4)-1];}
        int S_allele2(void){return genome[3*(genome.size()/4)-1];}
};



int random_gen_selected(gsl_rng *r){
	int x;
	if(gsl_rng_uniform(r)>0.5) x = 21;
	else  x = 23;
	return x;
}


int random_S_locus(gsl_rng *r, int S_alleles){
	int x = gsl_rng_uniform_int(r, S_alleles) + 100;
	return x;
}



void mutate_poisson_selected(gsl_rng *r,vector<int> &adn, double lambda) {
	int x = gsl_ran_poisson(r, lambda); // number of selected sites
	//cout << "Number of mutations " << x << " lambda "  << lambda << endl;
	int selected;
	double muback = 0.01;
	if(x>0){
		for(int i=0; i<x; i++){
		    selected = gsl_rng_uniform_int(r,adn.size()-1); //Selected sites
			if(adn[selected]==21) adn[selected] = 23;
			else if(adn[selected]==23){
				if(muback > gsl_rng_uniform(r))  adn[selected] = 21;
			}

		}    
	}
}	



void mutate_S_locus(gsl_rng *r,vector<int> &adn, double mu_slocus, vector<int> &Stype, int flag) {
	int L = adn.size(), diffs;
	int S_allele = adn[(L/2)-1];
	if(mu_slocus > gsl_rng_uniform(r) && S_allele>0){
		if(flag == 1){
		     adn[(L/2)-1] = -100;//S_allele;
		     //cout << " ------------------------------------------------------------------------------------> MUTATION S LOCUS " << adn[(L/2)-1] << endl;
		}
		//Update S alleles list
		diffs=0;
		for(vector<int>::iterator iter=Stype.begin(); iter!=Stype.end(); iter++){
			if(adn[(L/2)-1] == *iter)  diffs++;
		}
		if(diffs == 0){
			Stype.push_back(adn[(L/2)-1]);
			//cout << "Mutation S locus " << adn[(L/2)-1] << endl;
		}
	}
}	




void mutations(gsl_rng *r, vector<int> &adn,  double lambda_sel, double mu_slocus, vector<int> &Stype, int flag){
	mutate_poisson_selected(r, adn, lambda_sel);
	mutate_S_locus(r,adn, mu_slocus, Stype, flag);
}

vector<int> genome_fun(gsl_rng *r,int selgenes, int S_alleles){
	vector<int> genome_vec;
	vector<int> genome_vec2;
	int L = selgenes;

	for(int gen=0; gen<L; gen++){
		

		if(gen==(L/2)-1)
		     genome_vec.push_back(random_S_locus(r,S_alleles));
		else
		     genome_vec.push_back(21); //random_gen_selected(r)
		
	}
	
	for(int gen=0; gen<L; gen++){
		
        if(gen==(L/2)-1)
		     genome_vec2.push_back(random_S_locus(r,S_alleles));
		else
		     genome_vec2.push_back(21); //random_gen_selected(r)
	}

	genome_vec.insert(genome_vec.end(), genome_vec2.begin(), genome_vec2.end());
	return genome_vec;
}

void show_vector(vector<int> &vec){
	cout << "Vector contains: ";
    for (vector<int>::iterator it=vec.begin(); it!=vec.end(); ++it)
         cout << ' ' << *it;
    cout << endl;
}


void show_sc(vector<indiv> *vec){
	int sc=0,n=0;
	cout << "Self-compatibles: ";
    for (vector<indiv>::iterator it=vec->begin(); it!=vec->end(); it++){
         if(it->S_allele1() < 0){
			  sc++;
			 // cout << "S allele 1 " << it->S_allele1() << endl;
			  //cout << "ind " << n << endl;
		 }
         if(it->S_allele2() < 0){
			  sc++;
			   //cout << "S allele 2 " << it->S_allele2() << endl;
		 }
		 n++;    
	 }
	
    cout << sc << endl;
}

void number_genes(vector<int> &vec){
	int sel=0, neutral=0, S=0;
	cout << "Genes: ";
    for (vector<int>::iterator it=vec.begin(); it!=vec.end(); ++it){
		if(*it==21 || *it==23)   sel++;
		else if(*it>=100 && *it<500)         S++;
	}
    cout <<  "Sel " << sel << " N " << neutral << " S " << S << endl;
}


void refresh_chosen(vector<indiv> *adults){
	for(vector<indiv>::iterator iter=adults->begin(); iter!=adults->end(); iter++){
		iter->chosen=-10;
	}
}

/// New gamete production and mating function
bool sorted(int i, int j){ return (i<j); }


vector<int> choose_dna(gsl_rng *r, vector<int> adn1, vector<int> adn2){
	vector<int> adn;
	if(gsl_rng_uniform(r)>0.5){
		adn = adn1;
	}else{
		adn = adn2;
	}
	return adn;
}



vector<int> create_haplotype(gsl_rng *r, vector<int> genome, double lambda, int genes){
	int i=0,x;
	vector<int> haplotype;
	vector<int> adn;
	vector<int> selected;
	vector<int> pool;
	/// First chromosome:
	vector<int>::const_iterator first_chrom1 = genome.begin();
    vector<int>::const_iterator last_chrom1 = genome.begin() + genome.size()/2;
    vector<int> adn1(first_chrom1,last_chrom1); // First chromosome
    /// Second chromosome:
    vector<int>::const_iterator first_chrom2 = genome.begin() + genome.size()/2;
    vector<int>::const_iterator last_chrom2 = genome.end();
    vector<int> adn2(first_chrom2,last_chrom2); // Second chromosome
	// Pick cross-over points:
    // number of selected crossover points
	do{
		x = gsl_ran_poisson(r, lambda); 
	}while(x>genes);
	if(x>0){
		 // Fisher-Yates algorithm
		 for(int j=0; j<genes; j++)
	          pool.push_back(j);
	     random_shuffle(pool.begin(), pool.end()); // Create shuffled array 
         selected.insert(selected.begin(), pool.begin(), pool.begin()+x); // Select x sites
	     sort(selected.begin(), selected.end(), sorted);
         for(vector<int>::iterator iter=selected.begin(); iter!=selected.end(); iter++){ // Recombination process
		        if(i%2==0) adn = adn1;
		        else       adn = adn2;
		        if(i==0){ 
					 haplotype.insert(haplotype.end(), adn.begin(), adn.begin() + *iter);
				}
				else if(i>0){
		             haplotype.insert(haplotype.end(), adn.begin() + *(iter-1), adn.begin() + *iter);
				}
		        i++;
		  }
		  if(i%2==0) adn = adn1;
		  else       adn = adn2;
		  haplotype.insert(haplotype.end(), adn.begin() + selected[x-1], adn.end());
	}else{
		  haplotype = choose_dna(r,adn1,adn2);
	}
	return haplotype;
}



vector<int> create_haplotype_freerec(gsl_rng *r, vector<int> genome, double lambda, int genes){
	vector<int> haplotype;
	/// First chromosome:
	vector<int>::const_iterator first_chrom1 = genome.begin();
    vector<int>::const_iterator last_chrom1 = genome.begin() + genome.size()/2;
    vector<int> adn1(first_chrom1,last_chrom1); // First chromosome
    /// Second chromosome:
    vector<int>::const_iterator first_chrom2 = genome.begin() + genome.size()/2;
    vector<int>::const_iterator last_chrom2 = genome.end();
    vector<int> adn2(first_chrom2,last_chrom2); // Second chromosome
	// Pick cross-over points:
    for(std::vector<int>::size_type i = 0; i != adn1.size(); i++) {
		if(gsl_rng_uniform(r) > lambda) haplotype.push_back(adn1[i]);
		else haplotype.push_back(adn2[i]);
		
	}
    
	return haplotype;
}



double calculatefitness(gsl_rng *r,vector<int> &gen, int genes, double mu, double h){
	int i;
	int L=genes, allele1, allele2;
	double x1,x2,x3,s = mu;
	vector<double> fitness;
	int count=0, countna=0, homoB=0, homob=0, hetero=0;
	for(i=0; i<genes; i++){
        allele1 = gen[i];
        allele2 = gen[i+L];
  
		if (allele1==21 || allele1 ==23 || allele2==21 || allele2==23){
			//cout << "Allele 1 " << allele1  <<  " Allele 2 " << allele2 << endl;
			if(allele1+allele2 == 42){
			    x1 = 1.0; 
			    homoB++;
			    fitness.push_back(x1);
			    //cout << "Non-del " << fitness[countna] << " " << x1 << endl;
		    }else if(allele2+allele1==46){
			    x2 = 1.0 - s;//*(mat_b + (countna)*3 + j);
			    homob++;
			    fitness.push_back(x2);
			    //cout << "Deleterious " << fitness[countna] << " " << allele1 << " " << allele2 << endl;
		    }else if(allele1+allele2==44){
			    x3 =  1.0 - h*s;//*(mat_b + (countna)*3 + j);
			    hetero++;
			    fitness.push_back(x3);
			    //cout << "Hetero "  << fitness[countna] << " " << allele1 << " " << allele2 << endl;
			}
			countna++;
		}
		count++;
			
	}
	//cout << "Fitness size vector " << fitness.size() << endl;
	// Product of fitness effects:
	double fitness_product = fitness[0];
	for(unsigned int j=1;j<fitness.size();j++){
		fitness_product *= fitness[j];
	}
	//cout << "Homo non-del: " << homoB  << " Hetero: "  << hetero << " Homo del: " << homob << " Fitness: " << fitness_product << endl;
	return fitness_product;
}



// Population initialization for fully populated space (metapopulation) or range expansion case 
void initializeplants(gsl_rng *r, vector<indiv> *p, int selgenes, int S_alleles, int subpop, int patchesx, double mu, double dominance, int *space)
{ 
	    int ii=0;
	    for(int i=0; i<patchesx; i++){
			if(i < 2){
				    for(int inds=0; inds<subpop; inds++){
						if(i == 0){
					        (*p)[inds + ii].genome = genome_fun(r, selgenes, S_alleles); // 100 % SI pop
						}else if(i == 1){
							(*p)[inds + ii].genome = genome_fun(r, selgenes, S_alleles); // 100 % SC pop
							(*p)[inds + ii].genome[((*p)[inds + ii].genome.size()/4)-1] = -100;
							(*p)[inds + ii].genome[3*((*p)[inds + ii].genome.size()/4)-1] = -100;
						}   
                       (*p)[inds + ii].patch_x = i;
                       (*p)[inds + ii].chosen = -10;
                       (*p)[inds + ii].fitness = calculatefitness(r, (*p)[inds+ii].genome , (*p)[inds+ii].genome.size()/2, mu, dominance); 
                       //(*p)[inds + ii].S_allele = s_alleles_fun( (*p)[inds+ii].genome );
                       //cout << "S alleles: " << (*p)[inds + ii].S_allele1() << " " << (*p)[inds + ii].S_allele2() << endl;
				    }
				space[i]=1;
				ii=ii+subpop;
			}
		}
							
}

void show_vector_pointer(vector<indiv> *vec){
    cout << "Vector length: " << vec->size() << endl;
    for (vector<indiv>::iterator it=vec->begin(); it!=vec->end(); ++it){
	    cout << "Vector contains: ";
		for(vector<int>::iterator gen=it->genome.begin(); gen!=it->genome.end(); gen++)
            cout << " " << *gen;
        cout << endl;
    }
}

void glass(vector<int> &vec){

    for (vector<int>::iterator it=vec.begin(); it!=vec.end(); ++it){
		if(*it == 0){
			show_vector(vec);
			cout << "FLAG 1" << endl;
			exit(EXIT_FAILURE);
		}else if(*it < -100){
			show_vector(vec);
			cout << "FLAG 2" << endl;
			exit(EXIT_FAILURE);
		}else if(*it > 100000){
			show_vector(vec);
			cout << "FLAG 3" << endl;
			exit(EXIT_FAILURE);
		}
	}

}


vector<int> fertilization2(vector<int> ovule, vector<int> pollen){
	vector<int> seed;
    seed.reserve(ovule.size() + pollen.size());
	seed.insert(seed.end(), ovule.begin(), ovule.end());
	seed.insert(seed.end(), pollen.begin(), pollen.end());
    return seed;
}



void build_subpop_vector_pointer2D(int patchx, vector<indiv> *vectors, vector<indiv> *subpop){
	int popsize = (*vectors).size();
	for(int i=0; i<popsize; i++){
			if((*vectors)[i].patch_x == patchx){
				subpop->push_back((*vectors)[i]);
			}
	}
}

void build_subpop_vector_pointer2D_indexes(int patchx, vector<indiv> *vectors, vector<int> *subpop){
	int popsize = (*vectors).size();
	for(int i=0; i<popsize; i++){
			if((*vectors)[i].patch_x == patchx){
				subpop->push_back(i);
			}
	}
}



void build_subpop_indexes(int patchx, vector<indiv> *seeds, vector<int> *indexes){
	int i=0;
	for(vector<indiv>::iterator it=seeds->begin(); it!=seeds->end(); it++){
		if(it->patch_x == patchx && it->chosen==-10){
			indexes->push_back(i);
		}
		i++;
	}
}


bool rankfunc (indiv i, indiv j) { return (i.fitness>j.fitness); }

double maxfit_indiv(vector<indiv> *subpop){
	sort (subpop->begin(), subpop->end(), rankfunc);
	return (*subpop)[0].fitness;
}

vector<int> find_mate_GSI(gsl_rng *r,vector<indiv> *subpop, vector<int> &ovule_S_allele, double rec, double lambda_sel, double mu_slocus, int genes, vector<int> &Stype, int flag){
	int k, L, fail = 0;
	double maxfit,siring_prob;
	vector<int> pollen_haplotype, temp;
	int pop_size = subpop->size();
	k = gsl_rng_uniform_int(r,pop_size);
	pollen_haplotype = create_haplotype_freerec(r, (*subpop)[k].genome, rec, genes);
	L = pollen_haplotype.size(); 
	//cout << " First Ovule S: " << ovule_S_allele[0] << " " << ovule_S_allele[1] << "  Pollen S: " << pollen_haplotype[(L/2)-1]  << endl;
	maxfit = maxfit_indiv(subpop);
	for(int counter = 0; counter < pop_size; counter++){
	/// Select SC pollen grains twice more than SI pollen	
	    if(pollen_haplotype[(L/2)-1] < 0){
	       do{
			  siring_prob = gsl_rng_uniform(r);
		   }while(siring_prob < 0.5);
	    }else
	       siring_prob = gsl_rng_uniform(r); 
     /// Functional S allele
			//cout << "Ovule S: " << ovule_S_allele[0] << " " << ovule_S_allele[1] << "  Pollen S: " << pollen_haplotype[(L/2)-1] << endl;
		    if(pollen_haplotype[(L/2)-1] != ovule_S_allele[0] &&  pollen_haplotype[(L/2)-1] != ovule_S_allele[1] && gsl_rng_uniform(r) < (*subpop)[k].fitness/maxfit && gsl_rng_uniform(r) < siring_prob){
				   temp = pollen_haplotype;
				   fail = 1;
			       break;
		    }
		k = gsl_rng_uniform_int(r,pop_size);
		pollen_haplotype = create_haplotype_freerec(r, (*subpop)[k].genome, rec, genes);
	}
	
	if(fail == 0){
		temp.push_back(-1);
		//cout << "NO compatible MATES !!" << endl;
	}else mutations(r, temp, lambda_sel, mu_slocus, Stype,flag);
	
	return temp;
}


vector<int> find_mate_GSI_SC(gsl_rng *r,vector<indiv> *subpop, vector<int> &ovule_S_allele, double rec, double lambda_sel, double mu_slocus, int genes, vector<int> &Stype, int flag){
	int counter=0,k, L;
	double maxfit, siring_prob;
	vector<int> pollen_haplotype, temp;
	int pop_size = subpop->size();
	k = gsl_rng_uniform_int(r,pop_size);
	pollen_haplotype = create_haplotype_freerec(r, (*subpop)[k].genome, rec, genes);
	L = pollen_haplotype.size(); 
	//cout << " First Ovule S: " << ovule_S_allele[0] << " " << ovule_S_allele[1] << "  Pollen S: " << pollen_haplotype[(L/2)-1]  << endl;
	maxfit = maxfit_indiv(subpop);	
	while(counter<pop_size){
	  /// Select SC pollen grains twice more than SI pollen	
	    if(pollen_haplotype[(L/2)-1] < 0){
	       do{
			  siring_prob = gsl_rng_uniform(r);
		   }while(siring_prob < 0.5);
	    }else
	       siring_prob = gsl_rng_uniform(r); 
     /// Functional S allele
		//cout << "Ovule S: " << ovule_S_allele[0] << " " << ovule_S_allele[1] << "  Pollen S: " << pollen_haplotype[(L/2)-1] << endl;
		if(pollen_haplotype[(L/2)-1] > 0){	
		    if(pollen_haplotype[(L/2)-1] != ovule_S_allele[0] &&  pollen_haplotype[(L/2)-1] != ovule_S_allele[1] && gsl_rng_uniform(r) < (*subpop)[k].fitness/maxfit && gsl_rng_uniform(r) < siring_prob){
				   temp = pollen_haplotype;
			       break;
		    }
        }else{ // SC pollen haplotype - functional S allele of ovule has to be different than pollen haplotype
			if(ovule_S_allele[0] > 0 &&  pollen_haplotype[(L/2)-1] != ovule_S_allele[0] && gsl_rng_uniform(r) < (*subpop)[k].fitness/maxfit && gsl_rng_uniform(r) < siring_prob){ 
				   temp = pollen_haplotype;
			       break;
		    }else if(ovule_S_allele[1] > 0 &&  pollen_haplotype[(L/2)-1] != ovule_S_allele[1] && gsl_rng_uniform(r) < (*subpop)[k].fitness/maxfit && gsl_rng_uniform(r) < siring_prob){
				   temp = pollen_haplotype;
			       break;
			}
		}
		k = gsl_rng_uniform_int(r,pop_size);
		pollen_haplotype = create_haplotype_freerec(r, (*subpop)[k].genome, rec, genes);
		counter++;
	}
	
	//cout << "Counter " << counter << endl;
	if(counter >= pop_size)   temp.push_back(-1);
	else mutations(r, temp, lambda_sel, mu_slocus, Stype,flag);
	
	return temp;
}

vector<int> find_pollen(gsl_rng *r,vector<indiv> *subpop, vector<int> &ovule_S_allele, double rec, double lambda_sel, double mu_slocus, int genes, vector<int> &Stype, int flag){
	int chosen_father;
	double maxfit, siring_prob;
	vector<int> pollen_haplotype;
	int pop_size = subpop->size();
    maxfit = maxfit_indiv(subpop);
	chosen_father = gsl_rng_uniform_int(r, pop_size);
	do{
	   chosen_father = gsl_rng_uniform_int(r, pop_size);
	   if((*subpop)[chosen_father].S_allele1() == -100 || (*subpop)[chosen_father].S_allele2() == -100){
	       
	       do{
			  siring_prob = gsl_rng_uniform(r);
		   }while(siring_prob < 0.5);
	   }else
	       siring_prob = gsl_rng_uniform(r);    
    }while( gsl_rng_uniform(r) > (*subpop)[chosen_father].fitness/maxfit && gsl_rng_uniform(r) > siring_prob);
	pollen_haplotype = create_haplotype_freerec(r, (*subpop)[chosen_father].genome, rec, genes);
	mutations(r, pollen_haplotype,lambda_sel, mu_slocus, Stype, flag);
	return pollen_haplotype;
}






int subpopsize(int patchx, vector<indiv> *pop){
	int inds=0;
	for(vector<indiv>::iterator iter=pop->begin(); iter!=pop->end(); iter++){
		if(iter->patch_x == patchx)   inds++;
	}
	return inds;
}





bool rank (indiv i, indiv j) { return (i.fitness>j.fitness); }



void experiment(gsl_rng *r, vector<indiv> *seeds, int sample1, int sample2, int *spatial_occupancy){
    vector<int> *subpopSI = new vector<int>;
    vector<int> *subpopSC = new vector<int>;
	build_subpop_indexes(0,seeds,subpopSI);
    build_subpop_indexes(1,seeds,subpopSC);
    random_shuffle ( subpopSI->begin(), subpopSI->end());
    for(int i=0; i<sample1; ++i){
		//cout << "From patch: " << (*seeds)[(*subpopSI)[i]].patch_x << " " << "Index: " << (*subpopSI)[i] << " S allele 1 " << (*seeds)[(*subpopSI)[i]].S_allele1() << endl;
		(*seeds)[(*subpopSI)[i]].patch_x = 2;
        (*seeds)[(*subpopSI)[i]].chosen = 10;
	}
	random_shuffle ( subpopSC->begin(), subpopSC->end());
	for(int i=0; i<sample2; ++i){
		//cout << "From patch: " << (*seeds)[(*subpopSC)[i]].patch_x << " " << "Index: " << (*subpopSC)[i] << " S allele 1 " << (*seeds)[(*subpopSC)[i]].S_allele1() << endl;
		(*seeds)[(*subpopSC)[i]].patch_x = 2;
        (*seeds)[(*subpopSC)[i]].chosen = 10;
	}
	spatial_occupancy[2] = 1;
	delete subpopSI;
	delete subpopSC;
}

void checkpopsize(int patchesx, int patchesy, vector<indiv> *seeds, int **spatial_occupancy){
	for(int i=0; i<patchesx; i++){
		for(int j=0; j<patchesy; j++){
			if(spatial_occupancy[i][j] == 1){
				vector<indiv> *subpop = new vector<indiv>;
				build_subpop_vector_pointer2D(i,seeds,subpop);
				//cout << "Subpop size patch: " << i << " " << j << " : " << subpop->size() << endl;
				delete subpop;
			}
		}
	} 
}


double mean_fitness(vector<indiv> *subpop){
	double sumfit=0;
	for(vector<indiv>::iterator iter=subpop->begin(); iter!=subpop->end(); iter++)
		 sumfit += iter->fitness;
    return sumfit/subpop->size();
}


void build_subpop_vector_pointer1D(int patchx, vector<indiv> *vectors, vector<indiv> *subpop){
	int popsize = (*vectors).size();
	for(int i=0; i<popsize; i++){
			if((*vectors)[i].patch_x == patchx){
				subpop->push_back((*vectors)[i]);
			}
	}
}
void gameteprod_mating(gsl_rng *r, vector<indiv> *parents, vector<indiv> *seeds, int genes, int patchesx, int N, double lambda_sel, double mu_slocus, double selection_exp_mu, double dominance, double rec, double selfing_prob, int *space, vector<int> &Stype, int flag, double **p2outputSelfrate, int expa, int k_si, int k_sc, int k_exp){
	int chosen_mother, popsize, self, stop;
	double sumfit, maxfit;
	vector<int> ovule_haplotype, pollen_haplotype;
	pair <vector<int>, int > pollen_pair;
    int ks[] = {k_si,k_sc,k_exp};
	for(int i=0; i<patchesx; i++){

			if(space[i]==1){
				

				   vector<indiv> *subpop = new  vector<indiv>;
				   build_subpop_vector_pointer1D(i,parents,subpop);
				   //meanfit = mean_fitness(subpop);
	               stop = 0, self = 0;				
				   vector<int> S_alleles_mother;
				   //cout << "subpop size: " << subpop->size() << endl;			
				   //cout << "mean fitness: " << meanfit << endl;	
				   maxfit = maxfit_indiv(subpop);	            
                   popsize = subpop->size();
				   sumfit = 0;
			       for(vector<indiv>::iterator iter=subpop->begin(); iter!=subpop->end(); iter++)
		              sumfit += iter->fitness;
				   for(int n=0; n < ks[i]; n++){
					  //cout << "Counter for loop: " << n << endl;
					  chosen_mother = gsl_rng_uniform_int(r, popsize);
					  do{
					    	chosen_mother = gsl_rng_uniform_int(r, popsize);
					  }while( gsl_rng_uniform(r) > ((*subpop)[chosen_mother].fitness)/maxfit);
					  //cout << "Fitness " << ((*subpop)[chosen_mother].fitness) << endl;
					
					  S_alleles_mother.push_back((*subpop)[chosen_mother].S_allele1());
					  S_alleles_mother.push_back((*subpop)[chosen_mother].S_allele2());
					  ovule_haplotype = create_haplotype_freerec(r, (*subpop)[chosen_mother].genome, rec, genes);
					  glass(ovule_haplotype);
					  //L = ovule_haplotype.size();
					  if (S_alleles_mother[0] < 0 || S_alleles_mother[1] < 0){
						
						//selfing
					    	if(selfing_prob > gsl_rng_uniform(r)){
						    	pollen_haplotype = create_haplotype_freerec(r, (*subpop)[chosen_mother].genome, rec, genes);
							    glass(pollen_haplotype);
							    mutations(r, pollen_haplotype, lambda_sel, mu_slocus, Stype, flag);
							    self++;
						    }else // outcrossing
						        if((S_alleles_mother[0] < 0 && S_alleles_mother[1] > 0)  || (S_alleles_mother[0] > 0 && S_alleles_mother[1] < 0)){
							       pollen_haplotype = find_mate_GSI_SC(r, subpop, S_alleles_mother, rec, lambda_sel, mu_slocus, genes, Stype, flag);
						        }else if(S_alleles_mother[0] < 0 && S_alleles_mother[1] < 0){
								   pollen_haplotype = find_pollen(r, subpop, S_alleles_mother, rec, lambda_sel, mu_slocus, genes, Stype, flag);
							    }
						
				 	  }else // GSI - outcrossing
					  {
						    pollen_haplotype = find_mate_GSI(r, subpop, S_alleles_mother, rec, lambda_sel, mu_slocus, genes, Stype, flag);
					  }
						
					
					  vector<int>().swap(S_alleles_mother);
					  // If there aren't compatible mates (Allee effect) or all with fitness ~= 0
					  if(pollen_haplotype[0] == -1)   stop++;	

					  if(stop >= N)   continue;
					  if(pollen_haplotype[0] != -1){
					     indiv seed;
					     seed.genome = fertilization2(ovule_haplotype, pollen_haplotype);
					     seed.patch_x = i;
					     seed.chosen = -10;
					     seed.fitness = calculatefitness(r, seed.genome, seed.genome.size()/2, selection_exp_mu, dominance);
					     seeds->push_back(seed);
				      } 
                      vector<int>().swap(pollen_haplotype);
                      vector<int>().swap(ovule_haplotype);
				  //}
				}

				//if(self>0)  cout << "-------------------- >>> SELFING RATE: " << (double) self/ks[i] << endl;
				if(expa>0)   p2outputSelfrate[expa][i] = (double) self/ks[i];

				
			    if(subpopsize(i,seeds) > 0) space[i] = 1;
			    else space[i] = 0;
				
			    delete subpop;
			   
			}
	}
	
}


void show_vector_pointer(vector<double> *vec){
	cout << "Vector contains: ";
    for (vector<double>::iterator it=vec->begin(); it!=vec->end(); ++it)
         cout << ' ' << *it;
    cout << endl;
}


double variance(double mean, double data[], int size){
	double sqsum=0;
	int n=0;
	while(n<size){
		sqsum += (data[n] - mean)*(data[n]-mean);
		n++;
	}
	return sqsum/size;
	
}


double inbreeding_depression_within_deme(gsl_rng *r, vector<indiv> *subpop, double selection_exp_mu, double rec, double dominance){
	vector<int> pollen_haplotype, ovule_haplotype;
	int chosen_father, chosen_mother, chosen_ind, N = subpop->size(), n=0, L = (*subpop)[0].genome.size()/2;
	double sumrand=0.0, suminb=0.0; 
	vector<double> offspring(N);
    vector<double> inbred(N);
    vector<int> seed_genomer, seed_genomei;
	while(n<N){
		// Random mating
		chosen_mother = gsl_rng_uniform_int(r,subpop->size());
		ovule_haplotype = create_haplotype_freerec(r, (*subpop)[chosen_mother].genome, rec, L);
		chosen_father = gsl_rng_uniform_int(r,subpop->size());
		pollen_haplotype = create_haplotype_freerec(r, (*subpop)[chosen_father].genome, rec, L);
		seed_genomer = fertilization2(ovule_haplotype, pollen_haplotype);
        offspring[n] = calculatefitness(r,seed_genomer, L, selection_exp_mu, dominance);
		// Inbreeding
		chosen_ind = gsl_rng_uniform_int(r,subpop->size());
		ovule_haplotype = create_haplotype_freerec(r, (*subpop)[chosen_ind].genome, rec, L);
		pollen_haplotype = create_haplotype_freerec(r, (*subpop)[chosen_ind].genome, rec, L);
		seed_genomei = fertilization2(ovule_haplotype, pollen_haplotype);
		inbred[n] = calculatefitness(r, seed_genomei, L, selection_exp_mu, dominance);
		n++;
	}
	//Calculate mean fitness of offspring and inbred
	for(vector<double>::iterator iter=offspring.begin(); iter!= offspring.end(); iter++){
	      sumrand += *iter;
	}
    for(vector<double>::iterator iter=inbred.begin(); iter!= inbred.end(); iter++){
	      suminb += *iter;
	} 
	return 1.0 - ( (suminb/N)/(sumrand/N) );
}



//---------> Create function to Calculate mate availability based on counting the number of compatible pollen donors for each plant available in the sample of the same population (Glemin et al 2008)
double mate_availability(gsl_rng *r, vector<indiv> *subpop)
{
	// Loop through all ovules in the subpop:
	int  N = subpop->size(), x, n=0, chosen_mother, chosen_father;
	double compat=0;
	vector<int> S_alleles_mother, S_alleles_father;
	while(n<N){
		chosen_mother = gsl_rng_uniform_int(r,subpop->size());
		S_alleles_mother.push_back((*subpop)[chosen_mother].S_allele1());
		S_alleles_mother.push_back((*subpop)[chosen_mother].S_allele2());
		chosen_father = gsl_rng_uniform_int(r,subpop->size());
	    S_alleles_father.push_back((*subpop)[chosen_father].S_allele1());
	    S_alleles_father.push_back((*subpop)[chosen_father].S_allele2());
	    if(0.5 > gsl_rng_uniform(r)) x = 0;
	    else x = 1;
	    if(S_alleles_mother[0] != S_alleles_father[x] && S_alleles_mother[1] != S_alleles_father[x])   compat++;
	    //if(S_alleles_mother[0] != S_alleles_father[1] && S_alleles_mother[1] != S_alleles_father[1])   comp++;
	    n++;
	}
	//cout << " ------------------------------------------------------------->>>>>> MATE AVAIL: " << compat << endl;
	return compat/N;
}


/// Function to record S alleles 

void s_alleles_record(int gen, vector<indiv> *population, int patchesx, int patchesy, int **space, double ***** p2_S_alleles, vector<int> &Stype){
	int L = (*population)[0].genome.size();
	unsigned int s_alleles =200;

	for(int k=0; k<patchesx; k++){
		for(int l=0; l<patchesy; l++){
			if(space[k][l]==1){
				int total_S_allele = 0;
				vector<indiv> *subpop = new vector<indiv>;
				vector<int> *s_alleles_count = new vector<int>;
				int *counts_S = new int[200];
				double *freqs_S = new double[200];
				build_subpop_vector_pointer2D(k,population,subpop);
				
				
				for(unsigned int s_al=0; s_al < s_alleles; s_al++){
					int delet = 0;
					for(std::vector<indiv>::iterator it = subpop->begin(); it != subpop->end(); it++){
						s_alleles_count->push_back(it->genome[(L/4) - 1]);
			            s_alleles_count->push_back(it->genome[(3*L/4) - 1]);
						if(Stype[s_al] ==  it->genome[(L/4) - 1]){
							// Counting deleterious alleles associated to one S allele
							for(int kk=0; kk<50; kk++){ 
							   if(it->genome[(L/4)+kk] == 23) delet++;
							   if(it->genome[(L/4)-kk] == 23) delet++;
						    }

					    }
					    if(Stype[s_al] ==  it->genome[(3*L/4) - 1]){
							// Counting deleterious alleles associated to one S allele
							for(int kk=0; kk<50; kk++){ 
							   if(it->genome[(3*L/4)+kk] == 23) delet++;
							   if(it->genome[(3*L/4)-kk] == 23) delet++;
						    }
					    }					    	            
					}
					//cout << "S allele " << Stype[s_al] << " index S-type " << s_al << " patch " << k << " " << l << " generation " << gen << endl;
					p2_S_alleles[gen][k][l][s_al][1] = (double) delet/(2*subpop->size());  // average number of deleterious mutations
					//cout << "S allele " << Stype[s_al] << " delet " << p2_S_alleles[gen][k][l][s_al][1] << endl;
				}
				for(unsigned int s_al=0; s_al < s_alleles; s_al++){
					counts_S[s_al] = std::count(s_alleles_count->begin(), s_alleles_count->end(), Stype[s_al]);
					total_S_allele += counts_S[s_al];
				}
				for(unsigned int s_al=0; s_al < s_alleles; s_al++){
					freqs_S[s_al] =  (double) counts_S[s_al]/total_S_allele;
					//cout << "Frequency of S allele: " << Stype[s_al] << "  "  << freqs_S[s_al] << endl;
					// Putting output to pointer array 
					p2_S_alleles[gen][k][l][s_al][0] = freqs_S[s_al];
				}
				
			    delete subpop;
			    delete s_alleles_count;
			    delete[] counts_S;
		        delete[] freqs_S;
			}
		}
	}
}



// Calculate allele freqs for all patches
void allele_frequencies_2(gsl_rng *r, vector<indiv> *population, double **stat_year_vec, int interval, int patchesx, int selgenes, int S_alleles, vector<int> &Stype, int year, int *space, double selection_exp_mu, double rec, double dominance, double ***p2output, int gen){
    double sumfit_total=0, sum_inbreeding=0, sum_Ae=0, sum_p=0,sum_p2=0, sum_q=0, freq_p_total; 
    int L, number_of_demes, demes=0, totselfcomp=0;
    vector<double> *sum_freq_S = new vector<double>;
    for(unsigned int s_al=0; s_al < Stype.size(); s_al++)  sum_freq_S->push_back(0.0);
    vector<double> *mean_subpop_allele_size_vec = new vector<double>;
    vector<double> *sum_sd_allele_size_vec = new vector<double>;


	vector<int> selfcomps;
	for(int k=0; k<patchesx; k++){
			if(space[k]==1){
				vector<indiv> *subpop = new vector<indiv>;
				build_subpop_vector_pointer2D(k,population,subpop);
				double allele_b, allele_B, hetero_b=0, deleterious_alleles=0, homo_b =0, homo_B=0, sumb=0, sumB=0, p_average=0, q_average=0, sum_quad_freq_S=0, inbdep;
				int total_S_allele=0, allele1, allele2;
				vector<double> *locus_b = new vector<double>;
				vector<double> *locus_B = new vector<double>;
				vector<int> *s_alleles_count = new vector<int>;
				vector<double> *freqs_S = new vector<double>;
				vector<int> *counts_S = new vector<int>;
				double sumfitness=0.0;
                //cout<< "subpop size: "<< subpop->size() << " patch: " << k << " " << l << endl;
                //if(k==0 && l==10){
					//cout << "Contains :";
					//for(vector<indiv>::iterator it=subpop->begin(); it!=subpop->end(); it++)
					  //   cout << it->S_allele1() << "-" << it->S_allele2() << " ";
					//cout << endl; 
				//}
			    //Calculate allele frequencies for each selected locus and neutral locus
			    L = (*subpop)[0].genome.size();
			    int delet = 0;
			    for(int gene=0; gene<L; gene++){
					allele_b=0;
					allele_B=0;
					// Selected and neutral loci
					for(vector<indiv>::iterator it=subpop->begin(); it!=subpop->end(); it++){
						if(it->genome[gene] == 21)  allele_b++;
						if(it->genome[gene] == 23)  allele_B++;
				    }
				    locus_b->push_back(allele_b/subpop->size());
				    locus_B->push_back(allele_B/subpop->size());
				    
				    delet += allele_B;
				    //cout << "Deleterious alleles per locus: " << delet << endl;
			    }

			    
			    number_genes((*subpop)[0].genome);
			    //cout << "Size genome " << (*subpop)[0].genome.size() << endl;
                int selected_genes = (L/2) - 1;
			    // Calculate number of heterozygotes
			    int selfcomp=0;		
			    for(std::vector<indiv>::iterator it = subpop->begin(); it != subpop->end(); it++){
			        for (int j=0; j < L/2; j++){
					   allele1 = it->genome[j];
					   allele2 = it->genome[j+L/2];
					   
					   if(allele1+allele2 == 44){
						   hetero_b++;
					   }else if(allele2+allele1 == 42){
						   homo_b++;
					   }else if(allele2+allele1 == 46){
						   homo_B++;
					   }
			        }
			        sumfitness += it->fitness;
			        s_alleles_count->push_back(it->genome[(L/4) - 1]);
			        s_alleles_count->push_back(it->genome[(3*L/4) - 1]);
			        if(it->genome[(L/4) - 1] < 0) selfcomp++;
			        if(it->genome[(3*L/4) - 1] < 0) selfcomp++;
			        
			        
			    }
			    selfcomps.push_back(selfcomp);
			    //cout << hetero_b/selected_genes << " " << homo_b/selected_genes << " " << homo_B/selected_genes << " " << endl;
			    

			    /// Calculate S allele frequencies:
			    for(unsigned int s_al=0; s_al < Stype.size(); s_al++){
					counts_S->push_back(std::count(s_alleles_count->begin(), s_alleles_count->end(), Stype[s_al]));
					total_S_allele += (*counts_S)[s_al];
				} 
				for(unsigned int s_al=0; s_al < Stype.size(); s_al++){
					freqs_S->push_back( (double) (*counts_S)[s_al]/total_S_allele);
					//cout << "Frequency of S allele: " << Stype[s_al] << "  "  << (*freqs_S)[s_al] << endl;
					sum_quad_freq_S += (*freqs_S)[s_al]*(*freqs_S)[s_al];
					(*sum_freq_S)[s_al] += (*freqs_S)[s_al];
				}
               // Calculate average allele frequecies:
                for(vector<double>::iterator iter=locus_b->begin(); iter!=locus_b->end(); iter++){
				   sumb += *iter;
	            }
	            p_average = ((homo_b/selected_genes)+0.5*(hetero_b/selected_genes))/subpop->size();
	            for(vector<double>::iterator iter2=locus_B->begin(); iter2!=locus_B->end(); iter2++){
				   sumB += *iter2;
			    }
			    q_average = ((homo_B/selected_genes)+0.5*(hetero_b/selected_genes))/subpop->size();  // <<--- Deleterious allele frequency
			    deleterious_alleles = ((homo_B)+0.5*(hetero_b))/subpop->size();
			    //p_onelocus = (*locus_B)[2];
                //Hexp_b = 1.0 - q_average*q_average - p_average*p_average;
                //Hs_S_d = 1.0 - sum_quad_freq_S;
                //Hobs_b = hetero_b/(selgenes*subpop->size());
                inbdep = inbreeding_depression_within_deme(r, subpop, selection_exp_mu, rec, dominance);

                //cout << "Selected locus freq: " << k << " " << p_onelocus << endl;
		        //cout << "Allele freq A in patch: " << k << " " << p_average << endl;
		        //cout << "Allele freq a in patch: " << k << " "  << q_average << endl;
	            //cout << "Mean fitness: " << sumfitness/subpop->size() << endl;
	            //cout << "Genetic load: " << 1 - sumfitness/subpop->size() << endl;
	            //if(k == 0) cout << "SI population " << endl;
	            //else if(k == 1) cout << "SC population  " << endl;
	            //cout << "Inbreeding depression:  " << inbdep << endl;
	            //cout << "Effective number of S alleles: " << 1/sum_quad_freq_S << endl;
	            //cout << "Mate availability: " << mate_availability(r,subpop) << endl;
	            //cout << "Deleterious allele number: " << deleterious_alleles << " allele B "  << delet/subpop->size() << endl;


                //sum_p_onelocus += p_onelocus;
                //sum_p2_onelocus += p_onelocus*p_onelocus;
	            sum_p+= p_average;
	            sum_p2+= p_average*p_average;
	            sum_q+= q_average;
	            //sum_genAA += homo_b/(selgenes*subpop->size());
	            sumfit_total += sumfitness/subpop->size();
	            sum_inbreeding += inbreeding_depression_within_deme(r, subpop, selection_exp_mu, rec, dominance);
	            sum_Ae += 1/sum_quad_freq_S;
	            /// OUTPUT DEME STATS

				p2output[interval][k][0] = inbdep;
				p2output[interval][k][1] = 1 - sumfitness/subpop->size();
				p2output[interval][k][2] = 1/sum_quad_freq_S;    
				p2output[interval][k][3] = mate_availability(r,subpop);
				p2output[interval][k][4] = deleterious_alleles;
				p2output[interval][k][5] = selfcomp;

 		        delete subpop;
		        delete locus_b;
		        delete locus_B;
		        delete s_alleles_count;
		        delete freqs_S;
		        delete counts_S;
		        demes++;
			}
	}
	/// METAPOP AVERAGE
	number_of_demes = demes;
	//cout << "number of demes occupied: " << number_of_demes << endl;
	// Calculate Fst at each intervals
	freq_p_total = sum_p/(number_of_demes);
    //freq_q_total = sum_q/(number_of_demes);
    //freq_p_onelocus = sum_p_onelocus/number_of_demes;
    //genAA = freq_p_total*freq_p_total; //Pooled HW
    //genAa = 2*freq_p_total*freq_q_total; //Pooled HW
    //genaa = freq_q_total*freq_q_total; //Pooled HW
    //average_genA = sum_p2_onelocus/number_of_demes;//sum_genAA/number_of_demes;
    //cout << "Homozygote A total freq: " << genAA << endl;
    //cout << "Homozygote a total freq: " << genaa << endl;
    //cout << "Heterozygote a total freq: " << genAa << endl;
    //cout << "Average genA " << average_genA << endl;

    
    /*for(int patch=0; patch<patchesx; patch++){
		cout << " SC " << selfcomps[patch];
		totselfcomp += selfcomps[patch];
	}
    cout << endl;*/

	stat_year_vec[interval][0] = (double) year;
	stat_year_vec[interval][1] = freq_p_total;  /// selected allele
    stat_year_vec[interval][2] = sumfit_total/number_of_demes; /// mean fitness
    stat_year_vec[interval][3] = 1.0 - sumfit_total/number_of_demes; /// mean genetic load
    stat_year_vec[interval][4] = sum_inbreeding/number_of_demes; /// Inbreeding depression mean
    stat_year_vec[interval][5] = sum_Ae/number_of_demes; /// Mean effective allele number
    stat_year_vec[interval][6] = totselfcomp; /// Self-compatibles total metapop
	delete sum_freq_S;
	delete mean_subpop_allele_size_vec;
	delete sum_sd_allele_size_vec;
}

void final_year_stat(double **stat_year_vec, double **p2output_sims, int sim, int interval, int generations)
{
	int final_year = (generations/interval) - 1;
	p2output_sims[sim][0] = (double) sim;
	p2output_sims[sim][1] = stat_year_vec[final_year][1]; // freq p
	p2output_sims[sim][2] = stat_year_vec[final_year][2]; // Mean W
	p2output_sims[sim][3] = stat_year_vec[final_year][3]; // Load
	p2output_sims[sim][4] = stat_year_vec[final_year][4]; // Inb Dep
	p2output_sims[sim][5] = stat_year_vec[final_year][5]; // AE number S
	p2output_sims[sim][6] = stat_year_vec[final_year][6]; // Total SC metapop

}


/// Output function to quantify the number of SC alleles in the metapopulation

void self_compatibles_output(gsl_rng *r, int generation, int sim, int patchesx, double ***p2outputSC, double **p2outputSelfrate, vector<indiv> *seeds, int *space, double selection_exp_mu, double rec, double dominance){
	int subpopsize;
	double sumfit=0, sc;
	for(int i=0; i<patchesx; i++){
			if(space[i] == 1){
			  vector<indiv> *subpop = new vector<indiv>;
			  build_subpop_vector_pointer2D(i,seeds,subpop);
			  //cout << "Subpop size " << subpop->size() << " " << i << " " << j << endl;
			  sc=0,sumfit=0;
			  for(vector<indiv>::iterator iter=subpop->begin(); iter != subpop->end(); iter++){
			  	if(iter->S_allele1() < 0) sc++;
				if(iter->S_allele2() < 0) sc++;
				sumfit += iter->fitness;
			  }
			  subpopsize = subpop->size();
			  p2outputSC[generation][i][0] = sc/(2*subpopsize);
		      p2outputSC[generation][i][1] = sumfit/subpopsize;
		      p2outputSC[generation][i][2] = inbreeding_depression_within_deme(r, subpop, selection_exp_mu, rec, dominance);
		      p2outputSC[generation][i][3] = p2outputSelfrate[generation][i];
		      delete subpop;
		  }
	}
	
}



/// Create output of one simulation -  statistics for every  X generations (intervals) - to show spatial stats

void makeoutput_between_years_spatial(double ***stat, int stat_length, int patchesx, int intervals, ofstream &outfile)
{
	    outfile << "Deme\tID\tLoad\tEf-S\tMA\tDelet\tSC\n";
	    
	    // output write
	    for(int gen=0; gen<intervals; gen++){
			outfile << "\t" << gen << endl;
			for(int i=0; i<patchesx; i++){
					outfile << i << "\t";
			        for(int k=0; k<stat_length; k++){
				       outfile << format("%+4.2f") % stat[gen][i][k] << "\t";
				    }
				    outfile << endl;
			}
			outfile << endl;
		}
        
}


/// Create output of one simulation -  statistics for every  X generations - to show spatial stats of SC during expansion

void makeoutput_between_years_spatial_sc(double ***stat, int stat_length, int patchesx, int generations, ofstream &outfile)
{
	    outfile <<  "Deme\tFreqSC\tFitness\tInb-Dep\tSelf_rate\n";
	    
	    // output write
	    for(int gen=0; gen<generations; gen++){
			outfile << "\t" << gen << endl;
			for(int i=0; i<patchesx; i++){
					outfile << i << "\t";
			        for(int k=0; k<stat_length; k++){
				       outfile << format("%4.2f") % stat[gen][i][k] << "\t";//setprecision(4) << stat[gen][i][j][k] << setfill(' ') << setw(10);
				    }
				    outfile << endl;
			}
			outfile << endl;
		}
        
}


/// Create output of one simulation -  statistics for every  generations- to show spatial stats of S alleles before and during expansion

void makeoutput_between_years_spatial_s_alleles(double *****stat, int stat_length, int patchesx, int generations, vector<int> &Stype, int *space, ofstream &outfile)
{
	    outfile << "Deme\tS allele\tFreq S allele\tDelet mutations";
	    // add for loop s alleles
	    // output write
	    for(int gen=0; gen<generations; gen++){
			outfile << "Generation: " << gen << endl;
			for(int i=0; i<patchesx; i++){
					if(space[i] == 1){
					   outfile << "Deme " << i << endl;
					   for(unsigned int s_al=1; s_al<200; s_al++){
						  outfile << Stype[s_al] << "\t";
						  for(int k=0; k<stat_length; k++){
				              outfile << format("%+4.2f") % stat[gen][i][s_al][k] << "\t";
				          }
				          outfile << endl;
					   }
				       outfile << endl;
				    }
			}
			outfile << endl;
		}
        
}


/// Create output of one simulation -  statistics for every  X generations (intervals) - to show the dynamics

void makeoutput_onesim(double **stat, int stat_length, int intervals, ofstream &outfile,  double mu_selection, double mu_slocus, int generations, double rec, double dominance, double mutrates)
{	    
        //first line: column titles
        outfile << " ** Selection mean : " << mu_selection << "  " << "Mutation S locus: " << mu_slocus << "  " << "Generations: " << generations << " " << "Recombination " << rec << " Dominance " << dominance << " Mut_sel " << mutrates << endl;
	    outfile << "Year\tFreqp\tMean_W\tLoad\tInb-Dep\tAE-S\tSC\n"; // Freq of p, mean Fis, Hs, Ht, Fst, Sdiv, Mate availability
	    outfile << endl;
	    
	    // output write
	    for(int gen=0; gen<intervals; gen++){
			for(int k=0; k<stat_length; k++){
				if(k==0){
					outfile << (int) stat[gen][k] << "\t";
				}else{
					outfile << format("%+4.2f") % stat[gen][k] << "\t";
				}
			 }
			 outfile << endl;
		}

}


/// Create output of several simulations -  statistics for every simulation run at equilibrium

void makeoutput(double **stat, int stat_length, int simulations, ofstream &outfile, double mu_selection, double mu_slocus, int generations, double rec, double selfing_prob, double dominance)
{	    
        //first line: column titles
                //first line: column titles
        outfile << " Parameters:  Selection mean : " << mu_selection << "  " << "Mutation S locus: " << mu_slocus << "  " << "Generations: " << generations << " " << "Recombination " << rec << "  " << "Selfng prob " << selfing_prob << " Dominance " << dominance << endl;
	    outfile << "Run\tp-freq\tMean W\tLoad\tID-within\tAE-S\tSC"; 
	    outfile << endl;
	    
	    // output write
	    for(int sim=0; sim<simulations; sim++){
			for(int k=0; k<stat_length; k++){
				if(k==0){
					outfile << (int) stat[sim][k] << "\t";
				}else{
					outfile << setprecision(4) << stat[sim][k] << "\t";
				}
			 }
			 outfile << endl;
		}

}




//////////////////////////////////////////////////////////////////

#define RUNS 500
#define GENS 3501
#define EXPANSION 3501
#define INTERS 122
#define INTERVALS 51
#define STATS 7
#define STATSD 6
#define STATSC 4
#define DIMX 3
#define S_AL 200

//////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
	srand (time(NULL));
	
	const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, rand()%100);
    
    //Make output file
    string outfilename_year = argv[1];            
	string outfilename = argv[2];            
	ofstream outfile (outfilename.c_str());
    string outfilename_space = argv[3];            
	string outfilename_sc = argv[4];
  
	double mu_selection = atof(argv[5]);
	double dominance = atof(argv[6]);
	double mutrates = atof(argv[7]);
	double pollenLim = atof(argv[8]);
	double rec = atof(argv[9]);
	double sample1 = atoi(argv[10]);
	double sample2 = atoi(argv[11]);
	      
    
    int subpop= 100, patchesx=DIMX;
	int simulations = RUNS;    // number of simulations
	int generations = GENS;
	int interval = 100;
	int popsize_init = subpop*2;
	int genes = 500;
	int S_alleles = S_AL;
	


	double lambda_sel = genes*mutrates;
	double mu_slocus = 0.0;
	int k_si = subpop;
	int k_sc = subpop;
	int k_exp = sample1 + sample2;
	int stat_length=STATS;
	int stat_lengthD = STATSD;
	int stat_lengthSC = STATSC;
	int intervalsD = INTERVALS;
    int expansion = EXPANSION;

    //clock_t t;
	

	// Vectors:
    //Vector of S alleles
    vector<int> Stype(S_alleles);  //vector of S alleles
    for( int s_al=0; s_al<S_alleles; s_al++){
		Stype[s_al] = 100 + s_al;
        //cout << "S type: " << Stype[s_al] << endl;
	}

    // Create dyn mem pointers to mult arrays for output
    double **p2stat_year_vec; /// pointers to stats per interval ---> single simulation
    double **p2output_sim;  /// pointers to stats per simulation
    int *p2space_occup;
    double ***p2output;  /// pointers to stats per interval and deme(x,y) ---> single simulation
    double ***p2outputSC;  /// pointers to SC freq and ID per generation during expansion and deme(x,y)
    double **p2outputSelfrate;  /// pointers to Selfing rate per generation during expansion and deme(x,y) 

    
    // Allocating memeory for simulation statistics pointer
 	p2output_sim = new double*[RUNS];
	for(int run=0;  run<RUNS; ++run){
		p2output_sim[run] = new double[STATS];
	}
    
    
    /// ------------ %Simulation tuns % -------------- //
    for( int sim=0; sim<simulations; sim++){
		vector<indiv> *parents = new vector<indiv>(popsize_init); 
        
        
        //Allocate memory deme-interval pointers
        p2output = new double**[INTERVALS];
        for(int i=0;  i<INTERVALS; ++i){
			p2output[i] = new double*[DIMX];
				for(int l=0; l<DIMX; ++l)
				    p2output[i][l] = new double[STATSD];
	
		}
        
        //Allocate memory deme-expansion SC pointers
        p2outputSC = new double**[EXPANSION];
        for(int i=0;  i<EXPANSION; ++i){
			p2outputSC[i] = new double*[DIMX];
				for(int l=0; l<DIMX; ++l)
				    p2outputSC[i][l] = new double[STATSC];
		}

		
		//Allocate memory deme-expansion SC pointers
        p2outputSelfrate = new double*[EXPANSION];
        for(int i=0;  i<EXPANSION; ++i){
			p2outputSelfrate[i] = new double[DIMX];
		}
		
		//Allocate memory interval pointers
		p2stat_year_vec = new double*[INTERS];
		for(int i=0;  i<INTERS; ++i){
			p2stat_year_vec[i] = new double[STATS];
		}
 

    // Allocate memory space occupancy  --- Fill spatial occupancy matrix with zeroes
         p2space_occup = new int[DIMX];
	     for(int i=0; i<patchesx; ++i){
				 p2space_occup[i]=0;
	     }


		// Initialize individuals
        initializeplants(r, parents, genes, S_alleles, subpop,  DIMX, mu_selection, dominance, p2space_occup);
        show_vector((*parents)[0].genome);
        int inter=0, flag=0;
        //t = clock();
	// ------------  % Generation loop % ----------- //
	    for(int gen=0; gen<generations ; gen++){
			//cout << "Generation: " << gen << " of " << generations << endl;
			// Allocate memory
		    vector<indiv> *seeds = new vector<indiv>;
            //vector<vector<int> > *pt_neilist = new vector<vector<int> >;
            //vector<vector<int> > *pt_listexp = new vector<vector<int> >;  
            vector<int> *pt_chosen = new vector<int>; 

			//int popsize = parents->size();
			//cout << "Popsize SIZE: " << popsize << endl;
			
			/// 1) Gamete production (recombination, mutation) and local mating
			//t1 = clock();

			gameteprod_mating(r, parents, seeds, genes, patchesx, subpop, lambda_sel, mu_slocus, mu_selection, dominance, rec, pollenLim, p2space_occup, Stype, flag, p2outputSelfrate, gen, k_si, k_sc, k_exp);
			//cout << "Time gametoprod: " << float(clock() - t1) << endl;
			
			/// Introduction of null mutation at the S locus (introduction of Self-Compatibility)
			if (gen == 2500){
				//cout << "Experiment starts: " << mu_slocus << endl;
				experiment(r, seeds, sample1, sample2, p2space_occup);
			}
			//cout << "S locus mutation rate " << mu_slocus << endl;
			self_compatibles_output(r,gen, sim, patchesx, p2outputSC, p2outputSelfrate, seeds, p2space_occup,mu_selection, rec, dominance);
			show_sc(seeds);

		    /// Calculate statistics:
	        if(gen%interval == 0){
				 allele_frequencies_2(r,seeds, p2stat_year_vec, inter, patchesx, genes, S_alleles, Stype, gen, p2space_occup, mu_selection, rec, dominance, p2output, gen);
				 //cout << "Generation: " << gen << endl;
				 //cout << "Size: " << seeds->size() << " Selection " << mu_selection << " Dominance " << dominance << endl;
				 inter++;
		    }


            /// Break when expansion reaches the boundary
			//if(check_expansion(seeds,patchesx,p2space_occup) == true)  break;
			
			parents->swap(*seeds);
			
			// Release memory
	        delete seeds;
            delete pt_chosen;
           
		} /// END OF GENERATION LOOP
		//cout << "It took one generation: " << clock() - t << endl; 
		delete parents;
		
       /// Output spatial stats 
		stringstream extension_file;//create a stringstream
		extension_file << sim;//add number to the stream
		string outfilespacefull1 = outfilename_space + extension_file.str();
		string outfilespacefull1_year = outfilename_year + extension_file.str();
		string outfilespacefull1_sc= outfilename_sc + extension_file.str();
		ofstream outfilespacefull2 (outfilespacefull1.c_str());
		ofstream outfilespacefull2_year (outfilespacefull1_year.c_str());
		ofstream outfilespacefull2_sc (outfilespacefull1_sc.c_str());
		makeoutput_between_years_spatial(p2output, stat_lengthD, patchesx, intervalsD, outfilespacefull2);
		makeoutput_onesim(p2stat_year_vec, stat_length, inter, outfilespacefull2_year, mu_selection, mu_slocus, generations, rec, dominance, mutrates);  
		makeoutput_between_years_spatial_sc(p2outputSC, stat_lengthSC, patchesx, expansion, outfilespacefull2_sc);
		//selection_coef_distribution(p2output_sel, outfilespacefull2_sel);
		
		// Put final year statistics (equilibrium) at double pointer p2output_sim
		final_year_stat(p2stat_year_vec, p2output_sim, sim, interval, generations);
		
		/// Deallocating memory:
        // Deallocating memory pointers stats-interval
        for(int i=0; i<INTERS; i++)
			delete [] p2stat_year_vec[i];
		delete [] p2stat_year_vec;
		
		//Deallocate memory pointers space-occupancy (DIMX)
        delete [] p2space_occup;
        
        
        // Deallocating memory pointers stats-deme-interval
		for (int i = 0; i < INTERVALS; ++i){
			for (int j = 0; j < DIMX; ++j){
				delete [] p2output[i][j];
		     }
			 delete [] p2output[i];
		}	    
        delete [] p2output;
        
        // Deallocating memory pointers generations(expansion)-demes
		for (int i=0; i < EXPANSION; ++i){
			for (int j=0; j < DIMX; ++j){
				delete [] p2outputSC[i][j];
		     }
			 delete [] p2outputSC[i];
		}	    
        delete [] p2outputSC;
        
        // Deallocating memory pointers (p2outputSelfrate) generations(expansion)-demes
		for (int i=0; i < EXPANSION; ++i){
			 delete [] p2outputSelfrate[i];
		}	    
        delete [] p2outputSelfrate;
        
      
	}/// END OF SIMULATION LOOP

    // Write output for all simulations
    makeoutput(p2output_sim, stat_length, simulations, outfile,  mu_selection, mu_slocus, generations, rec, pollenLim, dominance);
    
    // Deallocate memory pointers stats-simulation
    for(int run=0; run<RUNS; run++)
		delete [] p2output_sim[run];
	delete [] p2output_sim;



	gsl_rng_free (r);
	return 0;

}
