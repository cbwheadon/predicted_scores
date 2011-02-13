#!/usr/bin/r
# Function requires two parameters as matrices
# The first parameter is the probability of a correct response to each category with each item on a separate row
# The second parameter is a single column of categories per item
# i.e. A dichotomous item has two categories
# The function returns the predicted conditional number-correct score distribution 

wLord <- cxxfunction(signature(X = "numeric", Y = "numeric") ,
		"
		int *dimX, *dimY;
  		double *xptr, *yptr;
  		dimX = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP));
  		PROTECT(X = coerceVector(X, REALSXP));
  		xptr = REAL(X);
  		dimY = INTEGER(coerceVector(getAttrib(Y, R_DimSymbol), INTSXP));
  		PROTECT(Y = coerceVector(Y, REALSXP));
  		yptr = REAL(Y);
		
		/* Do some calculations */
  		int nrX, ncX, nrY, ncY;
  		nrX = dimX[0];   ncX = dimX[1];
  		nrY = dimY[0];   ncY = dimY[1];		

		int cats[nrY];
  		int k = 0;
  		for(int i=0;i<nrY;i++){
			cats[i]=yptr[i];
			k = k + cats[i];
  		}

		const int scrpts = k;
  		const int ret = k - nrX + 1;
		
		/* Create SEXP to hold the answer */
		SEXP ans;
  		double *ansptr;
  		PROTECT(ans = allocMatrix(REALSXP, ret, 1));
  		ansptr = REAL(ans);

		int min_c=1;
  		int max_c=0;
  		int prev_mx =0;
  		int prev_min=0;
  		double frx=0;
  		const size_t arr_size = nrX;

  		int *cat;
  		cat = &cats[0];
  		double *pia;
  		double *ppa;
  		double *pla;
  		double *ppia;
  		double *ppla;

		double pa[arr_size][scrpts];
 		double ia[nrX][ncX];

  		for(int i=0;i<nrX;i++){
			for(int j=0;j<ncX;j++){
				ia[i][j]=xptr[i+(j*nrX)];
			}
  		}

		//loop through probs
			    int itm = 0;
			    for (int *ibegin=cats, *iend=ibegin+arr_size; ibegin!=iend; ++ibegin){
			        if(itm==0){
			                   //loop through cats
			                   pia = &ia[0][0];
			                   ppa = &pa[0][0];
			                   for(int i=0; i!=*ibegin;i++){
			                           *ppa = *pia;
			                           pia++;
									   ppa++;
			                           max_c++;
			                   }
			        } else {
			                   pia = &ia[itm][0];
			                   pla = &pa[itm-1][0];
			                   ppa = &pa[itm][0];
		
				               //loop through potential score points
			                   //potential score max
			                  	prev_mx = max_c;
			                   	prev_min = min_c;
			                    max_c = max_c + *cat;
			                    min_c = min_c + 1;
				                int l = 0;
		
					            for (int j=min_c;j!=max_c+1;j++){
			                        frx = 0;
			                    	//loop through each category
			                       	for (int k = 0;k!=*cat;++k){
			                            //calculate x-Wjk
			                    		if(l-k >= 0 && l-k < 1 + prev_mx-prev_min) {
											ppia= &ia[itm][0];
											ppia= ppia + k;
			                                //p.last
											ppla= &pa[itm-1][0];
			                                ppla= ppla + (l-k);
			                                frx = frx + (*ppia * *ppla);
			         					}
		
			          				}
			                    	*ppa = frx;
			                       	ppa++;
			                    	l++;
			                   	}
			              }
			              itm++;
			              cat++;
		
			    }
		
				ppa = &pa[itm-1][0];
		
				for (int k=0;k!=1 + max_c - min_c;k++){
		        	ansptr[k] = *ppa;
		        	ppa++;
		    	}
		
		
		  /* Wrap up and return the result to R */

		UNPROTECT(3);
		return ans;
		",verbose=TRUE)